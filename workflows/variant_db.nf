include { GFFREAD } from '../modules/nf-core/gffread/main' 
include { sampleTissue; cohortTissue } from './tissue_metadata'

workflow GENERATE_VARIANT_DB {
    
    take:
    annotated_vcf   // channel: [ val(meta), path(vcf) ]  (meta.id is sample id)

    main:
    def ref_fasta = params.fasta ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].fasta : null)
    def ref_gtf   = params.gtf   ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].gtf   : null)
    
    def ch_protein_ref = Channel.value(
          file(params.protein_reference_db, checkIfExists: true)
        )
    
    if (!ref_fasta) {
        error("Missing required reference fasta: set --fasta or provide a genome entry with a fasta path.")
    }
    if (!ref_gtf) {
        error("Missing required reference GTF: set --gtf or provide a genome entry with a gtf path.")
    }

    def ch_ref_fasta = Channel.value( file(ref_fasta, checkIfExists:true) )

    gffread_in_ch = Channel.value( tuple([id: 'ref'], file(ref_gtf, checkIfExists: true)) )

    GFFREAD(
      gffread_in_ch,
      ch_ref_fasta)
    
    def ch_input_fa = GFFREAD.out.gffread_fasta
      .map { _meta, fa -> fa }

    // 1) generate variant peptide fasta
    var_peptides_ch = GENERATE_VAR_PEPTIDES(annotated_vcf, ch_input_fa, ref_gtf)

    // 2) add AA-change annotation / modify peptides
    mod_peptides_ch = ANNOTATE_VAR_PEPTIDES(var_peptides_ch, ch_protein_ref)

    // Variant calling remains run-level, but the manual workflow merges every
    // run-level FASTA into one cohort database before final DB assembly.
    def expected_runs = (params.cohort_expected_runs ?: 27) as int
    if (expected_runs < 1) {
      error("Invalid --cohort_expected_runs '${params.cohort_expected_runs}': expected a positive integer.")
    }

    def cohort_mod_peptides_ch = mod_peptides_ch
      .map { meta, fa -> tuple(params.cohort_db_id, meta, fa) }

    // groupTuple(size:) bounds each cohort group, while the count channel
    // prevents a complete group from being merged until the source closes and
    // confirms that no run is missing and no extra run arrived.
    def checked_cohort_count_ch = cohort_mod_peptides_ch
      .count()
      .map { observed_runs ->
        if (observed_runs != expected_runs) {
          error("Cohort '${params.cohort_db_id}' expected exactly ${expected_runs} run-level variant FASTAs, but received ${observed_runs}.")
        }
        tuple(params.cohort_db_id, observed_runs)
      }

    def grouped_mod_peptides_ch = cohort_mod_peptides_ch
      .groupTuple(size: expected_runs)
      .combine(checked_cohort_count_ch)
      .map { db_id, metas, fas, checked_db_id, observed_runs ->
        if (db_id != checked_db_id || metas.size() != expected_runs || fas.size() != expected_runs) {
          error("Cohort aggregation mismatch for '${db_id}': expected ${expected_runs} entries, received ${metas.size()} metadata records and ${fas.size()} FASTAs.")
        }

        def run_ids = metas.collect { run_meta -> run_meta.rnaseq_run ?: run_meta.id }
        def duplicate_run_ids = run_ids
          .countBy { run_id -> run_id }
          .findAll { _run_id, count -> count > 1 }
          .keySet()
          .sort()
        def unique_run_ids = run_ids.unique()
        if (unique_run_ids.size() != expected_runs) {
          error("Cohort '${db_id}' expected ${expected_runs} unique RNA-seq runs, but received ${unique_run_ids.size()}. Duplicate run IDs: ${duplicate_run_ids ?: 'none'}.")
        }

        def first_meta = metas[0]
        def dia_raw = metas.collect { run_meta -> run_meta.dia_raw }.find { raw_files -> raw_files }
        def meta = first_meta + [
          id           : db_id,
          db_id        : db_id,
          rnaseq_runs  : unique_run_ids,
          tissue       : cohortTissue(metas),
          expected_runs: expected_runs,
          observed_runs: observed_runs,
          dia_raw      : dia_raw ?: (first_meta.dia_raw ?: [])
        ]
        tuple(meta, fas)
      }

    MERGE_VARIANT_PEPTIDES(grouped_mod_peptides_ch)

    // Keep the full deduplicated FASTA as an intermediate, but use compact
    // headers for the variant DB passed into final protein DB assembly. Full
    // run/tissue provenance is retained in a companion TSV keyed by protein_id.
    COMPACT_VARIANT_FASTA(MERGE_VARIANT_PEPTIDES.out.merged_variant_db)

    emit:
    variant_db         = COMPACT_VARIANT_FASTA.out.compact_variant_db
    variant_provenance = COMPACT_VARIANT_FASTA.out.provenance
}

process MERGE_VARIANT_PEPTIDES {
  tag "${meta.id}"

  input:
    tuple val(meta), path(var_modified_peptides)

  output:
    tuple val(meta), path("${meta.id}_var_modified_peptides.fa"), emit: merged_variant_db

  script:
    def input_files = var_modified_peptides instanceof List ? var_modified_peptides : [var_modified_peptides]
    def expected_runs = (params.cohort_expected_runs ?: 27) as int
    if (input_files.size() != expected_runs) {
      error("Cohort '${meta.id}' merge expected ${expected_runs} FASTAs, but received ${input_files.size()}.")
    }
    def inputs = input_files.join(' ')
    def min_samples = params.variant_peptide_min_runs ?: 1
    """
    set -euo pipefail

    output_fasta="${meta.id}_var_modified_peptides.fa"

    python3 ${projectDir}/bin/deduplicate_fasta.py \\
      --output "\$output_fasta" \\
      --dedup-by sequence \\
      --min-samples ${min_samples} \\
      $inputs

    if [[ ! -s "\$output_fasta" ]] || ! grep -q '^>' "\$output_fasta"; then
      echo "ERROR: Cohort '${meta.id}' produced an empty or invalid variant FASTA: \$output_fasta" >&2
      exit 1
    fi
    """
}

process COMPACT_VARIANT_FASTA {
  tag "${meta.id}"

  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/python:3.8.3' :
      'quay.io/biocontainers/python:3.8.3' }"

  input:
    tuple val(meta), path(var_modified_peptides)

  output:
    tuple val(meta), path("${meta.id}_var_modified_peptides.compact.fa"), emit: compact_variant_db
    tuple val(meta), path("${meta.id}_var_modified_peptides.provenance.tsv"), emit: provenance

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/compact_variant_fasta.py \
      "${var_modified_peptides}" \
      "${meta.id}_var_modified_peptides.compact.fa" \
      "${meta.id}_var_modified_peptides.provenance.tsv"

    if [[ ! -s "${meta.id}_var_modified_peptides.compact.fa" ]] || \
       ! grep -q '^>' "${meta.id}_var_modified_peptides.compact.fa"; then
      echo "ERROR: Compact variant FASTA is empty or invalid for '${meta.id}'." >&2
      exit 1
    fi

    if [[ ! -s "${meta.id}_var_modified_peptides.provenance.tsv" ]]; then
      echo "ERROR: Variant provenance TSV is empty for '${meta.id}'." >&2
      exit 1
    fi
    """
}


process GENERATE_VAR_PEPTIDES {
  tag "${meta.id}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pypgatk:0.0.24--pyhdfd78af_0' :
      'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0' }"

  input:
    tuple val(meta), path(annotated_vcf)
    path ref_fasta // transcript_fasta
    path ref_gtf

  output:
    tuple val(meta), path(annotated_vcf), path("${meta.id}_var_peptides.fa")

  script:
    def using_snpeff = !params.skip_variantannotation && params.tools && params.tools.split(',').contains('snpeff')
    def use_consequence_filter = !params.skip_variantannotation && params.tools && params.tools != 'none'
    def consequence_args = use_consequence_filter ? '--include_consequences missense_variant,frameshift_variant,inframe_insertion,inframe_deletion' : ''
    def pass_filter = params.variant_peptide_require_pass ? """
    awk 'BEGIN{FS=OFS="\\t"} /^#/ || \$7=="PASS" {print}' "${meta.id}.pypgatk_input.vcf" > "${meta.id}.pypgatk_input.pass.vcf"
    mv "${meta.id}.pypgatk_input.pass.vcf" "${meta.id}.pypgatk_input.vcf"
    """ : ''
    def prepare_vcf = using_snpeff ? """
    python3 ${projectDir}/bin/snpeff_ann_to_vep_csq.py \\
      --input "${annotated_vcf}" \\
      --output "${meta.id}.pypgatk_input.vcf"
    """ : """
    ln -sf "${annotated_vcf}" "${meta.id}.pypgatk_input.vcf"
    """
    """
    set -euo pipefail

    ${prepare_vcf}
    ${pass_filter}

    python3 -m pypgatk.pypgatk_cli vcf-to-proteindb \
      --vcf "${meta.id}.pypgatk_input.vcf" \
      --input_fasta "${ref_fasta}" \
      --gene_annotations_gtf "${ref_gtf}" \
      ${consequence_args} \
      --af_field AF \
      --af_threshold 0.05 \
      -o ${meta.id}_var_peptides.fa
    """
}

 
process ANNOTATE_VAR_PEPTIDES {
  tag "${meta.id}"
  
  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"


  input:
    tuple val(meta), path(annotated_vcf), path(var_peptides)
    path protein_reference_db

  output:
    tuple val(meta), path("${meta.id}_var_modified_peptides.fa")

  script:
    def tissue = sampleTissue(meta)
    """
    set -euo pipefail

    python3 ${projectDir}/bin/get_var_aa_change.py \
      "${annotated_vcf}" \
      "${var_peptides}" \
      "${protein_reference_db}" \
      "${meta.id}_var_modified_peptides.fa" \
      --sample-id "${meta.id}" \
      --tissue "${tissue}" \
      --organism "${params.organism}" \
      --taxonomy-id "${params.taxonomy_id}"
    """
}
