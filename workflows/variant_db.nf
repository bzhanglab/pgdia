workflow GENERATE_VARIANT_DB {
    
    take:
    annotated_vcf   // channel: [ val(meta), path(vcf) ]  (meta.id is sample id)

    main:
    def ref_fasta = params.fasta ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].fasta : null)
    def ref_gtf   = params.gtf   ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].gtf   : null)

    if (!ref_fasta) {
        error("Missing required reference fasta: set --fasta or provide a genome entry with a fasta path.")
    }
    if (!ref_gtf) {
        error("Missing required reference GTF: set --gtf or provide a genome entry with a gtf path.")
    }

    // 1) generate variant peptide fasta
    var_peptides_ch = gen_var_db(annotated_vcf, ref_fasta, ref_gtf)

    // 2) add AA-change annotation / modify peptides
    mod_peptides_ch = mod_var_peptides(var_peptides_ch)

    emit:
    variant_db = mod_peptides_ch
}


process gen_var_db {
  tag {meta.id}
  container 'python:3.11-slim'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pypgatk:0.0.24--pyhdfd78af_0' :
      'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0' }"

  input:
    tuple val(meta), path(annotated_vcf)
    path ref_fasta
    path ref_gtf

  output:
    tuple val(meta), path(annotated_vcf), path("${meta.id}_var_peptides.fa")

  script:
    """
    set -euo pipefail

    python3 -m pypgatk.pypgatk_cli vcf-to-proteindb \
      --vcf "${annotated_vcf}" \
      --input_fasta "${ref_fasta}" \
      --gene_annotations_gtf "${ref_gtf}" \
      --include_consequences missense_variant,frameshift_variant,inframe_insertion,inframe_deletion \
      --af_field AF \
      --af_threshold 0.05 \
      -o ${meta.id}_var_peptides.fa
    """
}

 
process mod_var_peptides {
  tag {meta.id}
  container 'python:3.11-slim'

  input:
    tuple val(meta), path(annotated_vcf), path(var_peptides)

  output:
    tuple val(meta), path("${meta.id}_var_modified_peptides.fa")

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/get_var_aa_change.py \
      "${annotated_vcf}" \
      "${var_peptides}" \
      "${meta.id}_var_modified_peptides.fa"
    """
}

