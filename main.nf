nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.dict              = getGenomeAttribute('dict')
params.gtf               = getGenomeAttribute('gtf')
params.gff               = getGenomeAttribute('gff')
params.exon_bed          = getGenomeAttribute('exon_bed')
params.star_index        = getGenomeAttribute('star')
params.dbsnp             = getGenomeAttribute('dbsnp')
params.dbsnp_tbi         = getGenomeAttribute('dbsnp_tbi')
params.known_indels      = getGenomeAttribute('known_indels')
params.known_indels_tbi  = getGenomeAttribute('known_indels_tbi')
params.snpeff_db         = getGenomeAttribute('snpeff_db') ?: params.snpeff_db
params.snpeff_genome_name = getGenomeAttribute('snpeff_genome_name') ?: params.snpeff_genome_name ?: params.snpeff_db
params.vep_cache         = getGenomeAttribute('vep_cache') ?: params.vep_cache
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')


include { RNAVAR_MINI                          } from './workflows/rnavar_mini'
include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { PREPARE_GENOME                  } from './subworkflows/local/prepare_genome'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'
include { ANNOTATION_CACHE_INITIALISATION } from './subworkflows/local/annotation_cache_initialisation'


include { RUN_STRINGTIE } from './workflows/stringtie/stringtie'
include { GENERATE_VARIANT_DB         } from './workflows/variant_db'
include { GENERATE_NOVEL_ISOFORM_DB   } from './workflows/novel_isoform_db'
include { COMBINE_PROTEIN_DBS         } from './workflows/combine_db'


process BUILD_SNPEFF_REFSEQ_CACHE {
    tag "$snpeff_db"
    label 'process_high_memory'

    conda "${projectDir}/modules/nf-core/snpeff/snpeff/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
        'biocontainers/snpeff:5.1--hdfd78af_2' }"

    input:
    tuple val(fasta_meta), path(fasta)
    tuple val(gtf_meta), path(gtf)
    val snpeff_db
    val genome_name

    output:
    tuple val([ id: snpeff_db ]), path("snpeff_cache"), emit: cache
    path "versions.yml", emit: versions

    script:
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 6144
    """
    set -euo pipefail

    mkdir -p "snpeff_cache/${snpeff_db}"
    cp -L "${fasta}" "snpeff_cache/${snpeff_db}/sequences.fa"
    cp -L "${gtf}" "snpeff_cache/${snpeff_db}/genes.gtf"

    cat > snpeff_cache/snpEff.config <<-END_SNPEFF_CONFIG
    data.dir = \${PWD}/snpeff_cache
    ${snpeff_db}.genome : ${genome_name}
    END_SNPEFF_CONFIG

    SNPEFF_JAR="\$(find /usr/local/share /opt/conda/share -path '*/snpEff.jar' -print -quit)"
    if [[ -z "\${SNPEFF_JAR}" ]]; then
        echo "Could not locate snpEff.jar" >&2
        exit 1
    fi

    java \\
        -Xmx${avail_mem}M \\
        -jar "\${SNPEFF_JAR}" \\
        build \\
        -c "\${PWD}/snpeff_cache/snpEff.config" \\
        -dataDir "\${PWD}/snpeff_cache" \\
        -gtf22 \\
        -noCheckCds \\
        -noCheckProtein \\
        -v \\
        -noLog \\
        "${snpeff_db}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(java -jar "\${SNPEFF_JAR}" -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}

process CHECK_READ_LENGTH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads, stageAs: 'reads/*')

    output:
    tuple val(meta), path("${meta.id}.read_length.txt"), path("${meta.id}.read_length.tsv"), emit: result
    path "versions.yml", emit: versions

    script:
    def expected_arg = params.read_length ? "--expected-length ${params.read_length}" : ''
    def read_inputs = reads instanceof List ? reads.join(' ') : reads
    """
    set -euo pipefail

    python3 ${projectDir}/bin/check_read_length.py \\
      --max-reads ${params.read_length_check_reads} \\
      ${expected_arg} \\
      --length-output "${meta.id}.read_length.txt" \\
      --report "${meta.id}.read_length.tsv" \\
      ${read_inputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}



workflow NFCORE_RNAVAR {
  take:
    samplesheet
    align

  main:
    ch_versions = Channel.empty()

    // Inspect the first available FASTQ before STAR. Its observed read length
    // is used for the cohort STAR index and every alignment. If --read_length
    // is supplied, CHECK_READ_LENGTH also verifies it against this FASTQ.
    def ch_read_length_input = samplesheet
      .map { row_id, meta, fastq_1, fastq_2, bam, bai, cram, crai, vcf, tbi, dia_raw ->
        tuple(meta, [fastq_1].findAll { it })
      }
      .filter { meta, reads -> !reads.isEmpty() }
      .first()

    CHECK_READ_LENGTH(ch_read_length_input)
    ch_versions = ch_versions.mix(CHECK_READ_LENGTH.out.versions)

    def ch_read_length = CHECK_READ_LENGTH.out.result
      .map { meta, length_file, report -> length_file }
      .collect()
      .map { length_files ->
        if (!length_files) {
          return params.read_length ?: 0
        }
        def lengths = length_files.collect { it.text.trim().toInteger() }
        def detected = lengths.max()
        log.info("Detected FASTQ read length: ${detected} bp (checked the first fastq_1 entry in the input samplesheet)")
        return detected
      }

    /*
     * Basic parameter checks (keep only what RNAVAR really needs)
     */
    if( params.gtf && params.gff ) {
      error("Using both --gtf and --gff is not supported. Please use only one of these parameters")
    } else if( !params.gtf && !params.gff ) {
      error("Missing required parameters: --gtf or --gff")
    }

    if( params.extract_umi && !params.umitools_bc_pattern ) {
      error("Expected --umitools_bc_pattern when --extract_umi is specified.")
    }

    if( !params.skip_baserecalibration && !params.dbsnp && !params.known_indels ) {
      error("Known sites are required for performing base recalibration. Supply them with either --dbsnp and/or --known_indels or disable base recalibration with --skip_baserecalibration")
    }

    /*
     * Reference channels (raw)
     * Keep these in the same "collect into a singleton list" style as nf-core
     */
    ch_fasta_raw      = params.fasta      ? Channel.fromPath(params.fasta).map{ it -> [ [id: it.baseName], it ] }.collect() : Channel.empty()
    ch_dict_raw       = params.dict       ? Channel.fromPath(params.dict).map{ it -> [ [id: it.baseName], it ] }.collect()  : Channel.empty()
    ch_fai_raw        = params.fasta_fai  ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id: it.baseName], it ] }.collect() : Channel.empty()

    ch_dbsnp_raw      = params.dbsnp      ? Channel.fromPath(params.dbsnp).map{ f -> [ [id: f.baseName], f ] }.collect() : Channel.value([])
    ch_known_indels_raw     = params.known_indels     ? Channel.fromPath(params.known_indels)     : Channel.empty()
    ch_known_indels_tbi_raw = params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi) : Channel.empty()

    ch_gff_raw        = params.gff        ? Channel.fromPath(params.gff).map{ it -> [ [id: it.baseName], it ] }.collect() : Channel.empty()
    ch_gtf_raw        = params.gtf        ? Channel.fromPath(params.gtf).map{ it -> [ [id: it.baseName], it ] }.collect() : Channel.empty()

    ch_star_index_raw = params.star_index ? Channel.fromPath(params.star_index).map{ it -> [ [id: it.baseName], it ] }.collect()
                                         : Channel.value([[],[]])

    ch_exon_bed_raw   = params.exon_bed   ? Channel.fromPath(params.exon_bed).map{ it -> [ [id: it.baseName], it ] }.collect()
                                         : Channel.empty()

    /*
     * VEP (VEP-only or VEP+others depends on your RNAVAR implementation)
     * Keep as plain Groovy values where RNAVAR expects values (like nf-core does).
     */
    seq_platform = params.seq_platform ?: []
    seq_center   = params.seq_center   ?: []

    vep_extra_files = []
    if (params.dbnsfp && params.dbnsfp_tbi) {
      vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
      vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }
    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
      vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
      vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
      vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
      vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    /*
     * Prepare genome (creates the correctly-shaped channels required by the modules)
     */
    PREPARE_GENOME(
      ch_fasta_raw,
      ch_dict_raw,
      ch_fai_raw,
      ch_star_index_raw,
      ch_gff_raw,
      ch_gtf_raw,
      ch_exon_bed_raw,
      ch_dbsnp_raw,
      ch_known_indels_raw,
      ch_known_indels_tbi_raw,
      params.feature_type,
      align,
      ch_read_length
    )

    ch_fasta            = PREPARE_GENOME.out.fasta
    ch_star_index       = PREPARE_GENOME.out.star_index
    ch_gtf              = PREPARE_GENOME.out.gtf
    ch_dict             = PREPARE_GENOME.out.dict
    ch_fasta_fai        = PREPARE_GENOME.out.fasta_fai
    ch_exon_bed         = PREPARE_GENOME.out.exon_bed
    ch_dbsnp            = params.dbsnp && params.dbsnp.endsWith(".gz") ? ch_dbsnp_raw : PREPARE_GENOME.out.dbsnp

    ch_dbsnp_tbi        = (params.dbsnp?.toString()?.endsWith(".gz") && params.dbsnp_tbi)
                          ? Channel.fromPath(params.dbsnp_tbi).map{ f -> [ [id: f.baseName], f ] }.collect()
                          : PREPARE_GENOME.out.dbsnp_tbi

    ch_known_indels     = params.known_indels ? PREPARE_GENOME.out.known_indels     : Channel.value([])
    ch_known_indels_tbi = params.known_indels ? PREPARE_GENOME.out.known_indels_tbi : Channel.value([])

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    /*
     * Cache init (keep your existing logic if you need it; otherwise, pass empty cache channels/values)
     */
    def snpeff_cache = Channel.empty()
    def vep_cache    = Channel.empty()

    def snpeff_requested = params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))
    def selected_genome_attrs = (params.genomes && params.genome && params.genomes.containsKey(params.genome)) ? params.genomes[params.genome] : [:]
    def snpeff_db_value = params.snpeff_db ?: selected_genome_attrs.snpeff_db
    def snpeff_genome_name_value = params.snpeff_genome_name ?: selected_genome_attrs.snpeff_genome_name ?: snpeff_db_value

    if (params.build_snpeff_refseq_cache && !params.skip_variantannotation && snpeff_requested) {
      if (!snpeff_db_value) {
        error("Missing required parameter for local snpEff cache build: --snpeff_db")
      }
      if (!params.gtf) {
        error("Missing required parameter for local snpEff cache build: --gtf or genome GTF")
      }
      if (!params.fasta) {
        error("Missing required parameter for local snpEff cache build: --fasta or genome FASTA")
      }

      def snpeff_fasta_for_build = Channel.value(tuple([ id: 'snpeff_fasta' ], file(params.fasta, checkIfExists: true)))
      def snpeff_gtf_for_build   = Channel.value(tuple([ id: 'snpeff_gtf' ], file(params.gtf, checkIfExists: true)))

      BUILD_SNPEFF_REFSEQ_CACHE(
        snpeff_fasta_for_build,
        snpeff_gtf_for_build,
        snpeff_db_value,
        snpeff_genome_name_value
      )

      snpeff_cache = BUILD_SNPEFF_REFSEQ_CACHE.out.cache.first()
      ch_versions  = ch_versions.mix(BUILD_SNPEFF_REFSEQ_CACHE.out.versions)

    } else if (params.download_cache && !params.skip_variantannotation) {
      ensemblvep_info = Channel.of([ [ id: "${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
      snpeff_info     = Channel.of([ [ id: "${snpeff_db_value}" ], snpeff_db_value ])

      DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)

      snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
      vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map { meta, cache -> [ cache ] }

      ch_versions  = ch_versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    } else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (!params.skip_variantannotation && params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            snpeff_db_value,
            (!params.skip_variantannotation && params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_custom_args,
            "Please refer to https://nf-co.re/rnavar/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

            snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
            vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    /*
     * Call your modified RNAVAR workflow
     * (Your RNAVAR should emit markdup_bams and annotated_vcf.)
     */
    RNAVAR_MINI(
      samplesheet,
      ch_dbsnp,
      ch_dbsnp_tbi,
      ch_dict,
      ch_exon_bed,
      ch_fasta,
      ch_fasta_fai,
      ch_gtf,
      ch_known_indels,
      ch_known_indels_tbi,
      ch_star_index,
      snpeff_cache,
      snpeff_db_value,
      params.vep_genome,
      params.vep_species,
      params.vep_cache_version,
      params.vep_include_fasta,
      vep_cache,
      vep_extra_files,
      seq_center,
      seq_platform,
      ch_read_length
    )

    ch_versions = ch_versions.mix(RNAVAR_MINI.out.versions)

  emit:
    markdup_bams = RNAVAR_MINI.out.markdup_bams
    annotated_vcf = RNAVAR_MINI.out.annotated_vcf

    multiqc_report = RNAVAR_MINI.out.multiqc_report
    versions       = ch_versions
}



workflow PGDIA {

    take:
        samplesheet
        align
        gene_annotation_gtf

    main:

      // 1. Run nf-core/rnavar
      def rnavar_run = NFCORE_RNAVAR(
          samplesheet,
          align
      )

      def ch_markdup_bams = rnavar_run.markdup_bams
      def ch_annotated_vcf = rnavar_run.annotated_vcf

      def ch_vcf_for_var = ch_annotated_vcf.map { item ->
          tuple(item[0], item[1])
      }

      GENERATE_VARIANT_DB(ch_vcf_for_var)

      variant_fasta = GENERATE_VARIANT_DB.out.variant_db
          .map { meta, fa -> tuple(meta, fa) }                // tuple(val(meta), var_modified_peptides.compact.fa)
      

      // 2. StringTie on markdup BAMs from RNAVAR
      RUN_STRINGTIE(
          ch_markdup_bams,
          gene_annotation_gtf
      )

      // 3. Generate novel isoform DBs from StringTie GTFs
      def isoform_results = GENERATE_NOVEL_ISOFORM_DB(
          RUN_STRINGTIE.out.stringtie_gtf   // emits: [ val(meta), path(gtf) ]
      )
      isoform_fasta = isoform_results.isoform_db
      isoform_tmap = isoform_results.isoform_tmap

      //def isoform_fasta_gated = isoform_fasta.after(vcf_done)
      
      //vcf_done.view { "VCF_DONE=$it ; class=${it.getClass()}" }
      //isoform_fasta_gated.view { "ISOFORM_GATED = $it ; class=${it.getClass()}" }


      // Step 3: Combine protein DBs
      // StringTie is run once on the merged RNA-seq BAM, so reuse the single
      // novel isoform DB for each per-sample variant peptide DB.
      combine_in_ch = variant_fasta
          .combine(isoform_fasta)
          .map { variant_meta, var_fa, isoform_meta, novel_fa ->
              tuple(variant_meta.subMap(['id','dia_raw']), var_fa, novel_fa)
          }
      
      def ch_protein_ref = Channel.value(
        file(params.protein_reference_db, checkIfExists: true)
      )

      def dbs = COMBINE_PROTEIN_DBS(combine_in_ch, ch_protein_ref)

      def combined_db_ch = dbs.combined_db
      def novel_db_ch = dbs.novel_db

    emit:
      combined_db = combined_db_ch
      novel_db = novel_db_ch
      isoform_tmap = isoform_tmap
}

workflow {
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir
    )
    PGDIA(PIPELINE_INITIALISATION.out.samplesheet, PIPELINE_INITIALISATION.out.align, file(params.gtf, checkIfExists: true))
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}
