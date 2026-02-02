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
params.snpeff_db         = getGenomeAttribute('snpeff_db')
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')


include { RNAVAR_MINI                          } from './workflows/rnavar_mini'
include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { PREPARE_GENOME                  } from './subworkflows/local/prepare_genome'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'


include { RUN_STRINGTIE } from './workflows/stringtie/stringtie'
include { GENERATE_VARIANT_DB         } from './workflows/variant_db'
include { GENERATE_NOVEL_ISOFORM_DB   } from './workflows/novel_isoform_db'
include { COMBINE_PROTEIN_DBS         } from './workflows/combine_db'



workflow NFCORE_RNAVAR {
  take:
    samplesheet
    align

  main:
    ch_versions = Channel.empty()

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
      align
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

    if (params.download_cache) {
      ensemblvep_info = Channel.of([ [ id: "${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
      snpeff_info     = Channel.of([ [ id: "${params.snpeff_db}" ], params.snpeff_db ])

      DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)

      snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
      vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map { meta, cache -> [ cache ] }

      ch_versions  = ch_versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
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
      params.snpeff_db,
      params.vep_genome,
      params.vep_species,
      params.vep_cache_version,
      params.vep_include_fasta,
      vep_cache,
      vep_extra_files,
      seq_center,
      seq_platform
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

        def ch_annotated_vcf = rnavar_run.annotated_vcf
        def ch_markdup_bams = rnavar_run.markdup_bams
        def vcf_done = ch_annotated_vcf.collect()

        GENERATE_VARIANT_DB(
            ch_annotated_vcf
        )

        variant_fasta = GENERATE_VARIANT_DB.out.variant_db
            .map { meta, fa -> tuple(meta.id, fa) }                // tuple(id, var_modified_peptides.fa)
        

        // 2. StringTie on markdup BAMs from RNAVAR
        RUN_STRINGTIE(
            ch_markdup_bams,
            gene_annotation_gtf
        )

        // 3. Generate novel isoform DBs from StringTie GTFs
        GENERATE_NOVEL_ISOFORM_DB(
            RUN_STRINGTIE.out.stringtie_gtf   // emits: [ val(meta), path(gtf) ]
        )
        isoform_fasta = GENERATE_NOVEL_ISOFORM_DB.out.isoform_db
        def isoform_fasta_gated = vcf_done
            .map { null }
            .concat(isoform_fasta)
            .filter { it != null }
         

        // Step 3: Combine protein DBs
        combine_in_ch = variant_fasta
            .combine(isoform_fasta_gated, by: 0)
            .map { id, var_fa, novel_fa ->
                tuple(id, var_fa, novel_fa)
            }

        COMBINE_PROTEIN_DBS(combine_in_ch)

    emit:
        combined_db_ch = COMBINE_PROTEIN_DBS.out.combined_db
        novel_db_ch    = COMBINE_PROTEIN_DBS.out.novel_db
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
