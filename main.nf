nextflow.enable.dsl=2

include { NFCORE_RNAVAR } from './submodules/rnavar'

include { RUN_STRINGTIE } from './workflows/stringtie/stringtie'
include { generate_variant_db         } from './workflows/variant_db'
include { generate_novel_isoform_db   } from './workflows/novel_isoform_db'
include { combine_protein_dbs         } from './workflows/combine_db'


workflow {

    samplesheet_path = file(params.input, checkIfExists: true)

    samplesheet_ch = Channel
    .fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        def meta = [
            id: row.sample,                  // sample_id
            strandedness: row.strandedness ?: "unstranded"
        ]
        tuple(meta, file(row.fastq_1), file(row.fastq_2))
    }
    gene_annotation_gtf = file(params.gtf, checkIfExists: true)


    // 1. Run nf-core/rnavar (assumes samplesheet CSV matches expected format)
    NFCORE_RNAVAR(
        Channel.fromPath(samplesheet_path),
        Channel.value(false)
    )

    // 2. StringTie on markdup BAMs from RNAVAR
    RUN_STRINGTIE(
        NFCORE_RNAVAR.out.markdup_bams,   // emits: [ val(meta), path(bam), path(bai) ]
        gene_annotation_gtf
    )

    // 3. Generate novel isoform DBs from StringTie GTFs
    generate_novel_isoform_db(
        RUN_STRINGTIE.out.stringtie_gtf   // emits: [ val(meta), path(gtf) ]
    )

    // 4. Variant DB from annotated VCFs
    generate_variant_db(
        samplesheet_ch.meta.id,                   // or meta.id if running for all samples
        NFCORE_RNAVAR.out.annotated_vcf
    )

    variant_fasta = generate_variant_db.out.variant_db

    // Step 2: Novel isoform-based DB generation
    isoform_fasta = generate_novel_isoform_db(samplesheet_ch.meta.id).out.isoform_db

    // Step 3: Combine protein DBs
    combine_protein_dbs(samplesheet_ch.meta.id, variant_fasta, isoform_fasta)
}


