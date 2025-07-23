nextflow.enable.dsl=2

include { NFCORE_RNAVAR } from './submodules/rnavar'
include { RUN_STRINGTIE } from './workflows/stringtie/stringtie'
include { generate_variant_db         } from './workflows/variant_db'
include { generate_novel_isoform_db   } from './workflows/novel_isoform_db'
include { combine_protein_dbs         } from './workflows/combine_db'


workflow {

    // Parse the samplesheet for all samples
    samplesheet_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def meta = [ id: row.sample, strandedness: row.strandedness ]
        [ meta, file(row.fastq_1), file(row.fastq_2) ]
    }


    // 1. Run nf-core/rnavar (assumes samplesheet CSV matches expected format)
    NFCORE_RNAVAR(
        Channel.fromPath(params.samplesheet),
        Channel.value(false)
    )

    // 2. StringTie on markdup BAMs from RNAVAR
    RUN_STRINGTIE(
        NFCORE_RNAVAR.out.markdup_bams,   // emits: [ val(meta), path(bam), path(bai) ]
        params.gtf
    )

    // 3. Generate novel isoform DBs from StringTie GTFs
    generate_novel_isoform_db(
        RUN_STRINGTIE.out.stringtie_gtf   // emits: [ val(meta), path(gtf) ]
    )

    // 4. Variant DB from annotated VCFs
    generate_variant_db(
        samplesheet_ch.meta.id,                   // or meta.id if running for all samples
        NFCORE_RNAVAR.out.annotated_vcfs
    )

    variant_fasta = generate_variant_db.out.out

    // Step 2: Novel isoform-based DB generation
    novel_output = generate_novel_isoform_db(params.sample_id)

    // Step 3: Combine protein DBs
    combine_protein_dbs(params.sample_id, variant_fasta, novel_output.out)
}


