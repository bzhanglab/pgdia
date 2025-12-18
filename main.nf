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
                id          : row.sample,                  // sample_id
                strandedness: row.strandedness ?: "unstranded"
            ]

            def columns = row.keySet()
            def hasFastqColumns = ['fastq_1', 'fastq_2'].any { col -> columns.contains(col) }
            def hasBamColumns   = ['bam', 'bai'].any { col -> columns.contains(col) }

            def useFastq = hasFastqColumns && row.fastq_1
            def useBam   = !useFastq && hasBamColumns && row.bam

            if (!useFastq && !useBam) {
                error("Samplesheet row for '${meta.id}' must provide either FASTQ (fastq_1/fastq_2) or BAM/BAI columns.")
            }

            def fastq1 = useFastq && row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : null
            def fastq2 = useFastq && row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            def bam    = useBam   && row.bam     ? file(row.bam,     checkIfExists: true) : null
            def bai    = useBam   && row.bai     ? file(row.bai,     checkIfExists: true) : null

            tuple(meta, fastq1, fastq2, bam, bai)
        }
        .view { row -> "Samplesheet entry: ${row[0].id}" }
    gene_annotation_gtf = file(params.gtf, checkIfExists: true)

    id_ch = samplesheet_ch.map { meta, fastq1, fastq2, bam, bai -> meta.id }

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
    def novel_isoform_run = generate_novel_isoform_db(
        RUN_STRINGTIE.out.stringtie_gtf   // emits: [ val(meta), path(gtf) ]
    )

    def variant_db_run = generate_variant_db(
        id_ch,                   // or meta.id if running for all samples
        NFCORE_RNAVAR.out.annotated_vcf
    )

    variant_fasta = variant_db_run.out.variant_db
    isoform_fasta = novel_isoform_run.out.isoform_db

    // Step 3: Combine protein DBs
    
    combine_protein_dbs(
        id_ch,
        variant_fasta,
        isoform_fasta
    )
}
