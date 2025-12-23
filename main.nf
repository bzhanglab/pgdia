nextflow.enable.dsl=2

include { NFCORE_RNAVAR } from './submodules/rnavar'

include { RUN_STRINGTIE } from './workflows/stringtie/stringtie'
include { generate_variant_db         } from './workflows/variant_db'
include { generate_novel_isoform_db   } from './workflows/novel_isoform_db'
include { combine_protein_dbs         } from './workflows/combine_db'


workflow MAIN {

    take:
        samplesheet_path
        gene_annotation_gtf

    main:
        parsed_samplesheet = Channel
            .fromPath(samplesheet_path, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def id = row.sample ?: row.id ?: row.sample_id
                if (!id) {
                    error("Samplesheet row is missing a 'sample' identifier.")
                }

                def hasFastq = row.fastq_1
                def hasBam   = !hasFastq && row.bam

                if (!hasFastq && !hasBam) {
                    error("Samplesheet row for '${id}' must provide either FASTQ (fastq_1/fastq_2) or BAM/BAI columns.")
                }
                if (hasFastq && hasBam) {
                    error("Samplesheet row for '${id}' must not mix FASTQ and BAM inputs.")
                }

                def fastq1 = hasFastq ? file(row.fastq_1, checkIfExists: true) : []
                def fastq2 = (hasFastq && row.fastq_2) ? file(row.fastq_2, checkIfExists: true) : []
                def bam    = hasBam   ? file(row.bam,   checkIfExists: true) : []
                def bai    = (hasBam && row.bai) ? file(row.bai, checkIfExists: true) : []

                def meta = [
                    id          : id,
                    strandedness: row.strandedness ?: "unstranded",
                    single_end  : hasFastq ? !row.fastq_2 : true
                ]

                tuple(meta, fastq1, fastq2, bam, bai)
            }
            .view { row -> "Samplesheet entry: ${row[0].id}" }

        rnavar_samplesheet = parsed_samplesheet.map { meta, fastq1, fastq2, bam, bai ->
            tuple(meta.id, meta, fastq1 ?: [], fastq2 ?: [], bam ?: [], bai ?: [], [], [], [], [])
        }

        // 1. Run nf-core/rnavar (assumes samplesheet CSV matches expected format)
        NFCORE_RNAVAR(
            rnavar_samplesheet,
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

        generate_variant_db(
            NFCORE_RNAVAR.out.annotated_vcf
        )

        variant_fasta = generate_variant_db.out.variant_db
            .map { meta, fa -> tuple(meta.id, fa) }                // tuple(id, var_modified_peptides.fa)
        isoform_fasta = generate_novel_isoform_db.out.isoform_db

        // Step 3: Combine protein DBs
        combine_in_ch = variant_fasta
            .join(isoform_fasta)
            .map { id, var_fa, novel_fa ->
                tuple(id, var_fa, novel_fa)
            }

        combine_protein_dbs(combine_in_ch)

    emit:
        combined_db_ch = combine_protein_dbs.out.combined_db
        novel_db_ch    = combine_protein_dbs.out.novel_db
}

workflow {
    samplesheet_path = file(params.input, checkIfExists: true)

    def resolved_gtf = params.gtf ?: (
        params.genomes && params.genome && params.genomes.containsKey(params.genome)
            ? params.genomes[params.genome].gtf
            : null
    )

    if (!resolved_gtf) {
        error("Missing required reference GTF: set --gtf or provide a genome entry with a gtf path.")
    }

    gene_annotation_gtf = file(resolved_gtf, checkIfExists: true)

    MAIN(samplesheet_path, gene_annotation_gtf)
}
