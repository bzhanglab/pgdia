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

        // 1. Run nf-core/rnavar (assumes samplesheet CSV matches expected format)
        NFCORE_RNAVAR(
            Channel.value(samplesheet_path),
            Channel.value(false)
        )

        // Debug + filter nulls before StringTie
        NFCORE_RNAVAR.out.markdup_bams.view { meta, bam, bai ->
            "markdup_bams: ${meta.id} bam=${bam} bai=${bai}"
        }       

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

    gene_annotation_gtf = file(params.gtf, checkIfExists: true)

    MAIN(samplesheet_path, gene_annotation_gtf)
}
