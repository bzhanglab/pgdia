nextflow.enable.dsl=2

workflow combine_protein_dbs {
    take:
        input_ch // channel: [ id, var_modified_peptides, novel_pep ]

    main:
        COMBINE_PROTEIN_DBS(input_ch)

    emit:
        combined_db = COMBINE_PROTEIN_DBS.out.combined_db
        novel_db    = COMBINE_PROTEIN_DBS.out.novel_db
}

process COMBINE_PROTEIN_DBS {
    tag "${sample_id}"
    container 'python:3.11-slim'

    input:
        tuple val(sample_id), path(var_modified_peptides), path(novel_pep)

    output:
        path "${sample_id}_combined_protein_db.fa", emit: combined_db
        path "r${sample_id}_novel_protein_db.fa"  , emit: novel_db

    script:
    """
    python3 assemble_protein_db.py \
        "${var_modified_peptides}" \
        "${novel_pep}" \
        "${sample_id}_combined_protein_db.fa" \
        "${sample_id}_novel_protein_db.fa"
    """
}
