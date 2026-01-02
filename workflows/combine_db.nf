nextflow.enable.dsl=2

process COMBINE_PROTEIN_DBS {
    tag "${sample_id}"
    container 'python:3.11-slim'

    publishDir "${params.outdir}/protein_db", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(var_modified_peptides), path(novel_pep)

    output:
    tuple val(sample_id), path("${sample_id}_combined_protein_db.fa"), emit: combined_db
    tuple val(sample_id), path("${sample_id}_novel_protein_db.fa"),    emit: novel_db

    script:
    """
    python3 ${projectDir}/bin/assemble_protein_db.py \
        "${var_modified_peptides}" \
        "${novel_pep}" \
        "${sample_id}_combined_protein_db.fa" \
        "${sample_id}_novel_protein_db.fa"
    """
}
