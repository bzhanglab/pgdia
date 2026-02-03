nextflow.enable.dsl=2

process COMBINE_PROTEIN_DBS {
    tag "${id}"
    container 'python:3.11-slim'

    publishDir "${params.outdir}/protein_db", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(var_modified_peptides), path(novel_pep)
        path protein_reference_db

    output:
    tuple val(id), path("${id}_combined_protein_db.fa"), emit: combined_db
    tuple val(id), path("${id}_novel_protein_db.fa"),    emit: novel_db

    script:
    """
    python3 ${projectDir}/bin/assemble_protein_db.py \
        "${var_modified_peptides}" \
        "${novel_pep}" \
        "${id}_combined_protein_db.fa" \
        "${id}_novel_protein_db.fa" \
        "${protein_reference_db}"
    """
}
