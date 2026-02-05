nextflow.enable.dsl=2

process COMBINE_PROTEIN_DBS {
    tag "${meta.id}"
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    publishDir "${params.outdir}/protein_db", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(var_modified_peptides), path(novel_pep)
        path protein_reference_db

    output:
    tuple val(meta), path("${meta.id}_combined_protein_db.fa"), emit: combined_db
    tuple val(meta), path("${meta.id}_novel_protein_db.fa"),    emit: novel_db

    script:
    """
    python3 ${projectDir}/bin/assemble_protein_db.py \
        "${var_modified_peptides}" \
        "${novel_pep}" \
        "${meta.id}_combined_protein_db.fa" \
        "${meta.id}_novel_protein_db.fa" \
        "${protein_reference_db}"
    """
}
