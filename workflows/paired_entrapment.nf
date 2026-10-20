nextflow.enable.dsl=2

include { RUN_DIANN as RUN_DIANN_ENTRAPMENT } from './diann_pipeline'

process BUILD_PAIRED_ENTRAPMENT_DB {
    tag { meta.id }
    label 'process_medium'
    conda (params.enable_conda ? 'conda-forge::python=3.8.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
    publishDir { "${params.outdir}/paired_entrapment/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(novel_fasta)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: database
    tuple val(meta), path("${meta.id}.pairs.tsv"), emit: pairs
    tuple val(meta), path("${meta.id}.manifest.json"), emit: manifest
    tuple val(meta), path("${meta.id}.excluded.tsv"), emit: excluded

    script:
    """
    set -euo pipefail
    python3 -m pip install --no-cache-dir --target .pydeps biopython
    export PYTHONPATH="\$PWD/.pydeps:\${PYTHONPATH:-}"
    python3 ${projectDir}/bin/build_paired_entrapment_db.py \\
      --reference-fasta "${reference_fasta}" \\
      --novel-fasta "${novel_fasta}" \\
      --sample "${meta.id}" \\
      --output-prefix "${meta.id}" \\
      --seed ${params.entrapment_seed} \\
      --missed-cleavages ${params.entrapment_missed_cleavages}
    """
}

process ASSESS_PAIRED_ENTRAPMENT {
    tag { meta.id }
    label 'process_medium'
    conda (params.enable_conda ? 'conda-forge::python=3.8.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
    publishDir { "${params.outdir}/paired_entrapment/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(report), path(pairs)

    output:
    tuple val(meta), path("${meta.id}_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.id}_method.json"), emit: method

    script:
    def recalculate = params.entrapment_recalculate_q ? '--recalculate-q-values' : ''
    """
    set -euo pipefail
    python3 -m pip install --no-cache-dir --target .pydeps pandas pyarrow
    export PYTHONPATH="\$PWD/.pydeps:\${PYTHONPATH:-}"
    python3 ${projectDir}/bin/calculate_paired_entrapment_fdp.py \\
      --report "${report}" \\
      --pairs "${pairs}" \\
      --output-prefix "${meta.id}" \\
      --q-column "${params.entrapment_q_column}" \\
      --cutoff ${params.entrapment_cutoff} \\
      ${recalculate}
    """
}

workflow PAIRED_ENTRAPMENT {
    take:
    samples_ch
    diann_image_name
    reference_fasta

    main:
    BUILD_PAIRED_ENTRAPMENT_DB(
        samples_ch.map { meta, db, raw, novel, annotation -> tuple(meta, novel) },
        reference_fasta
    )
    search_input = BUILD_PAIRED_ENTRAPMENT_DB.out.database
        .join(samples_ch.map { meta, db, raw, novel, annotation -> tuple(meta, raw) }, by: 0, failOnMismatch: true)
    RUN_DIANN_ENTRAPMENT(diann_image_name, search_input, true)
    assessment_input = RUN_DIANN_ENTRAPMENT.out.diann_report
        .join(BUILD_PAIRED_ENTRAPMENT_DB.out.pairs, by: 0, failOnMismatch: true)
    ASSESS_PAIRED_ENTRAPMENT(assessment_input)

    emit:
    summary = ASSESS_PAIRED_ENTRAPMENT.out.summary
    method = ASSESS_PAIRED_ENTRAPMENT.out.method
}
