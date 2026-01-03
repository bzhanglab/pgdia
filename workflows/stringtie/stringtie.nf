workflow RUN_STRINGTIE {
    take:
    aligned_bam        // channel: [ val(meta), bam, bai ]
    annotation_gtf     // reference annotation

    main:
    STRINGTIE_STRINGTIE(aligned_bam, annotation_gtf)

    emit:
    stringtie_gtf = STRINGTIE_STRINGTIE.out.coverage_gtf
    versions      = STRINGTIE_STRINGTIE.out.versions
}

process STRINGTIE_STRINGTIE {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf

    output:
    tuple val(meta), path("*.stringtie_output.gtf")   , emit: coverage_gtf
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }

    """
    stringtie \\
        $bam \\
        $strandedness \\
        -G $gtf \\
        -o ${prefix}.stringtie_output.gtf \\
        -p $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
