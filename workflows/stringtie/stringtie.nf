include { cohortTissue } from '../tissue_metadata'

workflow RUN_STRINGTIE {
    take:
    aligned_bam        // channel: [ val(meta), bam, bai ]
    annotation_gtf     // reference annotation

    main:
    ch_versions = Channel.empty()

    // Merge every run-level mark-duplicate BAM into one cohort BAM.
    def expected_runs = (params.cohort_expected_runs ?: 27) as int
    if (expected_runs < 1) {
        error("Invalid --cohort_expected_runs '${params.cohort_expected_runs}': expected a positive integer.")
    }

    def cohort_bams = aligned_bam
        .map { meta, bam, bai ->
            tuple(params.cohort_db_id, meta, bam)
        }

    def checked_cohort_count = cohort_bams
        .count()
        .map { observed_runs ->
            if (observed_runs != expected_runs) {
                error("StringTie cohort '${params.cohort_db_id}' expected exactly ${expected_runs} run-level BAMs, but received ${observed_runs}.")
            }
            tuple(params.cohort_db_id, observed_runs)
        }

    def stringtie_input = cohort_bams
        .groupTuple(size: expected_runs)
        .combine(checked_cohort_count)
        .map { db_id, metas, bams, checked_db_id, observed_runs ->
            if (db_id != checked_db_id || metas.size() != expected_runs || bams.size() != expected_runs) {
                error("StringTie cohort aggregation mismatch for '${db_id}': expected ${expected_runs}, received ${metas.size()} metadata records and ${bams.size()} BAMs.")
            }

            def run_ids = metas.collect { run_meta -> run_meta.rnaseq_run ?: run_meta.id }
            def duplicate_run_ids = run_ids
                .countBy { run_id -> run_id }
                .findAll { _run_id, count -> count > 1 }
                .keySet()
                .sort()
            def unique_run_ids = run_ids.unique()
            if (unique_run_ids.size() != expected_runs) {
                error("StringTie cohort '${db_id}' expected ${expected_runs} unique RNA-seq runs, but received ${unique_run_ids.size()}. Duplicate run IDs: ${duplicate_run_ids ?: 'none'}.")
            }

            def strandedness = metas.collect { it.strandedness ?: 'unstranded' }.unique()
            if (strandedness.size() > 1) {
                error("StringTie requires one strandedness setting for the merged RNA-seq BAM. Found ${strandedness.join(', ')}.")
            }
            def first_meta = metas[0]
            def dia_raw = metas.collect { it.dia_raw }.findAll { it }.flatten().unique()
            def meta = first_meta + [
                id          : db_id,
                db_id       : db_id,
                rnaseq_runs : unique_run_ids,
                tissue      : cohortTissue(metas),
                expected_runs: expected_runs,
                observed_runs: observed_runs,
                dia_raw     : dia_raw ?: []
            ]
            tuple(meta, bams)
        }

    MERGE_STRINGTIE_BAMS(stringtie_input)
    ch_versions = ch_versions.mix(MERGE_STRINGTIE_BAMS.out.versions)

    STRINGTIE_STRINGTIE(MERGE_STRINGTIE_BAMS.out.bam, annotation_gtf)
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

    // Gate all sample-level work on completion of the cohort-level StringTie
    // task and use that cohort GTF to guide every sample-level quantification.
    def cohort_stringtie_done = STRINGTIE_STRINGTIE.out.coverage_gtf.collect()
    def gated_sample_bams = aligned_bam
        .combine(cohort_stringtie_done)
        .map { meta, bam, bai, _cohort_meta, cohort_gtf ->
            tuple(meta, bam, cohort_gtf)
        }

    STRINGTIE_PER_SAMPLE(gated_sample_bams)
    ch_versions = ch_versions.mix(STRINGTIE_PER_SAMPLE.out.versions)

    emit:
    stringtie_gtf        = STRINGTIE_STRINGTIE.out.coverage_gtf
    sample_stringtie_gtf = STRINGTIE_PER_SAMPLE.out.coverage_gtf
    versions             = ch_versions
}

process MERGE_STRINGTIE_BAMS {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::samtools=1.21" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.stringtie_input.bam"), emit: bam
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = bams instanceof List ? bams : [bams]
    def expected_runs = (params.cohort_expected_runs ?: 27) as int
    if (input_files.size() != expected_runs) {
        error("StringTie cohort '${meta.id}' merge expected ${expected_runs} BAMs, but received ${input_files.size()}.")
    }
    def bam_inputs = input_files.join(' ')
    def output_bam = "${prefix}.stringtie_input.bam"

    """
    set -euo pipefail

    samtools quickcheck -v ${bam_inputs}
    samtools merge --threads ${task.cpus} -o "${output_bam}" ${bam_inputs}

    if [[ ! -s "${output_bam}" ]]; then
        echo "ERROR: StringTie cohort BAM is missing or empty: ${output_bam}" >&2
        exit 1
    fi
    samtools quickcheck -v "${output_bam}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process STRINGTIE_STRINGTIE {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    tuple val(meta), path(bam)
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

process STRINGTIE_PER_SAMPLE {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    tuple val(meta), path(bam), path(cohort_gtf)

    output:
    tuple val(meta), path("*.sample_stringtie_output.gtf"), emit: coverage_gtf
    path  "versions.yml"                                  , emit: versions

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
    stringtie \
        $bam \
        $strandedness \
        -G $cohort_gtf \
        -e \
        -o ${prefix}.sample_stringtie_output.gtf \
        -p $task.cpus \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
