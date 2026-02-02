process VCF_DECOMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? 'bioconda::htslib=1.20' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2' :
        'biocontainers/htslib:1.20--h5efdd21_2' }"

    input:
    tuple val(meta), path(vcf_gz), path(tbi)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    script:
    def out_name = vcf_gz.baseName.endsWith('.vcf') ? vcf_gz.baseName : "${vcf_gz.baseName}.vcf"
    """
    set -euo pipefail
    bgzip -dc ${vcf_gz} > ${out_name}
    """
}