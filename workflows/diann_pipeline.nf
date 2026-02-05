nextflow.enable.dsl=2

/*
 * DIA-NN (Docker image: diann-2.0.2)
 *
 * Your tested command:
 *   docker run diann-2.0.2 /diann-2.0.2/diann-linux -v --f <raw.d> ... --fasta <db.fa> ...
 *
 * This Nextflow process:
 * - takes the raw data path from meta.raw (staged as a path input)
 * - takes the per-sample FASTA from meta.fasta (staged as a path input)
 * - writes outputs as <id>_report.parquet and <id>_lib.parquet (and matrices parquet)
 * - runs inside the diann-2.0.2 container and calls /diann-2.0.2/diann-linux
 */

params.outdir      = params.outdir ?: "results"
params.diann_image = params.diann_image ?: "diann-2.0.2"
params.diann_bin   = params.diann_bin   ?: "/diann-2.0.2/diann-linux"
params.diann_cpus  = params.diann_cpus  ?: 25

process RUN_DIANN {
  tag { meta.id }
  container params.diann_image
  cpus params.diann_cpus

  publishDir { "${params.outdir}/diann_output/${meta.id}" }, mode: 'copy', overwrite: true

  input:
    tuple val(meta), path(protein_db_fa)

  output:
    tuple val(meta), path("${meta.id}_report.parquet"), emit: diann_report
          // path("${meta.id}_lib.parquet"),
          // path("${meta.id}_matrices.parquet")

  script:
    """
    set -euo pipefail

    ${params.diann_bin} -v \\
      --f "${meta.dia_raw}" \\
      --lib "" \\
      --threads ${task.cpus} \\
      --verbose 1 \\
      --out "${meta.id}_report.parquet" \\
      --out-lib "${meta.id}_lib.parquet" \\
      --qvalue 1 \\
      --matrices \\
      --out-matrices "${meta.id}_matrices.parquet" \\
      --predictor \\
      --fasta "${protein_db_fa}" \\
      --fasta-search \\
      --min-fr-mz 300 \\
      --max-fr-mz 1500 \\
      --met-excision \\
      --min-pep-len 7 \\
      --max-pep-len 30 \\
      --min-pr-mz 300 \\
      --max-pr-mz 1700 \\
      --min-pr-charge 1 \\
      --max-pr-charge 4 \\
      --cut K*,R* \\
      --missed-cleavages 1 \\
      --unimod4 \\
      --var-mods 1 \\
      --var-mod UniMod:35,15.994915,M \\
      --var-mod UniMod:1,42.010565,*n \\
      --peptidoforms \\
      --reanalyse \\
      --rt-profiling \\
      --report-decoys
    """
}

workflow DIANN_PIPELINE {
  take:
    // expected items: tuple(meta, path(raw_d), path(protein_db_fa))
    // meta.id is sample id
    samples_ch

  main:
    RUN_DIANN(samples_ch)

  emit:
    diann_out = RUN_DIANN.diann_report
}
