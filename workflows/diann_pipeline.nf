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
params.diann_tar_tag = params.diann_tar_tag ?: "diann-2.0.2"        // tag used if loading tar
params.diann_bin   = params.diann_bin   ?: "/diann-2.0.2/diann-linux"
params.diann_cpus  = params.diann_cpus  ?: 22


process LOAD_DIANN_IMAGE {
  tag "load_diann_image"
  // This process must run where Docker CLI is available and can talk to the daemon
  // Do NOT set a container here.

  input:
    val(diann_image)

  output:
    path("diann_image_name.txt")

  script:
    """
    set -euo pipefail

    img="${diann_image}"

    if [[ "\$img" == *.tar ]]; then
      if [[ ! -f "\$img" ]]; then
        echo "ERROR: DIA-NN tar not found: \$img" >&2
        exit 1
      fi

      # Load tar into Docker
      docker load -i "\$img" >/dev/null

      # Ensure the desired tag exists; if not, tag the most recent image id
      if ! docker image inspect "${params.diann_tar_tag}" >/dev/null 2>&1; then
        last_id=\$(docker images -q | head -n 1)
        if [[ -z "\$last_id" ]]; then
          echo "ERROR: docker load succeeded but cannot find an image id to tag" >&2
          exit 1
        fi
        docker tag "\$last_id" "${params.diann_tar_tag}" >/dev/null
      fi

      echo "${params.diann_tar_tag}" > diann_image_name.txt
    else
      echo "\$img" > diann_image_name.txt
    fi
    """
  // read the image name from file
  // Nextflow will capture this file, but we declared `val` output, so we need a small trick:
  // Use stdout as the value by printing it as the last line:
  // We therefore set diann_image_name as stdout.
  // (See below: use `shell:` to directly echo.)
}

process RUN_DIANN {
  label 'process_high_large_disk'
  tag { meta.id }
  cpus params.diann_cpus

  publishDir { "${params.outdir}/diann_output/${meta.id}" }, mode: 'copy', overwrite: true

  input:
    val(diann_image_name)
    tuple val(meta), path(protein_db_fa), path(dia_raw)

  output:
    tuple val(meta), path("${meta.id}_report.parquet"), emit: diann_report
    path("${meta.id}_lib.parquet"), optional: true
    path("${meta.id}_matrices.parquet"), optional: true

  container { diann_image_name }

  script:
    """
    set -euo pipefail

    ${params.diann_bin} -v \\
      --f "${dia_raw}" \\
      --lib "" \\
      --threads ${task.cpus} \\
      --verbose 1 \\
      --out "${meta.id}_report.parquet" \\
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
    // expected items: tuple(meta, path(protein_db_fa), path(dia_raw))
    // meta.id is sample id
    samples_ch

  main:

    def diann_name_file_ch = LOAD_DIANN_IMAGE( Channel.value(params.diann_image) )

    // convert diann_image_name.txt to a value
    def diann_image_name_ch = diann_name_file_ch.map { f -> f.text.trim() }

    diann_out_ch = RUN_DIANN(diann_image_name_ch, samples_ch).diann_report

  emit:
    diann_out = diann_out_ch
}
