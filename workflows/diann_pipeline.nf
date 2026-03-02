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
  cpus task.cpus

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

process POSTPROCESS_DIANN_REPORT {
  tag { meta.id }

  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/python:3.8.3' :
      'quay.io/biocontainers/python:3.8.3' }"

  publishDir { "${params.outdir}/diann_output/${meta.id}" }, mode: 'copy', overwrite: true

  input:
    tuple val(meta), path(diann_report), path(novel_fasta), path(isoform_annotation)
    path protein_reference_db

  output:
    tuple val(meta), path("${meta.id}_novel_matrix.tsv"), emit: processed_novel_matrix

  script:
    """
    set -euo pipefail

    WORKDIR="\${NXF_TASK_WORKDIR:-.}"
    DEPS_DIR="\$WORKDIR/.pydeps"

    python3 -m pip install --no-cache-dir --target "\$DEPS_DIR" pandas pyarrow biopython
    export PYTHONPATH="\$DEPS_DIR:\${PYTHONPATH:-}"

    python3 ${projectDir}/bin/process_parquet_report.py \\
      --report-parquet "${diann_report}" \\
      --output-prefix "${meta.id}" \\
      --novel-fasta "${novel_fasta}" \\
      --isoform-annotation "${isoform_annotation}" \\
      --reference-fasta "${protein_reference_db}"
    """
}

workflow DIANN_PIPELINE {
  take:
    // expected items: tuple(meta, path(protein_db_fa), path(dia_raw), path(novel_fasta), path(isoform_tmap))
    // meta.id is sample id
    samples_ch

  main:

    def diann_name_file_ch = LOAD_DIANN_IMAGE( Channel.value(params.diann_image) )
    def protein_reference_db_ch = Channel.value(file(params.protein_reference_db, checkIfExists: true))

    // convert diann_image_name.txt to a value
    def diann_image_name_ch = diann_name_file_ch.map { f -> f.text.trim() }

    def run_diann_input_ch = samples_ch
      .map { meta, protein_db_fa, dia_raw, _, _ ->
        tuple(meta, protein_db_fa, dia_raw)
      }

    def run_diann_out = RUN_DIANN(diann_image_name_ch, run_diann_input_ch)
    def diann_out_ch = run_diann_out.diann_report
    def postprocess_meta_ch = samples_ch
      .map { meta, _, _, novel_fasta, isoform_annotation ->
        tuple(meta, novel_fasta, isoform_annotation)
      }

    def postprocess_input_ch = diann_out_ch
      .join(postprocess_meta_ch, by: 0, failOnMismatch: true)
      .map { meta, diann_report, novel_fasta, isoform_annotation ->
        tuple(meta, diann_report, novel_fasta, isoform_annotation)
      }

    def postprocess_out = POSTPROCESS_DIANN_REPORT(postprocess_input_ch, protein_reference_db_ch)
    def postprocessed_novel_matrix_ch = postprocess_out.processed_novel_matrix

  emit:
    diann_out = diann_out_ch
    diann_postprocessed = postprocessed_novel_matrix_ch
}
