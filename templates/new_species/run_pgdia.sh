#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW="${NEXTFLOW:-/home/chenw/nextflow}"
PGDIA_PIPELINE="${PGDIA_PIPELINE:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
CUSTOM_CONFIG="${CUSTOM_CONFIG:-${SCRIPT_DIR}/species.config}"
INPUT_SAMPLESHEET="${INPUT_SAMPLESHEET:-${SCRIPT_DIR}/samplesheet.csv}"

# Use a separate run root for each cohort to isolate its cache and outputs.
RUN_ROOT="${RUN_ROOT:-${SCRIPT_DIR}/run}"
OUTDIR="${OUTDIR:-${RUN_ROOT}/results}"
WORK_DIR="${WORK_DIR:-${RUN_ROOT}/work}"
PROFILE="${PROFILE:-docker}"
RESUME_ID="${RESUME_ID:-}"
DRY_RUN="${DRY_RUN:-0}"

[[ -x "${NEXTFLOW}" ]] || { echo "Nextflow executable not found: ${NEXTFLOW}" >&2; exit 1; }
[[ -f "${PGDIA_PIPELINE}/main.nf" ]] || { echo "Pipeline not found: ${PGDIA_PIPELINE}" >&2; exit 1; }
[[ -s "${CUSTOM_CONFIG}" ]] || { echo "Config not found or empty: ${CUSTOM_CONFIG}" >&2; exit 1; }
[[ -s "${INPUT_SAMPLESHEET}" ]] || { echo "Samplesheet not found or empty: ${INPUT_SAMPLESHEET}" >&2; exit 1; }

# Resolve input paths before changing the launch directory.
PGDIA_PIPELINE="$(realpath "${PGDIA_PIPELINE}")"
CUSTOM_CONFIG="$(realpath "${CUSTOM_CONFIG}")"
INPUT_SAMPLESHEET="$(realpath "${INPUT_SAMPLESHEET}")"
RUN_ROOT="$(realpath -m "${RUN_ROOT}")"
OUTDIR="$(realpath -m "${OUTDIR}")"
WORK_DIR="$(realpath -m "${WORK_DIR}")"
NEXTFLOW="$(realpath "${NEXTFLOW}")"

if grep -n 'REPLACE_' "${CUSTOM_CONFIG}" "${INPUT_SAMPLESHEET}"; then
    echo "Fill in the template placeholders above before launching." >&2
    exit 1
fi

export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"
# Set NXF_OFFLINE=false if Nextflow needs to download a missing plugin.
export NXF_OFFLINE="${NXF_OFFLINE:-true}"

RESUME_ARGS=(-resume)
if [[ -n "${RESUME_ID}" ]]; then
    RESUME_ARGS=(-resume "${RESUME_ID}")
fi

# Species, references and typed parameters live in species.config.
COMMAND=("${NEXTFLOW}" -c "${CUSTOM_CONFIG}" run "${PGDIA_PIPELINE}"
    -profile "${PROFILE}" "${RESUME_ARGS[@]}"
    --input "${INPUT_SAMPLESHEET}" --outdir "${OUTDIR}"
    -work-dir "${WORK_DIR}")

if [[ "${DRY_RUN}" == 1 ]]; then
    printf 'Launch directory: %s\n' "${RUN_ROOT}"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "${RUN_ROOT}" "${OUTDIR}" "${WORK_DIR}"
cd "${RUN_ROOT}"
exec "${COMMAND[@]}"
