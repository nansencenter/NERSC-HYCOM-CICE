#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_ENV_DIR="${SCRIPT_DIR}/py-nhc-env"
PYTHON_POST_SCRIPT_FILE="${SCRIPT_DIR}/py-nhc-env-post.sh"
PYTHON_CONDA_ENV_FILE="${SCRIPT_DIR}/py-nhc-env.yml"
REQUIRED_FILES=("${PYTHON_POST_SCRIPT_FILE}" "${PYTHON_CONDA_ENV_FILE}")

for file in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Error: Required file '${file}' not found in ${SCRIPT_DIR}" >&2
        exit 1
    fi
done

if ! command -v conda-containerize &> /dev/null; then
    echo "Error: 'conda-containerize' command is not available in PATH." >&2
    echo "Please load the 'hpc-container-wrapper' module."
    exit 1
fi

conda-containerize new \
        --prefix "${PYTHON_ENV_DIR}" \
        --post-install "${PYTHON_POST_SCRIPT_FILE}" \
        "${PYTHON_CONDA_ENV_FILE}"
