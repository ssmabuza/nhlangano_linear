#!/usr/bin/env bash

set -euo pipefail

if [[ -z "${SOFTWARE:-}" ]]; then
  echo "ERROR: SOFTWARE is not set. Example: export SOFTWARE=\$HOME/Software" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cmake -S "${PROJECT_ROOT}" -B "${PROJECT_ROOT}/build/debug" \
  -DMPI_DIR="${SOFTWARE}/mpi" \
  -DTrilinos_ROOT="${SOFTWARE}/trilinos" \
  --preset debug

cmake --build --preset debug
