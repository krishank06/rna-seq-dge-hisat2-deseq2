#!/usr/bin/env bash
# ============================================================
# Script name: 00_build_hisat2_index.sh
#
# Description:
#   Build HISAT2 genome indices from a reference FASTA.
#
# ============================================================

set -euo pipefail

source config/reference_config.sh

mkdir -p "$(dirname ${HISAT2_INDEX_PREFIX})"

hisat2-build \
  -p ${THREADS} \
  ${GENOME_FA} \
  ${HISAT2_INDEX_PREFIX}
