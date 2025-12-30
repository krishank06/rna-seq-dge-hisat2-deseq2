#!/usr/bin/env bash
# ============================================================
# Script name: 03_featureCounts.sh
#
# Description:
#   Generate gene-level read counts from BAM files
#   using featureCounts.
#
# Tool:
#   Subread / featureCounts >= 2.1
#
# Author: Krishan Kumar
# License: MIT
# ============================================================

set -euo pipefail

# -------- User-configurable variables --------
PROJECT_DIR=/path/to/project
BAM_DIR=${PROJECT_DIR}/bam
GTF=${PROJECT_DIR}/annotation/gencode.v44.annotation.gtf
COUNT_DIR=${PROJECT_DIR}/counts
THREADS=8
# --------------------------------------------

mkdir -p "${COUNT_DIR}"

featureCounts \
  -T ${THREADS} \
  -p -B -C \
  -s 0 \
  -t exon \
  -g gene_id \
  -a ${GTF} \
  -o ${COUNT_DIR}/gene_counts.txt \
  ${BAM_DIR}/*.sorted.bam
