#!/usr/bin/env bash
# ============================================================
# Script name: 01_fastp_trimming.sh
#
# Description:
#   Adapter trimming and quality filtering of paired-end
#   RNA-seq FASTQ files using fastp.
#
# Tool:
#   fastp >= 0.23
#
# Author: Krishan Kumar
# License: MIT
# ============================================================

set -euo pipefail

# -------- User-configurable variables --------
PROJECT_DIR=/path/to/project
RAW_FASTQ_DIR=${PROJECT_DIR}/raw_fastq
TRIMMED_FASTQ_DIR=${PROJECT_DIR}/trimmed_fastq
QC_DIR=${PROJECT_DIR}/qc/fastp
THREADS=8
# --------------------------------------------

mkdir -p "${TRIMMED_FASTQ_DIR}" "${QC_DIR}"

for R1 in ${RAW_FASTQ_DIR}/*_R1.fastq.gz; do
    R2=${R1/_R1/_R2}
    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${TRIMMED_FASTQ_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
      -O "${TRIMMED_FASTQ_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 30 \
      --thread ${THREADS} \
      --html "${QC_DIR}/${SAMPLE}_fastp.html" \
      --json "${QC_DIR}/${SAMPLE}_fastp.json"
done
