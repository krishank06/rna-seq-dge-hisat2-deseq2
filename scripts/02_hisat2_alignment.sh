#!/usr/bin/env bash
# ============================================================
# Script name: 02_hisat2_alignment.sh
#
# Description:
#   Align trimmed paired-end RNA-seq reads to GRCh38
#   using HISAT2 and generate sorted BAM files.
#
# Tools:
#   HISAT2 >= 2.2.1
#   SAMtools >= 1.15
#
# Author: Krishan Kumar
# License: MIT
# ============================================================

set -euo pipefail

# -------- User-configurable variables --------
PROJECT_DIR=/path/to/project
TRIMMED_FASTQ_DIR=${PROJECT_DIR}/trimmed_fastq
HISAT2_INDEX=${PROJECT_DIR}/hisat2_index/grch38
BAM_DIR=${PROJECT_DIR}/bam
THREADS=8
# --------------------------------------------

mkdir -p "${BAM_DIR}"

for R1 in ${TRIMMED_FASTQ_DIR}/*_R1_trimmed.fastq.gz; do
    R2=${R1/_R1_trimmed/_R2_trimmed}
    SAMPLE=$(basename "$R1" _R1_trimmed.fastq.gz)

    hisat2 \
      -p ${THREADS} \
      -x ${HISAT2_INDEX} \
      -1 "$R1" \
      -2 "$R2" \
      --dta \
      2> "${BAM_DIR}/${SAMPLE}_hisat2.log" | \
    samtools sort -@ ${THREADS} -o "${BAM_DIR}/${SAMPLE}.sorted.bam"

    samtools index "${BAM_DIR}/${SAMPLE}.sorted.bam"
done
