# ============================================================
# Script name: 04_deseq2_analysis.R
#
# Description:
#   Differential gene expression analysis using DESeq2
#   with rRNA removal, low-count filtering, and
#   log2 fold-change shrinkage.
#
# Author: Krishan Kumar
# License: MIT
# ============================================================

# -------- User-configurable variables --------
project_dir <- "path/to/project"
count_file  <- file.path(project_dir, "counts/gene_counts.txt")
gtf_file    <- file.path(project_dir, "annotation/gencode.v44.annotation.gtf")
out_dir     <- file.path(project_dir, "results")
plot_dir    <- file.path(project_dir, "plots")
# --------------------------------------------

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

library(data.table)
library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(AnnotationDbi)

# -----------------------------
# Load featureCounts output
# -----------------------------
counts_raw <- fread(count_file, skip = "#")

# -----------------------------
# Remove rRNA genes
# -----------------------------
gtf <- fread(gtf_file, sep = "\t", header = FALSE)
colnames(gtf) <- c("chr","src","feat","start","end","score",
                   "strand","frame","attr")

gtf[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attr)]
gtf[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attr)]

rRNA_genes <- unique(gtf[gene_type == "rRNA", gene_id])

counts_no_rRNA <- counts_raw[!Geneid %in% rRNA_genes]

# -----------------------------
# Low-count filtering
# -----------------------------
keep <- rowSums(counts_no_rRNA[, -(1:6)] >= 10) >= 2
counts_filt <- counts_no_rRNA[keep]

count_matrix <- as.matrix(counts_filt[, -(1:6)])
rownames(count_matrix) <- counts_filt$Geneid

# -----------------------------
# Metadata
# -----------------------------
# NOTE:
# Replace 'GroupA' and 'GroupB' with your experimental groups.
# The first level (GroupB) is used as the reference for log2FC.

coldata <- data.frame(
  group = factor(
    c("GroupA", "GroupA", "GroupB", "GroupB"),
    levels = c("GroupB", "GroupA")
  ),
  row.names = colnames(count_matrix)
)

# Sanity check
print(colnames(count_matrix))
print(coldata)

# -----------------------------
# DESeq2 analysis
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ group
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "GroupA", "GroupB"))

# Log2FC shrinkage
resLFC <- lfcShrink(
  dds,
  coef = "group_GroupA_vs_GroupB",
  type = "apeglm"
)

# -----------------------------
# Annotate results
# -----------------------------
res_out <- as.data.frame(resLFC)
res_out$gene_id <- rownames(res_out)

# Remove Ensembl version for mapping
res_out$gene_id_clean <- sub("\\..*", "", res_out$gene_id)

res_out$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_out$gene_id_clean,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# Reorder columns
res_out <- res_out[, c(
  "gene_id",
  "symbol",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj"
)]

# Order by adjusted p-value
res_out <- res_out[order(res_out$padj), ]

# -----------------------------
# Save results
# -----------------------------
write.csv(
  res_out,
  file.path(out_dir, "dge_GroupA_vs_GroupB.csv"),
  row.names = FALSE
)

# Optional: save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
