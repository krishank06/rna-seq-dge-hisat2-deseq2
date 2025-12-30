# ============================================================
# Script name: 05_pathway_enrichment.R
#
# Description:
#   Pathway enrichment analysis using over-representation
#   analysis (ORA) for Gene Ontology Biological Process (GO BP)
#   and KEGG pathways. Upregulated and downregulated genes
#   are analyzed separately.
#
# Author: Krishan Kumar
# License: MIT
# ============================================================

# -------- User-configurable variables --------
project_dir <- "path/to/project"
dge_file    <- file.path(project_dir, "results/dge_GroupA_vs_GroupB.csv")
out_dir     <- file.path(project_dir, "results")
plot_dir    <- file.path(project_dir, "plots")
# --------------------------------------------

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# -----------------------------
# Load DGE results
# -----------------------------
res <- read.csv(dge_file)

# Remove rows without gene symbols
res <- res %>% filter(!is.na(symbol))

# -----------------------------
# Define background universe
# -----------------------------
background_genes <- unique(res$symbol)

# -----------------------------
# Define DE gene sets
# -----------------------------
# NOTE:
# GroupA is compared against GroupB.
# Positive log2FC = higher expression in GroupA
# Negative log2FC = lower expression in GroupA

genes_up <- res %>%
  filter(padj < 0.05, log2FoldChange >= 1) %>%
  pull(symbol) %>%
  unique()

genes_down <- res %>%
  filter(padj < 0.05, log2FoldChange <= -1) %>%
  pull(symbol) %>%
  unique()

# -----------------------------
# GO Biological Process enrichment
# -----------------------------
ego_up <- enrichGO(
  gene          = genes_up,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

ego_down <- enrichGO(
  gene          = genes_down,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# -----------------------------
# KEGG enrichment
# -----------------------------
# Convert SYMBOL → ENTREZ
gene_map <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

up_entrez <- gene_map$ENTREZID[gene_map$SYMBOL %in% genes_up]
down_entrez <- gene_map$ENTREZID[gene_map$SYMBOL %in% genes_down]

ekegg_up <- enrichKEGG(
  gene         = up_entrez,
  universe     = gene_map$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_down <- enrichKEGG(
  gene         = down_entrez,
  universe     = gene_map$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# -----------------------------
# Save enrichment tables
# -----------------------------
write.csv(as.data.frame(ego_up),
          file.path(out_dir, "GO_BP_UP_GroupA.csv"),
          row.names = FALSE)

write.csv(as.data.frame(ego_down),
          file.path(out_dir, "GO_BP_DOWN_GroupA.csv"),
          row.names = FALSE)

write.csv(as.data.frame(ekegg_up),
          file.path(out_dir, "KEGG_UP_GroupA.csv"),
          row.names = FALSE)

write.csv(as.data.frame(ekegg_down),
          file.path(out_dir, "KEGG_DOWN_GroupA.csv"),
          row.names = FALSE)

# -----------------------------
# Visualization: top 10 pathways
# -----------------------------
png(file.path(plot_dir, "GO_BP_UP_dotplot.png"),
    width = 1800, height = 1400, res = 300)
dotplot(ego_up, showCategory = 10) +
  ggtitle("GO BP enrichment (Upregulated in GroupA)")
dev.off()

png(file.path(plot_dir, "GO_BP_DOWN_dotplot.png"),
    width = 1800, height = 1400, res = 300)
dotplot(ego_down, showCategory = 10) +
  ggtitle("GO BP enrichment (Downregulated in GroupA)")
dev.off()

png(file.path(plot_dir, "KEGG_UP_dotplot.png"),
    width = 1800, height = 1400, res = 300)
dotplot(ekegg_up, showCategory = 10) +
  ggtitle("KEGG enrichment (Upregulated in GroupA)")
dev.off()

png(file.path(plot_dir, "KEGG_DOWN_dotplot.png"),
    width = 1800, height = 1400, res = 300)
dotplot(ekegg_down, showCategory = 10) +
  ggtitle("KEGG enrichment (Downregulated in GroupA)")
dev.off()
