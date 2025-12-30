# RNA-seq Differential Gene Expression Workflow

This repository contains a transparent, step-by-step RNA-seq differential gene
expression (DGE) analysis workflow for human paired-end RNA-seq data. The
workflow performs gene-level analysis using fastp, HISAT2, featureCounts, and
DESeq2, with downstream visualization and pathway enrichment.

The pipeline is designed for **exploratory gene-level differential expression**
analysis and emphasizes reproducibility, interpretability, and explicit
parameter choices.

---

## Overview of the workflow

1. Adapter trimming and quality filtering (fastp)
2. Read-level QC (FastQC)
3. Splice-aware alignment to the human genome (HISAT2)
4. BAM processing and indexing (SAMtools)
5. Gene-level quantification (featureCounts)
6. Removal of rRNA genes and low-count filtering
7. Differential expression analysis (DESeq2)
8. Log2 fold-change shrinkage (apeglm)
9. Visualization (PCA, volcano plots, heatmaps)
10. Pathway enrichment analysis (GO BP and KEGG)

---

## Reference genome and annotation

- Genome: GRCh38 (primary assembly)
- Annotation: GENCODE v44
- Species: Homo sapiens

---

## Software and tools

| Step | Tool |
|----|----|
| Trimming & QC | fastp |
| Read QC | FastQC |
| Alignment | HISAT2 |
| BAM processing | SAMtools |
| Quantification | featureCounts (Subread) |
| Differential expression | DESeq2 |
| Log2FC shrinkage | apeglm |
| Visualization | ggplot2, pheatmap, EnhancedVolcano |
| Pathway analysis | clusterProfiler |
| Annotation | org.Hs.eg.db |

---

## Filtering and normalization

- rRNA genes were removed post-quantification using GENCODE annotations
- Low-count genes were filtered (≥10 reads in ≥2 samples)
- DESeq2 median-of-ratios normalization was applied
- Log2 fold-change shrinkage was performed using apeglm

---

## Differential expression analysis

- Experimental design: GroupA vs GroupB
- Statistical testing: DESeq2 Wald test
- Multiple testing correction: Benjamini–Hochberg
- Significance threshold: FDR < 0.1 (exploratory)
- Effect size emphasis: |log2FC| ≥ 1

---

## Pathway enrichment analysis

Over-representation analysis was performed separately for upregulated and
downregulated genes using Gene Ontology Biological Process and KEGG pathways.
Differentially expressed genes were defined using FDR < 0.05 and |log2FC| ≥ 1.
All expressed genes were used as the background universe. The top 10 pathways
(rank-ordered by adjusted p-value) were visualized for each analysis.

---

## Notes and limitations

- This workflow is optimized for **gene-level** analysis
- Transcript-level quantification and alternative splicing are not addressed
- The example dataset uses a small sample size (n = 2 per group); results
  should be interpreted as exploratory and validated independently

---

## Citation

If you use this workflow, please cite the relevant software tools and consider
citing this repository (see `CITATION.cff`).

---

## License

MIT License
