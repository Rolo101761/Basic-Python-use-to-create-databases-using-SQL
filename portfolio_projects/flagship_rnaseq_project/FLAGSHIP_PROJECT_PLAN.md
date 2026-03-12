# Flagship Project: RNA-seq Differential Expression Analysis Pipeline

**Target roles:** Industry (pharma/biotech) | Clinical bioinformatics | Data scientist (bio focus)

---

## Project overview

An end-to-end RNA-seq analysis pipeline spanning Python, R, SQL, and Power BI/Tableau.
Uses publicly available TCGA or GEO data so it is fully reproducible by any employer.

---

## Recommended dataset

**GEO accession GSE62944** — TCGA RNA-seq data (tumour vs normal) across 9,264 samples.
Or use a smaller starter dataset: **GSE183947** (breast cancer, ~20 samples, quick to run).

Download with R:
```r
# BiocManager::install("GEOquery")
library(GEOquery)
gse <- getGEO("GSE183947", GSEMatrix = TRUE)
```

---

## Pipeline architecture

```
Raw counts (GEO/TCGA)
        │
        ▼
[Step 1 — Python]  Data ingestion + QC
  - Download via GEOparse or recount3
  - Check library sizes, outlier detection
  - Output: cleaned_counts.tsv, sample_metadata.tsv

        │
        ▼
[Step 2 — R]  Differential expression (DESeq2)
  - Normalisation (VST / rlog)
  - PCA plot
  - DESeq2 model fitting
  - Volcano plot (ggplot2)
  - Heatmap (pheatmap / ComplexHeatmap)
  - Output: DESeq2_results.tsv

        │
        ▼
[Step 3 — R]  Pathway enrichment
  - Gene Ontology (clusterProfiler)
  - KEGG pathway analysis
  - Dot plots, enrichment maps
  - Output: enrichment_results.tsv

        │
        ▼
[Step 4 — Python + SQL]  Results database
  - Load DESeq2 results into SQLite (bioinformatics_db.py)
  - Store enrichment results, sample metadata
  - Query: top DE genes per pathway, per condition

        │
        ▼
[Step 5 — Power BI / Tableau]  Interactive dashboard
  - Connect to SQLite via ODBC or export to CSV
  - Dashboard pages:
      • Sample overview (PCA, library sizes)
      • DE gene table (filterable by log2FC, padj)
      • Pathway enrichment summary
      • Variant-expression correlation (if TCGA somatic data included)
```

---

## File structure (suggested repo layout)

```
rnaseq-de-pipeline/
├── README.md
├── data/
│   ├── raw/                    # raw count matrices (gitignored if large)
│   └── processed/              # cleaned_counts.tsv, sample_metadata.tsv
├── scripts/
│   ├── 01_data_download.py     # Python: GEO download + QC
│   ├── 02_deseq2_analysis.R    # R: DESeq2 DE analysis
│   ├── 03_pathway_enrichment.R # R: clusterProfiler enrichment
│   ├── 04_load_to_db.py        # Python: load results into SQLite
│   └── 05_db_queries.py        # Python: example queries
├── figures/
│   ├── pca_plot.png
│   ├── volcano_plot.png
│   ├── heatmap.png
│   └── enrichment_dotplot.png
├── database/
│   └── rnaseq_results.db       # SQLite results database
├── dashboard/
│   └── rnaseq_dashboard.pbix   # Power BI dashboard (or .twbx for Tableau)
└── requirements.txt
```

---

## R packages needed

```r
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "ggrepel"))
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "GEOquery", "ComplexHeatmap"))
```

## Python packages needed

```
pip install pandas numpy matplotlib seaborn GEOparse requests
```

---

## Key R code snippets

### DESeq2 analysis (02_deseq2_analysis.R)

```r
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load data
counts   <- read.delim("data/processed/cleaned_counts.tsv", row.names=1)
metadata <- read.delim("data/processed/sample_metadata.tsv", row.names=1)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumour", "normal"))
res <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res), "data/processed/DESeq2_results.csv")

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene     <- rownames(res_df)
res_df$sig      <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                          "Significant", "NS")
res_df$label    <- ifelse(res_df$sig == "Significant" &
                          abs(res_df$log2FoldChange) > 3, res_df$gene, NA)

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), colour=sig)) +
  geom_point(alpha=0.6, size=1.5) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=20) +
  scale_colour_manual(values=c("Significant"="red", "NS"="grey60")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  labs(title="Tumour vs Normal — Differential Expression",
       x="log2 Fold Change", y="-log10(adjusted p-value)") +
  theme_bw()
ggsave("figures/volcano_plot.png", width=8, height=6, dpi=300)
```

### Pathway enrichment (03_pathway_enrichment.R)

```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

sig_genes <- rownames(res_df[res_df$sig == "Significant", ])

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO enrichment
go_res <- enrichGO(
  gene         = gene_ids$ENTREZID,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",          # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

dotplot(go_res, showCategory=20) + ggtitle("GO Biological Process Enrichment")
ggsave("figures/enrichment_dotplot.png", width=10, height=8, dpi=300)

# KEGG pathway enrichment
kegg_res <- enrichKEGG(gene=gene_ids$ENTREZID, organism='hsa', pvalueCutoff=0.05)
write.csv(as.data.frame(kegg_res), "data/processed/KEGG_enrichment.csv")
```

---

## Power BI / Tableau dashboard pages

### Page 1: Sample QC Overview
- Bar chart: library sizes per sample (coloured by condition)
- Scatter plot: PCA plot (PC1 vs PC2, coloured by condition)
- KPI cards: total samples, conditions, sequencing depth range

### Page 2: Differential Expression
- Interactive volcano plot (filter by log2FC and padj sliders)
- Table: top 50 DE genes with symbol, log2FC, padj, biotype
- Bar chart: upregulated vs downregulated gene counts

### Page 3: Pathway Enrichment
- Horizontal bar chart: top 20 GO terms by -log10(padj)
- KEGG pathway table with gene ratio and count
- Word cloud (Tableau) or treemap (Power BI) of pathway categories

### Page 4: Gene Deep-Dive
- Slicer/filter: select any gene
- Bar chart: expression (TPM) across all samples
- Annotation card: gene description, biotype, chromosome

---

## Why this project impresses employers

| Employer type | What they see |
|--------------|---------------|
| **Pharma/biotech** | End-to-end NGS pipeline, DESeq2, pathway analysis, reproducibility |
| **Clinical bioinformatics** | Variant/expression integration, clinically relevant genes (BRCA1, TP53) |
| **Data scientist** | Python + R + SQL + BI tool stack, clean code, visualisations |
| **All** | Public data (reproducible), GitHub-ready, README with figures |
