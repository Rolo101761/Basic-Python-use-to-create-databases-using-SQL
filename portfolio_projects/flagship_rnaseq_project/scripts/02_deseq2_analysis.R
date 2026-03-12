# =============================================================================
# Step 2 — Differential Expression Analysis (DESeq2)
# RNA-seq Differential Expression Pipeline
# Dataset: GSE183947 (breast cancer tumour vs normal)
# =============================================================================
#
# Reads cleaned count matrix and sample metadata produced by 01_data_download.py,
# runs DESeq2 differential expression analysis, and produces:
#
#   data/processed/DESeq2_results.csv   full results table (all genes)
#   data/processed/DESeq2_sig.csv       significant genes (padj<0.05, |log2FC|>1)
#   data/processed/vst_counts.csv       variance-stabilised counts for PCA/heatmap
#   figures/pca_plot.png                PCA of VST-normalised counts
#   figures/volcano_plot.png            volcano plot coloured by significance
#   figures/heatmap.png                 heatmap of top 50 DE genes
#
# Usage:
#   Rscript scripts/02_deseq2_analysis.R
#   Rscript scripts/02_deseq2_analysis.R --indir data/processed --figdir figures
#
# Requirements:
#   BiocManager::install(c("DESeq2", "pheatmap", "ComplexHeatmap"))
#   install.packages(c("ggplot2", "ggrepel", "RColorBrewer"))
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# Parse command line arguments
# ---------------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
in_dir    <- if ("--indir"  %in% args) args[which(args == "--indir")  + 1] else "data/processed"
fig_dir   <- if ("--figdir" %in% args) args[which(args == "--figdir") + 1] else "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("=============================================================\n")
cat("DESeq2 Differential Expression — GSE183947\n")
cat("=============================================================\n\n")
cat("[1/6] Loading count matrix and metadata...\n")

counts_path   <- file.path(in_dir, "cleaned_counts.tsv")
metadata_path <- file.path(in_dir, "sample_metadata.tsv")

if (!file.exists(counts_path)) stop(sprintf("Count matrix not found: %s\nRun 01_data_download.py first.", counts_path))
if (!file.exists(metadata_path)) stop(sprintf("Metadata not found: %s\nRun 01_data_download.py first.", metadata_path))

counts   <- read.delim(counts_path,   row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)

# Align samples
shared   <- intersect(colnames(counts), rownames(metadata))
counts   <- counts[, shared]
metadata <- metadata[shared, , drop = FALSE]

# Ensure condition is a factor with tumour as the comparison level
metadata$condition <- factor(metadata$condition, levels = c("normal", "tumour"))

cat(sprintf("      Genes   : %d\n", nrow(counts)))
cat(sprintf("      Samples : %d  (%s)\n", ncol(counts),
            paste(table(metadata$condition), names(table(metadata$condition)), collapse = " / ")))

# ---------------------------------------------------------------------------
# Build DESeq2 object
# ---------------------------------------------------------------------------
cat("\n[2/6] Building DESeq2 object and fitting model...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Pre-filter: keep genes with >= 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat(sprintf("      Genes after pre-filter: %d\n", nrow(dds)))

# Run DESeq2
dds <- DESeq(dds, quiet = TRUE)

# Results: tumour vs normal
res <- results(dds, contrast = c("condition", "tumour", "normal"),
               alpha = 0.05)
res <- res[order(res$padj, na.last = TRUE), ]

res_df       <- as.data.frame(res)
res_df$gene  <- rownames(res_df)

n_sig <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE)
n_up  <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
n_dn  <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)

cat(sprintf("      Significant DE genes (padj<0.05, |log2FC|>1): %d\n", n_sig))
cat(sprintf("        Upregulated in tumour  : %d\n", n_up))
cat(sprintf("        Downregulated in tumour: %d\n", n_dn))

# Save full results
write.csv(res_df, file.path(in_dir, "DESeq2_results.csv"), row.names = FALSE)

# Save significant genes only
sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.csv(sig_df, file.path(in_dir, "DESeq2_sig.csv"), row.names = FALSE)

cat(sprintf("\n      Saved: %s\n", file.path(in_dir, "DESeq2_results.csv")))
cat(sprintf("      Saved: %s\n", file.path(in_dir, "DESeq2_sig.csv")))

# ---------------------------------------------------------------------------
# VST normalisation (for PCA and heatmap)
# ---------------------------------------------------------------------------
cat("\n[3/6] Variance-stabilising transformation (VST)...\n")
vst <- vst(dds, blind = FALSE)
vst_mat <- assay(vst)
write.csv(as.data.frame(vst_mat), file.path(in_dir, "vst_counts.csv"))
cat(sprintf("      Saved: %s\n", file.path(in_dir, "vst_counts.csv")))

# ---------------------------------------------------------------------------
# PCA plot
# ---------------------------------------------------------------------------
cat("\n[4/6] PCA plot...\n")

pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition, label = name)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(size = 2.8, max.overlaps = 15) +
  scale_colour_manual(values = c("normal" = "#2166ac", "tumour" = "#d6604d")) +
  labs(
    title    = "PCA — VST-normalised counts (GSE183947)",
    subtitle = "Breast cancer: tumour vs normal",
    x        = sprintf("PC1 (%d%% variance)", pct_var[1]),
    y        = sprintf("PC2 (%d%% variance)", pct_var[2]),
    colour   = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "pca_plot.png"), p_pca, width = 7, height = 5, dpi = 300)
cat(sprintf("      Saved: %s\n", file.path(fig_dir, "pca_plot.png")))

# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------
cat("\n[5/6] Volcano plot...\n")

plot_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
plot_df$sig <- "NS"
plot_df$sig[plot_df$padj < 0.05 & plot_df$log2FoldChange >  1] <- "Up"
plot_df$sig[plot_df$padj < 0.05 & plot_df$log2FoldChange < -1] <- "Down"
plot_df$sig <- factor(plot_df$sig, levels = c("Up", "Down", "NS"))

# Label top 15 most significant genes
plot_df$label <- NA
top_genes      <- head(plot_df[plot_df$sig != "NS", "gene"], 15)
plot_df$label[plot_df$gene %in% top_genes] <- plot_df$gene[plot_df$gene %in% top_genes]

p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  box.padding = 0.4, fontface = "italic") +
  scale_colour_manual(values = c("Up" = "#d6604d", "Down" = "#2166ac", "NS" = "#cccccc")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", alpha = 0.7) +
  annotate("text", x = max(plot_df$log2FoldChange) * 0.85,
           y = max(-log10(plot_df$padj), na.rm = TRUE) * 0.97,
           label = sprintf("Up: %d", n_up), colour = "#d6604d", size = 4) +
  annotate("text", x = min(plot_df$log2FoldChange) * 0.85,
           y = max(-log10(plot_df$padj), na.rm = TRUE) * 0.97,
           label = sprintf("Down: %d", n_dn), colour = "#2166ac", size = 4) +
  labs(
    title    = "Tumour vs Normal — Differential Expression (DESeq2)",
    subtitle = "GSE183947 breast cancer | padj < 0.05 & |log2FC| > 1",
    x        = "log2 Fold Change (tumour / normal)",
    y        = "-log10(adjusted p-value)",
    colour   = NULL
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(fig_dir, "volcano_plot.png"), p_vol, width = 8, height = 6, dpi = 300)
cat(sprintf("      Saved: %s\n", file.path(fig_dir, "volcano_plot.png")))

# ---------------------------------------------------------------------------
# Heatmap — top 50 DE genes
# ---------------------------------------------------------------------------
cat("\n[6/6] Heatmap of top 50 DE genes...\n")

top50 <- head(sig_df$gene, 50)
top50 <- top50[top50 %in% rownames(vst_mat)]

if (length(top50) >= 2) {
  mat        <- vst_mat[top50, ]
  mat_scaled <- t(scale(t(mat)))   # z-score per gene

  ann_col <- data.frame(
    condition = metadata$condition,
    row.names = rownames(metadata)
  )
  ann_colours <- list(condition = c(normal = "#2166ac", tumour = "#d6604d"))

  pheatmap(
    mat_scaled,
    annotation_col  = ann_col,
    annotation_colors = ann_colours,
    color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    show_rownames   = length(top50) <= 50,
    show_colnames   = TRUE,
    fontsize_row    = 7,
    fontsize_col    = 8,
    main            = "Top 50 DE Genes — z-scored VST counts\nGSE183947 Breast Cancer",
    filename        = file.path(fig_dir, "heatmap.png"),
    width = 10, height = 10
  )
  cat(sprintf("      Saved: %s\n", file.path(fig_dir, "heatmap.png")))
} else {
  cat("      Not enough significant genes for heatmap (need >= 2)\n")
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 55), "\n", sep = "")
cat("DESeq2 analysis complete\n")
cat(strrep("=", 55), "\n", sep = "")
cat(sprintf("  Total genes tested  : %d\n", nrow(res_df)))
cat(sprintf("  Significant (DE)    : %d\n", n_sig))
cat(sprintf("  Upregulated         : %d\n", n_up))
cat(sprintf("  Downregulated       : %d\n", n_dn))
cat("\nNext step:\n")
cat("  Rscript scripts/03_pathway_enrichment.R\n")
