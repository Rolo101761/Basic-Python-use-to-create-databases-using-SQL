# =============================================================================
# Step 3 — Pathway Enrichment Analysis (clusterProfiler)
# RNA-seq Differential Expression Pipeline
# Dataset: GSE183947 (breast cancer tumour vs normal)
# =============================================================================
#
# Reads DESeq2 significant gene list and runs:
#   - Gene Ontology (GO) enrichment — Biological Process, Molecular Function
#   - KEGG pathway enrichment
#
# Outputs:
#   data/processed/GO_BP_enrichment.csv    GO Biological Process results
#   data/processed/GO_MF_enrichment.csv    GO Molecular Function results
#   data/processed/KEGG_enrichment.csv     KEGG pathway results
#   data/processed/enrichment_results.tsv  combined table for SQLite/dashboard
#   figures/go_dotplot.png                 GO BP top 20 dot plot
#   figures/kegg_dotplot.png               KEGG top 20 dot plot
#   figures/enrichment_barplot.png         top GO terms coloured by direction
#
# Usage:
#   Rscript scripts/03_pathway_enrichment.R
#   Rscript scripts/03_pathway_enrichment.R --indir data/processed --figdir figures
#
# Requirements:
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if ("--indir"  %in% args) args[which(args == "--indir")  + 1] else "data/processed"
fig_dir <- if ("--figdir" %in% args) args[which(args == "--figdir") + 1] else "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load DE results
# ---------------------------------------------------------------------------
cat("=============================================================\n")
cat("Pathway Enrichment Analysis — GSE183947\n")
cat("=============================================================\n\n")
cat("[1/5] Loading DESeq2 results...\n")

sig_path <- file.path(in_dir, "DESeq2_sig.csv")
all_path <- file.path(in_dir, "DESeq2_results.csv")

if (!file.exists(sig_path)) stop(sprintf("Significant genes file not found: %s\nRun 02_deseq2_analysis.R first.", sig_path))

sig_df <- read.csv(sig_path)
all_df <- read.csv(all_path)

cat(sprintf("      Significant DE genes: %d\n", nrow(sig_df)))

# Separate up and down for directional enrichment
up_genes   <- sig_df$gene[sig_df$log2FoldChange > 0]
down_genes <- sig_df$gene[sig_df$log2FoldChange < 0]
cat(sprintf("      Upregulated  : %d\n", length(up_genes)))
cat(sprintf("      Downregulated: %d\n", length(down_genes)))

# ---------------------------------------------------------------------------
# Gene ID conversion: symbol → Entrez
# ---------------------------------------------------------------------------
cat("\n[2/5] Converting gene symbols to Entrez IDs...\n")

convert_ids <- function(genes) {
  tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) {
      message("  Warning: bitr conversion failed: ", e$message)
      data.frame(SYMBOL = character(), ENTREZID = character())
    }
  )
}

sig_ids  <- convert_ids(sig_df$gene)
all_ids  <- convert_ids(all_df$gene[!is.na(all_df$gene)])

cat(sprintf("      Converted: %d / %d significant genes\n", nrow(sig_ids), nrow(sig_df)))

# Background gene universe (all tested genes)
universe <- all_ids$ENTREZID

# ---------------------------------------------------------------------------
# GO enrichment — Biological Process
# ---------------------------------------------------------------------------
cat("\n[3/5] GO Biological Process enrichment...\n")

run_go <- function(gene_ids, universe, ont, label) {
  if (length(gene_ids) < 5) {
    cat(sprintf("      Skipping %s — fewer than 5 genes\n", label))
    return(NULL)
  }
  res <- tryCatch(
    enrichGO(
      gene          = gene_ids,
      universe      = universe,
      OrgDb         = org.Hs.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    ),
    error = function(e) { message("  Warning: GO failed: ", e$message); NULL }
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat(sprintf("      No significant %s terms\n", label))
    return(NULL)
  }
  cat(sprintf("      Significant %s terms: %d\n", label, nrow(as.data.frame(res))))
  res
}

go_bp <- run_go(sig_ids$ENTREZID, universe, "BP", "GO-BP")
go_mf <- run_go(sig_ids$ENTREZID, universe, "MF", "GO-MF")

# Save GO results
if (!is.null(go_bp)) {
  write.csv(as.data.frame(go_bp), file.path(in_dir, "GO_BP_enrichment.csv"), row.names = FALSE)
  cat(sprintf("      Saved: %s\n", file.path(in_dir, "GO_BP_enrichment.csv")))
}
if (!is.null(go_mf)) {
  write.csv(as.data.frame(go_mf), file.path(in_dir, "GO_MF_enrichment.csv"), row.names = FALSE)
  cat(sprintf("      Saved: %s\n", file.path(in_dir, "GO_MF_enrichment.csv")))
}

# ---------------------------------------------------------------------------
# KEGG pathway enrichment
# ---------------------------------------------------------------------------
cat("\n[4/5] KEGG pathway enrichment...\n")

kegg_res <- tryCatch(
  enrichKEGG(
    gene         = sig_ids$ENTREZID,
    universe     = universe,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  ),
  error = function(e) { message("  Warning: KEGG failed: ", e$message); NULL }
)

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_df <- as.data.frame(kegg_res)
  # Convert Entrez IDs back to gene symbols in the geneID column
  kegg_df$geneSymbols <- sapply(kegg_df$geneID, function(ids) {
    entrez_ids <- strsplit(ids, "/")[[1]]
    symbols    <- mapIds(org.Hs.eg.db, keys = entrez_ids,
                         column = "SYMBOL", keytype = "ENTREZID",
                         multiVals = "first")
    paste(symbols, collapse = "/")
  })
  write.csv(kegg_df, file.path(in_dir, "KEGG_enrichment.csv"), row.names = FALSE)
  cat(sprintf("      Significant KEGG pathways: %d\n", nrow(kegg_df)))
  cat(sprintf("      Saved: %s\n", file.path(in_dir, "KEGG_enrichment.csv")))
} else {
  cat("      No significant KEGG pathways\n")
  kegg_df <- data.frame()
}

# ---------------------------------------------------------------------------
# Combined enrichment table for SQL/dashboard
# ---------------------------------------------------------------------------
combined_rows <- list()

if (!is.null(go_bp)) {
  bp_df <- as.data.frame(go_bp)[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID")]
  bp_df$source <- "GO_BP"
  combined_rows[["go_bp"]] <- bp_df
}
if (!is.null(go_mf)) {
  mf_df <- as.data.frame(go_mf)[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID")]
  mf_df$source <- "GO_MF"
  combined_rows[["go_mf"]] <- mf_df
}
if (nrow(kegg_df) > 0) {
  kk_df <- kegg_df[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID")]
  kk_df$source <- "KEGG"
  combined_rows[["kegg"]] <- kk_df
}

if (length(combined_rows) > 0) {
  combined_df <- do.call(rbind, combined_rows)
  write.table(combined_df, file.path(in_dir, "enrichment_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\n      Saved combined: %s (%d terms)\n",
              file.path(in_dir, "enrichment_results.tsv"), nrow(combined_df)))
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
cat("\n[5/5] Generating enrichment figures...\n")

# GO BP dot plot
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p_go <- dotplot(go_bp, showCategory = 20) +
    labs(title = "GO Biological Process Enrichment",
         subtitle = "GSE183947 breast cancer DE genes (tumour vs normal)") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(fig_dir, "go_dotplot.png"), p_go, width = 10, height = 8, dpi = 300)
  cat(sprintf("      Saved: %s\n", file.path(fig_dir, "go_dotplot.png")))
}

# KEGG dot plot
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  p_kegg <- dotplot(kegg_res, showCategory = 20) +
    labs(title = "KEGG Pathway Enrichment",
         subtitle = "GSE183947 breast cancer DE genes (tumour vs normal)") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(fig_dir, "kegg_dotplot.png"), p_kegg, width = 10, height = 7, dpi = 300)
  cat(sprintf("      Saved: %s\n", file.path(fig_dir, "kegg_dotplot.png")))
}

# Horizontal bar chart — top 20 GO-BP terms by -log10(padj)
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) >= 5) {
  bp_bar <- head(as.data.frame(go_bp), 20)
  bp_bar$Description <- factor(bp_bar$Description,
                                levels = rev(bp_bar$Description))
  p_bar <- ggplot(bp_bar, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
    geom_col() +
    scale_fill_gradient(low = "#c6dbef", high = "#2166ac") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "red", alpha = 0.6) +
    labs(title    = "Top 20 GO Biological Process Terms",
         subtitle = "GSE183947 | dashed line = padj 0.05",
         x        = "-log10(adjusted p-value)",
         y        = NULL,
         fill     = "Gene\ncount") +
    theme_bw(base_size = 11)
  ggsave(file.path(fig_dir, "enrichment_barplot.png"), p_bar, width = 11, height = 7, dpi = 300)
  cat(sprintf("      Saved: %s\n", file.path(fig_dir, "enrichment_barplot.png")))
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 55), "\n", sep = "")
cat("Pathway enrichment complete\n")
cat(strrep("=", 55), "\n", sep = "")
if (!is.null(go_bp))  cat(sprintf("  GO-BP significant terms : %d\n", nrow(as.data.frame(go_bp))))
if (!is.null(go_mf))  cat(sprintf("  GO-MF significant terms : %d\n", nrow(as.data.frame(go_mf))))
if (nrow(kegg_df) > 0) cat(sprintf("  KEGG significant pathways: %d\n", nrow(kegg_df)))
cat("\nNext step:\n")
cat("  python scripts/04_load_to_db.py\n")
