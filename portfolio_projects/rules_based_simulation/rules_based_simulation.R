# =============================================================================
# Rules-Based Writer/Eraser Simulation of Histone Methylation Cross-Talk
# Using GenomicLayers (Gerrard, BMC Bioinformatics 2025)
# =============================================================================
#
# Biological focus
# ================
# Models the antagonistic cross-talk between H3K4me3 (active) and H3K27me3
# (repressive Polycomb) marks using three mechanistically grounded rules:
#
#   Rule 1  Writer is blocked by repressive mark (mutual exclusivity)
#           → H3K4me3 cannot be deposited where H3K27me3 already exists
#           → Implements the observed genome-wide anti-correlation of these marks
#
#   Rule 2  PRC2 is inhibited by H3K36me3 (in-cis protection)
#           → H3K27me3 cannot spread into actively transcribed gene bodies
#           → Based on: Molecular Cell 2010, Nature Struct Mol Biol 2025
#           → This is a key challenge for genome-scale modelling
#
#   Rule 3  EED reader-writer feedback (positive feedback for H3K27me3)
#           → PRC2 reads existing H3K27me3 via EED WD40 domain
#           → Allosteric activation ~7-fold (Margueron et al. 2009)
#           → Creates bistability: once established, repression is self-sustaining
#
#   Rule 4  Set2-coupled H3K36me3 deposition (transcription feedback)
#           → Active mark (H3K4me3) recruits transcription, which recruits Set2
#           → Set2 deposits H3K36me3 in gene bodies → reinforces Rule 2
#
#   Rule 5  KDM6A/B (UTX) erases H3K27me3 at active loci
#           → Active promoter → recruits KDM6 demethylase → removes H3K27me3
#           → Reinforces mutual exclusivity
#
# Layers
# ======
#   LAYER.1  Active mark        (H3K4me3-like)
#   LAYER.2  Repressive mark    (H3K27me3-like)
#   LAYER.3  Protection mark    (H3K36me3-like, transcription-coupled)
#
# Novel prediction (biologically challenging)
# ==========================================
# By running from DNA sequence alone (no expression data), the model attempts
# to predict which gene promoters will be:
#   (a) Active    — TATA-box near CpG-free context, H3K36me3-protected
#   (b) Repressed — CpG-rich seed, no TATA-box, H3K36me3 absent
#   (c) Bivalent  — CpG-rich + TATA-box context, no H3K36me3 present
#                   (poised developmental gene — the unsolved prediction)
#   (d) Unmarked  — Neither writer nor repressor recruited
#
# Bivalent domain prediction is one of the open challenges in epigenomics:
# sequence features alone are insufficient in current models. This script
# tests whether the bistability rules can delineate bivalent loci from
# sequence context.
#
# Cross-talk score
# ================
# Quantifies how strongly the three interaction rules enforce mutual
# exclusivity of L1 and L2 across the genome. Without the rules (naive model)
# active and repressive marks co-localise randomly. With rules, the Pearson
# correlation of binned L1/L2 should approach -1.
#
# Installation (run once):
#   install.packages("devtools")
#   BiocManager::install(c("BSgenome", "GenomicRanges", "Biostrings",
#                          "GenomicFeatures", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"))
#   BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
#   devtools::install_github("davetgerrard/GenomicLayers", build_vignettes = TRUE)
#
# Usage:
#   Rscript rules_based_simulation.R             # default: full cross-talk model
#   Rscript rules_based_simulation.R --naive     # naive model (no interaction rules)
#   Rscript rules_based_simulation.R --compare   # both models, head-to-head
#   Rscript rules_based_simulation.R --sweep     # parameter sweep
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicLayers)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
params <- list(
  chromosome     = "chrI",    # S. cerevisiae chrI (230 kb — tractable)
  n_layers       = 3L,
  cycles         = 25L,
  mods_per_cycle = 600L,
  bin_size       = 1000L,     # 1 kb bins for cross-talk analysis
  promoter_up    = 500L,      # bp upstream of TSS for promoter classification
  rng_seed       = 42L
)

run_naive   <- "--naive"   %in% args || "--compare" %in% args
run_crosstalk <- !("--naive" %in% args)

set.seed(params$rng_seed)

# -----------------------------------------------------------------------------
# Helper: summarise layer state at current cycle
# -----------------------------------------------------------------------------
summarise_layers <- function(layerSet, chr) {
  lyr <- layerSet[[chr]]$layers
  data.frame(
    frac_active    = mean(lyr[, 1]),
    frac_repressed = mean(lyr[, 2]),
    frac_protected = mean(lyr[, 3]),
    frac_bivalent  = mean(lyr[, 1] & lyr[, 2]),
    frac_unmarked  = mean(!lyr[, 1] & !lyr[, 2])
  )
}

# -----------------------------------------------------------------------------
# Helper: run simulation and collect timeseries
# -----------------------------------------------------------------------------
run_simulation <- function(layerSet, factorSet, params, label) {
  cat(sprintf("  Running %s model: %d cycles x %d events...\n",
              label, params$cycles, params$mods_per_cycle))
  ts <- vector("list", params$cycles)
  for (cyc in seq_len(params$cycles)) {
    layerSet <- runLayerBinding(
      layerList  = layerSet,
      factorSet  = factorSet,
      iterations = params$mods_per_cycle,
      verbose    = FALSE
    )
    ts[[cyc]] <- cbind(data.frame(cycle = cyc, model = label),
                       summarise_layers(layerSet, params$chromosome))
    if (cyc %% 5 == 0) {
      s <- ts[[cyc]]
      cat(sprintf("    Cycle %2d/%d  active=%.3f  repressed=%.3f  protected=%.3f  bivalent=%.4f\n",
                  cyc, params$cycles, s$frac_active, s$frac_repressed,
                  s$frac_protected, s$frac_bivalent))
    }
  }
  list(layerSet = layerSet, ts = do.call(rbind, ts))
}

# =============================================================================
cat(strrep("=", 65), "\n", sep = "")
cat("GenomicLayers: Histone Mark Cross-Talk Simulation\n")
cat(strrep("=", 65), "\n\n", sep = "")

# =============================================================================
# SECTION 1: Build genome layer set
# =============================================================================
cat("[1/7] Building LayerSet — S. cerevisiae sacCer3, ", params$chromosome, "\n", sep = "")

genome   <- BSgenome.Scerevisiae.UCSC.sacCer3
chr_len  <- length(genome[[params$chromosome]])
cat(sprintf("      Chromosome length : %d bp\n\n", chr_len))

make_layerset <- function() {
  createLayerSet.BSgenome(
    genome   = genome,
    seqnames = params$chromosome,
    n.layers = params$n_layers,
    verbose  = FALSE
  )
}

# =============================================================================
# SECTION 2: Define binding factors
# =============================================================================
cat("[2/7] Defining binding factors\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# CROSS-TALK FACTOR SET (Rules 1–5 active)
# ─────────────────────────────────────────────────────────────────────────────

# Rule 1 + Rule 3 combined:
# Writer (Set1/COMPASS, H3K4me3) — deposits L1=1 at TATA-box motifs,
# but ONLY where L2 = 0 (blocked by existing H3K27me3)
# profile.layers/profile.marks implement the conditional binding rule.
bf_writer <- createBindingFactor.DNA_motif(
  name           = "Writer_H3K4me3",
  patternString  = "TATAAA",
  profile.layers = "LAYER.2",      # Rule 1: check L2 state before binding
  profile.marks  = 0L,             # only bind where L2 = 0 (not repressed)
  mod.layers     = "LAYER.1",
  mod.marks      = 1L,
  stateWidth     = 500L
)

# Rule 2:
# Repressive writer (PRC2, H3K27me3) — deposits L2=1 at CpG-rich seeds,
# but ONLY where L3 = 0 (blocked by H3K36me3 — the in-cis protection rule)
# Mechanistic basis: H3K36me3 inhibits PRC2 catalytic activity in cis
# (Nature Struct Mol Biol 2025; Mol Cell 2010)
bf_repwriter <- createBindingFactor.DNA_motif(
  name           = "RepWriter_PRC2",
  patternString  = "GCGCGC",
  profile.layers = "LAYER.3",      # Rule 2: check L3 state before binding
  profile.marks  = 0L,             # only bind where L3 = 0 (no H3K36me3)
  mod.layers     = "LAYER.2",
  mod.marks      = 1L,
  stateWidth     = 1000L
)

# PRC2 also removes the active mark at CpG seeds (eraser activity)
# Same conditional on L3 = 0
bf_eraser <- createBindingFactor.DNA_motif(
  name           = "Eraser_PRC2",
  patternString  = "GCGCGC",
  profile.layers = "LAYER.3",
  profile.marks  = 0L,
  mod.layers     = "LAYER.1",
  mod.marks      = 0L,
  stateWidth     = 1000L
)

# Rule 3:
# EED reader-writer feedback — PRC2 reads H3K27me3 via EED WD40 domain
# and deposits more H3K27me3 (allosteric activation ~7x, Margueron 2009)
# Implemented as layer_region factor: reads L2=1, writes L2=1 in wider window
# This creates the positive feedback loop that stabilises repressed domains
bf_eed <- createBindingFactor.layer_region(
  name           = "EED_ReaderWriter",
  profile.layers = "LAYER.2",      # reads existing H3K27me3
  profile.marks  = 1L,
  mod.layers     = "LAYER.2",      # deposits more H3K27me3
  mod.marks      = 1L,
  stateWidth     = 2000L           # broader spreading window
)

# Rule 4:
# Set2 (H3K36me3 writer) — transcription-coupled mark deposition
# Active mark → transcription → Set2 recruited to elongating RNAPII
# → H3K36me3 deposited in gene bodies → reinforces protection (Rule 2)
bf_set2 <- createBindingFactor.layer_region(
  name           = "Set2_H3K36me3",
  profile.layers = "LAYER.1",      # reads active mark (transcription proxy)
  profile.marks  = 1L,
  mod.layers     = "LAYER.3",      # deposits H3K36me3-like protection mark
  mod.marks      = 1L,
  stateWidth     = 800L
)

# Rule 5:
# KDM6A/B (UTX/JMJD3) — H3K27me3 demethylase recruited to active loci
# Active mark → KDM6 recruitment → erases H3K27me3
# This is the demethylase that resolves bivalent domains upon activation
bf_kdm6 <- createBindingFactor.layer_region(
  name           = "KDM6_Demethylase",
  profile.layers = "LAYER.1",      # reads active mark
  profile.marks  = 1L,
  mod.layers     = "LAYER.2",      # erases H3K27me3
  mod.marks      = 0L,
  stateWidth     = 400L
)

factorSet_crosstalk <- list(
  Writer    = bf_writer,
  RepWriter = bf_repwriter,
  Eraser    = bf_eraser,
  EED       = bf_eed,
  Set2      = bf_set2,
  KDM6      = bf_kdm6
)

# ─────────────────────────────────────────────────────────────────────────────
# NAIVE FACTOR SET (no interaction rules — marks deposited independently)
# ─────────────────────────────────────────────────────────────────────────────
bf_writer_naive <- createBindingFactor.DNA_motif(
  name = "Writer_naive", patternString = "TATAAA",
  mod.layers = "LAYER.1", mod.marks = 1L, stateWidth = 500L
)
bf_repwriter_naive <- createBindingFactor.DNA_motif(
  name = "RepWriter_naive", patternString = "GCGCGC",
  mod.layers = "LAYER.2", mod.marks = 1L, stateWidth = 1000L
)
bf_eraser_naive <- createBindingFactor.DNA_motif(
  name = "Eraser_naive", patternString = "GCGCGC",
  mod.layers = "LAYER.1", mod.marks = 0L, stateWidth = 1000L
)

factorSet_naive <- list(
  Writer    = bf_writer_naive,
  RepWriter = bf_repwriter_naive,
  Eraser    = bf_eraser_naive
)

# Print factor summary
cat("  Cross-talk model factors:\n")
cat("  ┌─────────────────────┬───────────────┬──────────┬───────┬──────────┐\n")
cat("  │ Factor              │ Type          │ Reads    │ Writes│ Window   │\n")
cat("  ├─────────────────────┼───────────────┼──────────┼───────┼──────────┤\n")
cat("  │ Writer_H3K4me3      │ DNA TATAAA    │ L2=0 req │ L1→1  │  500 bp  │\n")
cat("  │ RepWriter_PRC2      │ DNA GCGCGC    │ L3=0 req │ L2→1  │ 1000 bp  │\n")
cat("  │ Eraser_PRC2         │ DNA GCGCGC    │ L3=0 req │ L1→0  │ 1000 bp  │\n")
cat("  │ EED_ReaderWriter    │ layer_region  │ L2=1     │ L2→1  │ 2000 bp  │\n")
cat("  │ Set2_H3K36me3       │ layer_region  │ L1=1     │ L3→1  │  800 bp  │\n")
cat("  │ KDM6_Demethylase    │ layer_region  │ L1=1     │ L2→0  │  400 bp  │\n")
cat("  └─────────────────────┴───────────────┴──────────┴───────┴──────────┘\n\n")

# =============================================================================
# SECTION 3: Run simulation(s)
# =============================================================================
cat("[3/7] Running simulation(s)\n\n")

ts_all <- NULL

if (run_crosstalk) {
  set.seed(params$rng_seed)
  result_ct <- run_simulation(make_layerset(), factorSet_crosstalk, params, "cross-talk")
  layerSet_ct <- result_ct$layerSet
  ts_all      <- result_ct$ts
}

if (run_naive) {
  set.seed(params$rng_seed)
  result_nv <- run_simulation(make_layerset(), factorSet_naive, params, "naive")
  layerSet_nv <- result_nv$layerSet
  ts_all      <- rbind(ts_all, result_nv$ts)
}

cat("\n")

# =============================================================================
# SECTION 4: Cross-talk analysis
# =============================================================================
cat("[4/7] Cross-talk analysis\n\n")

compute_crosstalk_stats <- function(layerSet, chr, bin_size, label) {
  lyr   <- layerSet[[chr]]$layers
  n_bp  <- nrow(lyr)
  n_bin <- floor(n_bp / bin_size)

  bins <- data.frame(
    bin_mid    = seq(bin_size / 2, by = bin_size, length.out = n_bin),
    l1         = sapply(seq_len(n_bin), function(i) mean(lyr[((i-1)*bin_size+1):(i*bin_size), 1])),
    l2         = sapply(seq_len(n_bin), function(i) mean(lyr[((i-1)*bin_size+1):(i*bin_size), 2])),
    l3         = sapply(seq_len(n_bin), function(i) mean(lyr[((i-1)*bin_size+1):(i*bin_size), 3])),
    model      = label
  )

  r_l1_l2 <- cor(bins$l1, bins$l2)
  r_l3_l2 <- cor(bins$l3, bins$l2)  # H3K36me3 should anti-correlate with H3K27me3

  # Mutual exclusivity score: fraction of bins where at most one mark > 0.3
  me_score <- mean(!(bins$l1 > 0.3 & bins$l2 > 0.3))

  # Bivalent bins (both marks > 0.3) — the poised/bivalent prediction
  bivalent_bins <- sum(bins$l1 > 0.3 & bins$l2 > 0.3)

  cat(sprintf("  [%s model]\n", label))
  cat(sprintf("    Pearson r(L1, L2)     = %+.3f  (mutual exclusivity; target: -1)\n", r_l1_l2))
  cat(sprintf("    Pearson r(L3, L2)     = %+.3f  (H3K36me3 protection; target: -1)\n", r_l3_l2))
  cat(sprintf("    Mutual exclusivity    = %.3f  (fraction of bins with ≤1 active mark)\n", me_score))
  cat(sprintf("    Bivalent bins (>0.3)  = %d / %d  (%.1f%%)\n",
              bivalent_bins, n_bin, 100 * bivalent_bins / n_bin))

  list(bins = bins, r_l1_l2 = r_l1_l2, r_l3_l2 = r_l3_l2,
       me_score = me_score, bivalent_bins = bivalent_bins)
}

if (run_crosstalk) {
  ct_stats <- compute_crosstalk_stats(layerSet_ct, params$chromosome,
                                      params$bin_size, "cross-talk")
}
if (run_naive) {
  nv_stats <- compute_crosstalk_stats(layerSet_nv, params$chromosome,
                                      params$bin_size, "naive")
}
cat("\n")

# =============================================================================
# SECTION 5: Gene-level prediction
# =============================================================================
cat("[5/7] Gene-level promoter state prediction\n\n")

classify_promoters <- function(layerSet, chr, txdb, params, label) {
  lyr   <- layerSet[[chr]]$layers
  n_bp  <- nrow(lyr)

  # Get genes on this chromosome
  all_genes  <- suppressWarnings(genes(txdb))
  chr_genes  <- all_genes[seqnames(all_genes) == chr]

  if (length(chr_genes) == 0) {
    cat(sprintf("  No genes found on %s in TxDb\n", chr))
    return(NULL)
  }

  # Promoter = upstream window from TSS
  prom <- suppressWarnings(
    promoters(chr_genes, upstream = params$promoter_up, downstream = 0L)
  )
  prom <- trim(prom)   # clip to chromosome boundaries

  # Extract mean layer signals per promoter
  gene_df <- data.frame(
    gene_id  = names(prom),
    tss      = ifelse(strand(chr_genes) == "+",
                      start(chr_genes), end(chr_genes)),
    strand   = as.character(strand(chr_genes)),
    prom_l1  = sapply(seq_along(prom), function(i) {
      idx <- max(1, start(prom[i])):min(n_bp, end(prom[i]))
      if (length(idx) == 0) return(NA_real_)
      mean(lyr[idx, 1])
    }),
    prom_l2  = sapply(seq_along(prom), function(i) {
      idx <- max(1, start(prom[i])):min(n_bp, end(prom[i]))
      if (length(idx) == 0) return(NA_real_)
      mean(lyr[idx, 2])
    }),
    prom_l3  = sapply(seq_along(prom), function(i) {
      idx <- max(1, start(prom[i])):min(n_bp, end(prom[i]))
      if (length(idx) == 0) return(NA_real_)
      mean(lyr[idx, 3])
    })
  )
  gene_df <- gene_df[complete.cases(gene_df), ]

  # Classify each promoter state
  gene_df$state <- with(gene_df, ifelse(
    prom_l1 > 0.3 & prom_l2 > 0.3, "Bivalent",
    ifelse(prom_l1 > 0.3 & prom_l2 <= 0.3, "Active",
    ifelse(prom_l1 <= 0.3 & prom_l2 > 0.3, "Repressed",
    "Unmarked"))
  ))

  # Telomere-proximal test:
  # In S. cerevisiae SIR silencing represses subtelomeric genes.
  # Here chrI is ~230 kb; subtelomeric = within 20 kb of either end.
  gene_df$subtelomeric <- gene_df$tss < 20000 | gene_df$tss > (chr_len - 20000)

  n_sub_rep  <- sum(gene_df$subtelomeric & gene_df$state == "Repressed",  na.rm = TRUE)
  n_sub_tot  <- sum(gene_df$subtelomeric, na.rm = TRUE)
  n_int_act  <- sum(!gene_df$subtelomeric & gene_df$state == "Active",    na.rm = TRUE)
  n_int_tot  <- sum(!gene_df$subtelomeric, na.rm = TRUE)

  gene_df$model <- label

  cat(sprintf("  [%s model]  %d genes on %s\n", label, nrow(gene_df), chr))
  cat(sprintf("    Active     : %d  (%.1f%%)\n",
              sum(gene_df$state == "Active"),
              100 * mean(gene_df$state == "Active")))
  cat(sprintf("    Repressed  : %d  (%.1f%%)\n",
              sum(gene_df$state == "Repressed"),
              100 * mean(gene_df$state == "Repressed")))
  cat(sprintf("    Bivalent   : %d  (%.1f%%)\n",
              sum(gene_df$state == "Bivalent"),
              100 * mean(gene_df$state == "Bivalent")))
  cat(sprintf("    Unmarked   : %d  (%.1f%%)\n",
              sum(gene_df$state == "Unmarked"),
              100 * mean(gene_df$state == "Unmarked")))
  cat(sprintf("\n    Telomere-proximity test (SIR silencing proxy)\n"))
  cat(sprintf("      Subtelomeric genes (%d total): %.1f%% predicted repressed\n",
              n_sub_tot, 100 * n_sub_rep / max(1, n_sub_tot)))
  cat(sprintf("      Internal genes (%d total):     %.1f%% predicted active\n",
              n_int_tot, 100 * n_int_act / max(1, n_int_tot)))

  gene_df
}

txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

gene_pred_ct <- NULL
gene_pred_nv <- NULL

if (run_crosstalk) {
  gene_pred_ct <- classify_promoters(layerSet_ct, params$chromosome, txdb, params, "cross-talk")
  cat("\n")
}
if (run_naive) {
  gene_pred_nv <- classify_promoters(layerSet_nv, params$chromosome, txdb, params, "naive")
  cat("\n")
}

# =============================================================================
# SECTION 6: Bivalent domain characterisation
# =============================================================================
cat("[6/7] Bivalent domain characterisation\n\n")

# Bivalent domains are predicted at loci where the cross-talk rules produce
# co-occupancy of L1 and L2 — a stable poised state analogous to ESC bivalency.
# These are hypothesised to occur at CpG-rich promoters that also have a TATA-box
# nearby, but lack H3K36me3 protection (not actively transcribed).

if (run_crosstalk && !is.null(gene_pred_ct)) {
  biv_genes <- gene_pred_ct[gene_pred_ct$state == "Bivalent", ]
  cat(sprintf("  Cross-talk model: %d bivalent gene promoters\n", nrow(biv_genes)))

  if (nrow(biv_genes) > 0) {
    # Extract sequence around bivalent TSS and check GC content
    biv_seqs <- tryCatch({
      ranges_biv <- GRanges(
        seqnames = params$chromosome,
        ranges   = IRanges(
          start = pmax(1, biv_genes$tss - 500),
          end   = pmin(chr_len, biv_genes$tss + 500)
        )
      )
      as.character(getSeq(genome, ranges_biv))
    }, error = function(e) NULL)

    if (!is.null(biv_seqs)) {
      gc_content <- sapply(biv_seqs, function(s) {
        n <- nchar(s)
        if (n == 0) return(NA)
        (nchar(gsub("[^GC]", "", s))) / n
      })
      biv_genes$gc_500bp <- gc_content

      cat(sprintf("    Mean GC content (±500 bp from TSS): %.1f%%\n",
                  100 * mean(gc_content, na.rm = TRUE)))
      cat(sprintf("    (active genes for comparison: %.1f%%)\n",
                  100 * mean(
                    gene_pred_ct$gc_500bp[gene_pred_ct$state == "Active"],
                    na.rm = TRUE)))
    }

    cat("\n  Prediction: bivalent loci should have higher GC (CpG-like) content")
    cat("\n  than active-only loci — this would recapitulate the known CpG island")
    cat("\n  enrichment at bivalent domains in ESCs (Bernstein et al. 2006, Cell).\n\n")
  }
}

# =============================================================================
# SECTION 7: Visualisations
# =============================================================================
cat("[7/7] Generating visualisations\n\n")

COLOURS <- c(
  "Active (H3K4me3-like)"     = "#2166ac",
  "Repressive (H3K27me3-like)"= "#d6604d",
  "Bivalent"                  = "#4dac26",
  "Protected (H3K36me3-like)" = "#762a83",
  "Unmarked"                  = "#bdbdbd",
  "cross-talk"                = "#1b7837",
  "naive"                     = "#d73027"
)

# ── Plot 1: Timeseries of layer occupancy ─────────────────────────────────
if (!is.null(ts_all)) {
  ts_long <- reshape2::melt(ts_all,
    id.vars      = c("cycle", "model"),
    measure.vars = c("frac_active", "frac_repressed", "frac_protected", "frac_bivalent"),
    variable.name = "mark", value.name = "fraction")
  ts_long$mark <- factor(ts_long$mark,
    levels = c("frac_active", "frac_repressed", "frac_protected", "frac_bivalent"),
    labels = c("Active (H3K4me3-like)", "Repressive (H3K27me3-like)",
               "Protected (H3K36me3-like)", "Bivalent"))
  ts_long$model <- factor(ts_long$model, levels = c("cross-talk", "naive"))

  p1 <- ggplot(ts_long, aes(x = cycle, y = fraction, colour = mark,
                             linetype = model)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = COLOURS) +
    scale_linetype_manual(values = c("cross-talk" = "solid", "naive" = "dashed")) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    labs(title    = "Chromatin Mark Occupancy Over Simulation Cycles",
         subtitle = "Solid = cross-talk model (Rules 1–5)  |  Dashed = naive model",
         x = "Cycle", y = "Fraction of nucleotides",
         colour = "Histone mark", linetype = "Model") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave("layer_occupancy_timeseries.png", p1, width = 9, height = 5, dpi = 300)
  cat("  Saved: layer_occupancy_timeseries.png\n")
}

# ── Plot 2: Cross-talk scatter — L1 vs L2 per bin ─────────────────────────
# Expected: anti-correlated scatter (L-shaped cloud), more so in cross-talk model
all_bins <- NULL
if (run_crosstalk) all_bins <- rbind(all_bins, ct_stats$bins)
if (run_naive)     all_bins <- rbind(all_bins, nv_stats$bins)

if (!is.null(all_bins)) {
  # Label quadrants
  all_bins$domain_class <- with(all_bins,
    ifelse(l1 > 0.3 & l2 > 0.3, "Bivalent",
    ifelse(l1 > 0.3 & l2 <= 0.3, "Active",
    ifelse(l1 <= 0.3 & l2 > 0.3, "Repressed", "Unmarked"))))

  r_labels <- all_bins |>
    (\(d) tapply(seq_len(nrow(d)), d$model, function(i) {
      r <- cor(d$l1[i], d$l2[i])
      data.frame(model = d$model[i[1]], label = sprintf("r = %.3f", r),
                 x = 0.7, y = 0.85)
    }))() |> do.call(what = rbind)

  p2 <- ggplot(all_bins, aes(x = l1, y = l2, colour = domain_class)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text(data = r_labels, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_colour_manual(values = c(
      "Active"    = "#2166ac", "Repressed" = "#d6604d",
      "Bivalent"  = "#4dac26", "Unmarked"  = "#999999")) +
    facet_wrap(~ model, labeller = labeller(model = c(
      "cross-talk" = "Cross-talk model (Rules 1–5)",
      "naive"      = "Naive model (no interaction rules)"))) +
    labs(title    = "Active vs. Repressive Mark Cross-Talk Per Genomic Bin",
         subtitle = sprintf("1 kb bins on %s — anti-correlation = mutual exclusivity",
                            params$chromosome),
         x = "Active mark occupancy (L1)", y = "Repressive mark occupancy (L2)",
         colour = "Domain class") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave("crosstalk_scatter.png", p2, width = 9, height = 5, dpi = 300)
  cat("  Saved: crosstalk_scatter.png\n")
}

# ── Plot 3: Chromosome-wide mark profile ──────────────────────────────────
if (run_crosstalk) {
  bins_ct_long <- reshape2::melt(ct_stats$bins,
    id.vars = c("bin_mid", "model"), measure.vars = c("l1", "l2", "l3"),
    variable.name = "layer", value.name = "occupancy")
  bins_ct_long$layer <- factor(bins_ct_long$layer,
    levels = c("l1", "l2", "l3"),
    labels = c("Active (H3K4me3-like)", "Repressive (H3K27me3-like)",
               "Protected (H3K36me3-like)"))

  p3 <- ggplot(bins_ct_long, aes(x = bin_mid / 1000, y = occupancy, fill = layer)) +
    geom_area(alpha = 0.75) +
    scale_fill_manual(values = c(
      "Active (H3K4me3-like)"     = "#2166ac",
      "Repressive (H3K27me3-like)"= "#d6604d",
      "Protected (H3K36me3-like)" = "#762a83")) +
    facet_wrap(~ layer, ncol = 1) +
    labs(title    = sprintf("Chromatin Mark Distribution — %s (cross-talk model)",
                            params$chromosome),
         subtitle = "1 kb resolution | Three-layer cross-talk model",
         x = "Position (kb)", y = "Occupancy", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))

  ggsave("chromosome_mark_profile.png", p3, width = 10, height = 7, dpi = 300)
  cat("  Saved: chromosome_mark_profile.png\n")
}

# ── Plot 4: Gene promoter state classification ────────────────────────────
gene_pred_all <- NULL
if (!is.null(gene_pred_ct)) gene_pred_all <- rbind(gene_pred_all, gene_pred_ct)
if (!is.null(gene_pred_nv)) gene_pred_all <- rbind(gene_pred_all, gene_pred_nv)

if (!is.null(gene_pred_all)) {
  state_counts <- aggregate(gene_id ~ state + model, data = gene_pred_all, FUN = length)
  names(state_counts)[3] <- "count"
  total_by_model <- aggregate(count ~ model, data = state_counts, FUN = sum)
  state_counts <- merge(state_counts, total_by_model, by = "model", suffixes = c("", "_total"))
  state_counts$pct <- 100 * state_counts$count / state_counts$count_total
  state_counts$state <- factor(state_counts$state,
    levels = c("Active", "Bivalent", "Repressed", "Unmarked"))

  p4 <- ggplot(state_counts, aes(x = state, y = pct, fill = state)) +
    geom_col(width = 0.6, colour = "white") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c(
      "Active"    = "#2166ac", "Repressed" = "#d6604d",
      "Bivalent"  = "#4dac26", "Unmarked"  = "#bdbdbd")) +
    facet_wrap(~ model, labeller = labeller(model = c(
      "cross-talk" = "Cross-talk model", "naive" = "Naive model"))) +
    labs(title    = "Predicted Gene Promoter States",
         subtitle = sprintf("%s | Promoter = %d bp upstream of TSS",
                            params$chromosome, params$promoter_up),
         x = NULL, y = "% of genes", fill = "State") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  ggsave("gene_promoter_states.png", p4, width = 8, height = 5, dpi = 300)
  cat("  Saved: gene_promoter_states.png\n")
}

# ── Plot 5: Telomere-proximity enrichment (SIR silencing test) ────────────
if (!is.null(gene_pred_all)) {
  # Test: subtelomeric vs. internal gene repression rate
  gene_pred_all$region <- ifelse(gene_pred_all$subtelomeric,
                                  "Subtelomeric\n(<20 kb from end)", "Internal")

  telo_df <- aggregate(gene_id ~ state + region + model, data = gene_pred_all, FUN = length)
  names(telo_df)[4] <- "count"
  total_r <- aggregate(count ~ region + model, data = telo_df, FUN = sum)
  telo_df  <- merge(telo_df, total_r, by = c("region", "model"), suffixes = c("", "_total"))
  telo_df$pct   <- 100 * telo_df$count / telo_df$count_total
  telo_df$state <- factor(telo_df$state,
    levels = c("Active", "Bivalent", "Repressed", "Unmarked"))

  p5 <- ggplot(telo_df, aes(x = state, y = pct, fill = state)) +
    geom_col(width = 0.6, colour = "white") +
    scale_fill_manual(values = c(
      "Active"    = "#2166ac", "Repressed" = "#d6604d",
      "Bivalent"  = "#4dac26", "Unmarked"  = "#bdbdbd")) +
    facet_grid(model ~ region, labeller = labeller(model = c(
      "cross-talk" = "Cross-talk", "naive" = "Naive"))) +
    labs(title    = "Promoter State by Chromosomal Region",
         subtitle = "Test: are subtelomeric genes (SIR-silenced in yeast) predicted repressed?",
         x = NULL, y = "% of genes in region", fill = "State") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))

  ggsave("telomere_proximity_test.png", p5, width = 9, height = 6, dpi = 300)
  cat("  Saved: telomere_proximity_test.png\n")
}

# ── Plot 6: Final state breakdown bar chart ───────────────────────────────
if (run_crosstalk) {
  fl <- layerSet_ct[[params$chromosome]]$layers
  state_df <- data.frame(
    state = c("Active only", "Repressive only", "Bivalent", "Protected only", "Unmarked"),
    pct   = 100 * c(
      mean( fl[,1] & !fl[,2] & !fl[,3]),
      mean(!fl[,1] &  fl[,2] & !fl[,3]),
      mean( fl[,1] &  fl[,2]),
      mean(!fl[,1] & !fl[,2] &  fl[,3]),
      mean(!fl[,1] & !fl[,2] & !fl[,3])
    )
  )
  state_df$state <- factor(state_df$state,
    levels = c("Active only", "Bivalent", "Protected only", "Repressive only", "Unmarked"))

  p6 <- ggplot(state_df, aes(x = state, y = pct, fill = state)) +
    geom_col(width = 0.6, colour = "white") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c(
      "Active only"     = "#2166ac",
      "Repressive only" = "#d6604d",
      "Bivalent"        = "#4dac26",
      "Protected only"  = "#762a83",
      "Unmarked"        = "#bdbdbd")) +
    labs(title    = "Final Chromatin State Breakdown (cross-talk model)",
         subtitle = sprintf("%s after %d cycles", params$chromosome, params$cycles),
         x = NULL, y = "% of chromosome") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave("final_state_breakdown.png", p6, width = 7, height = 5, dpi = 300)
  cat("  Saved: final_state_breakdown.png\n")
}

# ── CSVs ──────────────────────────────────────────────────────────────────
if (!is.null(ts_all))         write.csv(ts_all,         "simulation_timeseries.csv",    row.names = FALSE)
if (!is.null(all_bins))       write.csv(all_bins,       "binned_mark_profiles.csv",     row.names = FALSE)
if (!is.null(gene_pred_all))  write.csv(gene_pred_all,  "gene_promoter_predictions.csv", row.names = FALSE)
cat("  Saved: simulation_timeseries.csv\n")
cat("  Saved: binned_mark_profiles.csv\n")
cat("  Saved: gene_promoter_predictions.csv\n")

# =============================================================================
# OPTIONAL: Parameter sweep
# =============================================================================
if ("--sweep" %in% args) {
  cat("\n[SWEEP] Parameter sweep: writer stateWidth vs EED spreading window\n\n")

  writer_widths <- c(200L, 500L, 1000L)
  eed_widths    <- c(500L, 1000L, 2000L, 4000L)
  sweep_results <- NULL

  for (ww in writer_widths) {
    for (ew in eed_widths) {
      set.seed(params$rng_seed)
      ls_s <- make_layerset()

      bf_w_s <- createBindingFactor.DNA_motif("Writer", "TATAAA",
        profile.layers = "LAYER.2", profile.marks = 0L,
        mod.layers = "LAYER.1", mod.marks = 1L, stateWidth = ww)
      bf_r_s <- createBindingFactor.DNA_motif("RepWriter", "GCGCGC",
        profile.layers = "LAYER.3", profile.marks = 0L,
        mod.layers = "LAYER.2", mod.marks = 1L, stateWidth = ww * 2L)
      bf_eed_s <- createBindingFactor.layer_region("EED", profile.layers = "LAYER.2",
        profile.marks = 1L, mod.layers = "LAYER.2", mod.marks = 1L, stateWidth = ew)
      bf_s2_s <- createBindingFactor.layer_region("Set2", profile.layers = "LAYER.1",
        profile.marks = 1L, mod.layers = "LAYER.3", mod.marks = 1L, stateWidth = ww)

      fs_s <- list(W = bf_w_s, R = bf_r_s, EED = bf_eed_s, S2 = bf_s2_s)

      for (cyc in seq_len(params$cycles)) {
        ls_s <- runLayerBinding(ls_s, fs_s, iterations = params$mods_per_cycle, verbose = FALSE)
      }

      fl_s <- ls_s[[params$chromosome]]$layers
      r_val <- cor(
        sapply(seq_len(floor(nrow(fl_s)/1000)), function(i) mean(fl_s[((i-1)*1000+1):(i*1000), 1])),
        sapply(seq_len(floor(nrow(fl_s)/1000)), function(i) mean(fl_s[((i-1)*1000+1):(i*1000), 2]))
      )

      sweep_results <- rbind(sweep_results, data.frame(
        writer_width  = ww,
        eed_width     = ew,
        frac_active   = round(mean(fl_s[, 1]), 4),
        frac_repressed = round(mean(fl_s[, 2]), 4),
        frac_bivalent  = round(mean(fl_s[, 1] & fl_s[, 2]), 4),
        r_l1_l2        = round(r_val, 4)
      ))
      cat(sprintf("  writer=%4d  EED=%4d  active=%.3f  repressed=%.3f  bivalent=%.4f  r=%.3f\n",
                  ww, ew, tail(sweep_results$frac_active, 1),
                  tail(sweep_results$frac_repressed, 1),
                  tail(sweep_results$frac_bivalent, 1),
                  tail(sweep_results$r_l1_l2, 1)))
    }
  }

  write.csv(sweep_results, "parameter_sweep_results.csv", row.names = FALSE)

  p_sw <- ggplot(sweep_results,
                 aes(x = factor(ew), y = factor(writer_width), fill = r_l1_l2)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r_l1_l2)), size = 3.5) +
    scale_fill_gradient2(low = "#d6604d", mid = "#f7f7f7", high = "#2166ac",
                         midpoint = 0, name = "r(L1,L2)\nmutual\nexclusivity") +
    labs(title = "Sweep: L1/L2 Anti-Correlation (mutual exclusivity) by stateWidth parameters",
         subtitle = "Target r = -1; more negative = better mutual exclusivity from interaction rules",
         x = "EED spreading window (bp)", y = "Writer stateWidth (bp)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  ggsave("parameter_sweep_heatmap.png", p_sw, width = 8, height = 5, dpi = 300)
  cat("  Saved: parameter_sweep_heatmap.png\n")
  cat("  Saved: parameter_sweep_results.csv\n")
}

cat("\nSimulation complete.\n")
