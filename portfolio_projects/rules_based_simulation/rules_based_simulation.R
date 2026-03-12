# =============================================================================
# Rules-Based Writer/Eraser Simulation of Histone Methylation
# Using GenomicLayers (Gerrard, BMC Bioinformatics 2025)
# =============================================================================
# Investigates how well a sequence-informed writer/eraser model explains
# epigenome-wide patterns of successive histone methylation.
#
# GenomicLayers models epigenetic state changes genome-wide by simulating
# the accumulation of changes to binary "layers" by user-specified binding
# factors that recognise both DNA sequence motifs and existing layer states.
#
# This script demonstrates the model on S. cerevisiae (yeast), matching the
# package vignette, then extends to a writer/eraser analysis of H3K27me3-like
# spreading from CpG-rich seed regions.
#
# Installation (run once):
#   install.packages("devtools")
#   BiocManager::install(c("BSgenome", "GenomicRanges", "Biostrings"))
#   BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
#   devtools::install_github("davetgerrard/GenomicLayers", build_vignettes = TRUE)
#
# Usage:
#   Rscript rules_based_simulation.R
#   Rscript rules_based_simulation.R --sweep
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicLayers)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(GenomicRanges)
  library(ggplot2)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
params <- list(
  chromosome    = "chrI",          # yeast chromosome I (230 kb — tractable)
  n_layers      = 3L,              # Layer 1 = Active mark, Layer 2 = Repressive mark, Layer 3 = Seed
  cycles        = 20L,             # simulation cycles
  mods_per_cycle = 500L,           # binding events per cycle
  rng_seed      = 42L
)

set.seed(params$rng_seed)

# -----------------------------------------------------------------------------
# 1. Build the LayerSet from yeast genome
# -----------------------------------------------------------------------------
cat("=============================================================\n")
cat("GenomicLayers: Writer/Eraser Histone Methylation Simulation\n")
cat("=============================================================\n\n")
cat(sprintf("[1/5] Building LayerSet for %s (S. cerevisiae sacCer3)...\n",
            params$chromosome))

genome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Create layer set for a single chromosome
layerSet <- createLayerSet.BSgenome(
  genome   = genome,
  seqnames = params$chromosome,
  n.layers = params$n_layers,
  verbose  = FALSE
)

cat(sprintf("      Chromosome length : %d bp\n",
            length(layerSet[[params$chromosome]]$sequence)))
cat(sprintf("      Layers created    : %d\n\n", params$n_layers))


# -----------------------------------------------------------------------------
# 2. Define binding factors (writer and eraser)
# -----------------------------------------------------------------------------
cat("[2/5] Defining binding factors...\n")

# WRITER — adds active mark (Layer 1) near AT-rich motifs
# Mimics H3K4me3 writer (Set1/COMPASS) recruited to open chromatin
bf_writer <- createBindingFactor.DNA_motif(
  name          = "Writer_H3K4me3",
  patternString = "TATAAA",          # TATA-box like — marks active promoters
  mod.layers    = "LAYER.1",
  mod.marks     = 1L,                # set Layer 1 = 1 (active)
  stateWidth    = 500L               # influence window (bp)
)

# ERASER — removes active mark (Layer 1), dependent on Layer 2 being set
# Mimics PRC2 / H3K27me3 spreading (erases active, writes repressive)
bf_eraser <- createBindingFactor.DNA_motif(
  name          = "Eraser_PRC2",
  patternString = "GCGCGC",          # CpG-rich seed — Polycomb recruitment
  mod.layers    = "LAYER.1",
  mod.marks     = 0L,                # set Layer 1 = 0 (removes active mark)
  stateWidth    = 1000L              # broader spreading window
)

# REPRESSIVE WRITER — deposits repressive mark (Layer 2) at same CpG sites
bf_repwriter <- createBindingFactor.DNA_motif(
  name          = "RepWriter_H3K27me3",
  patternString = "GCGCGC",
  mod.layers    = "LAYER.2",
  mod.marks     = 1L,                # set Layer 2 = 1 (repressive)
  stateWidth    = 1000L
)

factorSet <- list(
  Writer    = bf_writer,
  Eraser    = bf_eraser,
  RepWriter = bf_repwriter
)

for (nm in names(factorSet)) {
  cat(sprintf("      %-20s pattern=%-10s  layer=%s  mark=%d  window=%d bp\n",
              nm,
              factorSet[[nm]]$patternString,
              factorSet[[nm]]$mod.layers,
              factorSet[[nm]]$mod.marks,
              factorSet[[nm]]$stateWidth))
}
cat("\n")


# -----------------------------------------------------------------------------
# 3. Run the simulation
# -----------------------------------------------------------------------------
cat(sprintf("[3/5] Running simulation: %d cycles x %d modifications...\n",
            params$cycles, params$mods_per_cycle))

# Record layer state at each cycle for timeseries analysis
timeseries <- vector("list", params$cycles)

for (cycle in seq_len(params$cycles)) {
  layerSet <- runLayerBinding(
    layerList  = layerSet,
    factorSet  = factorSet,
    iterations = params$mods_per_cycle,
    verbose    = FALSE
  )

  # Summarise layer occupancy for this cycle
  chr_layers <- layerSet[[params$chromosome]]$layers
  timeseries[[cycle]] <- data.frame(
    cycle         = cycle,
    frac_active   = mean(chr_layers[, 1]),   # Layer 1: active mark
    frac_repressed = mean(chr_layers[, 2]),  # Layer 2: repressive mark
    frac_both     = mean(chr_layers[, 1] & chr_layers[, 2])  # bivalent
  )

  if (cycle %% 5 == 0) {
    cat(sprintf("      Cycle %2d/%d — active: %.3f  repressive: %.3f  bivalent: %.3f\n",
                cycle, params$cycles,
                timeseries[[cycle]]$frac_active,
                timeseries[[cycle]]$frac_repressed,
                timeseries[[cycle]]$frac_both))
  }
}

ts_df <- do.call(rbind, timeseries)
cat("\n")


# -----------------------------------------------------------------------------
# 4. Analyse final state
# -----------------------------------------------------------------------------
cat("[4/5] Analysing final chromatin state...\n\n")

final_layers <- layerSet[[params$chromosome]]$layers
n_bp         <- nrow(final_layers)

frac_active    <- mean(final_layers[, 1])
frac_repressed <- mean(final_layers[, 2])
frac_bivalent  <- mean(final_layers[, 1] & final_layers[, 2])
frac_unmarked  <- mean(!final_layers[, 1] & !final_layers[, 2])

cat(strrep("=", 55), "\n", sep = "")
cat("Final Chromatin State Summary\n")
cat(strrep("=", 55), "\n", sep = "")
cat(sprintf("  Chromosome     : %s (%d bp)\n", params$chromosome, n_bp))
cat(sprintf("  Cycles run     : %d\n\n", params$cycles))
cat(sprintf("  Active (L1=1)    : %5.1f%%  [H3K4me3-like]\n",    frac_active    * 100))
cat(sprintf("  Repressive (L2=1): %5.1f%%  [H3K27me3-like]\n",  frac_repressed * 100))
cat(sprintf("  Bivalent (both)  : %5.1f%%  [poised state]\n",   frac_bivalent  * 100))
cat(sprintf("  Unmarked         : %5.1f%%  [silent/default]\n", frac_unmarked  * 100))

# Pattern classification
if (frac_repressed > 0.6) {
  pattern <- "Predominantly REPRESSED (heterochromatin-like)"
} else if (frac_active > 0.6) {
  pattern <- "Predominantly ACTIVE (euchromatin-like)"
} else if (frac_bivalent > 0.1) {
  pattern <- "BIVALENT — poised for activation or silencing"
} else {
  pattern <- "MIXED / INTERMEDIATE state"
}
cat(sprintf("\n  Classification : %s\n", pattern))
cat(strrep("=", 55), "\n\n", sep = "")


# -----------------------------------------------------------------------------
# 5. Visualisations
# -----------------------------------------------------------------------------
cat("[5/5] Generating visualisations...\n")

# --- Plot 1: Timeseries of layer occupancy ---
ts_long <- reshape2::melt(ts_df, id.vars = "cycle",
                          variable.name = "mark", value.name = "fraction")
ts_long$mark <- factor(ts_long$mark,
                       levels = c("frac_active", "frac_repressed", "frac_both"),
                       labels = c("Active (H3K4me3-like)",
                                  "Repressive (H3K27me3-like)",
                                  "Bivalent"))

p1 <- ggplot(ts_long, aes(x = cycle, y = fraction, colour = mark)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    "Active (H3K4me3-like)"     = "#2166ac",
    "Repressive (H3K27me3-like)"= "#d6604d",
    "Bivalent"                  = "#4dac26"
  )) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
  labs(
    title    = "Chromatin Layer Occupancy Over Simulation Cycles",
    subtitle = sprintf("GenomicLayers writer/eraser model — %s, S. cerevisiae",
                       params$chromosome),
    x        = "Simulation cycle",
    y        = "Fraction of nucleotides",
    colour   = "Histone mark"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("layer_occupancy_timeseries.png", p1, width = 8, height = 5, dpi = 300)
cat("  Saved: layer_occupancy_timeseries.png\n")


# --- Plot 2: Final state profile along chromosome (binned) ---
bin_size  <- 5000L  # 5 kb bins
n_bins    <- floor(n_bp / bin_size)
bin_df    <- data.frame(
  bin_start  = seq(1, by = bin_size, length.out = n_bins),
  bin_mid    = seq(bin_size / 2, by = bin_size, length.out = n_bins),
  active     = sapply(seq_len(n_bins), function(i) {
    idx <- ((i - 1) * bin_size + 1):(i * bin_size)
    mean(final_layers[idx, 1])
  }),
  repressive = sapply(seq_len(n_bins), function(i) {
    idx <- ((i - 1) * bin_size + 1):(i * bin_size)
    mean(final_layers[idx, 2])
  })
)

bin_long <- reshape2::melt(bin_df, id.vars = c("bin_start", "bin_mid"),
                            variable.name = "mark", value.name = "fraction")
bin_long$mark <- factor(bin_long$mark,
                        levels = c("active", "repressive"),
                        labels = c("Active (H3K4me3-like)", "Repressive (H3K27me3-like)"))

p2 <- ggplot(bin_long, aes(x = bin_mid / 1000, y = fraction, fill = mark)) +
  geom_area(alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c(
    "Active (H3K4me3-like)"      = "#2166ac",
    "Repressive (H3K27me3-like)" = "#d6604d"
  )) +
  facet_wrap(~ mark, ncol = 1) +
  labs(
    title    = sprintf("Chromatin Mark Distribution Along %s", params$chromosome),
    subtitle = sprintf("Binned at %d kb resolution", bin_size / 1000),
    x        = "Chromosomal position (kb)",
    y        = "Fraction marked per bin",
    fill     = "Mark"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave("chromosome_mark_profile.png", p2, width = 10, height = 5, dpi = 300)
cat("  Saved: chromosome_mark_profile.png\n")


# --- Plot 3: Bar chart of final state breakdown ---
state_df <- data.frame(
  state    = c("Active only", "Repressive only", "Bivalent", "Unmarked"),
  fraction = c(
    mean( final_layers[, 1] & !final_layers[, 2]),
    mean(!final_layers[, 1] &  final_layers[, 2]),
    frac_bivalent,
    frac_unmarked
  )
)
state_df$state <- factor(state_df$state,
                         levels = c("Active only", "Bivalent",
                                    "Repressive only", "Unmarked"))

p3 <- ggplot(state_df, aes(x = state, y = fraction * 100, fill = state)) +
  geom_col(width = 0.6, colour = "white") +
  scale_fill_manual(values = c(
    "Active only"     = "#2166ac",
    "Repressive only" = "#d6604d",
    "Bivalent"        = "#4dac26",
    "Unmarked"        = "#bdbdbd"
  )) +
  geom_text(aes(label = sprintf("%.1f%%", fraction * 100)),
            vjust = -0.4, size = 4.5) +
  labs(
    title    = "Final Chromatin State Breakdown",
    subtitle = sprintf("%s after %d simulation cycles", params$chromosome, params$cycles),
    x        = NULL,
    y        = "% of chromosome"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("final_state_breakdown.png", p3, width = 6, height = 5, dpi = 300)
cat("  Saved: final_state_breakdown.png\n")


# --- Save CSVs ---
write.csv(ts_df,    "simulation_timeseries.csv",      row.names = FALSE)
write.csv(bin_df,   "chromosome_mark_profile.csv",    row.names = FALSE)
write.csv(state_df, "final_state_breakdown.csv",      row.names = FALSE)
cat("  Saved: simulation_timeseries.csv\n")
cat("  Saved: chromosome_mark_profile.csv\n")
cat("  Saved: final_state_breakdown.csv\n")


# -----------------------------------------------------------------------------
# Parameter sweep (optional: Rscript rules_based_simulation.R --sweep)
# -----------------------------------------------------------------------------
if ("--sweep" %in% args) {
  cat("\n[SWEEP] Running writer/eraser parameter sweep...\n")

  state_widths  <- c(200L, 500L, 1000L, 2000L)
  cycle_counts  <- c(5L, 10L, 20L, 40L)
  sweep_results <- vector("list", length(state_widths) * length(cycle_counts))
  idx           <- 1L

  for (sw in state_widths) {
    for (nc in cycle_counts) {
      set.seed(params$rng_seed)
      ls_tmp <- createLayerSet.BSgenome(genome, params$chromosome,
                                        n.layers = params$n_layers, verbose = FALSE)
      bf_w_tmp <- createBindingFactor.DNA_motif(
        "Writer", "TATAAA", mod.layers = "LAYER.1", mod.marks = 1L, stateWidth = sw)
      bf_e_tmp <- createBindingFactor.DNA_motif(
        "Eraser", "GCGCGC", mod.layers = "LAYER.1", mod.marks = 0L, stateWidth = sw * 2L)
      bf_r_tmp <- createBindingFactor.DNA_motif(
        "RepWriter", "GCGCGC", mod.layers = "LAYER.2", mod.marks = 1L, stateWidth = sw * 2L)

      for (cyc in seq_len(nc)) {
        ls_tmp <- runLayerBinding(ls_tmp,
                                  list(w = bf_w_tmp, e = bf_e_tmp, r = bf_r_tmp),
                                  iterations = params$mods_per_cycle, verbose = FALSE)
      }
      fl <- ls_tmp[[params$chromosome]]$layers
      sweep_results[[idx]] <- data.frame(
        state_width    = sw,
        n_cycles       = nc,
        frac_active    = round(mean(fl[, 1]), 4),
        frac_repressed = round(mean(fl[, 2]), 4),
        frac_bivalent  = round(mean(fl[, 1] & fl[, 2]), 4)
      )
      cat(sprintf("  stateWidth=%4d  cycles=%2d  active=%.3f  repressed=%.3f\n",
                  sw, nc, sweep_results[[idx]]$frac_active,
                  sweep_results[[idx]]$frac_repressed))
      idx <- idx + 1L
    }
  }

  sweep_df <- do.call(rbind, sweep_results)
  write.csv(sweep_df, "parameter_sweep_results.csv", row.names = FALSE)

  # Heatmap: repressive fraction
  p_sweep <- ggplot(sweep_df, aes(x = factor(n_cycles), y = factor(state_width),
                                   fill = frac_repressed)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", frac_repressed)), size = 3.5) +
    scale_fill_gradient2(low = "#4393c3", mid = "#f7f7f7", high = "#d6604d",
                         midpoint = 0.3, name = "Fraction\nrepressed") +
    labs(
      title    = "Parameter Sweep: Repressive Mark Fraction at Steady State",
      subtitle = "GenomicLayers writer/eraser model — varying stateWidth and cycles",
      x        = "Number of cycles",
      y        = "Binding factor stateWidth (bp)"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())

  ggsave("parameter_sweep_heatmap.png", p_sweep, width = 8, height = 6, dpi = 300)
  cat("  Saved: parameter_sweep_heatmap.png\n")
  cat("  Saved: parameter_sweep_results.csv\n")
}

cat("\nSimulation complete.\n")
