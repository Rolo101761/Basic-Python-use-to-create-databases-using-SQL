# =============================================================================
# Rules-Based Writer/Eraser Simulation of Histone Methylation
# =============================================================================
# Investigates how well a simple writer/eraser model explains epigenome-wide
# patterns of successive histone methylation (H3K4, H3K9, H3K27, H3K36).
#
# Model:
#   - Nucleosomes sit on a 1D chromatin fibre (periodic boundary conditions)
#   - Each nucleosome carries a methylation state: 0 = me0, 1 = me1, 2 = me2, 3 = me3
#   - WRITER: adds +1 methylation with probability p_write * (neighbour_state / 3)
#             (cooperative spreading — higher neighbour state = more likely)
#   - ERASER: removes -1 methylation with probability p_erase (constitutive)
#   - NOISE:  autonomous +1 addition with probability p_auto (basal activity)
#
# Usage:
#   Rscript rules_based_simulation.R                  # default parameters
#   Rscript rules_based_simulation.R --sweep          # parameter sweep
#
# Or interactively in RStudio: source("rules_based_simulation.R")
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
params <- list(
  n_nucleosomes = 200,     # number of nucleosomes on the fibre
  n_steps       = 50000,   # Monte Carlo steps
  p_write       = 0.08,    # writer probability (neighbour-dependent)
  p_erase       = 0.02,    # eraser probability (constitutive)
  p_auto        = 0.001,   # autonomous (noise) write probability
  record_every  = 500,     # record snapshot every N steps
  seed_pos      = NULL,    # seed position (NULL = centre)
  seed_state    = 3L,      # initial methylation state at seed
  rng_seed      = 42       # reproducibility
)


# -----------------------------------------------------------------------------
# Initialise chromatin fibre
# -----------------------------------------------------------------------------
init_fibre <- function(n, seed_pos = NULL, seed_state = 3L) {
  states <- integer(n)
  seed   <- if (is.null(seed_pos)) as.integer(n / 2) else seed_pos
  states[seed] <- seed_state
  states
}


# -----------------------------------------------------------------------------
# Monte Carlo step
# -----------------------------------------------------------------------------
mc_step <- function(states, p_write, p_erase, p_auto) {
  n       <- length(states)
  indices <- sample.int(n)  # random update order

  for (i in indices) {
    left  <- states[((i - 2) %% n) + 1]
    right <- states[(i %% n) + 1]
    nb_max <- max(left, right)

    # Writer (cooperative spreading)
    if (states[i] < 3L && nb_max > 0L) {
      p <- p_write * (nb_max / 3.0)
      if (runif(1) < p) {
        states[i] <- states[i] + 1L
      }
    }

    # Eraser (constitutive)
    if (states[i] > 0L && runif(1) < p_erase) {
      states[i] <- states[i] - 1L
    }

    # Autonomous write (noise)
    if (states[i] < 3L && runif(1) < p_auto) {
      states[i] <- states[i] + 1L
    }
  }
  states
}


# -----------------------------------------------------------------------------
# Run simulation
# -----------------------------------------------------------------------------
run_simulation <- function(p) {
  set.seed(p$rng_seed)

  states  <- init_fibre(p$n_nucleosomes, p$seed_pos, p$seed_state)
  history <- vector("list", length = ceiling(p$n_steps / p$record_every))
  snap_i  <- 1L

  cat(sprintf("Running %d MC steps on %d nucleosomes...\n", p$n_steps, p$n_nucleosomes))

  for (step in seq_len(p$n_steps)) {
    states <- mc_step(states, p$p_write, p$p_erase, p$p_auto)

    if (step %% p$record_every == 0L) {
      history[[snap_i]] <- list(
        step               = step,
        mean_state         = mean(states),
        frac_methylated    = mean(states > 0),
        me0 = sum(states == 0L),
        me1 = sum(states == 1L),
        me2 = sum(states == 2L),
        me3 = sum(states == 3L)
      )
      snap_i <- snap_i + 1L
    }
  }

  history_df <- do.call(rbind, lapply(history[seq_len(snap_i - 1)], as.data.frame))
  list(final_states = states, history = history_df)
}


# -----------------------------------------------------------------------------
# Analysis helpers
# -----------------------------------------------------------------------------
spreading_length <- function(states) {
  # Half-width at half-maximum around the peak methylation position
  if (all(states == 0)) return(0L)
  half_max   <- max(states) / 2.0
  above_half <- which(states >= half_max)
  if (length(above_half) == 0) return(0L)
  as.integer(max(above_half) - min(above_half) + 1L)
}

classify_pattern <- function(frac_met, spreading, n) {
  if      (frac_met > 0.8)              "HYPERMETHYLATED (heterochromatin-like)"
  else if (frac_met < 0.1)              "HYPOMETHYLATED (euchromatin-like)"
  else if (spreading < n / 5)           "FOCAL domain (sharp boundary)"
  else                                  "BROAD spreading domain"
}

print_final_analysis <- function(states, p) {
  n        <- length(states)
  frac_met <- mean(states > 0)
  mean_st  <- mean(states)
  spread   <- spreading_length(states)
  pattern  <- classify_pattern(frac_met, spread, n)
  counts   <- table(factor(states, levels = 0:3))

  cat("\n", strrep("=", 55), "\n", sep = "")
  cat("Final State Analysis\n")
  cat(strrep("=", 55), "\n", sep = "")
  cat(sprintf("  Nucleosomes  : %d\n",   n))
  cat(sprintf("  MC steps     : %d\n",   p$n_steps))
  cat(sprintf("  p_write      : %.4f\n", p$p_write))
  cat(sprintf("  p_erase      : %.4f\n", p$p_erase))
  cat(sprintf("  p_auto       : %.4f\n", p$p_auto))
  cat("\n  Final nucleosome states:\n")
  for (s in 0:3) {
    cnt <- as.integer(counts[s + 1])
    bar <- strrep("#", round(cnt * 40 / n))
    cat(sprintf("    me%d: %4d (%5.1f%%)  %s\n", s, cnt, cnt / n * 100, bar))
  }
  cat(sprintf("\n  Fraction methylated (me1+me2+me3): %.3f\n", frac_met))
  cat(sprintf("  Mean methylation state           : %.3f\n", mean_st))
  cat(sprintf("  Spreading length (HWHM)          : %d nucleosomes\n", spread))
  cat(sprintf("  Pattern classification           : %s\n", pattern))
  cat(strrep("=", 55), "\n", sep = "")
}


# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------
plot_final_profile <- function(states, output_dir = ".") {
  df <- data.frame(
    position = seq_along(states),
    state    = as.factor(states)
  )
  p <- ggplot(df, aes(x = position, y = 1, fill = state)) +
    geom_tile(height = 0.9) +
    scale_fill_manual(
      values = c("0" = "#d9d9d9", "1" = "#fee08b", "2" = "#f46d43", "3" = "#a50026"),
      labels = c("me0", "me1", "me2", "me3"),
      name   = "Methylation"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      title    = "Final Methylation Profile",
      subtitle = sprintf("p_write=%.3f  p_erase=%.3f  p_auto=%.4f",
                         params$p_write, params$p_erase, params$p_auto),
      x = "Nucleosome position",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank()
    )
  ggsave(file.path(output_dir, "final_methylation_profile.png"),
         p, width = 10, height = 2, dpi = 300)
  cat("  Saved: final_methylation_profile.png\n")
  p
}

plot_timeseries <- function(history, output_dir = ".") {
  # Fraction methylated over time
  p1 <- ggplot(history, aes(x = step, y = frac_methylated)) +
    geom_line(colour = "#a50026", linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    labs(title = "Fraction of Methylated Nucleosomes Over Time",
         x = "MC step", y = "Fraction methylated") +
    theme_bw()
  ggsave(file.path(output_dir, "timeseries_frac_methylated.png"),
         p1, width = 8, height = 4, dpi = 300)

  # State composition stacked area
  long_df <- reshape2::melt(
    history[, c("step", "me0", "me1", "me2", "me3")],
    id.vars = "step", variable.name = "state", value.name = "count"
  )
  p2 <- ggplot(long_df, aes(x = step, y = count, fill = state)) +
    geom_area(position = "fill", alpha = 0.85) +
    scale_fill_manual(
      values = c(me0 = "#d9d9d9", me1 = "#fee08b", me2 = "#f46d43", me3 = "#a50026")
    ) +
    labs(title = "Methylation State Composition Over Time",
         x = "MC step", y = "Proportion", fill = "State") +
    theme_bw()
  ggsave(file.path(output_dir, "timeseries_state_composition.png"),
         p2, width = 8, height = 4, dpi = 300)

  cat("  Saved: timeseries_frac_methylated.png\n")
  cat("  Saved: timeseries_state_composition.png\n")
  list(p1, p2)
}

plot_sweep_heatmap <- function(sweep_results, output_dir = ".") {
  p <- ggplot(sweep_results, aes(x = factor(p_write), y = factor(p_erase),
                                  fill = frac_methylated)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", frac_methylated)), size = 3.5) +
    scale_fill_gradient2(
      low = "#4393c3", mid = "#f7f7f7", high = "#d6604d",
      midpoint = 0.5, name = "Fraction\nmethylated"
    ) +
    labs(
      title    = "Parameter Sweep: Fraction Methylated at Steady State",
      subtitle = sprintf("%d nucleosomes, %d MC steps", params$n_nucleosomes, 20000),
      x        = "p_write",
      y        = "p_erase"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
  ggsave(file.path(output_dir, "parameter_sweep_heatmap.png"),
         p, width = 8, height = 6, dpi = 300)
  cat("  Saved: parameter_sweep_heatmap.png\n")
  p
}


# -----------------------------------------------------------------------------
# Parameter sweep
# -----------------------------------------------------------------------------
run_sweep <- function(n_nucleosomes = 100, n_steps = 20000, output_dir = ".") {
  p_write_vals <- c(0.02, 0.05, 0.08, 0.12, 0.20)
  p_erase_vals <- c(0.01, 0.02, 0.05, 0.10, 0.20)
  total        <- length(p_write_vals) * length(p_erase_vals)

  cat(sprintf("\nParameter sweep: %d combinations (%d nucleosomes, %d steps each)\n",
              total, n_nucleosomes, n_steps))

  results <- vector("list", total)
  idx     <- 1L
  for (pw in p_write_vals) {
    for (pe in p_erase_vals) {
      p_sweep <- modifyList(params, list(
        n_nucleosomes = n_nucleosomes,
        n_steps       = n_steps,
        p_write       = pw,
        p_erase       = pe,
        record_every  = n_steps  # only record final state
      ))
      sim    <- run_simulation(p_sweep)
      states <- sim$final_states
      results[[idx]] <- data.frame(
        p_write          = pw,
        p_erase          = pe,
        frac_methylated  = round(mean(states > 0), 4),
        mean_state       = round(mean(states), 4),
        spreading_length = spreading_length(states)
      )
      cat(sprintf("  p_write=%.2f  p_erase=%.2f  frac_met=%.3f\n",
                  pw, pe, results[[idx]]$frac_methylated))
      idx <- idx + 1L
    }
  }

  sweep_df <- do.call(rbind, results)
  write.csv(sweep_df, file.path(output_dir, "parameter_sweep_results.csv"), row.names = FALSE)
  cat(sprintf("\n  Sweep results saved to: %s\n",
              file.path(output_dir, "parameter_sweep_results.csv")))

  plot_sweep_heatmap(sweep_df, output_dir)
  sweep_df
}


# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
save_timeseries_csv <- function(history, output_dir = ".") {
  path <- file.path(output_dir, "simulation_timeseries.csv")
  write.csv(history, path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", path))
}

save_profile_csv <- function(states, output_dir = ".") {
  df   <- data.frame(position = seq_along(states), methylation_state = states)
  path <- file.path(output_dir, "final_methylation_profile.csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", path))
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

output_dir <- "."
dir.create(output_dir, showWarnings = FALSE)

if ("--sweep" %in% args) {
  run_sweep(n_nucleosomes = 100, n_steps = 20000, output_dir = output_dir)
} else {
  result <- run_simulation(params)
  print_final_analysis(result$final_states, params)

  cat("\nSaving outputs...\n")
  save_timeseries_csv(result$history,   output_dir)
  save_profile_csv(result$final_states, output_dir)
  plot_final_profile(result$final_states, output_dir)
  plot_timeseries(result$history,         output_dir)
  cat("\nDone.\n")
}
