# Rules-Based Epigenome Simulation

Sequence-informed writer/eraser simulation of histone methylation spreading,
using the [GenomicLayers](https://github.com/davetgerrard/GenomicLayers) package
(Gerrard, *BMC Bioinformatics* 2025).

## Biological background

Histone methylation patterns are maintained by opposing writer and eraser enzymes
that recognise both DNA sequence motifs and existing epigenetic states. This script
models that process on chromosome I of *S. cerevisiae* using three binding factors:

| Factor | Motif | Layer | Biological analogue |
|--------|-------|-------|---------------------|
| Writer | `TATAAA` (TATA-box) | Layer 1 → 1 | H3K4me3 writer (Set1/COMPASS) |
| Eraser | `GCGCGC` (CpG-rich) | Layer 1 → 0 | PRC2 — removes active mark |
| RepWriter | `GCGCGC` (CpG-rich) | Layer 2 → 1 | PRC2 — deposits H3K27me3 |

## Installation

```r
install.packages("devtools")
BiocManager::install(c("BSgenome", "GenomicRanges", "Biostrings"))
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
devtools::install_github("davetgerrard/GenomicLayers", build_vignettes = TRUE)
```

## Usage

```bash
# Default simulation (20 cycles on chrI)
Rscript rules_based_simulation.R

# Parameter sweep: stateWidth x n_cycles heatmap
Rscript rules_based_simulation.R --sweep
```

## Outputs

| File | Description |
|------|-------------|
| `layer_occupancy_timeseries.png` | Active / repressive / bivalent fraction per cycle |
| `chromosome_mark_profile.png` | Mark distribution along the chromosome (5 kb bins) |
| `final_state_breakdown.png` | Bar chart of final active / repressive / bivalent / unmarked |
| `parameter_sweep_heatmap.png` | Heatmap of repressive fraction across stateWidth × cycles |
| `simulation_timeseries.csv` | Numeric timeseries data |
| `chromosome_mark_profile.csv` | Binned chromosomal profile data |
| `parameter_sweep_results.csv` | Full sweep results |

## Skills demonstrated

- R, GenomicLayers, ggplot2, reshape2
- Epigenomics: histone modifications, Polycomb repression, bivalent chromatin
- Sequence-informed stochastic modelling
- BSgenome / Bioconductor ecosystem
- Parameter sweep and visualisation
