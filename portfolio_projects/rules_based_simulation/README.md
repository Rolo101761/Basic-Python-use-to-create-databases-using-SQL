# Rules-Based Epigenome Simulation with Histone Mark Cross-Talk

Sequence-informed writer/eraser simulation of histone methylation spreading,
using [GenomicLayers](https://github.com/davetgerrard/GenomicLayers)
(Gerrard, *BMC Bioinformatics* 2025) on *S. cerevisiae* chromosome I.

## Biological background

Histone modifications are deposited and removed by writer/eraser enzymes that
do not act independently — they sense the local chromatin environment. This
script models five mechanistically grounded interaction rules:

| Rule | Factor | Mechanism | Biological basis |
|------|--------|-----------|-----------------|
| 1 | Writer (H3K4me3) | Blocked if H3K27me3 already present | Mutual exclusivity of active/repressive marks |
| 2 | PRC2 (H3K27me3) | Blocked if H3K36me3 present (in *cis*) | H3K36me3 inhibits PRC2 catalytic activity; *Nat Struct Mol Biol* 2025 |
| 3 | EED reader-writer | Reads H3K27me3 → spreads H3K27me3 | Allosteric activation of PRC2 by EED ~7×; Margueron *et al.* 2009 |
| 4 | Set2 (H3K36me3) | Active mark → H3K36me3 deposition | Transcription-coupled Set2 recruitment to elongating RNAPII |
| 5 | KDM6A/B | Active mark → erases H3K27me3 | Demethylase resolves bivalency on activation |

Three layers represent the chromatin state at each nucleotide:

| Layer | Mark | Analogue |
|-------|------|---------|
| LAYER.1 | Active | H3K4me3 |
| LAYER.2 | Repressive | H3K27me3 |
| LAYER.3 | Protection | H3K36me3 |

## What makes this novel

Most histone modification simulation models do not:
- Use real DNA sequence to seed writer/eraser binding
- Implement all five interaction rules simultaneously
- Generate genome-scale predictions at single-nucleotide resolution
- Make testable gene-level predictions without expression data

This model attempts **bivalent domain prediction from sequence alone** — one of
the open challenges in computational epigenomics. Bivalent loci (H3K4me3 +
H3K27me3, observed at developmental gene promoters in ESCs) are predicted to
emerge at positions that have both a TATA-box motif (recruiting the writer) and
a CpG-rich motif (recruiting PRC2), but lack active transcription to establish
H3K36me3 protection.

The **cross-talk score** (Pearson r of binned L1 and L2 signals) quantifies how
strongly the interaction rules enforce mutual exclusivity. Without rules (naive
model) r ≈ 0. With Rules 1–5, r → −1.

A **telomere-proximity test** checks whether the model correctly predicts that
subtelomeric genes (known to be silenced by the SIR complex in *S. cerevisiae*)
are enriched for repressed promoter states — a biologically testable prediction
made purely from sequence and the five interaction rules.

## Installation

```r
install.packages("devtools")
BiocManager::install(c("BSgenome", "GenomicRanges", "Biostrings",
                       "GenomicFeatures",
                       "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"))
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
devtools::install_github("davetgerrard/GenomicLayers", build_vignettes = TRUE)
```

## Usage

```bash
# Full cross-talk model (Rules 1–5)
Rscript rules_based_simulation.R

# Naive model only (no interaction rules — baseline)
Rscript rules_based_simulation.R --naive

# Head-to-head comparison of both models
Rscript rules_based_simulation.R --compare

# Parameter sweep: writer stateWidth vs EED spreading window
Rscript rules_based_simulation.R --sweep
```

## Outputs

| File | Description |
|------|-------------|
| `layer_occupancy_timeseries.png` | Active / repressive / protected / bivalent fraction per cycle |
| `crosstalk_scatter.png` | L1 vs L2 scatter per 1 kb bin — mutual exclusivity test |
| `chromosome_mark_profile.png` | Three-layer mark distribution along the chromosome |
| `gene_promoter_states.png` | Predicted active / repressed / bivalent / unmarked per gene |
| `telomere_proximity_test.png` | Repression enrichment in subtelomeric vs. internal genes |
| `final_state_breakdown.png` | Final chromatin state proportions (cross-talk model) |
| `parameter_sweep_heatmap.png` | r(L1,L2) mutual exclusivity across writer × EED stateWidth |
| `simulation_timeseries.csv` | Numeric timeseries data |
| `binned_mark_profiles.csv` | 1 kb-binned chromosomal profiles |
| `gene_promoter_predictions.csv` | Per-gene promoter state predictions |
| `parameter_sweep_results.csv` | Full sweep results |

## Skills demonstrated

- R, GenomicLayers, ggplot2, GenomicRanges, TxDb, Bioconductor
- Epigenomics: histone modifications, Polycomb, reader-writer feedback, bivalency
- State-dependent (conditional) binding factor modelling
- Cross-talk score / mutual exclusivity quantification
- Gene-level prediction from sequence without expression data
- Telomere-proximity enrichment test

## Key references

- Gerrard (2025) *GenomicLayers: sequence-based simulation of epi-genomes.* BMC Bioinformatics. https://doi.org/10.1186/s12859-025-06224-y
- Margueron et al. (2009) *Role of the polycomb protein EED in the propagation of repressive histone marks.* Nature.
- Streubel et al. / Nat Struct Mol Biol (2025) *Structural basis for the inhibition of PRC2 by active transcription histone PTMs.*
- Finogenova et al. (2021) *Structural basis for PRC2 decoding of active histone methylation marks H3K36me2/3.* eLife.
- Blackledge et al. (2015) *Variant PRC1 complex-dependent H2A ubiquitylation drives PRC2 recruitment and polycomb domain formation.* Cell.
