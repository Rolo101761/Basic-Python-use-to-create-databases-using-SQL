# Enzyme Digestion Step-by-Step Pipeline

A pure-Python in-silico proteomics pipeline with **no external dependencies**.

## What it does

| Step | Module | Description |
|------|--------|-------------|
| 1 | ORF Finder | Extracts protein-coding ORFs from a bacterial genome FASTA (all 6 reading frames) |
| 2 | Protein Digester | Simulates proteolysis using Trypsin, Lys-C, Arg-C, or Glu-C with missed cleavage support |
| 3 | Mass Analyzer | Calculates monoisotopic or average peptide masses and converts to m/z (+1 charge) |
| 4 | Ion Statistics | Identifies which protease yields the most unique diagnostic ions in a target m/z window |

## Biological context

Mass spectrometry-based proteomics requires choosing the right protease for protein identification. This pipeline determines **which enzyme produces the highest proportion of unique peptides** in the 1000–1500 m/z range — maximising the chance of unambiguous protein identification.

## Usage

```bash
# Run all enzymes, compare statistics
python enzyme_digestion_pipeline.py --fasta genome.fasta --enzyme all

# Single enzyme, custom m/z window
python enzyme_digestion_pipeline.py --fasta genome.fasta --enzyme trypsin --mz_min 800 --mz_max 1200

# Average masses, 2 missed cleavages
python enzyme_digestion_pipeline.py --fasta genome.fasta --enzyme lysc --mode average --missed 2
```

## Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Input DNA FASTA file |
| `--enzyme` | `all` | `trypsin`, `lysc`, `argc`, `gluc`, or `all` |
| `--mode` | `monoisotopic` | Mass mode: `monoisotopic` or `average` |
| `--mz_min` | `1000` | Lower m/z window bound |
| `--mz_max` | `1500` | Upper m/z window bound |
| `--min_orf` | `100` | Minimum ORF length (amino acids) |
| `--missed` | `1` | Maximum missed cleavages allowed |

## Enzyme cleavage rules implemented

| Enzyme | Cleaves after | Not before | Notes |
|--------|--------------|------------|-------|
| Trypsin | K, R | P | Most common proteomics enzyme |
| Lys-C | K | — | Produces longer peptides |
| Arg-C | R | P | Complementary to Lys-C |
| Glu-C | D, E | P | Useful for acidic regions |

## Example output

```
[1/4] Parsing FASTA: genome.fasta
      → 1 sequences loaded

[2/4] Finding ORFs (min length=100 aa, all 6 frames)
      NC_000913: 287 ORFs found

[3/4] Digesting proteins with: trypsin, lysc, argc, gluc
      trypsin     :  18432 peptides  |   2341 in m/z 1000–1500
      lysc        :   9821 peptides  |   1876 in m/z 1000–1500
      argc        :  11203 peptides  |   2104 in m/z 1000–1500
      gluc        :  14567 peptides  |   1923 in m/z 1000–1500

[4/4] Running ion statistics

============================================================
Ion Statistics: m/z window 1000–1500 (monoisotopic)
============================================================
Enzyme         In window     Unique   % Unique
------------ ---------- ---------- ----------
trypsin            2341       1987      84.9%
lysc               1876       1654      88.2%
argc               2104       1743      82.8%
gluc               1923       1612      83.8%

  Best enzyme for unambiguous ID: lysc (88.2% unique peptides)
```

## Skills demonstrated

- Pure Python (no external libraries)
- Molecular biology: codon tables, 6-frame translation, enzyme kinetics
- Mass spectrometry: monoisotopic/average masses, m/z calculation
- Algorithm design: ORF detection, string parsing, set operations
- Command-line interface with argparse
