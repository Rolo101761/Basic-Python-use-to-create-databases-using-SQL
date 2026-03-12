# Enzyme Digestion Pipeline

A four-stage Python pipeline simulating the wet-lab workflow for mass spectrometry-based proteomics.

## Pipeline overview

```
DNA/genome FASTA
      │
      ▼  task1_orf_finder.py
Protein ORF FASTA  (6-frame translation, deduplicated)
      │
      ▼  task2_protein_digestion.py
Peptide FASTA  (in silico enzyme digestion, ≤1 missed cleavage)
      │
      ▼  task3_mass_to_charge.py
m/z TSV  (monoisotopic or average mass, user-defined charge state)
      │
      ▼  task4_peptide_analysis.py
Histograms / unique protein counts  (4 analysis modes)
```

## Scripts

### Task 1 — ORF Finder (`task1_orf_finder.py`)

Finds all open reading frames in all 6 reading frames (3 forward, 3 reverse complement) of a DNA FASTA file.

```bash
python task1_orf_finder.py input.fasta orfs.fasta --min_size 50
```

Output FASTA header: `>ORGID_F{frame}_{index} frame={frame} length={aa} start={pos}`

### Task 2 — Protein Digestion (`task2_protein_digestion.py`)
*Amin Miah*

In silico digestion of protein sequences using four proteases, with up to 1 missed cleavage.

```bash
python task2_protein_digestion.py orfs.fasta -e t -m 1 -min 5 -o peptides.fasta
```

| Enzyme flag | Enzyme | Cuts after |
|-------------|--------|-----------|
| `t` | Trypsin | K, R (not before P) |
| `lys-c` | Lys-C | K |
| `arg-c` | Arg-C | R |
| `glu-c` | Glu-C | E |

### Task 3 — Mass-to-Charge Calculator (`task3_mass_to_charge.py`)

Calculates m/z values for each peptide using monoisotopic or average amino acid masses.

```bash
python task3_mass_to_charge.py --mode m --charge 2 -i peptides.fasta -o mz.tsv
```

Formula: `m/z = (peptide_mass + z × proton_mass) / z`

Output columns: `protein  peptide_num  m/z  z  missed_cleavages  sequence`

### Task 4 — Peptide Analysis (`task4_peptide_analysis.py`)

Four analysis modes on the m/z output:

```bash
python task4_peptide_analysis.py mz.tsv --mode 1 -s 500 -e 3000          # count in range
python task4_peptide_analysis.py mz.tsv --mode 2 -b 1.0 -o hist.tsv      # fixed bin histogram
python task4_peptide_analysis.py mz.tsv --mode 3 --window 10 --step 1.0  # sliding window
python task4_peptide_analysis.py mz.tsv --mode 4 -b 0.01                  # unique protein count
```

| Mode | Description |
|------|-------------|
| 1 | Count peptides within m/z range |
| 2 | Fixed-bin m/z histogram (TSV output) |
| 3 | Sliding window density across m/z range |
| 4 | Count uniquely identifiable proteins by m/z |

## Skills demonstrated

- Python, argparse, pandas, csv, collections
- Molecular biology: ORF finding, 6-frame translation, proteomics, mass spectrometry
- FASTA parsing, sequence manipulation, codon tables
- Enzyme cleavage rules, missed cleavage modelling
- Mass spectrometry: m/z calculation, monoisotopic vs average mass
- Histogram and sliding window analysis
