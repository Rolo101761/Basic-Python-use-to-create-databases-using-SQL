"""
Enzyme Digestion Step-by-Step Pipeline
=======================================
In-silico proteomics pipeline:
  1. Find ORFs from a bacterial genome FASTA file
  2. Digest proteins with multiple enzymes (Trypsin, Lys-C, Arg-C, Glu-C)
  3. Calculate peptide m/z values (+1 charge, monoisotopic/average mass)
  4. Perform ion statistics to identify best protease for unambiguous identification

Usage:
    python enzyme_digestion_pipeline.py --fasta genome.fasta --enzyme trypsin --mode monoisotopic
    python enzyme_digestion_pipeline.py --fasta genome.fasta --enzyme all --mz_min 1000 --mz_max 1500
"""

import argparse
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Amino acid masses (monoisotopic and average)
# ---------------------------------------------------------------------------
MONOISOTOPIC = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}
AVERAGE = {
    'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
    'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
    'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
    'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
}
WATER_MONO = 18.01056
WATER_AVG  = 18.0153
PROTON     = 1.00728


# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------
def parse_fasta(filepath):
    """Parse a FASTA file, returning {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    return sequences


# ---------------------------------------------------------------------------
# ORF Finder
# ---------------------------------------------------------------------------
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def reverse_complement(seq):
    comp = str.maketrans('ACGT', 'TGCA')
    return seq.translate(comp)[::-1]


def translate(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def find_orfs(dna_seq, min_length=100):
    """
    Find all ORFs in all 6 reading frames.
    Returns list of (protein_seq, frame, start, strand).
    """
    orfs = []
    for strand, seq in [('+', dna_seq), ('-', reverse_complement(dna_seq))]:
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                codon = seq[i:i+3]
                if codon == 'ATG':
                    # Scan forward for stop codon
                    j = i + 3
                    while j < len(seq) - 2:
                        stop = seq[j:j+3]
                        if stop in STOP_CODONS:
                            orf_dna = seq[i:j+3]
                            if len(orf_dna) >= min_length * 3:
                                protein = translate(orf_dna)
                                if len(protein) >= min_length:
                                    orfs.append((protein, frame + 1, i, strand))
                            break
                        j += 3
                    i = j + 3
                else:
                    i += 3
    return orfs


# ---------------------------------------------------------------------------
# Enzyme cleavage rules
# ---------------------------------------------------------------------------
ENZYMES = {
    'trypsin': {
        'cleave_after': set('KR'),
        'not_before':   set('P'),
        'description':  'Cleaves after K/R, not before P'
    },
    'lysc': {
        'cleave_after': set('K'),
        'not_before':   set(),
        'description':  'Cleaves after K only'
    },
    'argc': {
        'cleave_after': set('R'),
        'not_before':   set('P'),
        'description':  'Cleaves after R, not before P'
    },
    'gluc': {
        'cleave_after': set('DE'),
        'not_before':   set('P'),
        'description':  'Cleaves after D/E, not before P'
    },
}


def digest_protein(protein, enzyme_name, missed_cleavages=1):
    """
    Digest a protein sequence using the specified enzyme.
    Returns a list of peptide strings.
    """
    rules = ENZYMES[enzyme_name.lower()]
    cleave_after = rules['cleave_after']
    not_before   = rules['not_before']

    # Find cleavage sites
    sites = [0]
    for i, aa in enumerate(protein[:-1]):
        next_aa = protein[i + 1]
        if aa in cleave_after and next_aa not in not_before:
            sites.append(i + 1)
    sites.append(len(protein))

    # Generate peptides with up to N missed cleavages
    peptides = []
    for i in range(len(sites) - 1):
        for mc in range(missed_cleavages + 1):
            if i + mc + 1 < len(sites):
                pep = protein[sites[i]:sites[i + mc + 1]]
                if len(pep) >= 6:  # discard very short peptides
                    peptides.append(pep)
    return peptides


# ---------------------------------------------------------------------------
# Mass calculator
# ---------------------------------------------------------------------------
def peptide_mass(peptide, mode='monoisotopic'):
    masses = MONOISOTOPIC if mode == 'monoisotopic' else AVERAGE
    water  = WATER_MONO   if mode == 'monoisotopic' else WATER_AVG
    mass = water
    for aa in peptide:
        mass += masses.get(aa, 0.0)
    return mass


def peptide_mz(peptide, charge=1, mode='monoisotopic'):
    mass = peptide_mass(peptide, mode)
    return (mass + charge * PROTON) / charge


# ---------------------------------------------------------------------------
# Ion statistics
# ---------------------------------------------------------------------------
def ion_statistics(all_enzyme_peptides, mz_min=1000, mz_max=1500, mode='monoisotopic'):
    """
    For each enzyme, count peptides whose m/z falls in [mz_min, mz_max],
    and calculate the fraction of those that are unique across ALL enzymes.
    """
    # Collect all peptides per enzyme in the m/z window
    enzyme_window_peptides = {}
    for enzyme, peptides in all_enzyme_peptides.items():
        in_window = set()
        for pep in peptides:
            mz = peptide_mz(pep, charge=1, mode=mode)
            if mz_min <= mz <= mz_max:
                in_window.add(pep)
        enzyme_window_peptides[enzyme] = in_window

    # Count how many enzymes each peptide appears in
    peptide_counts = defaultdict(int)
    for peptides in enzyme_window_peptides.values():
        for pep in peptides:
            peptide_counts[pep] += 1

    # Report
    print(f"\n{'='*60}")
    print(f"Ion Statistics: m/z window {mz_min}–{mz_max} ({mode})")
    print(f"{'='*60}")
    print(f"{'Enzyme':<12} {'In window':>10} {'Unique':>10} {'% Unique':>10}")
    print(f"{'-'*12} {'-'*10} {'-'*10} {'-'*10}")

    best_enzyme = None
    best_score  = -1
    for enzyme, peptides in enzyme_window_peptides.items():
        unique = sum(1 for p in peptides if peptide_counts[p] == 1)
        pct    = (unique / len(peptides) * 100) if peptides else 0
        print(f"{enzyme:<12} {len(peptides):>10} {unique:>10} {pct:>9.1f}%")
        if pct > best_score:
            best_score  = pct
            best_enzyme = enzyme

    print(f"\n  Best enzyme for unambiguous ID: {best_enzyme} ({best_score:.1f}% unique peptides)")
    return enzyme_window_peptides


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pipeline(fasta_path, enzyme_choice='all', mz_min=1000, mz_max=1500,
                 mode='monoisotopic', min_orf=100, missed=1):

    print(f"[1/4] Parsing FASTA: {fasta_path}")
    sequences = parse_fasta(fasta_path)
    print(f"      → {len(sequences)} sequences loaded")

    print(f"\n[2/4] Finding ORFs (min length={min_orf} aa, all 6 frames)")
    all_proteins = []
    for header, dna in sequences.items():
        orfs = find_orfs(dna, min_length=min_orf)
        all_proteins.extend([orf[0] for orf in orfs])
        print(f"      {header}: {len(orfs)} ORFs found")
    print(f"      → Total proteins: {len(all_proteins)}")

    enzymes_to_use = list(ENZYMES.keys()) if enzyme_choice == 'all' else [enzyme_choice.lower()]

    print(f"\n[3/4] Digesting proteins with: {', '.join(enzymes_to_use)}")
    all_enzyme_peptides = {}
    for enzyme in enzymes_to_use:
        peptides = []
        for protein in all_proteins:
            peptides.extend(digest_protein(protein, enzyme, missed_cleavages=missed))
        all_enzyme_peptides[enzyme] = peptides
        mz_values = [peptide_mz(p, mode=mode) for p in peptides]
        in_window = sum(1 for m in mz_values if mz_min <= m <= mz_max)
        print(f"      {enzyme:<12}: {len(peptides):>6} peptides  |  {in_window:>5} in m/z {mz_min}–{mz_max}")

    print(f"\n[4/4] Running ion statistics")
    ion_statistics(all_enzyme_peptides, mz_min=mz_min, mz_max=mz_max, mode=mode)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Enzyme Digestion Pipeline — in-silico proteomics'
    )
    parser.add_argument('--fasta',   required=True,          help='Input FASTA file (DNA)')
    parser.add_argument('--enzyme',  default='all',
                        choices=['trypsin', 'lysc', 'argc', 'gluc', 'all'],
                        help='Enzyme to use (default: all)')
    parser.add_argument('--mode',    default='monoisotopic',
                        choices=['monoisotopic', 'average'],
                        help='Mass mode (default: monoisotopic)')
    parser.add_argument('--mz_min',  type=float, default=1000, help='m/z window lower bound')
    parser.add_argument('--mz_max',  type=float, default=1500, help='m/z window upper bound')
    parser.add_argument('--min_orf', type=int,   default=100,  help='Minimum ORF length in aa')
    parser.add_argument('--missed',  type=int,   default=1,    help='Max missed cleavages')
    args = parser.parse_args()

    run_pipeline(
        fasta_path=args.fasta,
        enzyme_choice=args.enzyme,
        mz_min=args.mz_min,
        mz_max=args.mz_max,
        mode=args.mode,
        min_orf=args.min_orf,
        missed=args.missed,
    )


if __name__ == '__main__':
    main()
