"""
miRNA QC Pipeline
=================
Quality control and analysis tool for miRNA–gene–disease association datasets.

Pipeline steps:
  1. Load miRNA–gene–disease TSV dataset
  2. Apply confidence-score filtering (user-defined threshold)
  3. Remove entries with missing / incomplete fields
  4. Report summary statistics
  5. Identify top organisms by association count
  6. Output a cleaned TSV ready for downstream analysis

Compatible with miRTarBase, HMDD, and similar TSV-format databases.

Usage:
    python miRNA_QC_pipeline.py --input associations.tsv --confidence 0.7
    python miRNA_QC_pipeline.py --input associations.tsv --confidence 0.8 --output cleaned.tsv --top 10
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# Constants — expected column names (case-insensitive matching applied)
# ---------------------------------------------------------------------------
REQUIRED_FIELDS = ['mirna', 'gene', 'disease']
OPTIONAL_FIELDS = ['organism', 'confidence', 'pmid', 'experiment']


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def detect_delimiter(filepath):
    """Detect whether file uses tab or comma as delimiter."""
    with open(filepath, 'r', encoding='utf-8') as fh:
        sample = fh.read(2048)
    tabs   = sample.count('\t')
    commas = sample.count(',')
    return '\t' if tabs >= commas else ','


def load_dataset(filepath):
    """Load a TSV/CSV miRNA association file. Returns (headers, rows)."""
    if not os.path.exists(filepath):
        print(f"[ERROR] File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    delim = detect_delimiter(filepath)
    rows  = []
    with open(filepath, 'r', encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        headers = [h.strip().lower() for h in (reader.fieldnames or [])]
        for row in reader:
            normalised = {k.strip().lower(): v.strip() for k, v in row.items()}
            rows.append(normalised)

    print(f"[1/5] Loaded {len(rows):,} entries from '{filepath}'")
    print(f"      Columns detected: {', '.join(headers)}")
    return headers, rows


# ---------------------------------------------------------------------------
# QC filters
# ---------------------------------------------------------------------------
def check_required_columns(headers):
    """Verify required columns are present (flexible partial matching)."""
    missing = []
    for req in REQUIRED_FIELDS:
        if not any(req in h for h in headers):
            missing.append(req)
    if missing:
        print(f"[WARNING] Could not find expected columns: {missing}")
        print(f"          Available columns: {headers}")
        print(f"          Proceeding with available data.\n")


def remove_incomplete_entries(rows, headers):
    """Remove rows with empty or null values in any column."""
    before = len(rows)
    clean  = []
    null_tokens = {'', 'na', 'n/a', 'none', 'null', '-', 'unknown'}

    for row in rows:
        incomplete = any(
            str(v).strip().lower() in null_tokens
            for v in row.values()
        )
        if not incomplete:
            clean.append(row)

    removed = before - len(clean)
    print(f"[2/5] Removed {removed:,} incomplete entries  ({before:,} → {len(clean):,})")
    return clean


def filter_by_confidence(rows, threshold, conf_col='confidence'):
    """Keep only rows where confidence score >= threshold."""
    # Check if a confidence column exists
    sample = rows[0] if rows else {}
    conf_key = next((k for k in sample if conf_col in k), None)

    if conf_key is None:
        print(f"[3/5] No confidence column found — skipping confidence filter")
        return rows

    before = len(rows)
    clean  = []
    skipped_non_numeric = 0

    for row in rows:
        raw = row.get(conf_key, '')
        try:
            score = float(raw)
            if score >= threshold:
                clean.append(row)
        except ValueError:
            skipped_non_numeric += 1

    removed = before - len(clean)
    print(f"[3/5] Confidence filter (>= {threshold}): {removed:,} removed  ({before:,} → {len(clean):,})")
    if skipped_non_numeric:
        print(f"      (Skipped {skipped_non_numeric} non-numeric confidence values)")
    return clean


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def compute_summary_statistics(rows, headers):
    """Compute and print summary statistics."""
    print(f"\n[4/5] Summary Statistics")
    print(f"{'='*50}")
    print(f"  Total associations (after QC): {len(rows):,}")

    # Column-level stats
    for col in headers:
        values = [r[col] for r in rows if col in r]
        non_empty = [v for v in values if v.strip()]
        unique    = len(set(non_empty))
        print(f"  {col:<25} {len(non_empty):>7,} entries  |  {unique:>6,} unique values")

    # Confidence distribution (if present)
    conf_key = next((k for k in headers if 'confidence' in k), None)
    if conf_key:
        scores = []
        for row in rows:
            try:
                scores.append(float(row[conf_key]))
            except (ValueError, KeyError):
                pass
        if scores:
            scores.sort()
            n = len(scores)
            mean = sum(scores) / n
            median = scores[n // 2]
            print(f"\n  Confidence scores:")
            print(f"    Min:    {min(scores):.4f}")
            print(f"    Max:    {max(scores):.4f}")
            print(f"    Mean:   {mean:.4f}")
            print(f"    Median: {median:.4f}")


def top_organisms(rows, top_n=10):
    """Identify and display top organisms by association count."""
    org_key = next(
        (k for k in (rows[0] if rows else {}) if 'organism' in k or 'species' in k),
        None
    )
    if org_key is None:
        print(f"\n  No organism/species column found — skipping organism analysis")
        return {}

    counts = Counter(row[org_key] for row in rows if row.get(org_key))

    print(f"\n  Top {top_n} organisms by association count:")
    print(f"  {'Organism':<35} {'Count':>8}")
    print(f"  {'-'*35} {'-'*8}")
    for organism, count in counts.most_common(top_n):
        print(f"  {organism:<35} {count:>8,}")

    return dict(counts.most_common(top_n))


def mirna_disease_summary(rows, top_n=10):
    """Show top miRNAs and diseases by association count."""
    mirna_key   = next((k for k in (rows[0] if rows else {}) if 'mirna' in k or 'mir' in k), None)
    disease_key = next((k for k in (rows[0] if rows else {}) if 'disease' in k), None)
    gene_key    = next((k for k in (rows[0] if rows else {}) if 'gene' in k), None)

    if mirna_key:
        mirna_counts = Counter(row[mirna_key] for row in rows if row.get(mirna_key))
        print(f"\n  Top {top_n} miRNAs:")
        for mirna, count in mirna_counts.most_common(top_n):
            print(f"    {mirna:<25} {count:>6,} associations")

    if disease_key:
        disease_counts = Counter(row[disease_key] for row in rows if row.get(disease_key))
        print(f"\n  Top {top_n} diseases:")
        for disease, count in disease_counts.most_common(top_n):
            print(f"    {disease:<35} {count:>6,} associations")

    if gene_key:
        gene_counts = Counter(row[gene_key] for row in rows if row.get(gene_key))
        print(f"\n  Top {top_n} target genes:")
        for gene, count in gene_counts.most_common(top_n):
            print(f"    {gene:<20} {count:>6,} associations")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_cleaned_tsv(rows, output_path):
    """Write cleaned dataset to a TSV file."""
    if not rows:
        print(f"\n[5/5] No data to write.")
        return
    headers = list(rows[0].keys())
    with open(output_path, 'w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n[5/5] Cleaned dataset written to: {output_path}  ({len(rows):,} rows)")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pipeline(input_path, output_path, confidence_threshold=0.7, top_n=10):
    # Load
    headers, rows = load_dataset(input_path)
    check_required_columns(headers)

    initial_count = len(rows)

    # QC steps
    rows = remove_incomplete_entries(rows, headers)
    rows = filter_by_confidence(rows, confidence_threshold)

    # Statistics
    compute_summary_statistics(rows, headers)
    top_organisms(rows, top_n=top_n)
    mirna_disease_summary(rows, top_n=top_n)

    # Output
    write_cleaned_tsv(rows, output_path)

    # Final summary
    removed_total = initial_count - len(rows)
    pct_retained  = len(rows) / initial_count * 100 if initial_count else 0
    print(f"\n{'='*50}")
    print(f"QC complete: {removed_total:,} entries removed  |  {len(rows):,} retained ({pct_retained:.1f}%)")
    print(f"{'='*50}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='miRNA QC Pipeline — quality control for miRNA–gene–disease datasets'
    )
    parser.add_argument('--input',      required=True,       help='Input TSV/CSV file')
    parser.add_argument('--output',     default='cleaned_miRNA_associations.tsv',
                                                             help='Output cleaned TSV file')
    parser.add_argument('--confidence', type=float, default=0.7,
                                                             help='Minimum confidence score (0–1)')
    parser.add_argument('--top',        type=int,   default=10,
                                                             help='Number of top entries to display')
    args = parser.parse_args()

    run_pipeline(
        input_path=args.input,
        output_path=args.output,
        confidence_threshold=args.confidence,
        top_n=args.top,
    )


if __name__ == '__main__':
    main()
