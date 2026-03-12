"""
Step 1 — Data ingestion and quality control
RNA-seq Differential Expression Pipeline
Dataset: GSE183947 (breast cancer tumour vs normal, GEO)

Downloads raw count matrix and sample metadata from GEO, performs
library-size QC, flags outlier samples, and writes cleaned outputs
ready for DESeq2.

Outputs (written to data/processed/):
    cleaned_counts.tsv      integer count matrix (genes x samples)
    sample_metadata.tsv     sample annotations with condition column
    qc_report.tsv           per-sample QC metrics
    library_sizes.png       bar chart of library sizes coloured by condition

Usage:
    python scripts/01_data_download.py
    python scripts/01_data_download.py --geo GSE183947 --outdir data/processed
    python scripts/01_data_download.py --geo GSE183947 --outdir data/processed --min-library-size 500000
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    import GEOparse
except ImportError:
    sys.exit("GEOparse not found. Run: pip install GEOparse")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GSE183947 from GEO and run library-size QC."
    )
    parser.add_argument(
        "--geo", default="GSE183947",
        help="GEO accession to download (default: GSE183947)"
    )
    parser.add_argument(
        "--outdir", default="data/processed",
        help="Directory for processed output files (default: data/processed)"
    )
    parser.add_argument(
        "--rawdir", default="data/raw",
        help="Directory for raw GEO download cache (default: data/raw)"
    )
    parser.add_argument(
        "--min-library-size", type=int, default=500_000,
        dest="min_lib",
        help="Minimum total read count per sample (default: 500000)"
    )
    parser.add_argument(
        "--max-zero-fraction", type=float, default=0.95,
        dest="max_zeros",
        help="Remove samples where more than this fraction of genes are zero (default: 0.95)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Download from GEO
# ---------------------------------------------------------------------------
def download_geo(accession, raw_dir):
    os.makedirs(raw_dir, exist_ok=True)
    print(f"[1/5] Downloading {accession} from GEO...")
    gse = GEOparse.get_GEO(geo=accession, destdir=raw_dir, silent=True)
    print(f"      {len(gse.gsms)} samples found in {accession}")
    return gse


# ---------------------------------------------------------------------------
# Extract count matrix and metadata
# ---------------------------------------------------------------------------
def extract_counts_and_metadata(gse):
    print("[2/5] Extracting count matrix and sample metadata...")

    sample_data = {}
    metadata_rows = []

    for gsm_name, gsm in gse.gsms.items():
        # Extract count table (stored in gsm.table)
        if gsm.table is None or gsm.table.empty:
            print(f"      WARNING: no count table for {gsm_name}, skipping")
            continue

        tbl = gsm.table.copy()

        # GEO tables typically have 'ID_REF' and 'VALUE' columns
        # Try to find gene identifier and count columns
        id_col    = tbl.columns[0]
        count_col = tbl.columns[1] if len(tbl.columns) > 1 else tbl.columns[0]

        tbl = tbl[[id_col, count_col]].copy()
        tbl.columns = ["gene", "count"]
        tbl["count"] = pd.to_numeric(tbl["count"], errors="coerce")
        tbl = tbl.dropna(subset=["count"])
        tbl["count"] = tbl["count"].astype(int)

        sample_data[gsm_name] = tbl.set_index("gene")["count"]

        # Metadata from gsm.metadata
        meta = gsm.metadata
        condition = "unknown"
        title = meta.get("title", [""])[0].lower()
        source = meta.get("source_name_ch1", [""])[0].lower()
        char  = " ".join(meta.get("characteristics_ch1", [])).lower()

        for term in ["tumor", "tumour", "cancer", "carcinoma"]:
            if term in title or term in source or term in char:
                condition = "tumour"
                break
        if condition == "unknown":
            for term in ["normal", "adjacent", "healthy"]:
                if term in title or term in source or term in char:
                    condition = "normal"
                    break

        metadata_rows.append({
            "sample_id": gsm_name,
            "title":     meta.get("title", [""])[0],
            "condition": condition,
            "organism":  meta.get("organism_ch1", [""])[0],
            "source":    meta.get("source_name_ch1", [""])[0],
        })

    # Build count matrix
    counts_df = pd.DataFrame(sample_data).fillna(0).astype(int)
    metadata_df = pd.DataFrame(metadata_rows).set_index("sample_id")

    # Align metadata to count matrix columns
    metadata_df = metadata_df.loc[counts_df.columns]

    print(f"      Count matrix: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
    cond_counts = metadata_df["condition"].value_counts().to_dict()
    print(f"      Conditions: {cond_counts}")

    return counts_df, metadata_df


# ---------------------------------------------------------------------------
# QC filtering
# ---------------------------------------------------------------------------
def run_qc(counts_df, metadata_df, min_lib, max_zeros):
    print(f"[3/5] Running QC (min library size={min_lib:,}, max zero fraction={max_zeros})...")

    n_start = counts_df.shape[1]

    # Per-sample metrics
    lib_sizes    = counts_df.sum(axis=0)
    zero_fracs   = (counts_df == 0).mean(axis=0)
    median_count = counts_df.median(axis=0)

    qc_df = pd.DataFrame({
        "library_size":   lib_sizes,
        "zero_fraction":  zero_fracs,
        "median_count":   median_count,
        "condition":      metadata_df["condition"],
    })

    # Flag failing samples
    qc_df["pass_lib"]   = qc_df["library_size"]  >= min_lib
    qc_df["pass_zeros"] = qc_df["zero_fraction"]  <= max_zeros
    qc_df["pass_qc"]    = qc_df["pass_lib"] & qc_df["pass_zeros"]

    n_fail = (~qc_df["pass_qc"]).sum()
    if n_fail > 0:
        failing = qc_df[~qc_df["pass_qc"]].index.tolist()
        print(f"      Removing {n_fail} sample(s) that failed QC: {failing}")
    else:
        print(f"      All {n_start} samples passed QC")

    passing = qc_df[qc_df["pass_qc"]].index
    counts_clean   = counts_df[passing]
    metadata_clean = metadata_df.loc[passing]

    # Filter low-count genes: keep genes with >= 10 counts in at least 2 samples
    keep_genes = (counts_clean >= 10).sum(axis=1) >= 2
    counts_clean = counts_clean[keep_genes]
    print(f"      Genes retained after low-count filter: {keep_genes.sum():,} / {len(keep_genes):,}")

    return counts_clean, metadata_clean, qc_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_library_sizes(qc_df, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    condition_colours = {"tumour": "#d6604d", "normal": "#2166ac", "unknown": "#999999"}
    colours = [condition_colours.get(c, "#999999") for c in qc_df["condition"]]
    pass_markers = ["solid" if p else "dashed" for p in qc_df["pass_qc"]]

    # Library size bar chart
    ax = axes[0]
    bars = ax.bar(range(len(qc_df)), qc_df["library_size"] / 1e6,
                  color=colours, edgecolor="white", linewidth=0.5)
    # Mark failing samples
    for i, (bar, passed) in enumerate(zip(bars, qc_df["pass_qc"])):
        if not passed:
            bar.set_edgecolor("black")
            bar.set_linewidth(2)
            bar.set_linestyle("--")

    ax.set_xticks(range(len(qc_df)))
    ax.set_xticklabels(qc_df.index, rotation=90, fontsize=7)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Library size (millions of reads)")
    ax.set_title("Library Sizes by Sample")
    legend_patches = [
        mpatches.Patch(color="#d6604d", label="Tumour"),
        mpatches.Patch(color="#2166ac", label="Normal"),
    ]
    ax.legend(handles=legend_patches, loc="upper right")

    # Zero fraction scatter
    ax2 = axes[1]
    for cond, col in condition_colours.items():
        mask = qc_df["condition"] == cond
        if mask.any():
            ax2.scatter(
                qc_df.loc[mask, "library_size"] / 1e6,
                qc_df.loc[mask, "zero_fraction"] * 100,
                c=col, label=cond, alpha=0.8, s=60, edgecolors="white"
            )
    ax2.set_xlabel("Library size (M reads)")
    ax2.set_ylabel("Zero fraction (%)")
    ax2.set_title("Library Size vs Zero Fraction")
    ax2.legend()

    plt.suptitle("GSE183947 — Sample QC Metrics", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"      Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.rawdir, exist_ok=True)

    # Step 1: Download
    gse = download_geo(args.geo, args.rawdir)

    # Step 2: Extract
    counts_df, metadata_df = extract_counts_and_metadata(gse)

    # Step 3: QC
    counts_clean, metadata_clean, qc_df = run_qc(
        counts_df, metadata_df, args.min_lib, args.max_zeros
    )

    # Step 4: Save outputs
    print("[4/5] Writing outputs...")

    counts_path   = os.path.join(args.outdir, "cleaned_counts.tsv")
    metadata_path = os.path.join(args.outdir, "sample_metadata.tsv")
    qc_path       = os.path.join(args.outdir, "qc_report.tsv")
    fig_path      = os.path.join(args.outdir, "library_sizes.png")

    counts_clean.to_csv(counts_path,   sep="\t")
    metadata_clean.to_csv(metadata_path, sep="\t")
    qc_df.to_csv(qc_path,         sep="\t")

    print(f"      Saved: {counts_path}")
    print(f"      Saved: {metadata_path}")
    print(f"      Saved: {qc_path}")

    # Step 5: Plot
    print("[5/5] Generating QC figures...")
    plot_library_sizes(qc_df, fig_path)

    # Summary
    print("\n" + "=" * 55)
    print("Data ingestion complete")
    print("=" * 55)
    print(f"  Samples (post-QC) : {counts_clean.shape[1]}")
    print(f"  Genes retained    : {counts_clean.shape[0]:,}")
    cond = metadata_clean["condition"].value_counts().to_dict()
    for k, v in cond.items():
        print(f"  {k:18s}: {v} samples")
    print(f"\nNext step:")
    print(f"  Rscript scripts/02_deseq2_analysis.R")


if __name__ == "__main__":
    main()
