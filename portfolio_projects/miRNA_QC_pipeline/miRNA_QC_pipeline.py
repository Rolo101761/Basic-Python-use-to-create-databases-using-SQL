import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
# os function store figures in folder and allows python to work with folders/files
# define column names to easier names for easier coding
COL_MIRNA = "miRNA"
COL_GENE = "Gene"
COL_LOCUS = "Locus"
COL_DISEASE = "diseaseName"
COL_R2L = "rna2locus_conf_score"
COL_G2D = "gene2disease_conf_score"
COL_PMIDS = "NofPmids"
COL_SNPS = "Nofsnps"
COL_ORG = "Organism"

# argparse function. allows script to be run in terminal and forces user to provide
# required arguments.
def parse_args():
    parser = argparse.ArgumentParser(description="Reads a dataset of miRNA-gene-locus and gene-disease association, quality control perfomed under strict conditions prints basic statistics before and after QC with a figure produced. Thereafter output filtered TSV file.")
    # 4 key arguments
    parser.add_argument("--input", required=True, help="Path to the input CSV.")
    parser.add_argument("--output", required=True, help="Path to save the QC output TSV.")
    parser.add_argument("--rna2locus-threshold", type=float, required=True,
                        help="Keep rows where rna2locus_conf_score > 0.7.")
    parser.add_argument("--gene2disease-threshold", type=float, required=True,
                        help="Keep rows where gene2disease_conf_score > 0.7")
    # needed if user wants figures yes
    parser.add_argument("--fig-dir", default=None, help="Folders to save figures (optional).")
    return parser.parse_args()  # reads command lines arguments above and stores them for script usage

# print summary statistcis
def print_numeric_summary(df, cols):
    available = [c for c in cols if c in df.columns]  # ensure we get numeric columns of interest
    # convert numbers just in case of strings.
    numbers = df[available].apply(pd.to_numeric, errors="coerce")
    # calculate and display stats ddof means degree of freedom
    summary = pd.DataFrame({
        "mean": numbers.mean(),
        "std": numbers.std(ddof=1),
        "median": numbers.median(),
        "min": numbers.min(),
        "max": numbers.max()
    })
    # print clear table
    print(summary.to_string(float_format=lambda x: f"{x:.6g}"))
    # lambda creates small one line function in this script it helps format numbers
    # neatly for display

# print top 3 organisms
def print_top3_organisms(df):
    if COL_ORG not in df.columns:
        print("No organism column found.")
        return  # if org col is missing
    # count occurennces and print the top 3
    counts = df[COL_ORG].dropna().astype(str).value_counts().head(3)
    if counts.empty:
        print("Top 3 organisms: (no organism data)")
        return
    print("Top 3 organisms:")
    for org, count in counts.items():
        print(f"• {org} {int(count)}")

# quality control
def qc_filter(df, rna_thresh, gene_thresh):
    filtered_sections = []
    if all(c in df.columns for c in [COL_MIRNA, COL_GENE, COL_LOCUS, COL_R2L]):
        t = df.copy()
        t[COL_R2L] = pd.to_numeric(t[COL_R2L], errors="coerce")
        mask = (
            t[[COL_MIRNA, COL_GENE, COL_LOCUS, COL_R2L]].notna().all(axis=1) &
            (t[COL_R2L] > rna_thresh)
        )
        filtered_sections.append(t[mask])
    if all(c in df.columns for c in [COL_GENE, COL_DISEASE, COL_G2D]):
        t = df.copy()
        t[COL_G2D] = pd.to_numeric(t[COL_G2D], errors="coerce")
        mask = (
            t[[COL_GENE, COL_DISEASE, COL_G2D]].notna().all(axis=1) &
            (t[COL_G2D] > gene_thresh)
        )
        filtered_sections.append(t[mask])
    if not filtered_sections:
        return df.iloc[0:0].copy()  # returns empty dataframe if nothing quals
    return pd.concat(filtered_sections).drop_duplicates().reset_index(drop=True)

# figure generation
def generate_figures(df, fig_dir):
    os.makedirs(fig_dir, exist_ok=True)  # creates folders for figures

    def hist(col, title):
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        if values.empty:
            return
        plt.figure()
        plt.hist(values, bins=30, edgecolor="black")
        plt.axvline(values.mean(), color="red", label=f"Mean={values.mean():.3f}")
        plt.title(title)
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.legend()
        plt.savefig(f"{fig_dir}/{col}_hist.png")
        plt.close()

    for col in [COL_R2L, COL_G2D, COL_PMIDS, COL_SNPS]:
        if col in df.columns:
            hist(col, f"Distribution of {col}")

    counts = df[COL_ORG].dropna().astype(str).value_counts().head(3)
    if not counts.empty:
        plt.figure()
        counts.plot(kind="bar")
        plt.title("Top 3 organisms (post-QC)")
        plt.ylabel("count")
        plt.savefig(f"{fig_dir}/organisms_top3.png")
        plt.close()

# main function
def main():
    args = parse_args()
    df = pd.read_csv(args.input, sep=None, engine="python")
    print("\n=== BEFORE QC===")
    print_numeric_summary(df, [COL_R2L, COL_G2D, COL_PMIDS, COL_SNPS])
    print_top3_organisms(df)
    # apply QC
    filtered = qc_filter(df, args.rna2locus_threshold, args.gene2disease_threshold)
    print("\n=== AFTER QC ===")
    print_numeric_summary(filtered, [COL_R2L, COL_G2D, COL_PMIDS, COL_SNPS])
    print_top3_organisms(filtered)
    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved QC output to: {args.output}")
    if args.fig_dir:
        generate_figures(filtered, args.fig_dir)
        print(f"Figures saved to: {args.fig_dir}")

if __name__ == "__main__":
    main()  # starts the program.
