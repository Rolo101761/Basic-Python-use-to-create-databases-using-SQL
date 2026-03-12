"""
Step 4 — Load results into SQLite database
RNA-seq Differential Expression Pipeline
Dataset: GSE183947 (breast cancer tumour vs normal)

Reads the processed outputs from Steps 1–3 and loads them into a
SQLite database for querying and dashboard connection.

Schema:
    samples         — sample metadata (condition, organism, source)
    genes           — gene annotations (symbol, description)
    de_results      — full DESeq2 results per gene
    enrichment      — GO and KEGG enrichment results
    qc_metrics      — per-sample QC metrics from Step 1

Usage:
    python scripts/04_load_to_db.py
    python scripts/04_load_to_db.py --indir data/processed --db database/rnaseq_results.db
"""

import argparse
import os
import sys
import sqlite3

import pandas as pd


# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Load RNA-seq pipeline results into SQLite database."
    )
    parser.add_argument("--indir", default="data/processed",
                        help="Directory containing processed TSV/CSV files (default: data/processed)")
    parser.add_argument("--db", default="database/rnaseq_results.db",
                        help="Path to SQLite database (default: database/rnaseq_results.db)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and recreate all tables if database already exists")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_file(path, sep="\t", **kwargs):
    if not os.path.exists(path):
        print(f"  WARNING: file not found, skipping — {path}")
        return None
    df = pd.read_csv(path, sep=sep, **kwargs)
    print(f"  Loaded {len(df):,} rows from {os.path.basename(path)}")
    return df


def df_to_sqlite(df, table_name, conn, overwrite):
    if_exists = "replace" if overwrite else "append"
    df.to_sql(table_name, conn, if_exists=if_exists, index=False)
    print(f"  Wrote {len(df):,} rows → table '{table_name}'")


# ---------------------------------------------------------------------------
# Schema creation
# ---------------------------------------------------------------------------
CREATE_STATEMENTS = """
CREATE TABLE IF NOT EXISTS samples (
    sample_id   TEXT PRIMARY KEY,
    title       TEXT,
    condition   TEXT,
    organism    TEXT,
    source      TEXT
);

CREATE TABLE IF NOT EXISTS de_results (
    gene            TEXT PRIMARY KEY,
    baseMean        REAL,
    log2FoldChange  REAL,
    lfcSE           REAL,
    stat            REAL,
    pvalue          REAL,
    padj            REAL,
    significant     INTEGER   -- 1 if padj<0.05 and |log2FC|>1
);

CREATE TABLE IF NOT EXISTS enrichment (
    id          TEXT,
    description TEXT,
    source      TEXT,        -- GO_BP / GO_MF / KEGG
    gene_ratio  TEXT,
    bg_ratio    TEXT,
    pvalue      REAL,
    padj        REAL,
    gene_count  INTEGER,
    gene_ids    TEXT
);

CREATE TABLE IF NOT EXISTS qc_metrics (
    sample_id       TEXT PRIMARY KEY,
    library_size    INTEGER,
    zero_fraction   REAL,
    median_count    REAL,
    condition       TEXT,
    pass_qc         INTEGER
);
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.db) or ".", exist_ok=True)

    print("=" * 55)
    print("Loading RNA-seq results into SQLite")
    print("=" * 55)
    print(f"  Database : {args.db}")
    print(f"  Input dir: {args.indir}\n")

    conn = sqlite3.connect(args.db)

    if args.overwrite:
        for tbl in ["samples", "de_results", "enrichment", "qc_metrics"]:
            conn.execute(f"DROP TABLE IF EXISTS {tbl}")
        print("  Dropped existing tables (--overwrite)\n")

    conn.executescript(CREATE_STATEMENTS)
    conn.commit()

    # ── Samples ─────────────────────────────────────────────────────────
    print("[1/4] Loading sample metadata...")
    meta = load_file(os.path.join(args.indir, "sample_metadata.tsv"))
    if meta is not None:
        meta = meta.reset_index() if "sample_id" not in meta.columns else meta
        if meta.columns[0] != "sample_id":
            meta = meta.rename(columns={meta.columns[0]: "sample_id"})
        meta = meta[["sample_id", "title", "condition", "organism", "source"]
                    if all(c in meta.columns for c in ["title", "condition", "organism", "source"])
                    else meta.columns.tolist()]
        df_to_sqlite(meta, "samples", conn, overwrite=True)

    # ── DE results ───────────────────────────────────────────────────────
    print("\n[2/4] Loading DESeq2 results...")
    de = load_file(os.path.join(args.indir, "DESeq2_results.csv"), sep=",")
    if de is not None:
        de["significant"] = (
            (de["padj"] < 0.05) & (de["log2FoldChange"].abs() > 1)
        ).astype(int)
        keep_cols = [c for c in ["gene", "baseMean", "log2FoldChange", "lfcSE",
                                  "stat", "pvalue", "padj", "significant"]
                     if c in de.columns]
        df_to_sqlite(de[keep_cols], "de_results", conn, overwrite=True)

    # ── Enrichment ───────────────────────────────────────────────────────
    print("\n[3/4] Loading enrichment results...")
    enr = load_file(os.path.join(args.indir, "enrichment_results.tsv"))
    if enr is not None:
        col_map = {
            "ID": "id", "Description": "description", "source": "source",
            "GeneRatio": "gene_ratio", "BgRatio": "bg_ratio",
            "pvalue": "pvalue", "p.adjust": "padj",
            "Count": "gene_count", "geneID": "gene_ids"
        }
        enr = enr.rename(columns={k: v for k, v in col_map.items() if k in enr.columns})
        keep = [v for v in col_map.values() if v in enr.columns]
        df_to_sqlite(enr[keep], "enrichment", conn, overwrite=True)

    # ── QC metrics ───────────────────────────────────────────────────────
    print("\n[4/4] Loading QC metrics...")
    qc = load_file(os.path.join(args.indir, "qc_report.tsv"))
    if qc is not None:
        qc = qc.reset_index() if "sample_id" not in qc.columns else qc
        if qc.columns[0] != "sample_id":
            qc = qc.rename(columns={qc.columns[0]: "sample_id"})
        qc_cols = [c for c in ["sample_id", "library_size", "zero_fraction",
                                "median_count", "condition", "pass_qc"]
                   if c in qc.columns]
        df_to_sqlite(qc[qc_cols], "qc_metrics", conn, overwrite=True)

    conn.commit()

    # ── Verify with example queries ──────────────────────────────────────
    print("\n" + "=" * 55)
    print("Database verification queries")
    print("=" * 55)

    queries = {
        "Total DE genes":
            "SELECT COUNT(*) FROM de_results WHERE significant = 1",
        "Top 5 upregulated genes":
            "SELECT gene, ROUND(log2FoldChange,3) AS log2FC, ROUND(padj,6) AS padj "
            "FROM de_results WHERE significant=1 AND log2FoldChange>0 "
            "ORDER BY log2FoldChange DESC LIMIT 5",
        "Top 5 GO-BP terms":
            "SELECT description, gene_count, ROUND(padj,6) AS padj "
            "FROM enrichment WHERE source='GO_BP' ORDER BY padj LIMIT 5",
        "Sample counts by condition":
            "SELECT condition, COUNT(*) AS n FROM samples GROUP BY condition",
    }

    for label, sql in queries.items():
        print(f"\n  {label}:")
        try:
            result = pd.read_sql(sql, conn)
            print(result.to_string(index=False))
        except Exception as e:
            print(f"  (query skipped: {e})")

    conn.close()

    print(f"\nDatabase saved to: {args.db}")
    print("\nNext step:")
    print("  python scripts/05_db_queries.py")


if __name__ == "__main__":
    main()
