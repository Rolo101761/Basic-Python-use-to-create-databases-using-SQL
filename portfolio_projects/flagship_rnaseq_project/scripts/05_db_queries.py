"""
Step 5 — Example database queries
RNA-seq Differential Expression Pipeline

Demonstrates the kinds of queries an analyst or dashboard would run
against the SQLite results database. Exports CSVs for use in
Power BI / Tableau dashboards.

Exported CSVs (to data/processed/dashboard_exports/):
    top_de_genes.csv            top 50 DE genes for the gene table
    volcano_data.csv            all genes with fold change and significance
    go_terms.csv                GO enrichment for pathway bar chart
    kegg_pathways.csv           KEGG enrichment for pathway table
    sample_overview.csv         sample metadata + library sizes for QC page
    gene_expression_wide.csv    VST counts pivoted wide for gene deep-dive

Usage:
    python scripts/05_db_queries.py
    python scripts/05_db_queries.py --db database/rnaseq_results.db --outdir data/processed/dashboard_exports
"""

import argparse
import os
import sqlite3

import pandas as pd


# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Run example queries against RNA-seq SQLite database and export CSVs for dashboards."
    )
    parser.add_argument("--db", default="database/rnaseq_results.db",
                        help="SQLite database path (default: database/rnaseq_results.db)")
    parser.add_argument("--outdir", default="data/processed/dashboard_exports",
                        help="Output directory for CSV exports (default: data/processed/dashboard_exports)")
    parser.add_argument("--vst", default="data/processed/vst_counts.csv",
                        help="VST counts CSV from Step 2 (for gene expression export)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Query helper
# ---------------------------------------------------------------------------
def query(conn, sql, label):
    try:
        df = pd.read_sql(sql, conn)
        print(f"  {label}: {len(df)} rows")
        return df
    except Exception as e:
        print(f"  WARNING — {label} failed: {e}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    if not os.path.exists(args.db):
        raise FileNotFoundError(
            f"Database not found: {args.db}\nRun 04_load_to_db.py first."
        )

    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 55)
    print("Database queries and dashboard exports")
    print("=" * 55)
    print(f"  Database : {args.db}")
    print(f"  Exports  : {args.outdir}\n")

    conn = sqlite3.connect(args.db)

    # ── 1. Top DE genes (for gene table in dashboard) ───────────────────
    top_de = query(conn, """
        SELECT
            gene,
            ROUND(baseMean, 2)        AS baseMean,
            ROUND(log2FoldChange, 4)  AS log2FoldChange,
            ROUND(lfcSE, 4)           AS lfcSE,
            ROUND(pvalue, 8)          AS pvalue,
            ROUND(padj, 8)            AS padj,
            significant,
            CASE
                WHEN log2FoldChange > 0 THEN 'Up'
                WHEN log2FoldChange < 0 THEN 'Down'
                ELSE 'NC'
            END AS direction
        FROM de_results
        WHERE significant = 1
        ORDER BY ABS(log2FoldChange) DESC
        LIMIT 50
    """, "Top 50 DE genes")

    # ── 2. Volcano data (all genes with fold change) ─────────────────────
    volcano = query(conn, """
        SELECT
            gene,
            ROUND(log2FoldChange, 4) AS log2FoldChange,
            ROUND(padj, 8)           AS padj,
            ROUND(-LOG(padj) / LOG(10), 4) AS neg_log10_padj,
            significant,
            CASE
                WHEN padj < 0.05 AND log2FoldChange >  1 THEN 'Up'
                WHEN padj < 0.05 AND log2FoldChange < -1 THEN 'Down'
                ELSE 'NS'
            END AS category
        FROM de_results
        WHERE log2FoldChange IS NOT NULL AND padj IS NOT NULL
        ORDER BY padj
    """, "Volcano data")

    # ── 3. GO BP enrichment ───────────────────────────────────────────────
    go_terms = query(conn, """
        SELECT
            id,
            description,
            source,
            gene_ratio,
            gene_count,
            ROUND(pvalue, 8) AS pvalue,
            ROUND(padj,   8) AS padj,
            ROUND(-LOG(padj) / LOG(10), 4) AS neg_log10_padj,
            gene_ids
        FROM enrichment
        WHERE source = 'GO_BP'
        ORDER BY padj
        LIMIT 50
    """, "GO-BP terms")

    # ── 4. KEGG pathways ──────────────────────────────────────────────────
    kegg_paths = query(conn, """
        SELECT
            id,
            description,
            gene_ratio,
            gene_count,
            ROUND(pvalue, 8) AS pvalue,
            ROUND(padj,   8) AS padj,
            gene_ids
        FROM enrichment
        WHERE source = 'KEGG'
        ORDER BY padj
    """, "KEGG pathways")

    # ── 5. Sample overview (QC + metadata) ───────────────────────────────
    sample_overview = query(conn, """
        SELECT
            s.sample_id,
            s.title,
            s.condition,
            s.organism,
            q.library_size,
            ROUND(q.zero_fraction * 100, 2) AS zero_pct,
            ROUND(q.median_count, 1)         AS median_count,
            q.pass_qc
        FROM samples s
        LEFT JOIN qc_metrics q ON s.sample_id = q.sample_id
        ORDER BY s.condition, q.library_size DESC
    """, "Sample overview")

    # ── 6. DE gene counts summary ─────────────────────────────────────────
    de_summary = query(conn, """
        SELECT
            CASE
                WHEN log2FoldChange >  1 AND padj < 0.05 THEN 'Upregulated'
                WHEN log2FoldChange < -1 AND padj < 0.05 THEN 'Downregulated'
                ELSE 'Not significant'
            END AS category,
            COUNT(*) AS gene_count
        FROM de_results
        WHERE log2FoldChange IS NOT NULL AND padj IS NOT NULL
        GROUP BY category
        ORDER BY gene_count DESC
    """, "DE summary counts")

    # ── Export all CSVs ───────────────────────────────────────────────────
    print("\nExporting CSVs for dashboard...")
    exports = {
        "top_de_genes.csv":      top_de,
        "volcano_data.csv":      volcano,
        "go_terms.csv":          go_terms,
        "kegg_pathways.csv":     kegg_paths,
        "sample_overview.csv":   sample_overview,
        "de_summary.csv":        de_summary,
    }

    for fname, df in exports.items():
        if df is not None and not df.empty:
            out = os.path.join(args.outdir, fname)
            df.to_csv(out, index=False)
            print(f"  Saved: {out}")

    # ── VST counts wide export ────────────────────────────────────────────
    if os.path.exists(args.vst):
        print("\nExporting VST counts for gene deep-dive page...")
        vst = pd.read_csv(args.vst, index_col=0)
        # Keep only top DE genes to keep the file manageable
        if not top_de.empty:
            top_genes = top_de["gene"].tolist()
            vst_top   = vst.loc[vst.index.isin(top_genes)]
        else:
            vst_top = vst.head(200)
        vst_top.reset_index(inplace=True)
        vst_top.rename(columns={vst_top.columns[0]: "gene"}, inplace=True)
        out = os.path.join(args.outdir, "gene_expression_wide.csv")
        vst_top.to_csv(out, index=False)
        print(f"  Saved: {out}  ({len(vst_top)} genes × {vst_top.shape[1]-1} samples)")

    conn.close()

    print("\n" + "=" * 55)
    print("All exports complete")
    print("=" * 55)
    print("\nDashboard CSVs ready in:", args.outdir)
    print("\nNext step:")
    print("  Open dashboard/rnaseq_dashboard.pbix in Power BI")
    print("  Open dashboard/rnaseq_dashboard.twbx in Tableau")


if __name__ == "__main__":
    main()
