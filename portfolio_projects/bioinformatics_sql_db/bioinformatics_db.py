"""
Bioinformatics SQLite Database
================================
Relational database for storing and querying genomic variants,
gene annotations, and sample metadata.

Schema:
  genes        — gene annotations (symbol, chromosome, coordinates, biotype)
  samples      — biological samples (ID, condition, organism, tissue)
  variants     — genomic variants linked to genes and samples (VCF-style)
  pathways     — biological pathways (KEGG/Reactome)
  gene_pathways — many-to-many: genes <-> pathways
  expression   — gene expression values per sample (RNA-seq counts)

Usage:
    python bioinformatics_db.py --demo          # load demo data and run example queries
    python bioinformatics_db.py --db my.db      # use a custom database file
    python bioinformatics_db.py --vcf vars.vcf  # import a VCF file
"""

import argparse
import sqlite3
import csv
import os
import sys
from datetime import datetime


# ---------------------------------------------------------------------------
# Database setup
# ---------------------------------------------------------------------------
SCHEMA_SQL = """
-- Gene annotations
CREATE TABLE IF NOT EXISTS genes (
    gene_id     INTEGER PRIMARY KEY AUTOINCREMENT,
    symbol      TEXT    NOT NULL,
    chromosome  TEXT    NOT NULL,
    start_pos   INTEGER NOT NULL,
    end_pos     INTEGER NOT NULL,
    strand      TEXT    CHECK(strand IN ('+', '-', '.')),
    biotype     TEXT,
    description TEXT,
    UNIQUE(symbol, chromosome)
);

-- Biological samples
CREATE TABLE IF NOT EXISTS samples (
    sample_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_name TEXT    NOT NULL UNIQUE,
    condition   TEXT    NOT NULL,   -- e.g. 'tumour', 'normal', 'treated'
    organism    TEXT    NOT NULL,
    tissue      TEXT,
    sex         TEXT    CHECK(sex IN ('M', 'F', 'unknown', NULL)),
    age         INTEGER,
    created_at  TEXT    DEFAULT (datetime('now'))
);

-- Genomic variants (VCF-style)
CREATE TABLE IF NOT EXISTS variants (
    variant_id  INTEGER PRIMARY KEY AUTOINCREMENT,
    chromosome  TEXT    NOT NULL,
    position    INTEGER NOT NULL,
    ref_allele  TEXT    NOT NULL,
    alt_allele  TEXT    NOT NULL,
    variant_type TEXT,              -- SNV, INDEL, MNV
    effect      TEXT,               -- synonymous, missense, nonsense, splice_site
    gene_id     INTEGER REFERENCES genes(gene_id),
    sample_id   INTEGER REFERENCES samples(sample_id),
    qual_score  REAL,
    depth       INTEGER,
    allele_freq REAL,
    clinvar_sig TEXT                -- benign, pathogenic, VUS, etc.
);

-- Biological pathways
CREATE TABLE IF NOT EXISTS pathways (
    pathway_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    pathway_name TEXT    NOT NULL,
    database     TEXT    NOT NULL,  -- KEGG, Reactome, GO
    pathway_code TEXT    UNIQUE
);

-- Gene–pathway many-to-many
CREATE TABLE IF NOT EXISTS gene_pathways (
    gene_id    INTEGER NOT NULL REFERENCES genes(gene_id),
    pathway_id INTEGER NOT NULL REFERENCES pathways(pathway_id),
    PRIMARY KEY (gene_id, pathway_id)
);

-- Gene expression (RNA-seq raw counts or normalised values)
CREATE TABLE IF NOT EXISTS expression (
    expr_id     INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id     INTEGER NOT NULL REFERENCES genes(gene_id),
    sample_id   INTEGER NOT NULL REFERENCES samples(sample_id),
    raw_count   INTEGER,
    tpm         REAL,
    log2_fc     REAL,               -- log2 fold-change vs reference
    adj_pvalue  REAL,               -- adjusted p-value (DESeq2/edgeR)
    UNIQUE(gene_id, sample_id)
);

-- Indexes for common query patterns
CREATE INDEX IF NOT EXISTS idx_variants_gene    ON variants(gene_id);
CREATE INDEX IF NOT EXISTS idx_variants_sample  ON variants(sample_id);
CREATE INDEX IF NOT EXISTS idx_variants_chrom   ON variants(chromosome, position);
CREATE INDEX IF NOT EXISTS idx_expression_gene  ON expression(gene_id);
CREATE INDEX IF NOT EXISTS idx_expression_padj  ON expression(adj_pvalue);
"""


def connect(db_path=':memory:'):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def initialise_schema(conn):
    conn.executescript(SCHEMA_SQL)
    conn.commit()
    print(f"[DB] Schema initialised.")


# ---------------------------------------------------------------------------
# Demo data loader
# ---------------------------------------------------------------------------
DEMO_GENES = [
    ('BRCA1', '17', 43044295, 43125483, '-', 'protein_coding', 'Breast cancer type 1 susceptibility'),
    ('BRCA2', '13', 32315474, 32400266, '+', 'protein_coding', 'Breast cancer type 2 susceptibility'),
    ('TP53',  '17', 7565097,  7590856,  '-', 'protein_coding', 'Tumour protein p53'),
    ('KRAS',  '12', 25204789, 25250936, '-', 'protein_coding', 'KRAS proto-oncogene'),
    ('EGFR',  '7',  55019017, 55207337, '+', 'protein_coding', 'Epidermal growth factor receptor'),
    ('MIR21', '17', 57918627, 57918698, '+', 'miRNA',           'MicroRNA 21'),
    ('DNMT3A','2',  25234374, 25342811, '+', 'protein_coding', 'DNA methyltransferase 3A'),
    ('EZH2',  '7',  148504464,148581441,'-', 'protein_coding', 'Enhancer of zeste homolog 2'),
]

DEMO_SAMPLES = [
    ('SAMPLE_001', 'tumour',  'Homo sapiens', 'breast',  'F', 52),
    ('SAMPLE_002', 'normal',  'Homo sapiens', 'breast',  'F', 52),
    ('SAMPLE_003', 'tumour',  'Homo sapiens', 'lung',    'M', 67),
    ('SAMPLE_004', 'normal',  'Homo sapiens', 'lung',    'M', 67),
    ('SAMPLE_005', 'treated', 'Homo sapiens', 'blood',   'F', 44),
    ('SAMPLE_006', 'control', 'Homo sapiens', 'blood',   'F', 44),
]

DEMO_PATHWAYS = [
    ('DNA Damage Response',    'KEGG',     'hsa03430'),
    ('p53 Signalling Pathway', 'KEGG',     'hsa04115'),
    ('MAPK Signalling',        'KEGG',     'hsa04010'),
    ('Homologous Recombination','KEGG',    'hsa03440'),
    ('Epigenetic Regulation',  'Reactome', 'R-HSA-212165'),
]

DEMO_VARIANTS = [
    # (chrom, pos, ref, alt, type, effect, gene_symbol, sample_name, qual, depth, af, clinvar)
    ('17', 43071077, 'A', 'T', 'SNV', 'missense', 'BRCA1', 'SAMPLE_001', 60.0, 45, 0.48, 'pathogenic'),
    ('13', 32340300, 'C', 'T', 'SNV', 'missense', 'BRCA2', 'SAMPLE_001', 55.0, 38, 0.52, 'VUS'),
    ('17', 7577022,  'G', 'A', 'SNV', 'missense', 'TP53',  'SAMPLE_001', 80.0, 62, 0.45, 'pathogenic'),
    ('17', 7577022,  'G', 'A', 'SNV', 'missense', 'TP53',  'SAMPLE_003', 75.0, 58, 0.49, 'pathogenic'),
    ('12', 25245350, 'C', 'T', 'SNV', 'missense', 'KRAS',  'SAMPLE_003', 90.0, 71, 0.51, 'pathogenic'),
    ('7',  55174772, 'G', 'T', 'SNV', 'missense', 'EGFR',  'SAMPLE_003', 85.0, 55, 0.47, 'pathogenic'),
]

DEMO_EXPRESSION = [
    # (gene_symbol, sample_name, raw_count, tpm, log2_fc, adj_pvalue)
    ('BRCA1', 'SAMPLE_001', 120,  8.4,  -2.3, 0.001),
    ('BRCA1', 'SAMPLE_002', 980,  68.1,  0.0, 1.000),
    ('TP53',  'SAMPLE_001', 2400, 167.2,  1.8, 0.003),
    ('TP53',  'SAMPLE_002', 680,  47.3,   0.0, 1.000),
    ('KRAS',  'SAMPLE_003', 3100, 215.6,  2.1, 0.0008),
    ('KRAS',  'SAMPLE_004', 720,  50.1,   0.0, 1.000),
    ('EGFR',  'SAMPLE_003', 4200, 292.0,  2.6, 0.0001),
    ('EGFR',  'SAMPLE_004', 580,  40.3,   0.0, 1.000),
    ('EZH2',  'SAMPLE_005', 1800, 125.2,  1.5, 0.015),
    ('EZH2',  'SAMPLE_006', 760,  52.8,   0.0, 1.000),
]


def load_demo_data(conn):
    cur = conn.cursor()

    # Genes
    cur.executemany(
        "INSERT OR IGNORE INTO genes (symbol, chromosome, start_pos, end_pos, strand, biotype, description) VALUES (?,?,?,?,?,?,?)",
        DEMO_GENES
    )

    # Samples
    cur.executemany(
        "INSERT OR IGNORE INTO samples (sample_name, condition, organism, tissue, sex, age) VALUES (?,?,?,?,?,?)",
        DEMO_SAMPLES
    )

    # Pathways
    cur.executemany(
        "INSERT OR IGNORE INTO pathways (pathway_name, database, pathway_code) VALUES (?,?,?)",
        DEMO_PATHWAYS
    )

    # Gene–pathway associations
    gene_pathway_links = [
        ('BRCA1', 'DNA Damage Response'), ('BRCA1', 'Homologous Recombination'),
        ('BRCA2', 'DNA Damage Response'), ('BRCA2', 'Homologous Recombination'),
        ('TP53',  'p53 Signalling Pathway'), ('TP53', 'DNA Damage Response'),
        ('KRAS',  'MAPK Signalling'),
        ('EGFR',  'MAPK Signalling'),
        ('EZH2',  'Epigenetic Regulation'),
        ('DNMT3A','Epigenetic Regulation'),
    ]
    for gene_sym, pathway_name in gene_pathway_links:
        gene_id    = cur.execute("SELECT gene_id FROM genes WHERE symbol=?", (gene_sym,)).fetchone()
        pathway_id = cur.execute("SELECT pathway_id FROM pathways WHERE pathway_name=?", (pathway_name,)).fetchone()
        if gene_id and pathway_id:
            cur.execute("INSERT OR IGNORE INTO gene_pathways VALUES (?,?)", (gene_id[0], pathway_id[0]))

    # Variants
    for row in DEMO_VARIANTS:
        chrom, pos, ref, alt, vtype, effect, gene_sym, sample_name, qual, depth, af, clinvar = row
        gene_id   = cur.execute("SELECT gene_id   FROM genes   WHERE symbol=?",      (gene_sym,)).fetchone()
        sample_id = cur.execute("SELECT sample_id FROM samples WHERE sample_name=?", (sample_name,)).fetchone()
        if gene_id and sample_id:
            cur.execute(
                "INSERT INTO variants (chromosome, position, ref_allele, alt_allele, variant_type, effect, gene_id, sample_id, qual_score, depth, allele_freq, clinvar_sig) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (chrom, pos, ref, alt, vtype, effect, gene_id[0], sample_id[0], qual, depth, af, clinvar)
            )

    # Expression
    for row in DEMO_EXPRESSION:
        gene_sym, sample_name, raw, tpm, l2fc, padj = row
        gene_id   = cur.execute("SELECT gene_id   FROM genes   WHERE symbol=?",      (gene_sym,)).fetchone()
        sample_id = cur.execute("SELECT sample_id FROM samples WHERE sample_name=?", (sample_name,)).fetchone()
        if gene_id and sample_id:
            cur.execute(
                "INSERT OR IGNORE INTO expression (gene_id, sample_id, raw_count, tpm, log2_fc, adj_pvalue) VALUES (?,?,?,?,?,?)",
                (gene_id[0], sample_id[0], raw, tpm, l2fc, padj)
            )

    conn.commit()
    print(f"[DB] Demo data loaded.")


# ---------------------------------------------------------------------------
# Example queries
# ---------------------------------------------------------------------------
def run_example_queries(conn):
    print(f"\n{'='*60}")
    print("EXAMPLE QUERIES")
    print(f"{'='*60}")

    # Query 1: Pathogenic variants in tumour samples
    print("\n[Query 1] Pathogenic variants in tumour samples:")
    rows = conn.execute("""
        SELECT g.symbol, v.chromosome, v.position, v.ref_allele, v.alt_allele,
               v.effect, v.allele_freq, v.clinvar_sig, s.sample_name, s.tissue
        FROM   variants v
        JOIN   genes   g ON v.gene_id   = g.gene_id
        JOIN   samples s ON v.sample_id = s.sample_id
        WHERE  v.clinvar_sig = 'pathogenic'
          AND  s.condition   = 'tumour'
        ORDER  BY g.symbol, s.sample_name
    """).fetchall()
    _print_rows(rows)

    # Query 2: Differentially expressed genes (padj < 0.05)
    print("\n[Query 2] Significantly differentially expressed genes (adj p < 0.05):")
    rows = conn.execute("""
        SELECT g.symbol, g.biotype, s.condition, s.tissue,
               e.raw_count, e.tpm, e.log2_fc, e.adj_pvalue
        FROM   expression e
        JOIN   genes   g ON e.gene_id   = g.gene_id
        JOIN   samples s ON e.sample_id = s.sample_id
        WHERE  e.adj_pvalue < 0.05
          AND  ABS(e.log2_fc) > 1.5
        ORDER  BY e.adj_pvalue
    """).fetchall()
    _print_rows(rows)

    # Query 3: Genes in the DNA Damage Response pathway with variants
    print("\n[Query 3] Genes in DNA Damage Response pathway that carry variants:")
    rows = conn.execute("""
        SELECT DISTINCT g.symbol, p.pathway_name, v.effect, v.clinvar_sig
        FROM   gene_pathways gp
        JOIN   genes    g ON gp.gene_id    = g.gene_id
        JOIN   pathways p ON gp.pathway_id = p.pathway_id
        JOIN   variants v ON v.gene_id     = g.gene_id
        WHERE  p.pathway_name = 'DNA Damage Response'
        ORDER  BY g.symbol
    """).fetchall()
    _print_rows(rows)

    # Query 4: Summary variant counts per gene
    print("\n[Query 4] Variant count per gene:")
    rows = conn.execute("""
        SELECT g.symbol, COUNT(v.variant_id) AS n_variants,
               COUNT(CASE WHEN v.clinvar_sig='pathogenic' THEN 1 END) AS n_pathogenic
        FROM   genes g
        LEFT   JOIN variants v ON v.gene_id = g.gene_id
        GROUP  BY g.gene_id
        HAVING n_variants > 0
        ORDER  BY n_variants DESC
    """).fetchall()
    _print_rows(rows)

    # Query 5: Sample variant burden
    print("\n[Query 5] Variant burden per sample:")
    rows = conn.execute("""
        SELECT s.sample_name, s.condition, s.tissue,
               COUNT(v.variant_id)                                          AS total_variants,
               COUNT(CASE WHEN v.clinvar_sig='pathogenic' THEN 1 END)       AS pathogenic,
               ROUND(AVG(v.allele_freq), 3)                                 AS mean_AF
        FROM   samples s
        LEFT   JOIN variants v ON v.sample_id = s.sample_id
        GROUP  BY s.sample_id
        ORDER  BY total_variants DESC
    """).fetchall()
    _print_rows(rows)


def _print_rows(rows):
    if not rows:
        print("  (no results)")
        return
    keys = rows[0].keys()
    header = '  ' + '  '.join(f"{k:<20}" for k in keys)
    print(header)
    print('  ' + '-' * (len(header) - 2))
    for row in rows:
        print('  ' + '  '.join(f"{str(v):<20}" for v in row))


# ---------------------------------------------------------------------------
# VCF importer
# ---------------------------------------------------------------------------
def import_vcf(conn, vcf_path, sample_name, organism='Homo sapiens'):
    """Import variants from a standard VCF file into the database."""
    if not os.path.exists(vcf_path):
        print(f"[ERROR] VCF file not found: {vcf_path}", file=sys.stderr)
        return

    cur = conn.cursor()

    # Ensure sample exists
    cur.execute(
        "INSERT OR IGNORE INTO samples (sample_name, condition, organism) VALUES (?, 'unknown', ?)",
        (sample_name, organism)
    )
    sample_id = cur.execute("SELECT sample_id FROM samples WHERE sample_name=?", (sample_name,)).fetchone()[0]

    count = 0
    with open(vcf_path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            qual = float(parts[5]) if parts[5] not in ('.', '') else None
            cur.execute(
                "INSERT INTO variants (chromosome, position, ref_allele, alt_allele, sample_id, qual_score) VALUES (?,?,?,?,?,?)",
                (chrom, pos, ref, alt, sample_id, qual)
            )
            count += 1

    conn.commit()
    print(f"[VCF] Imported {count} variants from {vcf_path} (sample: {sample_name})")


# ---------------------------------------------------------------------------
# Database summary
# ---------------------------------------------------------------------------
def print_db_summary(conn):
    tables = ['genes', 'samples', 'variants', 'pathways', 'gene_pathways', 'expression']
    print(f"\n{'='*40}")
    print("Database Summary")
    print(f"{'='*40}")
    for table in tables:
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table:<20} {n:>6} rows")
    print(f"{'='*40}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Bioinformatics SQLite Database — genomic variants, expression, pathways'
    )
    parser.add_argument('--db',     default='bioinformatics.db', help='SQLite database file')
    parser.add_argument('--demo',   action='store_true',         help='Load demo data and run example queries')
    parser.add_argument('--vcf',    default=None,                help='VCF file to import')
    parser.add_argument('--sample', default='SAMPLE_VCF',        help='Sample name for VCF import')
    args = parser.parse_args()

    conn = connect(args.db)
    initialise_schema(conn)

    if args.demo:
        load_demo_data(conn)
        print_db_summary(conn)
        run_example_queries(conn)

    if args.vcf:
        import_vcf(conn, args.vcf, args.sample)
        print_db_summary(conn)

    if not args.demo and not args.vcf:
        print_db_summary(conn)
        print("\nTip: run with --demo to load example data and see queries.")

    conn.close()


if __name__ == '__main__':
    main()
