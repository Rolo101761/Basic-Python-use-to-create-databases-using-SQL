# Bioinformatics SQLite Database

A relational database for storing, querying, and analysing genomic variants,
gene expression data, and biological pathways — built with pure Python and SQLite.

## Schema

```
genes ─────────────────┐
  │                    │
  │ gene_pathways      │ expression
  │                    │
pathways            samples
                       │
                    variants
```

| Table | Contents |
|-------|----------|
| `genes` | Gene annotations: symbol, chromosome, coordinates, biotype |
| `samples` | Sample metadata: condition, organism, tissue, age |
| `variants` | Genomic variants (VCF-style): SNVs/INDELs with ClinVar significance |
| `pathways` | KEGG/Reactome pathway definitions |
| `gene_pathways` | Many-to-many: gene ↔ pathway membership |
| `expression` | RNA-seq values: raw counts, TPM, log2FC, adjusted p-value |

## Usage

```bash
# Load demo data and run example queries
python bioinformatics_db.py --demo

# Use a persistent database file
python bioinformatics_db.py --db my_project.db --demo

# Import variants from a VCF file
python bioinformatics_db.py --db my_project.db --vcf sample.vcf --sample PATIENT_001
```

## Example Queries Included

1. **Pathogenic variants in tumour samples** — joins variants + genes + samples
2. **Significantly differentially expressed genes** — filters by adj. p-value and log2FC
3. **Pathway-variant overlap** — genes in DNA Damage Response pathway carrying variants
4. **Variant count per gene** — summary statistics with pathogenic breakdown
5. **Sample variant burden** — total variants + mean allele frequency per sample

## Demo output

```
[DB] Schema initialised.
[DB] Demo data loaded.

========================================
Database Summary
========================================
  genes                    8 rows
  samples                  6 rows
  variants                 6 rows
  pathways                 5 rows
  gene_pathways           10 rows
  expression              10 rows
========================================

[Query 1] Pathogenic variants in tumour samples:
  symbol                chromosome            position    ...
  ----------------------------------------------------------
  BRCA1                 17                    43071077    ...
  EGFR                  7                     55174772    ...
  KRAS                  12                    25245350    ...
  TP53                  17                    7577022     ...
```

## Skills demonstrated

- Relational database design (normalisation, foreign keys, indexes)
- Complex SQL queries (multi-table JOINs, aggregations, conditional counts)
- Python + SQLite3 (no ORM — raw SQL for clarity and performance)
- Genomics domain knowledge: VCF format, ClinVar, RNA-seq data structures
- Biological pathway integration (KEGG/Reactome)
- Command-line interface with argparse
