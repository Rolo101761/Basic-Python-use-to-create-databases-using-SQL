# Bioinformatics Portfolio — Python, R, SQL, Power BI / Tableau

MSc Bioinformatics student portfolio. Projects span Python pipelines,
statistical analysis in R, relational databases (SQLite), and interactive dashboards.

---

## Projects

### 1. Bioinformatics SQLite Database
**`portfolio_projects/bioinformatics_sql_db/`**

Relational database for genomic variants, gene expression, and biological pathways.

- Normalised schema: genes, samples, variants (VCF-style), pathways, expression
- Complex SQL queries: multi-table JOINs, aggregations, pathway–variant overlap
- VCF importer, ClinVar significance tracking, RNA-seq integration
- **Skills:** Python, SQLite, relational database design, genomics

```bash
python portfolio_projects/bioinformatics_sql_db/bioinformatics_db.py --demo
```

---

### 2. Enzyme Digestion Pipeline
**`portfolio_projects/enzyme_digestion_pipeline/`**

Pure-Python in-silico proteomics pipeline (no external dependencies).

- 6-frame ORF detection from bacterial genome FASTA
- Enzyme digestion: Trypsin, Lys-C, Arg-C, Glu-C with missed cleavage support
- Monoisotopic/average mass calculation, m/z values (+1 charge)
- Ion statistics: identifies the best protease for unambiguous protein identification
- **Skills:** Python, mass spectrometry, proteomics, algorithm design

```bash
python portfolio_projects/enzyme_digestion_pipeline/enzyme_digestion_pipeline.py \
    --fasta genome.fasta --enzyme all --mz_min 1000 --mz_max 1500
```

---

### 3. miRNA QC Pipeline
**`portfolio_projects/miRNA_QC_pipeline/`**

Quality control and analysis tool for miRNA–gene–disease association datasets
(compatible with miRTarBase, HMDD, and similar TSV-format databases).

- Confidence-score filtering, incomplete entry removal
- Summary statistics: unique miRNAs, genes, diseases, confidence distribution
- Top organisms and top target genes by association count
- Cleaned TSV output ready for downstream analysis
- **Skills:** Python, data QC, non-coding RNA biology, CSV/TSV processing

```bash
python portfolio_projects/miRNA_QC_pipeline/miRNA_QC_pipeline.py \
    --input associations.tsv --confidence 0.7
```

---

### 4. Rules-Based Epigenome Simulation
**`portfolio_projects/rules_based_simulation/`**

Monte Carlo simulation of histone methylation spreading using a writer/eraser model.

- 1D chromatin fibre model with periodic boundary conditions
- Cooperative spreading rules (writer probability scales with neighbour methylation state)
- Constitutive erasure, autonomous writer noise
- Parameter sweep to explore phase space (hypermethylated vs euchromatin-like regimes)
- HWHM-based spreading length quantification
- **Skills:** Python, epigenomics, Monte Carlo simulation, stochastic modelling

```bash
python portfolio_projects/rules_based_simulation/rules_based_simulation.py \
    --nucleosomes 200 --steps 50000 --p_write 0.08 --p_erase 0.02

# Parameter sweep across p_write x p_erase space
python portfolio_projects/rules_based_simulation/rules_based_simulation.py --sweep
```

---

### 5. RNA-seq Differential Expression Pipeline (Flagship — in progress)
**`portfolio_projects/flagship_rnaseq_project/`**

End-to-end RNA-seq analysis using Python + R + SQL + Power BI / Tableau.

- Python: GEO data download, library-size QC, data cleaning
- R (DESeq2): normalisation, differential expression, PCA, volcano plots, heatmaps
- R (clusterProfiler): GO and KEGG pathway enrichment analysis
- SQLite: results database for variants, expression, and pathway data
- Power BI / Tableau: interactive 4-page dashboard for stakeholder reporting
- **Skills:** Python, R, DESeq2, clusterProfiler, SQL, Power BI, bioinformatics end-to-end

See [`FLAGSHIP_PROJECT_PLAN.md`](portfolio_projects/flagship_rnaseq_project/FLAGSHIP_PROJECT_PLAN.md)
for the full architecture, R code snippets, and dashboard design.

---

## Skills Overview

| Skill | Projects |
|-------|---------|
| Python | All projects |
| SQL / SQLite | Bioinformatics DB, RNA-seq pipeline |
| R (DESeq2, clusterProfiler, ggplot2) | Flagship RNA-seq |
| Power BI / Tableau | Flagship RNA-seq |
| Mass spectrometry / proteomics | Enzyme Digestion Pipeline |
| Epigenomics / stochastic modelling | Rules-Based Simulation |
| Non-coding RNA biology | miRNA QC Pipeline |
| Genomic variant analysis | Bioinformatics DB |

---

## Running the projects

All Python projects require **Python 3.8+** and **no external dependencies** unless stated.

```bash
# Clone and explore
git clone https://github.com/Rolo101761/Basic-Python-use-to-create-databases-using-SQL
cd Basic-Python-use-to-create-databases-using-SQL

# Run the database demo
python portfolio_projects/bioinformatics_sql_db/bioinformatics_db.py --demo

# Run the epigenome simulation
python portfolio_projects/rules_based_simulation/rules_based_simulation.py --nucleosomes 100 --steps 10000
```
