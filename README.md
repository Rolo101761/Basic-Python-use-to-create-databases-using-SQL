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

### 2. Enzyme Digestion Pipeline *(code coming soon)*
See: [Enzyme-Digestion-Step-by-Step-Pipeline](https://github.com/Rolo101761/Enzyme-Digestion-Step-by-Step-Pipeline)

---

### 3. miRNA QC Pipeline *(code coming soon)*
See: [miRNA-QC-Pipeline](https://github.com/Rolo101761/miRNA-QC-Pipeline)

---

### 4. Rules-Based Epigenome Simulation
**`portfolio_projects/rules_based_simulation/`**

Monte Carlo simulation of histone methylation spreading using a writer/eraser model.

- 1D chromatin fibre model with periodic boundary conditions
- Cooperative spreading rules (writer probability scales with neighbour methylation state)
- Constitutive erasure, autonomous writer noise
- Parameter sweep to explore phase space (hypermethylated vs euchromatin-like regimes)
- HWHM-based spreading length quantification
- ggplot2 visualisations: methylation profile, timeseries, parameter sweep heatmap
- **Skills:** R, ggplot2, epigenomics, Monte Carlo simulation, stochastic modelling

```bash
Rscript portfolio_projects/rules_based_simulation/rules_based_simulation.R

# Parameter sweep across p_write x p_erase space
Rscript portfolio_projects/rules_based_simulation/rules_based_simulation.R --sweep
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
| R (ggplot2, DESeq2, clusterProfiler) | Rules-Based Simulation, Flagship RNA-seq |
| Power BI / Tableau | Flagship RNA-seq |
| Mass spectrometry / proteomics | Enzyme Digestion Pipeline |
| Epigenomics / stochastic modelling | Rules-Based Simulation |
| Non-coding RNA biology | miRNA QC Pipeline |
| Genomic variant analysis | Bioinformatics DB |

---

## Running the projects

Python projects require **Python 3.8+**. R scripts require **R 4.0+** with `ggplot2` and `reshape2`.

```bash
# Clone and explore
git clone https://github.com/Rolo101761/Basic-Python-use-to-create-databases-using-SQL
cd Basic-Python-use-to-create-databases-using-SQL

# Run the database demo
python portfolio_projects/bioinformatics_sql_db/bioinformatics_db.py --demo

# Run the epigenome simulation (R)
Rscript portfolio_projects/rules_based_simulation/rules_based_simulation.R
```
