## Overview

```05_genotype_analysis``` follows Iman's snakemake pipeline for ONT reads pre-processing, mutations mapping, CB demultiplexing and genotype confidence assignment, and extends genotype analysis to the following aims:

- Filter genotyped variants and cells
- Tiering
- Assignment to scRNA metadata object
- Statistics

## Input
- *df_bc_all.csv*: result of CB demultiplexing after mutation analysis. It consists of a table of cells (=rows) and read counts for each genotyping class for each variant (=columns). For the time being, we refer to genes and variants interchangeably, but changing field names to "SYMBOL:CHR:START-END" may be warranted.
- *tg_wes_intersect.csv*: a table with the list of target variants
- ```/genotype_imputation_diff_appr_5%_err```: result of genotyping confidence assignment. This folder contains, for each variant, one file for mutated cells and one file for wt cells.

## Structure
Currently, ```05_genotype_analysis``` relies on the following folder structure:
```
genotype_analysis/
├── all_samples
│   ├── plots
│   └── tables
└── sample1
    ├── plots
    └── tables
```
