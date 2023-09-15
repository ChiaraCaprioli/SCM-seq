Currently, this part follows Iman's snakemake pipeline for ONT reads pre-processing, mutations mapping, CB demultiplexing and genotype confidence assignment.  

proposed structure for genotype analysis folder:
```
genotype_analysis/
├── data
│   └── sample1
│       ├── df_bc_all.csv
│       ├── genotype_imputation_diff_appr_5%_err
│       └── tg_wes_intersect.csv
└── results
    └── sample1
        ├── plots
        └── tables
```
