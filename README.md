# SCM-seq  
**SCM-seq** (**S**ingle-**C**ell and -**M**olecule sequencing) is an integrated platform that enables single-cell multiomic analysis by combining high-throughput scRNA-seq to Oxford Nanopore Technologies (ONT) single-molecule sequencing.

In particular, SCM-seq allows to explore and link distinct sources of intratumor heterogeneity on each individual cell:

- **Expressed somatic variants** 
- **Gene expression** 
- **Splicing isoforms**

This repo collects state-of-the-art scripts and pipelines for pre-processing and downstream analyses developed within the SRSF2-AML project.

For the time being, the workflow includes the following steps (order is mandatory):

1. **scRNA-seq**
2. **Cell barcode (CB) matching**
3. **Long-read alignment and transcript annotation**
4. **Target enrichment quality control**
5. **Mutation analysis and genotype assignment**
6. **Genotype analysis**
7. **Isoform analysis**

Steps 2-5 are enclosed into a dedicated [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html).
