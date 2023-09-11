# SCM-seq 
SCM-seq (Single-Cell and -Molecule sequencing) is an integrated experimental platform that enables single-cell multiomic analysis by combining high-throughput scRNA-seq to Oxford Nanopore Technologies (ONT) single-molecule sequencing.

In particular, SCM-seq allows the exploration and linking of distinct biological layers characterizing each single cell:

- ***Expressed somatic variants.*** By enriching for known 
we can 
 Established driver genes
Prognostic impact
Pharmacological targets
Mutation burden. Clonal evolution not supported at the moment
- ***Gene expression.*** 
It allows 
Tumor and immune subsets
Differentiation aberrancies
Functional subpopulations and gene modules.
- ***Splicing isoforms.***
Transcriptional regulation    and redundancy
Cancer-specific aberrancies
Proteome diversity (biomarkers? targets?)

we can reconstruct intratumoral heterogeneity

This repo collects scripts and pipelines for pre-processing and downstream analyses 
within the SRSF2 AML project.
We have applied this methodology to the analysis of three AML samples sharing a mutation in a spliceosome factor, with the aim to investigate how phenotypic heterogeneity is related to genetic complexity in both the malignant and immune compartments of a coherent AML subgroup.

For the time being, 
1. scRNA-seq
2. ONT basecalling and pre-processing
3. Cell barcode matching
4. Transcript alignment and annotation
5. Mutation analysis
6. Genotype analysis
7. Transcript analysis
