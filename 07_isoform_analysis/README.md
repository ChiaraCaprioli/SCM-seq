## Splicing Analysis
### Based on PSI-Sigma tool on the samples with SRSF2p95 comon mutation

In our research, we applied the FLAMES technique to match barcodes and quantify isoforms. Subsequently, we employed the iso-count matrix, derived from single-cell data analysis via 10x, to filter cells based on their specific lineage in the leukemia compartment, as well as their genotype determined through genotype imputation. Our method not only enables standard splicing analysis of ONT sequenced samples but also facilitates a comparative evaluation between the wild-type (WT) and mutant (MUT) cellular groups. Leveraging the SCM-seq strategy for single-cell mutational analysis, we delved into a comprehensive assessment.

We then proceeded to segregate the FASTQ files post barcode matching, followed by their realignment and mapping to the reference genome using the consolidated GTF file, which combines the gene and transcript annotations from all samples. Treating each resulting bam file as a pseudo-bulk sample offered several advantages. This approach not only amplified signals, improving the detection of subtle splicing events, but also aided in the identification of rare isoforms that might otherwise be obscured by individual cell variability. Our decision to merge the GTF files was based on the need to address discrepancies arising from the interpretation of transcript annotations in long read sequencing, as well as discrepancies at the intron-exon junctions during nanopore sequencing. This meticulous process ensures a more accurate depiction of the transcriptome, incorporating the internal exon structure and accounting for minor variations at intron-exon junctions within the nanopore sequencing data, ultimately enhancing the reliability of our findings.

To comprehensively investigate the interplay between mutations, phenotype, and splicing patterns in acute myeloid leukemia (AML), we meticulously conducted distinct splicing analyses on individual samples. This approach enabled a detailed examination of the unique splicing heterogeneity within each sample, crucial for understanding the underlying molecular mechanisms, especially in cases where the SRSF2p95 mutation appears across various AML samples. 

To enhance the statistical significance of identified splicing events, we devised a novel strategy involving the creation of synthetic internal replicates. This strategy involved the random division of mutated (MUT) cells into three non-redundant groups and wild-type (WT) cells into two groups due to their lower count. The splicing analysis was then performed between the MUT and WT groups, repeating the process 20 times to ensure robustness and consistency in the results. Events occurring frequently across iterations were selected for further analysis, while those occurring only once were excluded.

Utilizing the expectation-maximization approach (EM), we successfully identified key parameters and peaks in the bimodal distribution, helping us discern significant biological events from random fluctuations. Moreover, the Benjamini-Hochberg (BH) formula was employed to calculate the false discovery rate (FDR) for each run, and the Fisher formula was utilized to consolidate the FDRs of all runs, ensuring a comprehensive evaluation of the statistical significance. Additionally, the differential Percent Spliced-In (dPSI) values were averaged across the events from the 20 runs to provide a comprehensive understanding of the splicing dynamics.


This analysis is based on the PSI-Sigma tool's output.

- Different parts added here as it is listed
  - Splicing analysis of the all samples
  - Splicing analysis of All samples based on the sub-sampling method
  - pathway analysis
  - transcript expression analysis
  - GFF to GTF pipeline



