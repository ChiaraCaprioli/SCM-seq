#!/bin/bash
#PBS -q workq
#PBS -N fastq_split
#PBS -l select=2:ncpus=16:mem=32gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split
#PBS -m ae

# export PetaLinkMode="md5match"
# export PETASUITE_REFPATH=/hpcnfs/techunits/bioinformatics/software/petagene/petalink_1.3.15/species
# export LD_PRELOAD=/hpcnfs/techunits/bioinformatics/software/petagene/petalink_1.3.15/bin/petalink.so

source activate cpat
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split

seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xmut_AML4_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/AML4_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > AML4_mut.subset.fq.gz
seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xwt_AML4_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/AML4_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > AML4_wt.subset.fq.gz

seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xmut_AML5_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/AML5_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > AML5_mut.subset.fq.gz
seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xwt_AML5_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/AML5_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > AML5_wt.subset.fq.gz

seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xmut_sAML1_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/sAML1_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > sAML1_mut.subset.fq.gz
seqkit grep -r -f /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq/cell_ID_xwt_sAML1_HSC_genes.txt /hpcnfs/scratch/PGP/SCMseq/sAML1_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz > sAML1_wt.subset.fq.gz




# singularity exec -B /hpcnfs /hpcnfs/data/PGP/niman/images/SCM_rsudio.sif Rscript fastq_r.R
exit


