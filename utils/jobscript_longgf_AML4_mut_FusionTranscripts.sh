#!/bin/sh
#PBS -q workq
#PBS -N longgf_AML4
#PBS -l select=1:ncpus=16:mem=16g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF
#PBS -M iman.nazari@ieo.it
#PBS -m bae
PBS_O_WORKDIR="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF"
# source activate FLAMES
cd /hpcnfs/scratch/PGP/niman/Chiara/tools/LongGF/LongGF/bin
echo "The working directory is: 'PBS_O_WORKDIR'"
samtools sort -n /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/Realignment_bam_by_genotype_merged_gtf/AML4_mut_sorted.bam -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/Realignment_bam_by_genotype_merged_gtf/AML4_mut_sorted_sorted.bam
# source deactivate
LongGF /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/Realignment_bam_by_genotype_merged_gtf/AML4_mut_sorted_sorted.bam /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf 100 50 100 0 3 > /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF/longgf_AML4_mut_subsampled_srsf2.log
grep "SumGF" /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF/longgf_AML4_mut_subsampled_srsf2.log > /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/longGF/summary_AML4_mut.txt
