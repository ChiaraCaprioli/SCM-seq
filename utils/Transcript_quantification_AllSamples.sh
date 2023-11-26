#!/bin/bash
#PBS -q workq
#PBS -N salmon
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/salmon_quantification
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/salmon_quantification
#PBS -m ae

source activate salmon
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/salmon_quantification/result

# salmon index -t /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/GTF_FASTA/gffcomare_AML_12345_hBM_123.combined.fa -i transcriptome_index

AML2=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML2_ErlyLateHSC.subset.fq.gz
AML3=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML3_ErlyLateHSC.subset.fq.gz
AML4_mut=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML4_mut.subset.fq.gz
AML4_wt=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML4_wt.subset.fq.gz
AML5_mut=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML5_mut.subset.fq.gz
AML5_wt=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/AML5_wt.subset.fq.gz
hBM1=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/hBM1_ErlyLateHSC.subset.fq.gz
hBM2=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/hBM2_ErlyLateHSC.subset.fq.gz
hBM3=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/hBM3_ErlyLateHSC.subset.fq.gz
sAML1_mut=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/sAML1_mut.subset.fq.gz
sAML1_wt=/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/fastq_split/fastq_spli_allTR/sAML1_wt.subset.fq.gz


salmon quant -i transcriptome_index -l A -r $AML2 --validateMappings -o AML2
salmon quant -i transcriptome_index -l A -r $AML3 --validateMappings -o AML3
salmon quant -i transcriptome_index -l A -r $AML4_mut --validateMappings -o AML4_mut
salmon quant -i transcriptome_index -l A -r $AML4_wt --validateMappings -o AML4_wt
salmon quant -i transcriptome_index -l A -r $AML5_mut --validateMappings -o AML5_mut
salmon quant -i transcriptome_index -l A -r $AML5_wt --validateMappings -o AML5_wt
salmon quant -i transcriptome_index -l A -r $hBM1 --validateMappings -o hBM1
salmon quant -i transcriptome_index -l A -r $hBM2 --validateMappings -o hBM2
salmon quant -i transcriptome_index -l A -r $hBM3 --validateMappings -o hBM3
salmon quant -i transcriptome_index -l A -r $sAML1_mut --validateMappings -o sAML1_mut
salmon quant -i transcriptome_index -l A -r $sAML1_wt --validateMappings -o sAML1_wt
