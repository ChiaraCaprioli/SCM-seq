#!/bin/bash
#PBS -q workq
#PBS -N bam_split_salmon
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/bam_split
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/bam_split
#PBS -m ae

source activate FLAMES1
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split
python bam_split.py

samtools merge -f wt_merged.bam bam_wt_*.bam
samtools merge -f mut_merged.bam bam_mut_*.bam
samtools index wt_merged.bam
samtools index mut_merged.bam

source deactivate 

source activate salmon
mkdir /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/salmon_transcript_bam
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/salmon_transcript_bam

tr_as="/hpcnfs/scratch/PGP/SCMseq/sAML1_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xmut_sAML1_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xwt_sAML1_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut -o mut_sAML1
salmon quant -t $tr_as -l A -a $wt  -o wt_sAML1

tr_as="/hpcnfs/scratch/PGP/SCMseq/AML4_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xmut_AML4_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xwt_AML4_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut  -o mut_AML4
salmon quant -t $tr_as -l A -a $wt  -o wt_AML4

tr_as="/hpcnfs/scratch/PGP/SCMseq/AML5_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xmut_AML5_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/bam_xwt_AML5_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut  -o mut_AML5
salmon quant -t $tr_as -l A -a $wt  -o wt_AML5

