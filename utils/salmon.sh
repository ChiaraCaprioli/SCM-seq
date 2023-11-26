#!/bin/bash
#PBS -q workq
#PBS -N salmon
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/salmon_allTR_split
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/salmon_allTR_split
#PBS -m ae

source activate salmon
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/salmon_allTR_split

tr_as="/hpcnfs/scratch/PGP/SCMseq/sAML1_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xmut_sAML1_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xwt_sAML1_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut -o mut_sAML1
salmon quant -t $tr_as -l A -a $wt  -o wt_sAML1

tr_as="/hpcnfs/scratch/PGP/SCMseq/AML4_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xmut_AML4_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xwt_AML4_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut  -o mut_AML4
salmon quant -t $tr_as -l A -a $wt  -o wt_AML4

tr_as="/hpcnfs/scratch/PGP/SCMseq/AML5_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa"
mut="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xmut_AML5_HSC_genes_sorted.bam"
wt="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_split_bam_allTR/bam_xwt_AML5_HSC_genes_sorted.bam"

salmon quant -t $tr_as -l A -a $mut  -o mut_AML5
salmon quant -t $tr_as -l A -a $wt  -o wt_AML5