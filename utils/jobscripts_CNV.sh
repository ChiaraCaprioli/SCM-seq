#!/bin/sh
#PBS -q workq
#PBS -N CNVkit
#PBS -l select=1:ncpus=4:mem=32g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1
#PBS -e /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1
#PBS -M iman.nazari@ieo.it
#PBS -m bae

source activate cpat
PBS_O_WORKDIR="/hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1"
cnvkit.py access /hpcnfs/data/PGP/reference_genomes/UCSC/hg38/Indexes/bwa_0.7.8/hg38.fa -o /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1/access.hg38.bed
cnvkit.py batch /hpcnfs/data/PGP/lfiorenza/WES_AML_Caprioli/sAML1_D/realignment/sAML1_D_MarkDup_recal.bam --normal /hpcnfs/data/PGP/lfiorenza/WES_AML_Caprioli/sAML1_N/realignment/sAML1_N_MarkDup_recal.bam \
    --targets /hpcnfs/data/PGP/exome/referenceBed/hg19/ss_v7/hg38/ss_v7_regions_hg38.bed --annotate /hpcnfs/data/PGP/exome/referenceBed/hg19/ss_v7/hg38/S31285117_Covered.bed  \
   --short-names --drop-low-coverage  --fasta /hpcnfs/data/PGP/reference_genomes/UCSC/hg38/Indexes/bwa_0.7.8/hg38.fa --access /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1/access.hg38.bed \
    --output-reference /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1/my_reference.cnn --output-dir /hpcnfs/scratch/PGP/SCMseq/CNVkit/sAML1/ \
    --diagram --scatter

