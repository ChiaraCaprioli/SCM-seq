#!/bin/sh
#PBS -q workq
#PBS -N TRUST4
#PBS -l select=1:ncpus=4:mem=32g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/TRUST4_samples/saml1_promethion
#PBS -e /hpcnfs/scratch/PGP/SCMseq/TRUST4_samples/saml1_promethion
#PBS -M iman.nazari@ieo.it
#PBS -m bae

cd /hpcnfs/scratch/PGP/SCMseq/TRUST4/

run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -u /hpcnfs/scratch/PGP/SCMseq/sAML1_promethion/merged_non_enriched/bc_matching/merged_fastq/output_merged_samples.fastq.gz --barcode ../TRUST4_samples/saml1_promethion/barcodes.fastq --barcodeRange 0 15 + --od ../TRUST4_samples/saml1_promethion