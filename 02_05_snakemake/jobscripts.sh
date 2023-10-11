#!/bin/sh
#PBS -q workq
#PBS -N snakemake_flames_sAML1_new_panel
#PBS -l select=1:ncpus=16:mem=32g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/saml1_promethion
#PBS -e /hpcnfs/scratch/PGP/SCMseq/saml1_promethion
#PBS -M iman.nazari@ieo.it
#PBS -m bae

PBS_O_WORKDIR="/hpcnfs/scratch/PGP/SCMseq/saml1_promethion"

mkdir -p $PBS_O_WORKDIR/merged_enriched
mkdir -p $PBS_O_WORKDIR/merged_non_enriched



LD_PRELOAD=/hpcnfs/techunits/bioinformatics/software/petagene/petalink_1.3.15/bin/petalink.so
cd /hpcnfs/scratch/PGP/niman/snakemake_flames/saml1_promethion
source activate /hpcnfs/home/ieo5268/.conda/envs/snakemake_env
snakemake --resources mem_gb=16 -j 1 --latency-wait 10 --use-conda --use-singularity --singularity-args "--bind /hpcnfs:/hpcnfs" --unlock
snakemake --resources mem_gb=32 -j 1 --latency-wait 10 --use-conda --use-singularity --singularity-args "--bind /hpcnfs:/hpcnfs" 
snakemake --dag | dot -Tsvg > $PBS_O_WORKDIR/out.svg
snakemake --report $PBS_O_WORKDIR/report.html
