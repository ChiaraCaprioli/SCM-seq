#!/bin/sh
#PBS -q workq
#PBS -N Flame_sc_long
#PBS -l select=1:ncpus=16:mem=32g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts
#PBS -e /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts
#PBS -M iman.nazari@ieo.it
#PBS -m bae
PBS_O_WORKDIR="/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts"
cd /hpcnfs/home/ieo5268/FLAMES/python
source activate FLAMES
echo "The working directory is: 'PBS_O_WORKDIR'"
mkdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM1/FLAMES_out
mkdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM2/FLAMES_out
mkdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM3/FLAMES_out

sc_long_pipeline.py \
--gff3  /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf \
--infq /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM1/barcode_matching_output/out.fastq.gz \
--outdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM1/FLAMES_out \
--genomefa /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/GRCh38.primary_assembly.genome.fa \
--minimap2_dir /hpcnfs/scratch/PGP/niman/Chiara/tools/minimap2/minimap2 \
--config_file /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts/config_sclr_nanopore_non_enriched.json \
--downsample_ratio 1 \


sc_long_pipeline.py \
--gff3  /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf \
--infq /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM2/barcode_matching_output/out.fastq.gz \
--outdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM2/FLAMES_out \
--genomefa /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/GRCh38.primary_assembly.genome.fa \
--minimap2_dir /hpcnfs/scratch/PGP/niman/Chiara/tools/minimap2/minimap2 \
--config_file /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts/config_sclr_nanopore_non_enriched.json \
--downsample_ratio 1 \



sc_long_pipeline.py \
--gff3  /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf \
--infq /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM3/barcode_matching_output/out.fastq.gz \
--outdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/hBM3/FLAMES_out \
--genomefa /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/GRCh38.primary_assembly.genome.fa \
--minimap2_dir /hpcnfs/scratch/PGP/niman/Chiara/tools/minimap2/minimap2 \
--config_file /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts/config_sclr_nanopore_non_enriched.json \
--downsample_ratio 1 \
