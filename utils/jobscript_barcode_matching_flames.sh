#!/bin/sh
#PBS -q workq
#PBS -N barcode_matching
#PBS -l select=1:ncpus=16:mem=32g
#PBS -o /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts
#PBS -e /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/Scripts
#PBS -M iman.nazari@ieo.it
#PBS -m bae
PBS_O_WORKDIR="/hpcnfs/home/ieo5268/FLAMES/src"
cd /hpcnfs/home/ieo5268/FLAMES/src

mkdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/AML2/barcode_statistics
mkdir /hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/AML2/barcode_matching_output


echo "The working directory is: 'PBS_O_WORKDIR'"
match_cell_barcode  /hpcnfs/techunits/genomics/PublicData/TSSM/ccaprioli/FAST5/20230517_1307_P2S_00101-B_PAK60283_cf3e4d41/S50759_AML2/fastq_pass/ \
/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/AML2/barcode_statistics/stat_bc.txt \
/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/AML2/barcode_matching_output/out.fastq.gz \
/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/AML2/short_bardcodes.txt \
1 \
12
