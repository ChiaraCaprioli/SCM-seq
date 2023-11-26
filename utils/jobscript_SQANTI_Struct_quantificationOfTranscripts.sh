#!/bin/sh
#PBS -q workq
#PBS -N SQANTI
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -o /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/SQANTI
#PBS -e /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/SQANTI
#PBS -M iman.nazari@ieo.it
#PBS -m bae
PBS_O_WORKDIR="/hpcnfs/scratch/PGP/SCMseq/IsoSwitch/SQANTI"
cd /hpcnfs/scratch/PGP/niman/Chiara/tools/SQANTI/SQANTI3
source activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:/hpcnfs/scratch/PGP/niman/Chiara/tools/SQANTI/cDNA_Cupcake/sequence/
echo "The working directory is: 'PBS_O_WORKDIR'"

declare -a StringArray=("sAML1" "AML4" "AML5")



for val in ${StringArray[@]}; do
	python /hpcnfs/scratch/PGP/niman/Chiara/tools/SQANTI/SQANTI3/sqanti3_qc.py -g /hpcnfs/scratch/PGP/SCMseq/${val}_promethion/merged_non_enriched/out_flames_sc_non_enriched/isoform_annotated.filtered.gff3 /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/GRCh38.primary_assembly.genome.fa --cage_peak /hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/hg38.cage_peak_phase1and2combined_coord.bed -d /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/SQANTI/$val --force_id_ignore
done