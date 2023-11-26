#!/bin/bash
#PBS -q workq
#PBS -N CPAT
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/CPAT
#PBS -o /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/CPAT
#PBS -m ae

source ~/.bashrc
Source activate CPAT
cd /hpcnfs/scratch/PGP/SCMseq/IsoSwitch/CPAT
declare -a StringArray=("sAML1" "AML4" "AML5")



for val in ${StringArray[@]}; do
	cpat.py -x Human_Hexamer.tsv  -d  Human_logitModel.RData  --top-orf=100  --antisense -g /hpcnfs/scratch/PGP/SCMseq/${val}_promethion/merged_non_enriched/out_flames_sc_non_enriched/transcript_assembly.fa -o $val/cpat_
done
exit


