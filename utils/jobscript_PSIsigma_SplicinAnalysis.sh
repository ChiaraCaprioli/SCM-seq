#!/bin/bash
#PBS -q workq
#PBS -N PSIsigma
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma
#PBS -m ae

cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma

singularity exec --bind /hpcnfs:/hpcnfs psi_sigma_pipeline_3.9.sif perl /usr/local/bin/PSI-Sigma-2.1/dummyai.pl --gtf /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/GTF_FASTA/gffcomare.combined.gtf --name splicing_mut_wt_srsf2_HSC_myeloid --type 2 --nread 5 --output /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/Splicing_Results/Splicing_AllTR/results --groupa /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/Splicing_Results/Splicing_AllTR/groupa.txt --groupb /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/Splicing_Results/Splicing_AllTR/groupb.txt --threads 2 --fmode 2  --adjp 1

exit


