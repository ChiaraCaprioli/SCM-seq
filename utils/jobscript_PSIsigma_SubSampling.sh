#!/bin/bash
#PBS -q nocg_workq
#PBS -N PSIsigma
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -M iman.nazari@ieo.it
#PBS -e /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut
#PBS -o /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut
#PBS -m ae

declare -a sample_=("sAML1" "AML5" "AML4")

declare -a batch_=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")
declare -a n_=("1" "2" "3")
declare -a n_w=("1" "2")
cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut
d=()
for sample_x in ${sample_[@]}; do
	if [ "$sample_x" == "AML4" ]; then
    		d=(3) 
  	else
   		d=(10)
 	fi
	for batch_x in ${batch_[@]}; do
		mkdir "${sample_x}_${batch_x}"
		cd "${sample_x}_${batch_x}"
		for n_x in ${n_[@]}; do 
			find /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_bam_AML4_AML5_sAML1_GTF_hBM123_AML2345_sAML1 -iname "${sample_x}${batch_x}_${n_x}_HSC_genes_bam_sorted.bam" >> groupb.txt
		done
		for n_x in ${n_w[@]}; do 
			find /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_bam_AML4_AML5_sAML1_GTF_hBM123_AML2345_sAML1 -iname "${sample_x}${batch_x}_${n_x}_wt_HSC_genes_bam_sorted.bam" >> groupa.txt
		done
		singularity exec --bind /hpcnfs:/hpcnfs /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/psi_sigma_pipeline_3.9.sif \
 		perl /usr/local/bin/PSI-Sigma-2.1/dummyai.pl \
		--gtf /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/GTF_FASTA/gffcomare_AML_12345_hBM_123.combined.gtf \
		--name splicing_mut_wt_srsf2_HSC_myeloid --type 2 --nread $d \
		--output /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut/"${sample_x}_${batch_x}"/results \
		--groupa /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut/"${sample_x}_${batch_x}"/groupa.txt \
		--groupb /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut/"${sample_x}_${batch_x}"/groupb.txt \
		--irmode 0 --threads 4 --fmode 3 --adjp 1
		cd /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results/20epoch_2wt_3mut
		find  /hpcnfs/scratch/PGP/SCMseq/isoform_analysis/PSIsigma/subsampling_samplewise_barcode/subSampling_splicing_results -type d -name _wososatmp -exec rm -r {} +
	done
done

exit


