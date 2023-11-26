#! /usr/bin/env python

"""
extract_reads.py
Created by iman
"""


import pysam
import re

import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

a=os.listdir("/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq_filtered_geneByORF_and_MAD")
#cell_id_path=("/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/cell_ID_wt_AML4_HSC_genes.txt","jgh")
cell_id_path=find('cell_ID*.txt', "/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq_filtered_geneByORF_and_MAD")
samples=("AML4","AML5","sAML1")
output="/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Transcript_bam_split/"

#/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Tr_ID_mut_AML4_HSC_genes.txt
#s= fhand.readlines()
#x=[line.rstrip('\n') for line in s]
for i in cell_id_path:
    for ii in samples:
        if re.search(ii,i):
            if re.search("mut",i):
                fhand = open("/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq_filtered_geneByORF_and_MAD/Tr_ID_mut_"+ii+"_HSC_genes.txt","r")
                s= fhand.readlines()
                x=[line.rstrip('\n') for line in s]
            else:
                fhand = open("/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/isoswitch_result/Drim_seq_filtered_geneByORF_and_MAD/Tr_ID_wt_"+ii+"_HSC_genes.txt","r")
                s= fhand.readlines()
                x=[line.rstrip('\n') for line in s]
            infile= pysam.AlignmentFile("/hpcnfs/scratch/PGP/SCMseq/"+ii+"_promethion/merged_non_enriched/out_flames_sc_non_enriched/realign2transcript.bam")
            fq= open(i).readlines()
            fq= [o.strip() for o in fq]
            fq= set(fq)
            g=re.split(".txt",re.split("/cell_ID",i)[1])[0]
            bam_out=output+"bam"+g
            outfile= pysam.AlignmentFile(bam_out+".bam", template= infile, mode= 'wb')
            for aln in infile:
                w= re.split("_",aln.query_name)[0]
                if w in fq:
                    if aln.reference_name in x:
                    #print(w)
                        outfile.write(aln)
            outfile.close()
            pysam.sort("-o",bam_out+"_sorted.bam", bam_out+".bam")
            pysam.index(bam_out+"_sorted.bam")
            os.remove(bam_out+".bam")
