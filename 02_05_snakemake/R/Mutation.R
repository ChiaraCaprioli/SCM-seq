args = commandArgs(trailingOnly=TRUE)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)
library(biovizBase)
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(bamsignals)
library(GenomeInfoDb)
library(BSgenome)
library(Biostrings)
library(readxl)
library(tidyverse)
library(annotables)
library(motifStack)
library(data.table) 
library(dplyr)
#arg<- c("/hpcnfs/scratch/PGP/SCMseq/AML4/merged_enriched","/hpcnfs/scratch/PGP/niman/snakemake_flames/AML4")
#arg<- "/hpcnfs/scratch/PGP/niman/Chiara/FLAMES/sAML1_B/Long/enrichment_analysis/merged"
#arg[1]
#Rscript 

 args=c("/hpcnfs/scratch/PGP/SCMseq/aml5_promethion/merged_enriched/out_flames_sc_enriched",
        "/hpcnfs/scratch/PGP/niman/snakemake_flames/aml5_promethion",
       "/hpcnfs/scratch/PGP/niman/snakemake_flames",
        "/hpcnfs/scratch/PGP/SCMseq/aml5_promethion/merged_enriched/Mutation",
        "/hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/gencode.v24.annotation.gtf",
        "/hpcnfs/scratch/PGP/niman/Chiara/FLAMES/h38_cagepeak/GRCh38.primary_assembly.genome.fa")
#############
inputlist_<- readxl::read_xlsx(paste0(args[2],"/","AML.xlsx"))
sample_name<- inputlist_$Sample[1]
tg_genes<- read.csv(paste0(args[3],"/","new_panel.csv"),stringsAsFactors = F,sep = ";")

#inputlist<- inputlist_[which(inputlist_$Gene %in% tg_genes$Target),]
tg_genes<- tg_genes[which(tg_genes$Sample %in%  sample_name) ,]
inputlist<- inputlist_[which(inputlist_$Gene %in% tg_genes$Gene),]
genelist=c(inputlist$Gene[1:length(inputlist$Gene)])
dir.create(args[4])
write_csv(as.data.frame(genelist),paste0(args[4],"/","tg_wes_intersect.csv"))
gtfann<- args[5]
fag<- args[6]



#bampath=paste(arg[1],"/FLAMES_output_3FC/","align2genome.bam",sep = "")
#baminpath=paste(arg[1],"/FLAMES_output_3FC/","align2genome.bam.bai",sep = "")
#gffin=paste(arg[1],"/FLAMES_output_3FC/","isoform_annotated.filtered.gff3",sep = "")


bampath=paste(args[1],"/align2genome.bam",sep = "")
baminpath=paste(args[1],"/align2genome.bam.bai",sep = "")
gffin=paste(args[1],"/isoform_annotated.filtered.gff3",sep = "")
outpath=paste(args[4],"/",sep = "")
min_map_q=20 #### adjusted by minimap2 0 to 60
min_base_q_in= 10 ### qscore 0 to 40 , ascii to int from 33 to 74 ## ex: qscore 15 >asci 48 
min_base_q<- min_base_q_in + 33 ## change qscore to asscci
motifth = 10
options(ucscChromosomeNames=FALSE)
#mut_ref_count=function(gffin,outpath){
# min_base_q<- min_base_q_in + 33 ## change qscore to asscci
opts <- getOption("biovizBase")
opts$DNABasesNColor[1] <- "#006400"  ##A
opts$DNABasesNColor[2] <- "red"    ##T
opts$DNABasesNColor[3] <- "orange" ##G
opts$DNABasesNColor[4] <- "blue"   ##C
options(biovizBase = opts)
  
  chr_ide<- function(chrn,st,en,fp){
    #chrn="chr13"
    chr_ideogram <- Gviz::IdeogramTrack(genome = "hg38", chromosome = chrn)
    #Gviz::IdeogramTrack(genome = "hg38", chromosome = "chr13")
    pdf(file = paste(fp,"Ideogarm_",chrn,".pdf",sep = ""))
    Gviz::plotTracks(chr_ideogram,from =st,to =en,extend.left =0.5,extend.right = 0.5,chromosome =chrn)
    dev.off()
    return(chr_ideogram)
  }
  #############
  grtr<- function(chrn,st,en,sym,txdb,fp){
    seqlevels(txdb) <- chrn
    grtrack <- Gviz::GeneRegionTrack(txdb, start = st, end = en,showId = F,
                               name = "Gene Annotation",gene =sym
    )
    pdf(file = paste(fp,"GenTrack_",sym,".pdf",sep = ""))
    Gviz::plotTracks(list(grtrack),from = st,to = en,cex.main = 0.3)
    dev.off()
    return(grtrack)
  }
  
  ################ 
  hgtr<- function(chr_ideogram,grtrack,st,en,chrn,gnm,sym,fp){
    print(st)
    gtrack <- Gviz::GenomeAxisTrack(cex = 0.5)  # set the font size larger
    ht<-Gviz::HighlightTrack(trackList = list(gtrack,grtrack),start   =st,width =40 ,chromosome =chrn,genome =   gnm)
    pdf(file = paste(fp,"HL_chr_Gentr_GenLen_",sym,".pdf",sep = ""))
    Gviz::plotTracks(list(chr_ideogram,ht),from = st,to = en)
    dev.off()
    return(gtrack)
  }
  
 
  
  
  
  get_mutation <- function(i,alt,chr_edit,fp,altrack)
  {
    mut_ <- inputlist$Start_Position[i]
    ################## get the address of the gene plus 200 bases before and after and adding database
    st <- alt$start[1] -200
    en <- alt$end[1]+200
    txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
    ###################### create Ideogram based on the reffernce
    chr_ideogram<-chr_ide(inputlist$Chromosome[i],st,en,fp)
    ###################### create transcript based on the reffernce
    grtrack<- grtr(inputlist$Chromosome[i],st,en,inputlist$Gene[i],txdb,fp)
    ###################### create gene legth line with highlight based on the reffernce
    gtrack<-hgtr(chr_ideogram,grtrack,inputlist$Start_Position[i]-20,en,inputlist$Chromosome[i],"TxDb.Hsapiens.UCSC.hg38.knownGene",inputlist$Gene[i],fp)
    ########gff3 track
    #txdbFromGFF <- GenomicFeatures::makeTxDbFromGFF(file = gffin) 
    #s<- Gviz::GeneRegionTrack(txdbFromGFF,start =st,end = en, genome= "TxDb.Hsapiens.UCSC.hg38.knownGene",chromosome = inputlist$Chromosome[i],gene=inputlist$Gene[i] )
    ###################### Alignment track from bam with highlight long
    
    #ht<- Gviz::HighlightTrack(trackList = list(gtrack,altrack,grtrack,s),start   =mut_-20,width = 40 ,chromosome =inputlist$Chromosome[i],genome=  "TxDb.Hsapiens.UCSC.hg38.knownGene")
    # pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_",inputlist$Gene[i],".pdf",sep = ""))
    # Gviz::plotTracks(list(chr_ideogram,ht),from = st,to = en)
    
    # dev.off()
    
    
    ###################### Alignment track from bam with highlight , short
    ###   altrack_short <- AlignmentsTrack( 
    #   bampath_short, isPaired = T, col.mates =  "deeppink")
    #   ht<- HighlightTrack(trackList = list(gtrack,grtrack,altrack_short),start   =inputlist$Start_Position-20,width = 40 ,chromosome =inputlist$Chromosome[i],genome=  "TxDb.Hsapiens.UCSC.hg38.knownGene")
    # pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_short_",inputlist$Gene[i],".pdf"))
    #  plotTracks(list(altrack_short),chromosome = inputlist$Chromosome[i],from = st,to = en)
    
    # dev.off()
    ###
    ###################### sequence track from refference with highlight 
   # strack <- Gviz::SequenceTrack(
   #    BSgenome.Hsapiens.UCSC.hg38::Hsapiens,chromosome = inputlist$Chromosome[i],from =st,to = en,fontcolor=getBioColor("DNA_BASES_N"))
    #ht<- Gviz::HighlightTrack(trackList = list(gtrack,altrack,grtrack),start  =mut_-10,width = 20 ,chromosome =inputlist$Chromosome[i],genome =  "TxDb.Hsapiens.UCSC.hg38.knownGene")
    #turn on the pileup
    # pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_",inputlist$Gene[i],".pdf"))
    # plotTracks(list(chr_ideogram,ht),from = st,to = en,extend.left = 0.5, extend.right =  0.5,min.height = 0, coverageHeight = 0.08, minCoverageHeight = 0) 
    # dev.off() 
    
    ###################### turn off the pileup, zoom into region,sequence track from refference with highlight 
    #pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_1_",inputlist$Gene[i],".pdf",sep = ""))
    #Gviz::plotTracks(list(chr_ideogram,ht),from =mut_-1000,to =  mut_+1000,extend.left = 0.5, extend.right = 0.5,min.height = 0, 
    #            coverageHeight = 0.08, minCoverageHeight = 0)
    # dev.off()
    ####
    # plotTracks(altrack,from =mut_-1000,to =  mut_+1000,extend.left = 0.5, extend.right = 0.5,min.height = 0, 
    #          coverageHeight = 0.08, minCoverageHeight = 0)
    ######
    #pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_2",inputlist$Gene[i],".pdf",sep = ""))
    #Gviz::plotTracks(list(chr_ideogram,ht),from =mut_-50,to = mut_+150,extend.left = 0.5, extend.right = 0.5,min.height = 0, 
    #  coverageHeight = 0.08, minCoverageHeight = 0)
    #dev.off() 
    
    #ht<- Gviz::HighlightTrack(trackList = list(gtrack,altrack,grtrack,strack),start  =mut_, width = 0 ,chromosome =inputlist$Chromosome[i],genome = "TxDb.Hsapiens.UCSC.hg38.knownGene")
    
    ###zoom in
    # pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_3",inputlist$Gene[i],".pdf",sep = ""))
    #  Gviz::plotTracks(list(chr_ideogram,ht),from  =mut_-50,to =   mut_+150,extend.left = 0.5, extend.right = 0.5,min.height = 0,   coverageHeight = 0.08, minCoverageHeight = 0)
    #  dev.off() 
    
   # pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_4_",inputlist$Gene[i],".pdf",sep = ""))
   # Gviz::plotTracks(list(chr_ideogram,ht),from =mut_-10,to = mut_+10,extend.left = 0.5, extend.right = 0.5,add53 = TRUE,min.height = 3)
   # dev.off() 
    
    #pdf(file = paste(fp,"HL_chr_Gentr_GenLen_Bam_seq_5",inputlist$Gene[i],".pdf",sep = ""))
    # Gviz::plotTracks(list(chr_ideogram,ht),complement = TRUE,from  =mut_-10,to = mut_+10,extend.left = 0.5,  extend.right = 0.5,title.width = T,add53 = TRUE, min.height = 3)
    #dev.off() 
    
    ##################### Finding Mutation plus checking quality of mapping and base quality, also we extract the barcodes related to mutations
    stadd<- mut_- motifth
    endadd<-mut_+ motifth
    #muadd<- inputlist$Start_Position[i]
    chr<-inputlist$Chromosome[i]
    which <- GenomicRanges::GRanges(seqnames = chr, IRanges(stadd, endadd))
    p1 <- Rsamtools::ScanBamParam(which=which, what=scanBamWhat(),mapqFilter = 0)
    reads <- GenomicAlignments::readGAlignments(bampath, param=p1)
    write.csv(colSums(cigarOpTable(cigar(reads))),paste(fp,"cigar_table_reads_coverTheMutation.csv",sep = ""))
    
    bamcov<- bamsignals::bamCoverage(bampath, which, verbose=F,  mapqual = 0)
    write.csv(bamcov@signals[[1]][motifth+1],paste(fp,"bamcov_.csv",sep = ""))
    
    ################# Read refference to fins ref and allele seq
     genome <- BSgenome.Hsapiens.UCSC.hg38
    # read_ranges <- GenomicRanges::granges(reads)
    # strand(read_ranges) <- "+"
    # GenomeInfoDb::seqlevels(read_ranges)<- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
    # rseqs <- BSgenome::getSeq(genome, read_ranges)
    # ################ align the ref to our sample to extract the sequence in our sample it help  us to find the mutation position
    # rseqs2 <- GenomicAlignments::sequenceLayer(rseqs, cigar(reads), from="reference", to="query")
    # qseqs <- mcols(reads)$seq
    # ###############  query and ref should have the same length, lets check it here
    # identical(width(qseqs), width(rseqs2))  
    ###############  lets find the motif and plot it
    ############### check mapping quality
    filtered_bam<- GenomicAlignments::stackStringsFromBam(bampath, param=p1,use.names = T,what = "seq")
    all_reads<- filtered_bam%>% names(.) 
    data.table::fwrite(as.data.frame(all_reads),paste(fp,"allread_id_table_coverTheRegion.csv",sep = ""),row.names = T)
    
    rm("all_reads")
    rm("reads")

    #write.csv(all_reads,paste(fp,"allread_id_table_coverTheRegion.csv",sep = ""))
    bam_mapq_read_names <- filtered_bam[filtered_bam@elementMetadata$mapq >= min_map_q] %>% names(.)%>%BiocGenerics::unique()
    ############## check base quality , need to change asscii to number
    filtered_bam_q<- GenomicAlignments::stackStringsFromBam(bampath, param=p1,use.names = T,what = "qual")
    bam_bqual_read_names=list()
  
    
    
    
    q_<-as.data.frame(filtered_bam_q)
    q_list<- lapply(q_$x,function(z)utf8ToInt(z))
    s<- do.call(rbind, q_list)
    ls<- which(s[1:length(q_list),motifth+1]>= min_base_q)
    bam_bqual_read_names<- filtered_bam_q[ls]%>% names(.)%>%BiocGenerics::unique()
    
    #bam_bqual_read_names<- lapply(1:length(q_list),function(a) if (q_list[[a]][motifth+1] > min_base_q) {
   #   bam_bqual_read_names<- append(bam_bqual_read_names,(filtered_bam_q[a]%>% names(.)))
   # })
   # bam_bqual_read_names<- bam_bqual_read_names%>%unlist()%>%BiocGenerics::unique()
    # for (q in 1:length(filtered_bam_q)) {
  
    #}
    
    
    ############## intersect the read ids to get the list of the reads which passed QC
    final_list_read_id<- Reduce(intersect, list(bam_bqual_read_names,bam_mapq_read_names)) %>%BiocGenerics::unique()
    #write.csv(final_list_read_id,paste(fp,"read_id_table_coverTheRegion_","MP",min_map_q,"_BQ",min_base_q - 33,".csv",sep = ""))
    data.table::fwrite(as.data.frame(final_list_read_id),paste(fp,"read_id_table_coverTheRegion_","MP",min_map_q,"_BQ",min_base_q - 33,".csv",sep = ""),row.names = T)
    bam_qc_passed<- filtered_bam[names(filtered_bam) %in% final_list_read_id] 
    #bam_qc_passed_q<- filtered_bam_q[names(filtered_bam_q) %in% final_list_read_id]
  
    lfb<- length(filtered_bam)
    rm("filtered_bam_q")
    rm("filtered_bam")
    rm("q_list")
    ############## extract read ids support or not support mutation from the final list
    reads_support_mutation =list()
    reads_support_ref=list()
    reads_support_mismatch=list()
    

    
    ##########
    
    ########
    al2<- inputlist$Tumor_Seq_Allele2[i]
    al1<- inputlist$Reference_Allele[i]
    
    # case_<- function(x){
    #   case_when(
    #     x==al2 ~ "MUT",
    #     x==al1 ~ "REF",
    #     TRUE ~ "MIS" # Default value if none of previous conditional holds
    #   )
    # }
    # 
    bam_qc_passed<- bam_qc_passed[bam_qc_passed%>%names(.)%>%BiocGenerics::unique()]
    q_<-as.data.frame(bam_qc_passed)
    #####################
    a<- as.list(q_$x)
    a<-do.call(rbind, a)
    b<-do.call(rbind, strsplit(a,""))%>% as.matrix()
    mu_<- which(b[1:length(a),motifth+1]==al2)
    re_<- which(b[1:length(a),motifth+1]==al1)
    mis_<- which(b[1:length(a),motifth+1]!= al1 & b[1:length(a),motifth+1]!= al2)
    length(q_$x)
    length(mu_)+length(re_)+length(mis_)
    reads_support_mutation<- bam_qc_passed[mu_]%>% names(.) %>% BiocGenerics::unique()
    reads_support_ref<-  bam_qc_passed[re_]%>% names(.) %>% BiocGenerics::unique()
    reads_support_mismatch<-  bam_qc_passed[mis_]%>% names(.) %>% BiocGenerics::unique()
    length(reads_support_mutation)+length(reads_support_ref)+length(reads_support_mismatch)
    ######################
    # read_supp<-  lapply(1:length(q_$x), function(z){
    #   strsplit(q_$x[z], "")[[1]][motifth+1]%>%case_()
    # })
    # read_supp%>%unlist()
    # rm("q_")
    # reads_support_mutation<- bam_qc_passed[which(read_supp=="MUT")]%>% names(.) %>% BiocGenerics::unique()
    # reads_support_ref<-  bam_qc_passed[which(read_supp=="REF")]%>% names(.) %>% BiocGenerics::unique()
    # reads_support_mismatch<-  bam_qc_passed[which(read_supp=="MIS")]%>% names(.) %>% BiocGenerics::unique()
    
    #################
    
    
    vaf<- (length(reads_support_mutation)/(length(reads_support_mutation)+length(reads_support_ref)))
    write.csv(reads_support_mutation,paste(fp,"read_id_support_mutation_","MP",min_map_q,"_BQ",min_base_q-33,".csv",sep = ""))
    write.csv(reads_support_ref,paste(fp,"read_id_support_ref_","MP",min_map_q,"_BQ",min_base_q-33,".csv",sep = ""))
    write.csv(reads_support_mismatch,paste(fp,"read_id_mismatches_","MP",min_map_q,"_BQ",min_base_q-33,".csv",sep = ""))
    
    ###################Avarage base quality
    # av_bq<- matrix(1:125, nrow = 5, ncol = 21,data = 0)
    # rownames(av_bq)<- c("A","C","G","T","-")
    # 
    # x_<- as.data.frame(bam_qc_passed_q)
    # x__q<- lapply(x_$x,function(z)utf8ToInt(z)-33)
    # x__seq<- as.data.frame(bam_qc_passed)
    # x__seq_<- lapply(x__seq$x,function(z) strsplit(z, "")[[1]])
    # for (q in 1:length(bam_qc_passed)){
    #   for (l in 8:14)
    #     if (x__seq_[[q]][l]=="A"){ 
    #       av_bq[1,l]<-av_bq[1,l]+x__q[[q]][l]
    #     }
    #   else if (x__seq_[[q]][l]=="C"){ 
    #     av_bq[2,l]<-av_bq[2,l]+x__q[[q]][l]
    #   } 
    #   else if (x__seq_[[q]][l]=="G"){ 
    #     av_bq[3,l]<-av_bq[3,l]+x__q[[q]][l]
    #   } 
    #   else if (x__seq_[[q]][l]=="T"){ 
    #     av_bq[4,l]<-av_bq[4,l]+x__q[[q]][l]
    #   }
    #   else if (x__seq_[[q]][l]=="-"){ 
    #     av_bq[5,l]<-av_bq[5,l]+x__q[[q]][l]
    #   }
    # }
    # pcm<- Biostrings::consensusMatrix(bam_qc_passed)
    # #pcm<-t(as)%>% .[-5,]
    # pcm <-pcm %>% .[c(1,2,3,4,16),]
    # rownames(pcm) <- c("A","C","G","T","-") # must have rownames
    # colnames(av_bq)<- t(unlist(strsplit(as.character(BSgenome::getSeq(genome,which)[[1]]),""))) 
    # av_bq<-round((av_bq)/pcm,2)
    # write.csv((av_bq),paste(fp,"avarage_base_QC",".csv",sep = ""))
    ################ MOTIF batch 1 batch 2
    
    pcm<- Biostrings::consensusMatrix(bam_qc_passed)
    #pcm<-t(as)%>% .[-5,]
    pcm <-pcm %>% .[c(1,2,3,4),]
    rownames(pcm) <- c("A","C","G","T") # must have rownames
    markerRect <- new("marker", type="rect", start=10, stop=12, 
                      gp=gpar(lty=2, fill=NA, col="orange", lwd=3))
    motif <- new("pcm", mat=as.matrix(pcm),name=paste(inputlist$Gene[i],"_",inputlist$Protein_Change[i],"_",inputlist$Variant_Type[i],"_","Ref_",inputlist$Reference_Allele[i],"_","alt_",inputlist$Tumor_Seq_Allele2[i],"_",round(vaf,3),sep = ""))
    motif$markers <- list(markerRect)
    pdf(file = paste(fp,"MOTIF_",inputlist$Gene[i],"_MP",min_map_q[[1]],"_BQ",min_base_q[[1]]-33,".pdf",sep = ""),width = 12,height =4)
    plot(motif,ylab="probability")
    dev.off()
    ##################### probability of each base
    pcm_data<- Biostrings::consensusMatrix(bam_qc_passed) %>% .[c(1,2,3,4,16,18),] 
    colnames(pcm_data)<- t(unlist(strsplit(as.character(getSeq(genome,which)[[1]]),""))) 
    write.csv((pcm_data),paste(fp,"motif_table_","MP",min_map_q[[1]],"_BQ",min_base_q[[1]]-33,".csv",sep = ""))
    rm("pcm")
    rm("pcm_data")
    pcm <-Biostrings::consensusMatrix(bam_qc_passed) %>% .[c(1,2,3,4,16),] 
    rownames(pcm) <- c("A","C","G","T","-") # must have rownames
    pcm <- as.data.frame(t(pcm))
    pcm$num<- as.numeric(as.character(rownames(pcm)))
    pcm$ref<- t(t(unlist(strsplit(as.character(getSeq(genome,which)[[1]]),""))))
    
    t<-reshape2::melt(pcm,id.vars  =c("num","ref"),factorsAsStrings =F)%>%as.data.frame(.)%>% group_by(num) %>%
      mutate(frequency = as.numeric(value) / sum(as.numeric(value))) %>% ggplot(aes(as.factor(num), y = frequency, fill = variable))+
      geom_bar(stat = "identity", position = position_dodge(width=0.5), color = "black", width = 1) +labs(x =paste(inputlist$Gene[i],"_",inputlist$Protein_Change[i],"_",inputlist$Variant_Type[i],"_","Ref_",inputlist$Reference_Allele[i],"_","alt_",inputlist$Tumor_Seq_Allele2[i],"_","Allele_fre_",round(vaf,3),sep = ""), 
                                                                                                          y = "Frequency of each base from bam",fill = "Base") +scale_y_continuous(limits = c(0, 1)) +
      scale_fill_manual("legend", values = c("A" = "#006400", "T" = "red", "C" = "blue","G"="orange","-"="black"))+
      scale_x_discrete(label=pcm$ref)+ggplot2::annotate("rect", xmin=c(10.5), xmax=c(11.5), ymin=c(0) , ymax=c(0.9), alpha=0.2, color="black", fill="#16E2F5")+
      ggplot2::annotate("text", x = 11, y = c(1), label = "Mut") +theme_classic()
    pdf(file = paste(fp,"MOTIF_freq_",inputlist$Gene[i],"_MP",min_map_q[[1]],"_BQ",min_base_q[[1]]-33,".pdf",sep = ""),width = 12,height =4)
    print(t)
    dev.off()
    
    t<- reshape2::melt(pcm,id.vars  =c("num","ref"),factorsAsStrings =F)%>%as.data.frame(.)%>% group_by(num) %>%
      mutate(frequency = as.numeric(value) / sum(as.numeric(value))) %>% ggplot(aes(num,as.numeric(frequency),col=variable)) +
      geom_line(size=0.5, alpha=0.9, linetype=1) +
      geom_point()+
      labs(x =paste(inputlist$Gene[i],"_",inputlist$Protein_Change[i],"_",inputlist$Variant_Type[i],"_","Ref_",inputlist$Reference_Allele[i],"_","alt_",inputlist$Tumor_Seq_Allele2[i],"_","Allele_fre_",round(vaf,3),sep = ""), 
           y = "Frequency of each base from bam")  +facet_wrap(vars(variable), ncol = 1) +
      scale_color_manual(values =c("A" = "#006400", "T" = "red", "C" = "blue","G"="orange","-"="gray"))+
      ggplot2::annotate("rect", xmin=c(11), xmax=c(11), ymin=c(0) , ymax=c(0.9), alpha=0.0, color="black", fill="white")+
      scale_x_continuous(breaks = 1:21,
                         labels = pcm$ref)
    pdf(file = paste(fp,"MOTIF_freq_linegraph",inputlist$Gene[i],"_MP",min_map_q[[1]],"_BQ",min_base_q[[1]]-33,".pdf",sep = ""),width = 16,height =10)
    print(t) 
    dev.off() 
    ##############error rate
    t<-reshape2::melt(pcm,id.vars  =c("num","ref"),factorsAsStrings =F)%>%as.data.frame(.)%>% group_by(num) %>%
      mutate(frequency = as.numeric(value) / sum(as.numeric(value)))
    for (r in 1:length(t$num)){
      if (t$ref[r]==t$variable[r]) {t$value[r]=0
      t$frequency[r]=0}
    }
    err1<-reshape2::dcast(t,variable ~ num,value.var = "frequency")%>% replace(.,"11",0) %>% summarise(across(where(is.numeric), sum)) %>% rowMeans()
    ####################
    
    ####################### checking Mismatches according to quality control
    
    
    ############### adding table for number of the reads and etc ....
    table_out=list()
    table_out$Chr <- inputlist$Chromosome[i]
    table_out$Gene <- inputlist$Gene[i]
    table_out$G_start <- alt$start[1]
    table_out$G_end <- alt$end[1]
    table_out$Mutation_add <- inputlist$Start_Position[i]
    table_out$Ref <- inputlist$Reference_Allele[i]
    table_out$Alt <- inputlist$Tumor_Seq_Allele2[i]
    table_out$Description <- inputlist$Description[i]
    table_out$bamcov_mut_BQC<- bamcov@signals[[1]][motifth+1]
    table_out$All_reads<- lfb
    table_out$All_reads_filtered<- length(bam_qc_passed %>% names(.) %>%BiocGenerics::unique())
    table_out$Ref_reads<- length(reads_support_ref)  
    table_out$Alt_reads<- length(reads_support_mutation)
    table_out$Mism_reads<- length(reads_support_mismatch)
    table_out$QC1_base_map<- c(min_base_q[[1]]-33,min_map_q[[1]])
    table_out$Allele_freq<- vaf
    table_out$error_rate<- err1
    to<- data.frame(t(sapply(unlist(table_out),c)))
    write.csv(to,paste(fp,"output_summary.csv",sep = ""))
    rm("pcm")
    rm("pcm_data")
    ###########
    rm("bam_qc_passed")
    rm("bam_qc_passed")
    x_<- read.csv(paste(fp,"output_summary.csv",sep = ""),stringsAsFactors = T)%>%dplyr::select(All_reads,All_reads_filtered,Ref_reads,Alt_reads,Mism_reads)%>%t()%>%as.data.frame()%>%
      dplyr::rename(.,Number_of_reads=V1)%>%dplyr::mutate(Info_FLT3=factor(row.names(.),levels = c("All_reads","All_reads_filtered","Ref_reads","Alt_reads","Mism_reads")),GP=factor(c("All_bqc","All_aqc","ref_alt_mis","ref_alt_mis","ref_alt_mis"),levels =c("All_bqc","All_aqc","ref_alt_mis" )))
    a<-ggplot(x_, aes(x = GP, y = Number_of_reads, fill = Info_FLT3)) + 
      geom_bar(stat = 'identity', position = 'stack') +
      scale_fill_manual(paste("Info_",inputlist$Gene[i],sep = ""), values = c("All_reads" = "darkblue", "All_reads_filtered" = "orange", "Ref_reads" = "red", "Alt_reads" = "blue", "Mism_reads" = "gray"))+
      geom_text(aes(label=Number_of_reads), position=position_stack(vjust = 0.02),col="white")+ theme_dark()
    pdf(file = paste(fp,"barplot_1",inputlist$Gene[i],".pdf",sep = ""),width = 8,height =4)
    print(a)
    dev.off()
    ######################barplot of the table
  }  
  
  #####
  
  # inputlist= inputlist[6,]
  #lapply(inputlist$Gene, function(ii){
  for (i in 1:length(inputlist$Gene)){
    #ii= "WDR63"
    #i=1
    ii<- inputlist$Gene[i]
    i= which(inputlist$Gene==ii)
    altrack <- Gviz::AlignmentsTrack( 
      bampath, isPaired = F, col.mates =  "deeppink",showIndels=T,chromosome = inputlist$Chromosome[i])
    
    
    chr<-inputlist$Chromosome[i]
    which <- GenomicRanges::GRanges(seqnames = chr, IRanges(inputlist$Start_Position[i] - 10, inputlist$End_Position[i] + 10))
    bco<- bamsignals::bamCoverage(bampath = bampath,gr = which)
    if (inputlist$Chromosome[i]=="chrY" | bco@signals[[1]][11]<20){
      print(paste0("either ",inputlist$Gene[i]," is in ChrY or not having sufficient cov.","_cov= ",bco@signals[[1]][11]))
      next 
      
    }
    ###creating subfolders for results
    if (!is.na(inputlist$Gene[i])){ 
      fp= file.path(paste0(outpath,"Mutational_analysis/",inputlist$Gene[i],"_",inputlist$Chromosome[i],"_",inputlist$Start_Position[i],"/"))
      if (!dir.exists(fp)){
        dir.create(fp,recursive = T)
      }
    }
    
    
    ######################## finding gene and chr address
    chr_edit= sub("chr", "",inputlist$Chromosome[i])
    alt<- annotables::grch38 %>% filter(symbol==inputlist$Gene[i] & chr==chr_edit) %>% 
      dplyr::select(ensgene, symbol, chr, start, end, description) %>% dplyr::right_join(.,inputlist[i,1:12],by=c("symbol"="Gene"))
    write.csv(alt,paste0(fp,'gene_chr_description.csv'))
    print(paste0(i,fp))
    
    
    get_mutation(i,alt,chr_edit,fp,altrack)
    
    gc(full = T)
  }
  
  
  
  #####
  
   ###################
  