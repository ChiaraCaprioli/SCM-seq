args = commandArgs(trailingOnly=TRUE)
custom_colors <- c("#2d98da", "#a5b1c2", "#f7b731", "#4b6584", "#0fb9b1", "#fc5c65", "#20bf6b", "#fa8231")

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
library(Seurat)
#
args=c("/hpcnfs/scratch/PGP/SCMseq/aml5_promethion/merged_enriched/Mutation",
              "/hpcnfs/scratch/PGP/SCMseq/aml5_promethion/merged_enriched/BC_DEMUX",
            "/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML5.rds",
              "/hpcnfs/scratch/PGP/niman/snakemake_flames/aml5_promethion/BC.csv")
       ############
out_dir<- paste0(args[2],"/")
dir.create(args[2])
#out_dir<- paste0(args[2],"/BC_demux/")

raw_matrix_barcodes<- read_csv(args[4])
AML<- readRDS(args[3])
####
save = function(data,name){
  x<- as.data.frame(data) %>% tibble::rownames_to_column(.)
  write_csv(x,paste0(out_dir,name))
}
#target_genes <- 
r_<- read.csv(paste0(args[1],"/tg_wes_intersect.csv"))%>% as.data.frame()
target_genes<- r_$genelist %>% as.list()
names(target_genes) <- target_genes 
target_genes

g_list<-list.dirs(paste0(args[1],"/Mutational_analysis"),recursive = T)[-1]
ind_ <- lapply(target_genes,function(x){
  return(grep(x , g_list))
})
ind_<- ind_%>% keep((. != 0))
#target_genes<- target_genes[which(names(target_genes) %in% names(ind_))]
target_genes<- sort(as.data.frame(ind_)%>%t()%>%.[,1])%>%as.list()
target_genes<- names(target_genes) %>% as.list()
names(target_genes) <- target_genes 
##############
# All reads 
all_tier2_non_ed<- lapply(g_list,FUN = function(x){read.csv(paste0(x,"/allread_id_table_coverTheRegion.csv"))}) 
names(all_tier2_non_ed) <- target_genes

# reads supporting MUTATION
MUT_tier2_non_ed<- lapply(g_list,FUN = function(x){read.csv(paste0(x,"/read_id_support_mutation_MP20_BQ10.csv"))}) 
names(MUT_tier2_non_ed) <- target_genes

# reads supporting WT
WT_tier2_non_ed<- lapply(g_list,FUN = function(x){read.csv(paste0(x,"/read_id_support_ref_MP20_BQ10.csv"))}) 
names(WT_tier2_non_ed) <- target_genes

# mismatches
mis_tier2_non_ed<- lapply(g_list,FUN = function(x){read.csv(paste0(x,"/read_id_mismatches_MP20_BQ10.csv"))}) 
names(mis_tier2_non_ed) <- target_genes

for (i in names(target_genes)){
  if(length(MUT_tier2_non_ed[[i]]$x)<10){
    MUT_tier2_non_ed[[i]]=NULL
  }
  if(length(WT_tier2_non_ed[[i]]$x)<10){
    WT_tier2_non_ed[[i]]=NULL
  }
  if(length(mis_tier2_non_ed[[i]]$x)<10){
    mis_tier2_non_ed[[i]]=NULL
  }
}
names_<- intersect(names(MUT_tier2_non_ed),names(WT_tier2_non_ed))%>% intersect(.,names(mis_tier2_non_ed))
target_genes<- target_genes[which(names(target_genes) %in% names_)]
# MUT reads
MUT <- lapply(MUT_tier2_non_ed, separate, x, 
              into = c("barcode_UMI", "read"), sep = '#', remove = F)
MUT <- lapply(MUT, separate, barcode_UMI, into = c("barcode", "UMI"), sep = '_', remove = F)

# WT  reads
WT <- lapply(WT_tier2_non_ed, separate, x, 
             into = c("barcode_UMI", "read"), sep = '#', remove = F)
WT <- lapply(WT, separate, barcode_UMI, into = c("barcode", "UMI"), sep = '_', remove = F)

# mismatches
mis <- lapply(mis_tier2_non_ed, separate, x, 
              into = c("barcode_UMI", "read"), sep = '#', remove = F)
mis <- lapply(mis, separate, barcode_UMI, into = c("barcode", "UMI"), sep = '_', remove = F)

###########saving all barcodes without filetering
all_BC_bf <- rbind(MUT, WT, mis) %>% as.data.frame()
#bf_srsf2_bc_mut<- count(all_BC_bf$SRSF2$MUT,barcode)

#bf_srsf2_bc_wt<- count(all_BC_bf$SRSF2$WT,barcode)
#bf_srsf2_bc_mis<- count(all_BC_bf$SRSF2$mis,barcode)



#### Check if there are barcode-UMI combinations shared between MUT, WT or mismatch reads 
overlap <- lapply(target_genes, function(i){
  data.frame(
    "MUT_WT" = length(intersect(unique(MUT[[i]]$barcode_UMI), unique(WT[[i]]$barcode_UMI))),
    "MUT_mismatch" = length(intersect(unique(MUT[[i]]$barcode_UMI), unique(mis[[i]]$barcode_UMI))),
    "WT_mismatch" = length(intersect(unique(WT[[i]]$barcode_UMI), unique(mis[[i]]$barcode_UMI)))
  )  
}) %>% data.table::rbindlist() %>% as.data.frame()

rownames(overlap) <- target_genes
knitr::kable(overlap)
save(overlap,"overlap.csv")

#### Obtain list of barcode_UMI combinations to remove
bc_umi_overlap <- lapply(target_genes, function(i){
  list(
    "alt_ref" = intersect(unique(MUT[[i]]$barcode_UMI), unique(WT[[i]]$barcode_UMI)),
    "alt_mismatch" = intersect(unique(MUT[[i]]$barcode_UMI), unique(mis[[i]]$barcode_UMI)),
    "ref_mismatch" = intersect(unique(WT[[i]]$barcode_UMI), unique(mis[[i]]$barcode_UMI))
  )
}) %>% data.table::setattr('names', target_genes)

n_<- length(which(all_BC_bf$SRSF2$MUT$barcode_UMI %in% bc_umi_overlap$SRSF2$alt_ref))
length(which(all_BC_bf$SRSF2$MUT$barcode_UMI %in% bc_umi_overlap$SRSF2$alt_ref))
bf__<- length(MUT$SRSF2$X)
#### Filter out reads whose barcode-UMI combination is shared between MUT, WT and mismatch reads
# MUT
MUT<- lapply(target_genes, function(x){
  bcm<- bc_umi_overlap[[x]]
  MUT[[x]] %>%
    dplyr::filter(!barcode_UMI %in% c(bcm$alt_mismatch, bcm$alt_ref))
  })
names(MUT) <- target_genes
af__<- length(MUT$SRSF2$X)
# WT
WT<- lapply(target_genes, function(x){
  bcm<- bc_umi_overlap[[x]]
  WT[[x]] %>%
  dplyr::filter(!barcode_UMI %in% c(bcm$alt_mismatch, bcm$alt_ref)) })
names(WT) <- target_genes

#### Check if after removing barcode-UMI combination shared between MUT, WT and
#### mismatch reads there are still some overlaps between MUT and WT reads

too <- lapply(target_genes, function(i){
  data.frame(
    "MUT_WT" = length(intersect(unique(MUT[[i]]$barcode_UMI), unique(WT[[i]]$barcode_UMI))) 
  )
}) %>% data.table::rbindlist() %>% as.data.frame()
rownames(too) <- target_genes
knitr::kable(too) # no overlaps left

### Identify and filter out reads whose UMI is associated to multiple barcodes
#### Mutated reads

# identify UMIs associated to multiple barcodes
unique_UMI <- lapply(MUT, function(i){
  unique(i$UMI)
}) %>% data.table::setattr('names', target_genes)

umi_to_multiple_bc <- lapply(target_genes, function(gene){
  n_umi_to_bc <- lapply(unique_UMI[[gene]], function(i){
    x = as.list(MUT[[gene]] %>% dplyr::filter(UMI %in% i))
    y = length(unique(x$barcode))
  })
  unique_UMI[[gene]][which(n_umi_to_bc != 1)]
}) %>% data.table::setattr('names', target_genes)

# count UMIs associated to multiple barcodes
count = NULL
for (i in umi_to_multiple_bc) {
  count = rbind(count, data.frame(
    'UMI_to_multiple_BC' = length(i)
  ))
}
rownames(count) <- target_genes
knitr::kable(count)
save(count,"UMI_to_multiple_BC_MUT.csv")
# filter out
MUT <- lapply(target_genes, function(gene){
  MUT[[gene]] %>% dplyr::filter(!UMI %in% umi_to_multiple_bc[[gene]])
})

#### Wildtype reads

# identify UMIs associated to multiple barcodes
unique_UMI <- lapply(WT, function(i){
  unique(i$UMI)
}) %>% data.table::setattr('names', target_genes)

umi_to_multiple_bc <- lapply(target_genes, function(gene){
  n_umi_to_bc <- lapply(unique_UMI[[gene]], function(i){
    x = as.list(WT[[gene]] %>% dplyr::filter(UMI %in% i))
    y = length(unique(x$barcode))
  })
  unique_UMI[[gene]][which(n_umi_to_bc != 1)]
}) %>% data.table::setattr('names', target_genes)

# count UMIs associated to multiple barcodes
count = NULL
for (i in umi_to_multiple_bc) {
  count = rbind(count, data.frame(
    'UMI_to_multiple_BC' = length(i)
  ))
}
rownames(count) <- target_genes
knitr::kable(count)
save(count,"UMI_to_multiple_BC_WT.csv")
# filter out
WT <- lapply(target_genes, function(gene){
  WT[[gene]] %>% dplyr::filter(!UMI %in% umi_to_multiple_bc[[gene]])
})

#### Mismatch reads

# identify UMIs associated to multiple barcodes
unique_UMI <- lapply(mis, function(i){
  unique(i$UMI)
}) %>% data.table::setattr('names', target_genes)

umi_to_multiple_bc <- lapply(target_genes, function(gene){
  n_umi_to_bc <- lapply(unique_UMI[[gene]], function(i){
    x = as.list(mis[[gene]] %>% dplyr::filter(UMI %in% i))
    y = length(unique(x$barcode))
  })
  unique_UMI[[gene]][which(n_umi_to_bc != 1)]
}) %>% data.table::setattr('names', target_genes)

# count UMIs associated to multiple barcodes
count = NULL
for (i in umi_to_multiple_bc) {
  count = rbind(count, data.frame(
    'UMI_to_multiple_BC' = length(i)
  ))
}
rownames(count) <- target_genes
knitr::kable(count)
save(count,"UMI_to_multiple_BC_Mism.csv")
# filter out
mis <- lapply(target_genes, function(gene){
  mis[[gene]] %>% dplyr::filter(!UMI %in% umi_to_multiple_bc[[gene]])
})

# Demultiplexing by barcodes
## Derive list of barcodes in ONT data, uniquely associated to MUT, WT or mismatch reads
#MUT

MUT_BC <- lapply(target_genes, function(i){
  x<- as.data.frame(MUT[[i]]$barcode)
  x%>%mutate(gene= i)
})

MUT_BC <- data.table::rbindlist(MUT_BC) %>% as.data.frame() %>% 
  mutate(support = 'MUT') %>% rename(barcode = 'MUT[[i]]$barcode')

#REF
WT_BC <- lapply(target_genes, function(i){
  x<- as.data.frame(WT[[i]]$barcode)
  x%>%mutate(gene= i)
}) 
WT_BC <- data.table::rbindlist(WT_BC) %>% as.data.frame() %>% 
  mutate(support = 'WT') %>% rename(barcode = 'WT[[i]]$barcode')

#mismatch
mis_BC <- lapply(target_genes, function(i){
  x<- as.data.frame(mis[[i]]$barcode)
  x%>%mutate(gene= i)
}) 


mis_BC <- data.table::rbindlist(mis_BC) %>% as.data.frame() %>% 
  mutate(support = 'mis') %>% rename(barcode = 'mis[[i]]$barcode')

# unique list of barcodes, with read counts by gene-read type (colnames) and barcode (rownames)
all_BC <- rbind(MUT_BC, WT_BC, mis_BC) %>% as.data.frame()

length(unique(all_BC$barcode)) 

all_BC <- all_BC %>% pivot_wider(
  names_from = c(gene, support),
  values_from = gene,
  values_fill = 0,
  values_fn = length 
)
str(all_BC)
for (i in 1:length(target_genes)){
  all_BC<- relocate(all_BC,paste0(target_genes[i],"_WT"), .after = paste0(target_genes[i],"_MUT"))
}
all_BC$barcode <- str_remove_all(string =all_BC$barcode,pattern = '\"')
knitr::kable(head(all_BC))
length(all_BC$barcode) # 3031
save(all_BC,"all_BC.csv")




### Check correspondence between ONT barcodes and list from 10x raw unfiltered matrix

df_BC <- data.frame(
  "bc_10x" = length(raw_matrix_barcodes$x), 
  "bc_ONT" = length(all_BC$barcode), 
  "bc_ONT_and_10x" = length(intersect(raw_matrix_barcodes$x, all_BC$barcode)), 
  "bc_ONT_not_10x" = length(which(!all_BC$barcode %in% raw_matrix_barcodes$x)) 
) %>% t()

colnames(df_BC)<- c("Num_BC")
knitr::kable(df_BC)
save(df_BC,"overlapped_BC_ONT_10X.csv")

### Relationship between expression of 10x and ONT reads   
#### Load and check Seurat object

for (i in 1:length(target_genes)){
  if (target_genes[i] %in% colnames(AML@meta.data)){
    AML[[names(target_genes[i])]]<-NULL
  }
}


#### Coverage by 10x and ONT reads (estimate technical dropout)
ex_mat <- Seurat::GetAssayData(AML, slot = "counts") 
data_<- Seurat::GetAssayData(AML, slot = "data") 
dim(ex_mat)
#str(data_)
target_genes<- target_genes[which(target_genes %in% data_@Dimnames[[1]] )]
features <- unlist(target_genes)
plot_x<- function(x){
  Seurat::FeaturePlot(AML, dims = c(1, 2), reduction = "UMAP", features = x, cols = c("lightgrey", custom_colors[1]),
                      pt.size = 0.3, order = T, keep.scale = "all") & theme(legend.text = element_text(size = 10), 
                     axis.text = element_text(size = 10), axis.title = element_text(size = 10)) & ggtitle(x)  }
feature_plt<- lapply(target_genes, plot_x)

lapply(names(feature_plt), 
       function(x)ggsave(filename=paste(out_dir,"/",x,".jpeg",sep=""), plot=feature_plt[[x]]))

###########

# 10x read counts by cell and gene   #####number of reads or cells in the output table???

short_counts <- lapply(target_genes, function(i){
  ex_mat[(rownames(ex_mat) %in% i),]
}) %>% data.table::setattr('names', target_genes) %>% 
  as.data.frame() %>%
  rownames_to_column('barcode') %>%
  separate('barcode', into = c('barcode', 'sample'), remove = T) %>%
  dplyr::select(-sample) 
colnames(short_counts)[2:ncol(short_counts)] <- paste("short", target_genes, sep = '_') 

# 10x read counts across all cells by mutated gene 
SHORT <- lapply(target_genes, function(i){
  sum(ex_mat[(rownames(ex_mat) %in% i),])
}) %>% unlist() %>% data.table::setattr('names', target_genes) %>%
  rbind() %>% as.data.frame() %>% mutate('N cells' = length(colnames(ex_mat)))
rownames(SHORT) <- 'SHORT'

save(SHORT,"read_counts_across_all_cells_short.csv")



# ONT read counts by cell and gene
ONT_counts <- lapply(target_genes, function(i){
  mutate(all_BC, ONT = rowSums(dplyr::select(all_BC, contains(i)))) 
}) %>% as.data.frame()  %>% dplyr::select(contains(c('barcode', 'ONT'))) %>%
  dplyr::select(1,contains("ONT")) %>% dplyr::rename(barcode= paste0(target_genes[1],".barcode"))

# ONT read counts across all cells by mutated gene  
ONT_ex_mat <- ONT_counts %>%
  column_to_rownames(var = 'barcode' )%>% t()
rownames(ONT_ex_mat) <- target_genes


LONG_enriched <- lapply(target_genes, function(i){
  sum(ONT_ex_mat[(rownames(ONT_ex_mat) %in% i),])
}) %>% unlist() %>% data.table::setattr('names', target_genes) %>%
  rbind() %>% as.data.frame() %>% mutate('N cells' = length(ONT_counts$barcode))
rownames(LONG_enriched) <- 'LONG_enriched'
knitr::kable(rbind(SHORT, LONG_enriched))
x<- rbind(SHORT, LONG_enriched) %>% as.data.frame(.)
save(x,"short_long_bc_genes.csv")

####saving 
short_ONT_counts <- inner_join(short_counts, ONT_counts, by = "barcode") 
nrow(short_ONT_counts)

saveRDS(short_ONT_counts, paste0(out_dir,"/short_ONT_counts.rds"))



#### Filter by gene

# FLT3
# select barcodes with concordant expression
concordant_ <- lapply(target_genes,function(x){
  short_ONT_counts %>% 
  dplyr::select(barcode, contains(x)) %>% 
  dplyr::filter(across(contains("short"), ~ . >= 0) & across(contains("ONT"), ~ . > 0))
  })

# create ONT dataset with concordant expression
df_ <- lapply(target_genes,function(x){
  all_BC[which(all_BC$barcode %in% concordant_[[x]]$barcode), c('barcode', paste0(x,'_MUT'), paste0(x,'_WT'), paste0(x,'_mis'))]
} )
########
#### Plot distribution
thr<- lapply(target_genes,function(x){
  df_[[x]] %>% pivot_longer(cols = 2:ncol(df_[[x]]), names_to = "read", values_to = "N") %>% dplyr::filter(N != '0') 
})
#dim(df_$ABL1)
for (x in target_genes){
  if (dim(thr[[x]])[1]!=0){
  thr[[x]]$read <- str_replace(thr[[x]]$read,paste0(x,"_MUT"),  'mutated')
  thr[[x]]$read <- str_replace(thr[[x]]$read,paste0(x,'_WT'),  'wild_type')
  thr[[x]]$read <- str_replace(thr[[x]]$read,paste0(x,'_mis'),  'mismatch')
  thr[[x]]$read <- factor(thr[[x]]$read, levels =  c('mutated', 'wild_type', 'mismatch'))
  } else {thr[[x]]=NULL
  target_genes[[x]]=NULL}
}

ditribution_plot<- lapply(target_genes,function(x1){
  ggplot(thr[[x1]]) + 
    geom_density(aes(x = N, fill = read, color = read)) +
    geom_rug(aes(x = N, y = 0, color = read), position = position_jitter(height = 0)) +
    theme_bw() +
    xlim(min(thr[[x1]]$N), max(thr[[x1]]$N)) +
    labs(x = 'N reads', y = 'N cells (density)') +
    ggtitle(x1) +
    scale_fill_manual(values = c('darkred', 'navy', 'lightgrey')) +
    scale_color_manual(values = c('darkred', 'navy', 'lightgrey')) +
    theme(plot.title.position = "panel",
          plot.title = element_text(face = "bold", size = 20), 
          strip.text = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 13), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none') +
    facet_wrap(~read, scales = 'free_y')
})
lapply(names(ditribution_plot), 
       function(x)ggsave(filename=paste(out_dir,"/",x,"_distribution.jpeg",sep=""), plot=ditribution_plot[[x]]))
for (i in target_genes) {
  base::colnames(concordant_[[i]])<- c("barcode","Short","ONT")
}

#### Correlation between number of 10x reads and number of ONT reads by gene
concordant_<- concordant_[which(names(concordant_) %in% names(target_genes))]
cor_<- lapply(target_genes,function(x) {cor__<- cor(x = concordant_[[x]]$ONT, y = concordant_[[x]]$Short ,
                use = "pairwise.complete.obs", method = "pearson") 
              ggplot(concordant_[[x]], aes(x = ONT , y =Short)) +
                geom_point() +
                geom_smooth(method = 'lm') +
                theme_bw() +
                xlab('ONT, N reads') +
                ylab('10x, N reads') +
                labs(title = x,
                     subtitle = paste('R =', round(cor__, digits = 2), '(Pearson)', sep = " ")) +
                theme(plot.title = element_text(face = "bold", size = 20), 
                      plot.subtitle = element_text(face = "plain", size = 18),  
                      axis.text = element_text(size = 11),
                      axis.title = element_text(size = 13), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
              })

lapply(names(cor_), 
       function(x)ggsave(filename=paste(out_dir,"/",x,"_correlation.jpeg",sep=""), plot=cor_[[x]]))




df_bc_all<- data_frame(barcode="barcode")
for (i in target_genes){
  df_bc_all<- full_join(df_[[i]], df_bc_all, by = "barcode")
}
df_bc_all[is.na(df_bc_all)]=0


nrow(df_bc_all) 

saveRDS(df_bc_all,paste0(out_dir, 'df_bc_all.rds'))
write.csv(df_bc_all,paste0(out_dir, 'df_bc_all.csv'))
#write.csv(df_bc_all, paste0(out_dir,'df_bc_all.csv')) # input for Tirelli's algorithm

################### plot stack barplot
#reads_bf_demultiplexing<-read.csv('/hpcnfs/scratch/PGP/niman/Chiara/FLAMES/sAML1_B/Long/enrichment_analysis/merged/Mutation/tot_reads_bf_demultiplexing.csv',stringsAsFactors = F)
df_bc_all<- read.csv(paste0(out_dir,"df_bc_all.csv"))[-1]
reads_af_demultiplexing_ <- data.frame(N_reads=colSums(df_bc_all[-1])) %>% dplyr::mutate(Read_gp=str_split(row.names(.),pattern  = "_",simplify = T)[,2],Gene=str_split(row.names(.),pattern  = "_",simplify = T)[,1])
reads_af_demultiplexing_$Read_gp<- factor(reads_af_demultiplexing_$Read_gp,levels = c("MUT","WT","mis"))
#####reads just after demultiplexing

reads_af_demultiplexing_<-reads_af_demultiplexing_%>%
  group_by(Gene)%>%dplyr::mutate(pct = round(N_reads / sum(N_reads) * 100,2))
Tot<-reads_af_demultiplexing_%>%group_by(Gene)%>%summarize(total = sum(N_reads))%>%dplyr::mutate(pct=100)



a<-ggplot(reads_af_demultiplexing_, aes(x = Gene, y = pct,fill=Read_gp)) + 
  geom_bar(stat = 'identity', position = 'stack')+scale_fill_manual("Read_gp", values = c("MUT" = "darkblue", "WT" = "orange", "mis" = "azure"))+
  theme_dark()+theme(legend.position = "bottom",legend.title = element_blank()) +scale_y_continuous(name = 'Percentage [%]',labels = scales::percent_format(scale = 1))+
  geom_text(aes(Gene, label = total, fill = NULL), data = Tot,vjust=-0.1,col="white")
ggsave(filename=paste(out_dir,"N_reads_PCT.jpeg",sep=""), plot=a)



######### intersection of barcodes for srsf2 after and before filtering
# 
# 
# #######
# after_bc<- df_bc_all%>% dplyr::select(barcode,SRSF2_MUT) 
# after_bc<- after_bc[which(rowSums(after_bc[2])!=0),]
# len_bf<- length(bf_srsf2_bc_mut$barcode)
# len_af<-length(unique(after_bc$barcode))
# inter_sect_in<- which(bf_srsf2_bc_mut$barcod %in% intersect(bf_srsf2_bc_mut$barcode,after_bc$barcode))
# com_bc=bf_srsf2_bc_mut$barcode[inter_sect_in]
# inter_sect_out<- which(bf_srsf2_bc_mut$barcod %in% setdiff(bf_srsf2_bc_mut$barcod,after_bc$barcode))
# bc_bf<- length(unique(bf_srsf2_bc_mut$barcode))
# bc_af<- length(unique(after_bc$barcode))
# read_num_bf<- sum(bf_srsf2_bc_mut$n)
# read_num_af<-  sum(after_bc$SRSF2_MUT)
# read_num_diff<- sum(bf_srsf2_bc_mut$n[inter_sect_out])
# read_num_bf
# read_num_af
# read_num_diff
# bf_srsf2_bc_mut<- bf_srsf2_bc_mut%>% arrange(barcode)
# after_bc<- after_bc%>% arrange(barcode)
# bf_<- bf_srsf2_bc_mut[which(bf_srsf2_bc_mut$barcode %in% com_bc),]
# identical(bf_$barcode,after_bc$barcode)
# mut_df<- data_frame(reads_n_bf=bf_$n,reads_n_af=after_bc$SRSF2_MUT)
# row.names(mut_df) <- bf_$barcode
# mut_df<- rowid_to_column(mut_df)
# #############
# after_bc<- df_bc_all%>% dplyr::select(barcode,SRSF2_WT) 
# after_bc<- after_bc[which(rowSums(after_bc[2])!=0),]
# len_bf<- length(bf_srsf2_bc_wt$barcode)
# len_af<-length(unique(after_bc$barcode))
# inter_sect_in<- which(bf_srsf2_bc_wt$barcod %in% intersect(bf_srsf2_bc_wt$barcode,after_bc$barcode))
# com_bc=bf_srsf2_bc_wt$barcode[inter_sect_in]
# inter_sect_out<- which(bf_srsf2_bc_wt$barcod %in% setdiff(bf_srsf2_bc_wt$barcod,after_bc$barcode))
# bc_bf<- length(unique(bf_srsf2_bc_wt$barcode))
# bc_af<- length(unique(after_bc$barcode))
# read_num_bf<- sum(bf_srsf2_bc_wt$n)
# read_num_af<-  sum(after_bc$SRSF2_WT)
# read_num_diff<- sum(bf_srsf2_bc_wt$n[inter_sect_out])
# read_num_bf
# read_num_af
# read_num_diff
# bf_srsf2_bc_wt<- bf_srsf2_bc_wt%>% arrange(barcode)
# after_bc<- after_bc%>% arrange(barcode)
# bf_<- bf_srsf2_bc_wt[which(bf_srsf2_bc_wt$barcode %in% com_bc),]
# identical(bf_$barcode,after_bc$barcode)
# wt_df<- data_frame(reads_n_bf=bf_$n,reads_n_af=after_bc$SRSF2_WT)
# row.names(wt_df) <- bf_$barcode
# wt_df<- rowid_to_column(wt_df)
# #########
# after_bc<- df_bc_all%>% dplyr::select(barcode,SRSF2_mis) 
# after_bc<- after_bc[which(rowSums(after_bc[2])!=0),]
# len_bf<- length(bf_srsf2_bc_mis$barcode)
# len_bf
# len_af<-length(unique(after_bc$barcode))
# len_af
# inter_sect_in<- which(bf_srsf2_bc_mis$barcod %in% intersect(bf_srsf2_bc_mis$barcode,after_bc$barcode))
# com_bc=bf_srsf2_bc_mis$barcode[inter_sect_in]
# inter_sect_out<- which(bf_srsf2_bc_mis$barcod %in% setdiff(bf_srsf2_bc_mis$barcod,after_bc$barcode))
# bc_bf<- length(unique(bf_srsf2_bc_mis$barcode))
# bc_af<- length(unique(after_bc$barcode))
# read_num_bf<- sum(bf_srsf2_bc_mis$n)
# read_num_af<-  sum(after_bc$SRSF2_mis)
# read_num_diff<- sum(bf_srsf2_bc_mis$n[inter_sect_out])
# read_num_bf
# read_num_af
# read_num_diff
# bf_srsf2_bc_mis<- bf_srsf2_bc_mis%>% arrange(barcode)
# after_bc<- after_bc%>% arrange(barcode)
# bf_<- bf_srsf2_bc_mis[which(bf_srsf2_bc_mis$barcode %in% com_bc),]
# identical(bf_$barcode,after_bc$barcode)
# mis_df<- data_frame(reads_n_bf=bf_$n,reads_n_af=after_bc$SRSF2_mis)
# row.names(mis_df) <- bf_$barcode
# mis_df<- rowid_to_column(mis_df)
# #######################
# df <- reshape2::melt(mis_df ,  id.vars = 'rowid', variable.name = 'reads')
# x<- ggplot(df, aes(rowid, value)) +
#   geom_line(aes(colour = reads))
# 
# df <- reshape2::melt(wt_df ,  id.vars = 'rowid', variable.name = 'reads')
# x1<- ggplot(df, aes(rowid, value)) +
#   geom_line(aes(colour = reads))
# 
# df <- reshape2::melt(mut_df ,  id.vars = 'rowid', variable.name = 'reads')
# x2<- ggplot(df, aes(rowid, value)) +
#   geom_line(aes(colour = reads))
# x3<- x+x1+x2
# ggsave("/hpcnfs/scratch/PGP/SCMseq/aml5_promethion/merged_enriched/srsf2_comp.jpg",plot=x3,width = 15,height = 5)
