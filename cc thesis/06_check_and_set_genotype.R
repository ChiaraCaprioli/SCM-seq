# Check single-cell genotyping after smoothing 
# chiara caprioli
# Feb 01 2023

############################################
# SET UP

## Load packages
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(gmodels)
library(Seurat)

## Paths
path_main <- '/hpcnfs/scratch/PGP/SCMseq/sc_genotype/' 
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/')

## Define samples
samples <- c('AML4', 'AML5', 'sAML1')

## Set colors
colors <- list()

colors$sample <- setNames(
  c("#227093", "#ffb142", "#485460"),
  c("AML4", "AML5", "sAML1")
)

colors$genotype <- setNames(
  c('darkred', 'steelblue', 'lightgrey'),
  c('mutated', 'wild-type', 'mismatch'))

colors$col_keep <- setNames(
  c("#FEB24C", "#525252"),
  c('included', 'discarded'))

# EXPRESSION DROPOUT ############################################

## Percent dropout
### sample
df_sample <- data.frame()
for (s in samples) {
  dropout.matrix <- read_csv(paste0('/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/',s,'_dropout.matrix_modified.csv'))
  dropout.matrix <- dropout.matrix %>% column_to_rownames(var = '...1') %>% as.matrix()
  
  count_drop <- sum(matrixStats::rowCounts(dropout.matrix, value = 1)) / 
    (nrow(dropout.matrix)*ncol(dropout.matrix)) # note: '1' = dropout
  
  x <- data.frame(
    sample = s,
    dropout_percentage = count_drop
  )
  
  df_sample <- rbind(x, df_sample)
  
}

write_csv(df_sample, paste0(path_results, 'tables/df_sample_dropout.csv'))

### gene
df_gene <- data.frame()
for (s in samples) {
  dropout.matrix <- read_csv(paste0('/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/',s,'_dropout.matrix_modified.csv'))
  dropout.matrix <- dropout.matrix %>% column_to_rownames(var = '...1') %>% as.matrix()
  
  count_drop <- as.data.frame(rowSums(dropout.matrix) / ncol(dropout.matrix))
  colnames(count_drop) <- 'dropout_percentage'
  count_drop <- count_drop %>% rownames_to_column(var = 'gene')
  count_drop$sample <- s
    
  df_gene <- rbind(count_drop, df_gene)
  
}

write_csv(df_gene, paste0(path_results, 'tables/df_gene_dropout.csv'))

## Plot dropout 

### bar plot

p1 <- ggplot(df_sample, aes(x = sample, y = dropout_percentage)) +
  geom_col(width = 0.75) +
  ylab('Expression dropout (%)') +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0(path_results, 'plots/sample_dropout.png'), p1,
       width = 5, height = 4)

for (s in samples) {
  df_gene_sub <- df_gene %>% filter(sample == s) %>% arrange(desc(dropout_percentage))
  df_gene_sub$gene <- factor(df_gene_sub$gene, levels = df_gene_sub$gene)
  
  p2 <- ggplot(df_gene_sub, aes(x = gene, y = dropout_percentage)) +
    geom_col(width = 0.75) +
    ylab('Expression dropout (%)') +
    ylim(0,1) +
    ggtitle(s) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank())
  
  ggsave(paste0(path_results, 'plots/', s, '_gene_dropout.png'), p2,
         width = 5, height = 4)
}

### heatmap
for (s in samples) {
  dropout.matrix <- read_csv(paste0('/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/',s,'_dropout.matrix_modified.csv'))
  dropout.matrix <- dropout.matrix %>% column_to_rownames(var = '...1') %>% as.matrix()
  
  pheatmap(dropout.matrix, cluster_cols = F, cluster_rows = F, 
           show_colnames = F, show_rownames = T, scale = "none",
           legend_breaks = c(0, 1), 
           legend_labels = c('no dropout', 'dropout'), 
           color = c('#525252', '#F0F0F0'),
           main = paste0(s, ' - Expression dropout'),
           cellwidth = 0.1, cellheight = 10, 
           filename = paste0(path_results, 'plots/', s, '_drop_counts_heat.png'),
           width = 7, height = 2)
}

# MERGE SMOOTHED AND UNSMOOTHED DATA ############################################
first_geno_all <- read_csv(paste0(path_results, 'tables/', 'first_geno_all.csv'))

genotype_all <- data.frame()
for (s in samples) {
  first_geno_all_sample <- first_geno_all %>% filter(sample == s) 
  
  geno.matrix <- read_csv(paste0('/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/',s,'_geno.mat_modified.csv'))
  geno.matrix <- geno.matrix %>% rename('gene' = '...1') %>%
    pivot_longer(cols = 2:ncol(.), names_to = 'barcode', values_to = 'smoothed_genotype') 
  geno.matrix <- geno.matrix %>%  mutate(smoothed_genotype = recode(smoothed_genotype, 'MUT' = 'mutated', 'WT' = 'wild-type'))
  
  c.matrix <- read_csv(paste0('/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/',s,'_c.mat_modified.csv'))
  c.matrix <- c.matrix %>% rename('gene' = '...1') %>%
    pivot_longer(cols = 2:ncol(.), names_to = 'barcode', values_to = 'smoothed_confidence') 
  
  df_ <- full_join(geno.matrix, c.matrix)
  
  geno_all_sample <- left_join(first_geno_all_sample, df_, by = c('barcode', 'gene'))
  
  genotype_all <- rbind(geno_all_sample, genotype_all)
}

write.csv(genotype_all, paste0(path_results, 'tables/genotype_all.csv'))

# RE-CLASSIFICATION RATE ############################################

## confusion matrix to assess re-classification of cells imputed with high confidence with initial algorithm
for (s in samples) {
  expected <- genotype_all$genotype[which(genotype_all$confidence >= 0.75 & genotype_all$sample == s)]
  observed <- genotype_all$smoothed_genotype[which(genotype_all$confidence >= 0.75 & genotype_all$sample == s)]
  
  conf_mat <- CrossTable(expected, observed, dnn = c('no smoothing', 'smoothing'), chisq = T)
  
  pdf(paste0(path_results, 'plots/', s, '_reclass_highC_genotype_confusion_matrix.pdf'))
  ht = Heatmap(conf_mat$prop.tbl, 
          name = 'Proportion',
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", conf_mat$prop.tbl[i, j]), x, y, 
                      gp = gpar(font = 'bold', fontsize = 13, col = "red"))
          },
          cluster_rows = F,
          cluster_columns = F,
          row_title = 'No smoothing',
          row_title_gp = gpar(fontsize = 13.2),
          column_title = 'Smoothing',
          column_title_gp = gpar(fontsize = 13.2),
          column_title_rot = 0,
          col = viridisLite::mako(n = 100, direction = 1),
          heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1)),
          rect_gp = gpar(col = "white", lwd = 0.5),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          width = unit(8, "cm"), height = unit(8, "cm")
  )
  draw(ht)
  dev.off()
}

# COMPARE CONFIDENCE BEFORE AND AFTER SMOOTHING ############################################
for (s in samples) {
  df_ <- genotype_all %>% 
    filter(sample == s) %>%
    pivot_longer(cols = c('confidence', 'smoothed_confidence'), 
                 names_to = 'method', values_to = 'confidence') 
  df_ <- df_  %>% 
    mutate(method = as.factor(
      recode(method, 'confidence' = 'no smoothing', 'smoothed_confidence'='smoothing')) ) %>%
    arrange(desc(confidence))
  df_$barcode <- factor(df_$barcode, levels = unique(df_$barcode))
  
  p1 <- df_ %>% 
    filter(method == 'no smoothing') %>%
    ggplot(aes(x = barcode, y = confidence, fill = genotype)) +
    geom_col() +
    theme_bw() +
    ylim(0,1) +
    ggtitle('No smoothing') +
    scale_fill_manual(values = colors$genotype) +
    facet_wrap(~gene+genotype, ncol = 2) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(face = 'bold', size = 15),
          legend.position = 'none',
          plot.title = element_text(face = 'bold', size = 17, hjust = 0.5, vjust = 2))
  
  p2 <- df_ %>% 
    filter(method == 'smoothing') %>%
    ggplot(aes(x = barcode, y = confidence, fill = genotype)) +
    geom_col() +
    theme_bw() +
    ylim(0,1) +
    ggtitle('Smoothing') +
    scale_fill_manual(values = colors$genotype) +
    facet_wrap(~gene+genotype, ncol = 2) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(face = 'bold', size = 15),
          legend.position = 'none',
          plot.title = element_text(face = 'bold', size = 17, hjust = 0.5, vjust = 2))
  
  ggsave(paste0(path_results, 'plots/', s, '_confidence_density.png'), p1 + p2 +
           plot_layout(ncol = 2, widths = 10, heights = 25), width = 20, height = 25)
}

for (s in samples) {
  df_ <- genotype_all %>% 
    filter(sample == s) %>%
    group_by(gene) %>%
    summarise(mean_confidence_no_smoothing = mean(confidence),
              mean_confidence_smoothing = mean(smoothed_confidence)) 
  df_ <- df_ %>%
    pivot_longer(cols = c('mean_confidence_no_smoothing', 'mean_confidence_smoothing'), 
                 names_to = 'method', values_to = 'mean_confidence')
  
  p <- ggplot(df_, aes(gene, mean_confidence, fill = method)) +
    geom_col(position = 'dodge', width = 0.75) +
    theme_bw() +
    ggtitle(s) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          plot.title = element_text(face = 'bold', size = 17, hjust = 0.5, vjust = 2))
  
  ggsave(paste0(path_results, 'plots/', s, '_mean_confidence_method.png'), p,
           width = 6, height = 4.5)
  
}

# SUBSET BASED ON CONFIDENCE ############################################
genotype_all_final <- 
  genotype_all[which(genotype_all$confidence >= 0.75),]

write_csv(genotype_all_final, paste0(path_results, 'tables/genotype_all_final.csv'))

#a = nrow(genotype_all[which(genotype_all$confidence >= 0.75 & genotype_all$smoothed_confidence >= 0.75),]) #11361
#b = nrow(genotype_all[which(genotype_all$confidence >= 0.75 & genotype_all$smoothed_confidence < 0.75),]) #3637
#c = nrow(genotype_all[which(genotype_all$confidence < 0.75 & genotype_all$smoothed_confidence < 0.75),]) #3015
#d = nrow(genotype_all[which(genotype_all$confidence < 0.75 & genotype_all$smoothed_confidence >= 0.75),]) #1
#sum(a,b,c,d) == nrow(genotype_all) #TRUE

# Check density of known and missing, before and after smoothing ############################################

df_gain <- list()
for (s in samples) {
  data_unsmooth <- genotype_all %>%
    filter(sample == s & confidence >= 0.75) %>%
    dplyr::select(barcode, gene, genotype) %>%
    mutate(genotype = recode(genotype, 'mutated' = '1', 'wild-type' = '0'))
  data_unsmooth$genotype <- as.numeric(data_unsmooth$genotype)
  data_unsmooth <- data_unsmooth %>% 
    pivot_wider(names_from = gene, values_from = genotype, values_fn = mean) %>%
    column_to_rownames(var = 'barcode') %>% t()
  data_unsmooth[is.na(data_unsmooth)] <- '2'
  perc_unkn_unsmooth <- length(which(data_unsmooth == '2')) / (nrow(data_unsmooth)*ncol(data_unsmooth))
  
  data_smooth <- genotype_all_final %>%
    filter(sample == s & smoothed_confidence >= 0.75) %>%
    dplyr::select(barcode, gene, genotype) %>%
    mutate(genotype = recode(genotype, 'mutated' = '1', 'wild-type' = '0'))
  data_smooth$genotype <- as.numeric(data_smooth$genotype)
  data_smooth <- data_smooth %>% 
    pivot_wider(names_from = gene, values_from = genotype, values_fn = mean) %>%
    column_to_rownames(var = 'barcode') %>% t()
  data_smooth[is.na(data_smooth)] <- '2'
  perc_unkn_smooth <- length(which(data_smooth == '2')) / (nrow(data_smooth)*ncol(data_smooth))
  
  df_ <- data.frame(
    unknown_before_smoothing = perc_unkn_unsmooth,
    unknown_after_smoothing = perc_unkn_smooth
  ) 
  df_gain[[s]] <- df_
}

x <- do.call(rbind, df_gain)

write_csv(x, paste0(path_results, 'tables/unknown_before_after_smoothing.csv'))

# Plot genotype heatmap ############################################

for (s in samples) {
  
  dat_ <- genotype_all_final %>%
    filter(sample == s) %>%
    dplyr::select(barcode, gene, genotype) %>%
    mutate(genotype = recode(genotype, 'mutated' = '1', 'wild-type' = '0'))
  dat_$genotype <- as.numeric(dat_$genotype)
  
  dat_ <- dat_ %>% 
    pivot_wider(names_from = gene, values_from = genotype, values_fn = mean) %>%
    column_to_rownames(var = 'barcode') %>% t()
  dat_[is.na(dat_)] <- '2'
  
  dat_ <- dat_ %>% t() %>% as.data.frame() 
  dat_$n_unknown <- rowSums(dat_ == "2")
  
  dat_ <- dat_ %>% filter(n_unknown <= (ncol(dat_)-1)) %>%
    arrange(n_unknown) %>%
    dplyr::select(-n_unknown) %>% t() %>% as.data.frame()
  
  dat_$n_mut <- rowSums(dat_ == "1")
  dat_ <- dat_ %>%
    arrange(desc(n_mut)) %>%
    dplyr::select(-n_mut) %>% as.matrix()
  mode(dat_) <- 'numeric' 
  
  cols_h <- setNames(
    c('steelblue', 'darkred', '#F0F0F0'),
    c('0', '1', '2')
  )
  
  pheatmap::pheatmap(dat_, color = cols_h, cluster_rows = F, cluster_cols = F, 
                     show_colnames=F, fontsize = 6, legend = T,
                     legend_breaks = c(0,1,2), legend_labels=c('wild-type', 'mutated', 'unknown'),
                     cellheight = 10, cellwidth = 0.1, main = paste0(s, ', n cells = ', ncol(dat_)),
                     filename = paste0(path_results, 'plots/', s, '_genotype_heatmap.png'),
                     width = 7, height = 4)
  
}

# Correlation between VAF (from WES) and mutant cell fraction (from SCM-seq) ############################################
genotype_all_final <- genotype_all_final %>% 
  pivot_longer(cols = c(n_reads_mut, n_reads_wt, mis), names_to = 'read_type', values_to = 'n_reads_type')

df_vaf_scVAF <- list()
for (s in samples) {
  
  dat_ = genotype_all_final %>% filter(sample == s)
  gene = unique(dat_$gene)
  
  df <- data.frame()
  for (g in gene) {
    
    x = dat_ %>% 
      filter(gene == g & read_type != 'mis') %>% 
      dplyr::select(barcode, sample, gene, read_type, n_reads_type) %>%
      group_by(sample, gene, read_type) %>% 
      summarise(N_tot_reads = sum(n_reads_type))
    
    df = rbind(x, df)
    df_vaf_scVAF[[s]] <- df 
    
  }
}

df_vaf_scVAF <- do.call(rbind, df_vaf_scVAF)

df_vaf_scVAF <- df_vaf_scVAF %>% 
  mutate(read_type = recode(read_type, 'mutated'='n_reads_mut', 'wild-type' = 'wt')) %>%
  pivot_wider(names_from = read_type, values_from = N_tot_reads, values_fill = 0)

colnames(df_vaf_scVAF)[3] <- 'MUT_N_tot_reads'
colnames(df_vaf_scVAF)[4] <- 'WT_N_tot_reads'

#

df_vaf_mcf <- list()

for (s in samples) {
  
  data = genotype_all_final %>% filter(sample == s)
  gene = unique(data$gene)
  
  df <- data.frame()
  
  for (g in gene) {
    
    x = data %>% 
      filter(gene == g) %>% 
      distinct(barcode, .keep_all=T) %>%
      dplyr::select(barcode, sample, gene, genotype) %>%
      group_by(sample, gene, genotype) %>% 
      summarise(count = n())
    x = x %>% mutate(N_tot_cells = data %>% distinct(barcode) %>% nrow())
    
    df = rbind(x, df)
    df_vaf_mcf[[s]] <- df 
    
  }
}

df_vaf_mcf <- do.call(rbind, df_vaf_mcf)

df_vaf_mcf <- df_vaf_mcf %>% 
  pivot_wider(names_from = genotype, values_from = count, values_fill = 0)

colnames(df_vaf_mcf)[4] <- 'MUT_N_tot_cells'
colnames(df_vaf_mcf)[5] <- 'WT_N_tot_cells'

for (i in df_vaf_mcf$gene){
  N_cells_UNK = (df_vaf_mcf$N_tot_cells - (df_vaf_mcf$MUT_N_tot_cells + df_vaf_mcf$WT_N_tot_cells))
}

df_vaf_mcf$N_cells_UNK <- N_cells_UNK

df_vaf_mcf <- full_join(df_vaf_mcf, df_vaf_scVAF, by = c('sample', 'gene'))

#

vaf <- read_delim('/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/mutations_3_AML_VAF.csv', delim = ';')
vaf <- vaf %>% filter(Gene %in% df_vaf_mcf$gene) %>%
  dplyr::select(Sample, Gene, tumor_f)

df_vaf_mcf <- full_join(df_vaf_mcf, vaf, by = c("gene" = "Gene", 'sample' = 'Sample'))

# Mutant cell fraction (mcf) defined for each gene as the ratio between mutated cells and total number of genotyped cells
for (i in df_vaf_mcf$gene){
  mcf = df_vaf_mcf$MUT_N_tot_cells/(df_vaf_mcf$MUT_N_tot_cells + df_vaf_mcf$WT_N_tot_cells)
  scVAF = df_vaf_mcf$MUT_N_tot_reads/(df_vaf_mcf$MUT_N_tot_reads + df_vaf_mcf$WT_N_tot_reads) 
}

df_vaf_mcf <- df_vaf_mcf %>%
  cbind(mcf = mcf, scVAF = scVAF)


write_csv(df_vaf_mcf, paste0(path_results, 'tables/df_vaf_mcf.csv'))

# VAF vs MCF ############################################
cor_mcf <- cor.test(df_vaf_mcf$mcf, df_vaf_mcf$tumor_f, use = "everything", method = "pearson")
cor_scVAF <- cor.test(df_vaf_mcf$scVAF, df_vaf_mcf$tumor_f, use = "everything", method = "pearson")
tab_cor_all <- data.frame(
  sample = 'all',
  test = c('mcf_vaf', 'scVaf_vaf'),
  r = c(cor_mcf$estimate, cor_scVAF$estimate),
  p = c(cor_mcf$p.value, cor_scVAF$p.value)
)

tab_core_sample <- list()
for (s in samples) {
  data = df_vaf_mcf %>% filter(sample == s)
  cor_mcf <- cor.test(df_vaf_mcf$mcf, df_vaf_mcf$tumor_f, use = "everything", method = "pearson")
  cor_scVAF <- cor.test(df_vaf_mcf$scVAF, df_vaf_mcf$tumor_f, use = "everything", method = "pearson")
  tab_cor <- data.frame(
    sample = s,
    test = c('mcf_vaf', 'scVaf_vaf'),
    r = c(cor_mcf$estimate, cor_scVAF$estimate),
    p = c(cor_mcf$p.value, cor_scVAF$p.value)
  )
  tab_core_sample[[s]] <- tab_cor
}
tab_core_sample <- do.call(rbind, tab_core_sample)
tab_core <- rbind(tab_cor_all, tab_core_sample)
write.csv(tab_core, paste0(path_results, 'tables/cor_vaf.csv'))

## Plot correlations
p1 <- ggplot(df_vaf_mcf, aes(x = mcf, y = tumor_f)) +
  geom_point(size = 1.8) +
  geom_smooth(aes(mcf, tumor_f), method = 'lm', se = F, fullrange = T, size = 0.5)  +
  xlim(0,1) + ylim(0,1) +
  labs(x = 'Mutant cell fraction', 
       y = 'Variant allele frequency (WES)',
       subtitle = paste0(
         'r = ', round(tab_core$r[which(tab_core$sample == 'all' & tab_core$test == 'mcf_vaf')], digits = 2), 
         ', ', 
         'p = ', round(tab_core$p[which(tab_core$sample == 'all' & tab_core$test == 'mcf_vaf')], digits = 3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.subtitle = element_text(size = 17, hjust = 0.5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))

p2 <- ggplot(df_vaf_mcf, aes(x = scVAF, y = tumor_f)) +
  geom_point(size = 1.8) +
  geom_smooth(aes(mcf, tumor_f), method = 'lm', se = F, fullrange = T, size = 0.5)  +
  xlim(0,1) + ylim(0,1) +
  labs(x = 'Variant allele frequency (SCM-seq)', 
       y = 'Variant allele frequency (WES)',
       subtitle = paste0(
         'r = ', round(tab_core$r[which(tab_core$sample == 'all' & tab_core$test == 'scVaf_vaf')], digits = 2), 
         ', ', 
         'p = ', round(tab_core$p[which(tab_core$sample == 'all' & tab_core$test == 'scVaf_vaf')], digits = 3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.subtitle = element_text(size = 17, hjust = 0.5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))

ggsave(paste0(path_results, 'plots/cor_vaf_all.png'), p1 + p2 +
         plot_layout(ncol = 2, guides = 'collect', widths = 5, heights = 5),
       width = 10, height = 5)

# TRANSFER GENOTYPE INFORMATION TO SEURAT OBJECT ############################################
genotype_all_final <- read_csv("/hpcnfs/scratch/PGP/SCMseq/sc_genotype/results_and_plots/tables/genotype_all_final.csv")

geno_data_all <- list()
for (s in samples) {
  geno_data <- genotype_all_final %>% 
    filter(sample == s) %>% 
    mutate(barcode = paste0(barcode, '-', s)) %>%
    select(c(barcode, sample, gene, genotype)) %>%
    pivot_wider(names_from = gene, values_from = genotype) 
  
  geno_data[is.na(geno_data)] <- "unknown"
  
  geno_data_all[[s]] <- geno_data
}

geno_data_all <- Reduce(full_join, geno_data_all)

AML.combined.sct <- readRDS("/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML.combined.sct.rds")
AML.combined.sct$barcode <- colnames(AML.combined.sct)

AML.combined.sct@meta.data <- left_join(AML.combined.sct@meta.data, geno_data_all, by = c('barcode', 'sample'))
rownames(AML.combined.sct@meta.data) <- AML.combined.sct$barcode
saveRDS(AML.combined.sct, "/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML.combined.sct.rds")

sub_df <- AML.combined.sct@meta.data[,c('aggregated_lineage', 'barcode', 'cell_lineage')] %>%
  separate(barcode, into = c('barcode', 'sample'), sep = '-')

genotype_all_final <- left_join(genotype_all_final, sub_df, by = c('barcode', 'sample'))  

write_csv(genotype_all_final, paste0(path_results, 'tables/genotype_all_final.csv'))
