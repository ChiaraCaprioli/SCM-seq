#title: "07_genotype_mapping"
#author: "Chiara_Caprioli"
#date: "August 13th 2023"

## Set up
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(gmodels)
  library(Seurat)
  library(ggpubr)
})

load(paste0(path_data, "settings.RData"))

# set paths
path_main <- '/Users/ieo4874/Desktop/working/SCMseq/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/sc_genotype/')

# load Seurat object 
final_seurat <- readRDS(paste0(path_data, 'final_seurat.rds'))

# Expression dropout
## Percent dropout
### sample
df_sample <- data.frame()
for (s in samples_genotype) {
  dropout.matrix <- read_csv(paste0(path_results, 'tables/', s, '_dropout.matrix_modified.csv'))
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
for (s in samples_genotype) {
  dropout.matrix <- read_csv(paste0(path_results, 'tables/', s, '_dropout.matrix_modified.csv'))
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

for (s in samples_genotype) {
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
for (s in samples_genotype) {
  dropout.matrix <- read_csv(paste0(path_results, 'tables/', s, '_dropout.matrix_modified.csv'))
  dropout.matrix <- dropout.matrix %>% column_to_rownames(var = '...1') %>% as.matrix()
  
  dropout.matrix[which(dropout.matrix == 0)] <- "no"
  dropout.matrix[which(dropout.matrix == 1)] <- "yes"

  pdf(paste0(path_results, "plots/", s, "_drop_counts_heat.pdf"), width = 8, height = 8)
  h = Heatmap(
    dropout.matrix,
    col = c("navy", "#F0F0F0"),
    name = "Expression dropout",
    cluster_rows = T,
    cluster_columns = T,
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = T,
    row_names_side = "left",
    column_title = "cells",
    column_title_side = "top",
    show_column_names = F,
    width = unit(7, "cm"), height = unit(5, "cm")
  )
  
  print(h)
  dev.off()
  
}

# Merge smoothed and unsmoothed data
first_geno_all <- read_csv(paste0(path_results, 'tables/first_geno_all.csv'))

genotype_all <- data.frame()
for (s in samples_genotype) {
  first_geno_all_sample <- first_geno_all %>% filter(sample == s) 
  
  geno.matrix <- read_csv(paste0(path_results, 'tables/', s, '_geno.mat_modified.csv')) %>%
    pivot_longer(cols = 2:ncol(.), names_to = 'barcode', values_to = 'smoothed_genotype') %>%
    rename('...1' = "gene") 
  geno.matrix <- geno.matrix %>%  mutate(smoothed_genotype = recode(smoothed_genotype, 'MUT' = 'mutated', 'WT' = 'wild-type'))
  
  c.matrix <- read_csv(paste0(path_results, 'tables/', s, '_c.mat_modified.csv'))
  c.matrix <- c.matrix %>% rename('...1' = 'gene') %>%
    pivot_longer(cols = 2:ncol(.), names_to = 'barcode', values_to = 'smoothed_confidence') 
  
  df_ <- full_join(geno.matrix, c.matrix)
  
  geno_all_sample <- left_join(first_geno_all_sample, df_, by = c('barcode', 'gene'))
  
  genotype_all <- rbind(geno_all_sample, genotype_all)
  
}


# Re-classification rate
## confusion matrix to assess re-classification of cells imputed with high confidence with initial algorithm
for (s in samples_genotype) {
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

# Compare confidence before and after smoothing
for (s in samples_genotype) {
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

for (s in samples_genotype) {
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

# Select variants and cells based on coverage and confidence

genotype_all$barcode <- paste0(genotype_all$barcode, '-', genotype_all$sample)
genotype_all <- genotype_all[which(genotype_all$barcode %in% colnames(final_seurat)),]

genotype_all <- genotype_all %>% 
  separate(barcode, into = c("barcode", "stuff"), sep = "-") %>% 
  dplyr::select(-"stuff")

## We create different tiers based on the cutoff in coverage:
## Tier1: cell covered by ≥1 read
## Tier2: cell covered by ≥2 reads
## Tier3: cell covered by ≥3 reads

## Confidence by single variant

# table
L <- list()
for (i in seq_len(3)) {
  
  conf_table <- genotype_all %>%
    filter(n_reads_tot >= i) %>%
    group_by(sample, gene) %>%
    summarise(n_cells = n(),
              mean = mean(confidence),
              q1 = quantile(confidence)[1],
              q2 = quantile(confidence)[2],
              q3 = quantile(confidence)[3],
              q4 = quantile(confidence)[4],
              q5 = quantile(confidence)[5]
    ) %>%
    mutate(tier = paste0("Tier",i))
  L[[i]] <- conf_table
}

conf_table_tier <- do.call(rbind, L)

write_csv(conf_table_tier, paste0(path_results, "tables/stat_confidence_gene.csv"))

# plot
for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    temp_labels <- genotype_all %>%
      filter(sample == s, n_reads_tot >= i) %>%
      group_by(sample, gene) %>%
      summarise(n = n())
    
    p <- genotype_all %>%
      filter(sample == s, n_reads_tot >= i) %>%
      ggplot(aes(x = gene, y = confidence)) + 
      geom_jitter(width = 0.2, height = 0.001, cex = 0.1, color = "grey") +
      stat_summary(fun = "median", color = "red", geom = "point", shape = 18, size = 4) +
      ylab('Genotyping confidence') +
      ggtitle(paste0("Tier ", i)) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
      annotate(geom = "text",
               x = temp_labels$gene, y = -0.1, 
               label = paste0("n = ", temp_labels$n),
               size = 2.8) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none',
        plot.title = element_text(face = 'bold', vjust = 2, hjust = 0.5, size = 14)
      )
    
    ggsave(paste0(path_results, "plots/Tier",i, "_", s, "_gene_confidence_genotype.png"),
           p, width = 6, height = 5)
    
  }
}

var_to_keep <- conf_table_tier %>% 
  filter(n_cells >= 100 & q3 >= 0.75) %>%
  dplyr::select(sample, gene, tier) %>%
  mutate(tier = str_remove_all(tier, pattern = "Tier"))

## Confidence by single cell across all selected variants

df_mean_confidence_tier <- list()
df_mean_confidence_sample <- list()

for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    i_var_to_keep <- var_to_keep %>%
      filter(sample == s, tier == i) %>%
      pull(gene)
    
    x <- genotype_all %>%
      filter(sample == s, gene %in% i_var_to_keep, n_reads_tot >= i) %>%
      group_by(sample, barcode) %>%
      summarise(mean_confidence_cell = mean(confidence)) 
    x <- x %>% 
      arrange(desc(mean_confidence_cell)) %>%
      mutate(tier = paste0("Tier",i))
    x$barcode <- factor(x$barcode, levels = unique(x$barcode))
    
    df_mean_confidence_tier[[i]] <- x
    
    n_cells_above_thr <- x %>% filter(mean_confidence_cell >= 0.9) %>% nrow()
    
    p = ggplot(x, aes(x = barcode, y = mean_confidence_cell)) + 
      geom_point(size = 0.2) +
      geom_hline(yintercept = 0.9, color = "red", linetype = 2) +
      xlab("Cells") + ylab('Mean genotyping confidence') +
      scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,0.9,1)) +
      ggtitle(paste0("Tier ", i)) +
      annotate(geom = "text",
               x = x$barcode[60],
               hjust = 0,
               y = 0.93, 
               label = paste0("n cells = ", n_cells_above_thr),
               size = 3.5) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none',
        plot.title = element_text(face = 'bold', vjust = 2, hjust = 0.5, size = 14)
      )
    
    ggsave(paste0(path_results, "plots/Tier",i, "_", s, "_mean_cell_confidence_genotype.png"),
           p, width = 6, height = 5) 
    
  }
  df_mean_confidence_sample[[s]] <- do.call(rbind, df_mean_confidence_tier)
}

df_mean_confidence <- do.call(rbind, df_mean_confidence_sample)

df_mean_confidence <- df_mean_confidence %>%
  pivot_wider(values_from = mean_confidence_cell, 
              names_from = tier,
              names_glue = "{tier}_{.value}")

genotype_all <- full_join(genotype_all, df_mean_confidence,
                          by = c("sample", "barcode"))

# mark variants and cells to keep for each tier
df_sample <- list()
df_tier <- list()

for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    i_var_to_keep <- var_to_keep %>%
      filter(sample == s, tier == i) %>%
      pull(gene)
    
    subset_cells <- genotype_all %>% 
      filter(sample == s, n_reads_tot >=i) %>%
      mutate(tier = paste0("Tier",i))
    
    subset_cells$gene_keep <- ifelse(subset_cells$gene %in% i_var_to_keep, "yes", "no") # mark variants
    
    subset_cells$cell_keep =  ifelse(subset_cells[[paste0("Tier", i, "_mean_confidence_cell")]] >= 0.9, "yes", "no") # mark cells

    df_tier[[i]] <- subset_cells
        
  }
  
  df_sample[[s]] <- do.call(rbind, df_tier)
  
}

new_genotype_all <- do.call(rbind, df_sample)

write.csv(new_genotype_all, paste0(path_results, 'tables/genotype_all.csv'))

# Transfer genotype information to seurat object
## we only transfer tier 2

geno_data_all <- list()
for (s in samples_genotype) {
  geno_data <- new_genotype_all %>% 
    filter(sample == s,
           tier == "Tier2",
           gene_keep == "yes",
           cell_keep == "yes") %>% 
    mutate(barcode = paste0(barcode, '-', s)) %>%
    dplyr::select(c(barcode, sample, gene, genotype, cell_keep)) %>%
    pivot_wider(names_from = gene, values_from = genotype) 
  
  geno_data[is.na(geno_data)] <- "unknown"
  
  geno_data_all[[s]] <- geno_data
}

geno_data_all <- Reduce(full_join, geno_data_all)

final_seurat@meta.data <- left_join(final_seurat@meta.data, geno_data_all, by = c('barcode', 'sample'))
rownames(final_seurat@meta.data) <- final_seurat$barcode

saveRDS(final_seurat, paste0(path_data, "final_seurat.rds"))

sub_df <- final_seurat@meta.data[,c('aggregated_lineage2', 'barcode', 'lineage')] %>%
  separate(barcode, into = c('barcode', 'sample'), sep = '-')

new_genotype_all <- left_join(new_genotype_all, sub_df, by = c('barcode', 'sample'))  
new_genotype_all <- new_genotype_all %>%
  filter(!is.na(aggregated_lineage2))
write_csv(new_genotype_all, paste0(path_results, 'tables/genotype_all_final.csv'))

