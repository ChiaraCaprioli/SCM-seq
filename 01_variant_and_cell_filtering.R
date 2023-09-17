
#title: "Variant and cell filtering by genotyping confidence"
#author: "Chiara Caprioli"
#date: "September 15th 2023"

### Set up ###########################################

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(Seurat)
  library(ggpubr)
})

# set paths
path_main <- "/Users/ieo4874/Desktop/working/SCMseq/"
path_data <- "~/Desktop/temp/data/" 
path_results <- "~/Desktop/temp/genotype_analysis/"

# load custom settings
# including samples_genotype, a list of samples for which we currently have mutation analysis
load(paste0(path_main, "data/settings.RData"))  

# create folder structure
## all samples
if (!dir.exists(paste0(path_results, "all_samples/"))) {
  dir.create(paste0(path_results, "all_samples/plots"), recursive = T)
  dir.create(paste0(path_results, "all_samples/tables"))
}
  
## sample-specific folder with /plots and /tables sub-folders to store results
for (s in samples_genotype) {
  if (!dir.exists(paste0(path_results, s))) {
    dir.create(paste0(path_results, s, "/plots"), recursive = T)
    dir.create(paste0(path_results, s, "/tables"))
  }
}

# load Seurat object 
final_seurat <- readRDS(paste0(path_main, "data/final_seurat.rds"))

### Read counts by genotyping class  ###########################################

# For a given sample, before scoring cell genotyping confidence, we are interested in understanding
# if there is any imbalance in the coverage of genotyping classes across variants.

# Required input: result of CB demultiplexing after mutation analysis (i.e., df_bc_all.csv).
# It consists of a table of cells (=rows) and read counts for each genotyping class for each variant (=columns).

for (s in samples_genotype) {
  data_ <- read.csv(paste0(path_data, s, "/df_bc_all.csv"), row.names = 1)
  
  # correct error in df_bc_all.csv 
  data_ <- data_[1:(nrow(data_)-1),]
  
  # N reads
  data.1 <- data_ %>% 
    pivot_longer(cols = 2:ncol(data_), 
                 names_to = c('gene', 'genotype'), 
                 names_sep = '_', 
                 values_to = 'n_reads') %>%
    dplyr::select(-barcode) %>%
    mutate(genotype = factor(
      recode(genotype, 'MUT' = 'mutated', 'WT' = 'wild-type', 'mis' = 'mismatch'), 
      levels = c('mutated', 'wild-type', 'mismatch')))
  
  p1 <- ggplot(data.1, aes(x = gene, y = n_reads, fill = genotype)) +
    geom_col(width = 0.8) + 
    scale_fill_manual(name = 'Genotype', values = colors$genotype) +
    scale_y_continuous(labels = scales::comma) +
    ylab('N reads') +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      text = element_text(size = 13, color = 'black'),
      axis.text = element_text(color = 'black'),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 0, l = 10, unit = 'pt')
    )
  
  # Read percentage
  data.2 <- data_ %>%
    pivot_longer(cols = 2:ncol(data_), 
                 names_to = c('gene', 'genotype'), 
                 names_sep = '_', 
                 values_to = 'n_reads') %>%
    dplyr::select(-barcode) 
  
  data.2 <- data.2 %>%
    mutate(genotype = factor(
      recode(genotype, 'MUT' = 'mutated', 'WT' = 'wild-type', 'mis' = 'mismatch'), 
      levels = c('mutated', 'wild-type', 'mismatch'))) %>%
    group_by(gene, genotype) %>%
    summarize(count = sum(n_reads)) %>%
    spread(genotype, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_read_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('gene', 'total_read_count', everything())) %>%
    arrange(factor(gene, levels = levels(data.2$gene))) 
  
  labels <- data.2[,c(1,2)]
  
  data.2 <- data.2 %>%
    dplyr::select(-c('total_read_count')) %>%
    reshape2::melt(id.vars = 'gene') 
  colnames(data.2)[1] <- "variant"
  
  write_csv(data.2, paste0(path_results, s, "/tables/genotype_read_count.csv"))
  
  p2 <- ggplot(data.2, aes(variant, value)) +
    geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', width = 0.8) +
    scale_fill_manual(name = 'Genotype', values = colors$genotype) +
    scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      text = element_text(size = 13, color = 'black'),
      axis.text = element_text(color = 'black'),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 0, l = 10, unit = 'pt')
    )
  
  ggsave(paste0(path_results, s, "/plots/genotype_read_count.png"),
         p1 + p2 +
           patchwork::plot_layout(ncol = 2, guides = 'collect', heights = 5, width = 7) +
           patchwork::plot_annotation(
             title = s, 
             theme = theme(
               plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2))
           ),
         height = 5, width = 14)
  
}

### Merge genotyping confidence ###########################################

# store result for all variants and samples into a single file

# Required input: result of genotyping confidence assignment stored in /genotype_imputation_diff_appr_5%_err/.
# This folder contains, for each variant, one file for mutated cells and one file for wt cells. 
# To date, genotyping confidence assignment is obtained by:
# 1. defining candidate genotypes based on gene and cell error-based threshold (### mut = df[df.alt > df.cell_err], ### wt = df[df.alt <= df.cell_err] )
# 2. applying Tirelli's original algorithm, i.e. using different approaches to impute confidence of mutated and wt genotypes

df_mut_all = list()
df_wt_all = list()

for (s in samples_genotype) {
  genes = read.csv(paste0(path_data, s, "/tg_wes_intersect.csv"))
  path1 = paste0(path_data, s, '/genotype_imputation_diff_appr_5%_err/')
  
  df_mut = data.frame()
  df_wt = data.frame()
  
  for (g in genes$genelist) {
    
    if (file.exists(paste0(path1, g, '_wt_90acc_median.csv'))) {
      
      mut = read.csv(paste0(path1, g, '_mutated_NoFilter.csv'), row.names = 1)
      mut = mut %>% mutate(gene = g, sample = s, genotype = 'mutated')
      
      wt = read.csv(paste0(path1, g, '_wt_90acc_median.csv'), row.names = 1)
      wt = wt %>% mutate(gene = g, sample = s, genotype = 'wild-type')
      
    } else { 
      next
      
    }
    
    df_mut = rbind(mut, df_mut)
    df_mut_all[[s]] <- df_mut 
    
    df_wt = rbind(wt, df_wt)
    df_wt_all[[s]] <- df_wt 
    
  }
  
}

df_mut_all = do.call(rbind, df_mut_all)
df_wt_all = do.call(rbind, df_wt_all)
genotype_all <- rbind(df_mut_all, df_wt_all) 
rownames(genotype_all) <- NULL
colnames(genotype_all) <- c(
  "barcode", "n_reads_mut", "n_reads_wt", 
  "n_reads_mis", "n_reads_tot", "gene_error_rate", 
  "cell_error_rate", "accuracy", "confidence", 
  "variant", "sample", "genotype"
  )

write.csv(genotype_all, paste0(path_results, "all_samples/tables/genotype_all.csv")) 

### Expression dropout ###########################################

# We are interested in understanding how many cells/variants
# miss genotype information due to expression dropout of the target variant.

# Define function to get expression dropout matrix
make.dropout.matrix <- function(sample) {
  
  df_bc_all <- read.csv(paste0(path_data, s, "/df_bc_all.csv"), row.names = 1) %>%
    slice(1:(nrow(.)-1)) %>% #!! correct error in df_bc_all.csv !!
    column_to_rownames(var = 'barcode')
  
  ## get counts matrices
  M.ref <- df_bc_all %>% 
    select(contains('WT')) %>% 
    t() 
  rownames(M.ref) <- str_split_i(rownames(M.ref), pattern = "_", 1)
  
  M.alt <- df_bc_all %>% 
    select(contains('MUT')) %>% 
    t() 
  rownames(M.alt) <- str_split_i(rownames(M.alt), pattern = "_", 1)
  
  if (nrow(M.ref) == nrow(M.alt) &
      ncol(M.ref) == ncol(M.alt)) {
    
    m <- matrix(0, nrow(M.ref), ncol(M.ref))
    m[(M.ref == 0) & (M.alt == 0)] <- 1 # 1 = expression dropout, 0 = no expression dropout
    
    m <- as.data.frame(m)
    colnames(m)<- colnames(M.ref)
    m$variant <- rownames(M.ref)
    m <- m %>% column_to_rownames(var = 'variant')
    m
    
  } else {
    message("M.ref and M.alt have different dimensions!")
  }
  
}

L <- list()
for (s in samples_genotype) {
  
  # get and save dropout matrix
  dropout.matrix <- make.dropout.matrix(s)
  write.csv(dropout.matrix, paste0(path_results, s, "/tables/dropout.matrix.csv"))
  
  # dropout fraction by variant
  count_drop_var <- as.data.frame(rowSums(dropout.matrix) / ncol(dropout.matrix))
  colnames(count_drop_var) <- 'dropout_fraction'
  count_drop_var <- count_drop_var %>% 
    rownames_to_column(var = 'variant') %>% 
    arrange(desc(dropout_fraction))
  count_drop_var$variant <- factor(count_drop_var$variant, levels = count_drop_var$variant)
  write_csv(count_drop_var, paste0(path_results, s, '/tables/var_dropout.csv'))
  
  p1 <- ggplot(count_drop_var, aes(x = variant, y = dropout_fraction)) +
    geom_col(width = 0.75) +
    ylab('Expression dropout (%)') +
    ylim(0,1) +
    ggtitle(s) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
  
  ggsave(paste0(path_results, s, "/plots/var_dropout.png"), p1,
         width = 5, height = 4)
  
  # dropout fraction by sample overall
  count_drop_sample <- as.numeric(table(dropout.matrix == 1)[2]) # note: '1' = dropout
  total_count <- (nrow(dropout.matrix)*ncol(dropout.matrix)) 
  dropout_fraction <- round(count_drop_sample / total_count, digits = 2)
  x <- data.frame(
    sample = s,
    dropout_fraction = dropout_fraction
  )
  L[[s]] <- x
  
  # plot heatmap
  order_cells <- t(dropout.matrix) %>% as.data.frame()
  order_cells$order_cells <- rowSums(order_cells)
  order_cells <- order_cells %>%
    arrange(order_cells)
  
  dropout.matrix <- 
    dropout.matrix[match(rev(levels(count_drop_var$variant)), rownames(dropout.matrix)),
                   match(rownames(order_cells), colnames(dropout.matrix))] %>%
    as.matrix()

  dropout.matrix[which(dropout.matrix == 0)] <- "no"
  dropout.matrix[which(dropout.matrix == 1)] <- "yes"
  
  pdf(paste0(path_results, s, "/plots/exp_dropout_heat.pdf"), width = 8, height = 8)
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

df_sample <- do.call(rbind, L)
write_csv(df_sample, paste0(path_results, "all_samples/tables/sample_dropout.csv"))

p2 <- ggplot(df_sample, aes(x = sample, y = dropout_fraction)) +
  geom_col(width = 0.75) +
  ylab('Expression dropout (%)') +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(), 
        aspect.ratio = 1)

ggsave(paste0(path_results, "all_samples/plots/sample_dropout.png"), p2,
       width = 5, height = 4)

### Select variants and cells based on coverage and confidence ###########################################

# Ensure we only analyze high-quality cells included in the final seurat object
genotype_all$barcode <- paste0(genotype_all$barcode, '-', genotype_all$sample)
genotype_all <- genotype_all[which(genotype_all$barcode %in% colnames(final_seurat)),] 
genotype_all$barcode <- str_split_i(genotype_all$barcode, pattern = "-", 1)

# Remove low accuracy cells
genotype_all <- genotype_all[genotype_all$accuracy >= 0.9,] 

# Create different tiers based on coverage cutoff coverage:
## Tier1: cell covered by ≥1 read 
## Tier2: cell covered by ≥2 reads
## Tier3: cell covered by ≥3 reads

# Confidence by single variant
## For each variant within each tier, we check the distribution of confidence,
## then we filter variants based on median confidence and n cells

# table
L <- list()
for (i in seq_len(3)) {
  
  conf_table <- genotype_all %>%
    filter(n_reads_tot >= i) %>%
    group_by(sample, variant) %>%
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

write_csv(conf_table_tier, paste0(path_results, "all_samples/tables/stat_confidence_var.csv"))

# plot
for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    temp_labels <- genotype_all %>%
      filter(sample == s, n_reads_tot >= i) %>%
      group_by(sample, variant) %>%
      summarise(n = n())
    
    p <- genotype_all %>%
      filter(sample == s, n_reads_tot >= i) %>%
      ggplot(aes(x = variant, y = confidence)) + 
      geom_jitter(width = 0.2, height = 0.001, cex = 0.1, color = "grey") +
      stat_summary(fun = "median", color = "red", geom = "point", shape = 18, size = 4) +
      ylab('Genotyping confidence') +
      ggtitle(paste0("Tier ", i)) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
      annotate(geom = "text",
               x = temp_labels$variant, y = -0.1, 
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
    
    ggsave(paste0(path_results, s, "/plots/Tier", i, "_var_confidence_genotype.png"),
           p, width = 6, height = 5)
    
  }
}

# filter variants
var_to_keep <- conf_table_tier %>% 
  filter(n_cells >= 100 & q3 >= 0.75) %>%
  dplyr::select(sample, variant, tier) %>%
  mutate(tier = str_remove_all(tier, pattern = "Tier"))

df <- var_to_keep %>%
  mutate(tier = paste0("tier", tier)) %>%
  group_by(sample, tier) %>%
  summarise(n_variants = n()) %>%
  pivot_wider(names_from = tier, values_from = n_variants)

write_csv(df, paste0(path_results, "all_samples/tables/var_to_keep.csv"))

# Confidence by single cell across all selected variants
## For each cell within each tier, we check the mean confidence across all
## variants passing confidence threshold

df_mean_confidence_tier <- list()
df_mean_confidence_sample <- list()

for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    i_var_to_keep <- var_to_keep %>%
      filter(sample == s, tier == i) %>%
      pull(variant)
    
    x <- genotype_all %>%
      filter(sample == s, variant %in% i_var_to_keep, n_reads_tot >= i) %>%
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
    
    ggsave(paste0(path_results, s, "/plots/Tier", i, "_mean_cell_confidence_genotype.png"),
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

# Mark variants and cells to keep for each tier
df_sample <- list()
df_tier <- list()

for (s in samples_genotype) {
  for (i in seq_len(3)) {
    
    i_var_to_keep <- var_to_keep %>%
      filter(sample == s, tier == i) %>%
      pull(variant)
    
    subset_cells <- genotype_all %>% 
      filter(sample == s, n_reads_tot >=i) %>%
      mutate(tier = paste0("Tier",i))
    
    subset_cells$var_keep <- ifelse(subset_cells$variant %in% i_var_to_keep, "yes", "no") # mark variants
    
    subset_cells$cell_keep =  ifelse(subset_cells[[paste0("Tier", i, "_mean_confidence_cell")]] >= 0.9, "yes", "no") # mark cells
    
    df_tier[[i]] <- subset_cells
    
  }
  
  df_sample[[s]] <- do.call(rbind, df_tier)
  
}

new_genotype_all <- do.call(rbind, df_sample)

write.csv(new_genotype_all, paste0(path_results, "all_samples/tables/genotype_all.csv"))

### Transfer genotype information to seurat object ###########################################

# tier 2 seems to offer the best trade-off between confidence and number of variants/genotyped cells:
# we only transfer tier 2.

geno_data_all <- list()
for (s in samples_genotype) {
  geno_data <- new_genotype_all %>% 
    filter(sample == s,
           tier == "Tier2",
           var_keep == "yes",
           cell_keep == "yes") %>% 
    mutate(barcode = paste0(barcode, '-', s)) %>%
    dplyr::select(c(barcode, sample, variant, genotype, cell_keep)) %>%
    pivot_wider(names_from = variant, values_from = genotype) 
  
  geno_data[is.na(geno_data)] <- "unknown"
  
  geno_data_all[[s]] <- geno_data
}

geno_data_all <- Reduce(full_join, geno_data_all)

final_seurat@meta.data <- left_join(final_seurat@meta.data, geno_data_all, by = c('barcode', 'sample'))
rownames(final_seurat@meta.data) <- final_seurat$barcode

saveRDS(final_seurat, paste0(path_data, "final_seurat.rds"))

# we also transfer some metadata of interest to the genotype_all.csv, 
# to make later analyses easier
sub_df <- final_seurat@meta.data[,c('aggregated_lineage2', 'barcode', 'lineage')] %>%
  separate(barcode, into = c('barcode', 'sample'), sep = '-')

new_genotype_all <- left_join(new_genotype_all, sub_df, by = c('barcode', 'sample'))  

write_csv(new_genotype_all, paste0(path_results, "all_samples/tables/genotype_all.csv"))
