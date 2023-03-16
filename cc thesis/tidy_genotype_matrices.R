# Create input matrices to SCITE and SciFit
# chiara caprioli
# Feb 17 2023

# SET UP ############################################

## Load packages
library(tidyverse)
library(Seurat)

## Paths
path_main <- '/hpcnfs/scratch/PGP/SCMseq/sc_genotype/' 
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/')

## Define samples
samples <- c('AML4', 'AML5', 'sAML1')

# TIDY UP GENOTYPE MATRICES ############################################

df_meta_all <- list()
for (s in samples) {
  # define target mutated genes
  genotype_all_final <- read_delim(paste0(path_results, 'tables/genotype_all_final.csv'), delim = ',', show_col_types = FALSE) 
  
  target_genes <- genotype_all_final %>%
    filter(sample == s & !is.na(gene)) %>%
    pull('gene') %>% unique()
  
  # load and subset Seurat object
  seurat = readRDS('/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML.combined.sct.rds')
  subset = seurat[, colnames(seurat)[which(seurat$sample == s)] ]
  
  # recode matrix according to genotype
  geno_matrix = subset@meta.data %>% 
    dplyr::select(target_genes) %>%
    as.matrix()
  
  geno_matrix[geno_matrix == 'mutated'] <- 1
  geno_matrix[geno_matrix == 'wild-type'] <- 0
  geno_matrix[geno_matrix == 'unknown'] <- 2
  
  mode(geno_matrix) <- 'numeric'
  
  # count unknowns and n mutations per cell
  df <- geno_matrix %>% as.data.frame() 
  df <- df %>% mutate(
    n_unknown = apply(df == 2, 1, sum),
    n_mut_cell = apply(df == 1, 1, sum)
  ) %>%
    rownames_to_column(var = 'barcode')
  df$tot_mut_sample = length(df[2:ncol(df)])-2
  
  df <- df %>% 
    mutate(n_mut_cell_bin = ifelse(
      df$n_mut_cell == '0', '0', ifelse(
        df$n_mut_cell == '1', '1', ifelse(
          df$n_mut_cell == '2', '2', ifelse(
            df$n_mut_cell == '3', '3', '>3')
        )
      )))
  
  df <- df %>% 
    mutate(genotype = ifelse(df$n_unknown == length(target_genes), 'unknown', 'genotyped'))
  
  # pie chart n mutations by sample (for cells with at least one mutation)
  data <- df %>%
    dplyr::select(n_mut_cell_bin) %>%
    dplyr::filter(n_mut_cell_bin != 0) %>%
    table() %>% as.data.frame()
  data$n_mut_cell_bin <- factor(data$n_mut_cell_bin, levels = c('1', '2', '3', '>3'))
  cols <- c("#FFEDA0", "#FD8D3C", "#E31A1C", "#800026")
  
  pdf(paste0(path_results, 'plots/', s, '_pie_mutations.pdf'))
  pie(data$Freq, labels = c('1', '2', '3', '>3'), radius = 0.8, 
      clockwise = T, col = cols, border = T, cex = 2, main = s) 
  dev.off()
  
  # add metadata column to sample seurat object with n mutations per cell
  df <- df[match(colnames(subset), df$barcode),]
  subset$n_mut_cell <- df$n_mut_cell
  subset$n_mut_cell_bin <- df$n_mut_cell_bin
  subset$n_mut_cell_bin <- factor(subset$n_mut_cell_bin, levels = c('0', '1', '2', '3', '>3'))
  subset$tot_mut_sample <- df$tot_mut_sample
  subset$genotype <- df$genotype
  subset$genotype <- factor(subset$genotype, levels = c('genotyped', 'unknown'))
  
  df_meta <- data.frame(
    sample = subset$sample,
    barcode = colnames(subset),
    n_mut_cell = subset$n_mut_cell,
    n_mut_cell_bin = subset$n_mut_cell_bin,
    tot_mut_sample = subset$tot_mut_sample,
    genotype = subset$genotype
  )
  
  df_meta_all[[s]] <- df_meta 

  
  ##### SAVE MATRIX FOR WHOLE DATASET (tons of unknowns)
  matrix <- df %>% 
    dplyr::select(-c('n_mut_cell', 'n_mut_cell_bin', 'tot_mut_sample', 'genotype', 'n_unknown')) %>%
    column_to_rownames(var = 'barcode') %>%
    t() 
  
  # save matrix with rownames and colnames
  matrix_bc <- matrix 
  saveRDS(matrix, paste0(path_data, s, '_matrix_all_cells_row_colnames.rds'))
  
  # prepare SCITE input
  matrix[matrix == 2] <- 3 
  colnames(matrix) <- NULL 
  rownames(matrix) <- NULL 
  
  write_delim(as.data.frame(matrix), 
              paste0(path_data, s, '_SCITE_matrix_all_cells.csv'),
              col_names = F)
  
  # prepare SiFit input
  matrix_bc[matrix_bc == 2] <- 3
  matrix_bc <- matrix_bc %>%
    as.data.frame() %>%
    rownames_to_column(var = 'rowname')
  matrix_bc$rowname <- c(1:nrow(matrix_bc))
  
  write_delim(as.data.frame(matrix_bc), paste0(path_data, s, '_SiFit_matrix_all_cells.txt'),
              col_names = F, delim = " ")
  
  cellNames <- data.frame(
    colnames(matrix_bc)[2:ncol(matrix_bc)]
  ) %>% t()
  write_delim(as.data.frame(cellNames), paste0(path_data, s, '_cellNames_all.txt'), col_names = F, delim = " ")
  
  geneNames <- data.frame(
    rownames(matrix_bc)
  ) 
  write_delim(geneNames, paste0(path_data, s, '_geneNames_all.txt'),  col_names = F, delim = " ")
  
  ##### SAVE MATRIX FOR MOST INFORMATIVE DATA SUBSET 
  matrix <- df %>% 
    dplyr::select(-c('n_mut_cell', 'n_mut_cell_bin', 'tot_mut_sample', 'genotype', 'n_unknown')) %>%
    column_to_rownames(var = 'barcode') %>%
    t() 
  mode(matrix) <- 'numeric'
  matrix <- matrix[,!is.na(colSums(matrix))]

  dense_matrix <- phyloRNA::densest_subset(matrix, empty = 2, density = 0.6)
  dense_matrix <- dense_matrix$result
  
  # save dense matrix with rownames and colnames
  dense_matrix_bc <- dense_matrix
  saveRDS(dense_matrix_bc, paste0(path_data, s, '_matrix_dense_row_colnames.rds'))
  
  # prepare SCITE input
  dense_matrix[dense_matrix == 2] <- 3 
  colnames(dense_matrix) <- NULL 
  rownames(dense_matrix) <- NULL 
  
  write_delim(as.data.frame(dense_matrix), 
              paste0(path_data, s, '_SCITE_matrix_dense.csv'),
              col_names = F)
  
  # prepare SiFit input
  dense_matrix_bc[dense_matrix_bc == 2] <- 3
  dense_matrix_bc <- dense_matrix_bc %>%
    as.data.frame() %>%
    rownames_to_column(var = 'rowname')
  dense_matrix_bc$rowname <- c(1:nrow(dense_matrix_bc))
  
  write_delim(as.data.frame(dense_matrix_bc), paste0(path_data, s, '_SiFit_matrix_dense.txt'),
              col_names = F, delim = " ")
  
  cellNames <- data.frame(
    colnames(dense_matrix_bc)[2:ncol(dense_matrix_bc)]
  ) %>% t()
  write_delim(as.data.frame(cellNames), paste0(path_data, s, '_cellNames_dense.txt'), col_names = F, delim = " ")
  
  geneNames <- data.frame(
    rownames(dense_matrix_bc)
  ) 
  write_delim(geneNames, paste0(path_data, s, '_geneNames_dense.txt'),  col_names = F, delim = " ")
}

df_meta_all <- do.call(rbind, df_meta_all) 
AML.combined.sct <- readRDS("/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML.combined.sct.rds")

if(identical(df_meta_all$barcode, colnames(AML.combined.sct)) == TRUE) {
  AML.combined.sct$n_mut_cell <- df_meta_all$n_mut_cell
  AML.combined.sct$n_mut_cell_bin <- df_meta_all$n_mut_cell_bin
  AML.combined.sct$n_mut_cell_bin <- factor(AML.combined.sct$n_mut_cell_bin, levels = c('0', '1', '2', '3', '>3'))
  AML.combined.sct$tot_mut_sample <- df_meta_all$tot_mut_sample
  AML.combined.sct$genotype <- df_meta_all$genotype
  AML.combined.sct$genotype <- factor(AML.combined.sct$genotype, levels = c('genotyped', 'unknown'))
  saveRDS(AML.combined.sct, "/hpcnfs/scratch/PGP/SCMseq/short_read_analysis/data/AML.combined.sct.rds")
} else {
  message('adding metadata failed due to wrong barcode matching')
}
