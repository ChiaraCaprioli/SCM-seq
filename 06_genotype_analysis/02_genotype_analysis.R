
#title: "Genotype analysis"
#author: "Chiara Caprioli"
#date: "September 16th 2023"

### Set up ###########################################

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(Seurat)
  library(ggpubr)
  library(readxl)
})

# set paths
path_main <- "/Users/ieo4874/Desktop/working/SCMseq/"
path_data <- "~/Desktop/temp/data/" # path to final seurat object and result of genotype confidence assignment
path_results <- "~/Desktop/temp/genotype_analysis/"

# load custom settings
load(paste0(path_main, "data/settings.RData")) # including samples_genotype, a list of samples for which we currently have mutation analysis 

# load Seurat object 
final_seurat <- readRDS(paste0(path_main, "data/final_seurat.rds"))

# load filtered genotype table
genotype_all <- read_csv(paste0(path_results, "all_samples/tables/genotype_all.csv"))

### Correlation metrics ###########################################

# We obtain metrics to correlate WES and SCM-seq:
# - VAF (from WES)
# - pseudo-bulk VAF (pbVAF): ratio between number of mutated reads and total number of genotyped reads 
# (defined for each variant across all cells)
# - mutant cell fraction (MCF): ratio between number of mutated cells and total number of genotyped cells
# (defined for each variant and each genotyping tier)

for (t in unique(genotype_all$tier)) {
  
  # extract high-confidence variants and cells 
  tier_genotype_all <- 
    genotype_all %>%
    filter(var_keep == "yes",
           cell_keep == "yes",
           tier == t) %>% 
    pivot_longer(cols = c(n_reads_mut, n_reads_wt, n_reads_mis), 
                 names_to = 'read_type', values_to = 'n_reads_type')
  
  # pseudo-bulk VAF 
  df_vaf_pbVAF <- list()
  for (s in samples_genotype) {
    
    dat_ = tier_genotype_all %>% 
      filter(sample == s)
    variant = unique(dat_$variant)
    
    df <- data.frame()
    for (v in variant) {
      
      x = dat_ %>% 
        filter(variant == v & read_type != 'n_reads_mis') %>% 
        dplyr::select(barcode, sample, variant, read_type, n_reads_type) %>%
        group_by(sample, variant, read_type) %>% 
        summarise(N_tot_reads = sum(n_reads_type))
      
      df = rbind(x, df)
      df_vaf_pbVAF[[s]] <- df 
      
    }
  }
  
  df_vaf_pbVAF <- do.call(rbind, df_vaf_pbVAF)
  
  df_vaf_pbVAF <- 
    df_vaf_pbVAF %>% 
    mutate(read_type = recode(read_type, 'mutated'='n_reads_mut', 'wild-type' = 'wt')) %>%
    pivot_wider(names_from = read_type, values_from = N_tot_reads, values_fill = 0)
  
  colnames(df_vaf_pbVAF)[3] <- 'MUT_N_tot_reads'
  colnames(df_vaf_pbVAF)[4] <- 'WT_N_tot_reads'
  
  # Mutant cell fraction (MCF): defined for each variant as the ratio between mutated cells 
  # and total number of genotyped cells
  df_vaf_mcf <- list()
  for (s in samples_genotype) {
    
    data = tier_genotype_all %>% filter(sample == s)
    variant = unique(data$variant)
    
    df <- data.frame()
    
    for (v in variant) {
      
      x = data %>% 
        filter(variant == v) %>% 
        distinct(barcode, .keep_all=T) %>%
        dplyr::select(barcode, sample, variant, genotype) %>%
        group_by(sample, variant, genotype) %>% 
        summarise(count = n())
      x = x %>% mutate(N_tot_cells = data %>% distinct(barcode) %>% nrow())
      
      df = rbind(x, df)
      df_vaf_mcf[[s]] <- df 
      
    }
  }
  
  df_vaf_mcf <- do.call(rbind, df_vaf_mcf)
  df_vaf_mcf <- 
    df_vaf_mcf %>% 
    pivot_wider(names_from = genotype, values_from = count, values_fill = 0)
  
  colnames(df_vaf_mcf)[4] <- 'MUT_N_tot_cells'
  colnames(df_vaf_mcf)[5] <- 'WT_N_tot_cells'
  
  for (i in df_vaf_mcf$variant){
    N_cells_UNK = (df_vaf_mcf$N_tot_cells - (df_vaf_mcf$MUT_N_tot_cells + df_vaf_mcf$WT_N_tot_cells))
  }
  
  df_vaf_mcf$N_cells_UNK <- N_cells_UNK
  
  df_vaf_mcf <- full_join(df_vaf_mcf, df_vaf_pbVAF, by = c('sample', 'variant'))
  
  # associate WES VAF
  vaf <- read_delim(paste0(path_data, "mutations_AML_VAF.csv"), delim = ',')
  vaf <- vaf %>% filter(variant %in% df_vaf_mcf$variant) %>%
    dplyr::select(sample, variant, VAF)
  df_vaf_mcf <- full_join(df_vaf_mcf, vaf, by = c("variant", "sample"))
  
  for (i in df_vaf_mcf$variant){
    mcf = df_vaf_mcf$MUT_N_tot_cells/(df_vaf_mcf$MUT_N_tot_cells + df_vaf_mcf$WT_N_tot_cells)
    pbVAF = df_vaf_mcf$MUT_N_tot_reads/(df_vaf_mcf$MUT_N_tot_reads + df_vaf_mcf$WT_N_tot_reads) 
  }
  
  df_vaf_mcf <- df_vaf_mcf %>%
    cbind(mcf = mcf, pbVAF = pbVAF)
  
  write_csv(df_vaf_mcf, paste0(path_results, "all_samples/tables/", t, "_df_vaf_mcf.csv"))

}

### Pyclone ###########################################

# Rather than just correlating pbVAF and MCF to VAF from WES, we also want to know 
# the relationship between MCF and the cellular prevalence we can estimate from bulk analysis. 
# To this end, we run Pyclone on CNV-corrected WES data, only for high-confidence variants.

# Prepare input to pyclone
for (t in unique(genotype_all$tier)) {
  for (s in samples_genotype) {
    
    # get cnv data
    file <- list.files(paste0(path_main, 'clones/ascat/', s, "_D_vs_", s, "_N/"), 
                       pattern = paste0(s, "_D_vs_", s, "_N.cnvs.txt")) 
    CNV <- read_delim(paste0(path_main, "clones/ascat/", s, "_D_vs_", s, "_N/", file))
    CNV$normal_cn <- rep(2, nrow(CNV))
    
    # get sex information and change normal cn of chrX accordingly
    patient.metadata <- read_delim(paste0(path_main, "data/patient_metadata.tsv")) %>% filter(Sample == s)
    
    if (patient.metadata$Sex == "M") {
      CNV$normal_cn[which(CNV$chr == 'X')] <- 1
    } else {
      CNV <- CNV
    }
    
    # map variants from WES to CNV data 
    var_ <- read_csv(paste0(path_results, "all_samples/tables/genotype_all.csv")) %>%
      filter(sample == s, var_keep == "yes", tier == t) %>%
      pull(variant) %>% unique()
    
    file_var <- paste0(path_main, 'clones/WES/old/WES_to_SCMseq.xlsx')
    sheets <- excel_sheets(file_var)
    sample_sheets <- grep(s, sheets, value = TRUE, fixed = T)
    
    if (any(grep("reseq", sample_sheets))) {
      sample_sheets <- sample_sheets[grep("reseq", sample_sheets)]
    } else {
      sample_sheets <- sample_sheets
    }
    
    sample_WES <- do.call(rbind, lapply(sample_sheets, read_excel, path = file_var))
    sample_WES <- sample_WES[which(
      sample_WES$Hugo_Symbol %in% var_),
      c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "t_alt_count", "t_ref_count")]
    
    sample_WES <- sample_WES %>% 
      separate(Chromosome, into = c('stuff', 'chr'), sep = 'chr') %>% 
      dplyr::select(-'stuff')
    x <- left_join(sample_WES, CNV,  by = 'chr', multiple = "all")
    
    pyclone_tsv <- data.frame()
    for (i in rownames(x)) {
      row = x[i,]
      if (row$Start_Position >= row$startpos & row$End_Position <= row$endpos) {
        pyclone_tsv <- rbind(row, pyclone_tsv)
      } else {
        next
      }
    }
    
    pyclone_tsv <- pyclone_tsv %>%
      unite('mutation_id', 1:4, sep = ':', remove = T) %>%
      rename('minor_cn' = 'nMinor',
             'major_cn' = 'nMajor',
             'ref_counts' = 't_ref_count',
             'var_counts' = 't_alt_count') %>%
      dplyr::select(-c(startpos, endpos))
    
    write_tsv(pyclone_tsv, paste0(path_main, 'clones/WES/pyclone_expressed_var/', s, "_", t, '_pyclone_input.tsv'))
    
  } 
}

# Run pyclone

## note that we set tumour_contents = 1, to avoid the effect of blast percentage while
## computing cellular prevalence 
## This allows a more faithful comparison to SCM-seq results, 
## which for now does not consider tumor vs non tumor composition.

# PyClone run_analysis_pipeline --in_files /path/to/sample_tier_pyclone_input.tsv \
# --working_dir /path/to/sample/tier \
# --samples <sample> --tumour_contents 1 --num_iters 10000 --burnin 1000 

# Update VAF with CNV-corrected data for all mutations scored by SCMseq 

for (t in unique(genotype_all$tier)) {
  pyclone_out_all <- list()
  for (s in samples_genotype) {
    pyclone_out <- read_delim(paste0(path_main, 'clones/WES/pyclone_expressed_var/', s, "/", t, '/tables/loci.tsv')) %>%
      separate(mutation_id, into = c('variant', 'chr', 'stuff'), sep = ':', remove = F, extra = 'merge') %>%
      dplyr::select(-'stuff')
    pyclone_out$sample <- s
    pyclone_out$cnv_corrected_vaf <- ifelse(pyclone_out$chr != 'X', 
                                            pyclone_out$cellular_prevalence/2, pyclone_out$variant_allele_frequency)
    pyclone_out_all[[s]] <- pyclone_out
  }
  
  pyclone_out_all <- do.call(rbind, pyclone_out_all) 
  pyclone_out_all <- pyclone_out_all %>% 
    dplyr::select(variant, sample_id, cellular_prevalence, cnv_corrected_vaf, variant_allele_frequency)
  
  df_vaf_mcf <- read_csv(paste0(path_results, "all_samples/tables/", t, "_df_vaf_mcf.csv"))

  df_vaf_mcf <- full_join(df_vaf_mcf, pyclone_out_all, by = c('variant', 'sample'='sample_id') )
  
  write_csv(df_vaf_mcf, paste0(path_results, "all_samples/tables/", t, "_df_vaf_mcf.csv"))
  
}

### Correlation WES-to-SCMseq ###########################################

# For the metrics of interest between single-cell and bulk data, 
# we compute Spearman correlation and fit a linear model 
# using all variants from all samples.

for (t in unique(genotype_all$tier)) {
  
  df_vaf_mcf <- read_csv(paste0(path_results, "all_samples/tables/", t, "_df_vaf_mcf.csv"))
  
  # MCF vs cellular prevalence
  cor_test_cp <- cor.test(df_vaf_mcf$mcf, 
                          df_vaf_mcf$cellular_prevalence, 
                          use = "everything", 
                          method = "spearman")
  
  p1 <- ggplot(df_vaf_mcf, aes(mcf, cellular_prevalence)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, linewidth = 0.5) +
    xlim(0,1) + ylim(0,1) +
    labs(x = 'Mutant cell fraction (SCM-seq)', 
         y = 'Cellular prevalence (WES)') +
    annotate(geom = 'text',
             label = paste0('Spearman rho = ', round(cor_test_cp$estimate, digits = 2), 
                            ", p = ", round(cor_test_cp$p.value, digits = 3)), 
             x = 0.3, y = 0.9,
             size = 4) +
    theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.subtitle = element_text(size = 17, hjust = 0.5),
          plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))
  
  ggsave(paste0(path_results, "all_samples/plots/", t, "_cor_MCF_CP.png"), p1, height = 5, width = 5)
  
  # MCF vs VAF (non CNV-corrected)
  cor_test_vaf <- cor.test(df_vaf_mcf$mcf, 
                           df_vaf_mcf$VAF, 
                           use = "everything", 
                           method = "spearman")
  
  p2 <- ggplot(df_vaf_mcf, aes(mcf, VAF)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, linewidth = 0.5) +
    xlim(0,1) + ylim(0,1) +
    labs(x = 'Mutant cell fraction (SCM-seq)', 
         y = 'Variant allele frequency (WES)') +
    annotate(geom = 'text',
             label = paste0('Spearman rho = ', round(cor_test_vaf$estimate, digits = 2), 
                            ", p = ", round(cor_test_vaf$p.value, digits = 3)), 
             x = 0.3, y = 0.9,
             size = 4) +
    theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.subtitle = element_text(size = 17, hjust = 0.5),
          plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))
  
  ggsave(paste0(path_results, "all_samples/plots/", t, "cor_MCF_VAF.png"), p2, height = 5, width = 5)
  
  # MCF vs VAF (CNV-corrected)
  cor_test_vaf_corrected <- cor.test(df_vaf_mcf$mcf, 
                                     df_vaf_mcf$cnv_corrected_vaf, 
                                     use = "everything", 
                                     method = "spearman")
  
  p3 <- ggplot(df_vaf_mcf, aes(mcf, cnv_corrected_vaf)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, linewidth = 0.5) +
    xlim(0,1) + ylim(0,1) +
    labs(x = 'Mutant cell fraction (SCM-seq)', 
         y = 'Variant allele frequency, CNV-corrected (WES)') +
    annotate(geom = 'text',
             label = paste0('Spearman rho = ', round(cor_test_vaf_corrected$estimate, digits = 2), 
                            ", p = ", round(cor_test_vaf_corrected$p.value, digits = 3)), 
             x = 0.3, y = 0.9,
             size = 4) +
    theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.subtitle = element_text(size = 17, hjust = 0.5),
          plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))
  
  ggsave(paste0(path_results, "all_samples/plots/", t, "cor_MCF_VAF_corrected.png"), p3, height = 5, width = 5)
  
  # MCF vs pbVAF
  cor_test_vaf_pbvaf <- cor.test(df_vaf_mcf$mcf, 
                                     df_vaf_mcf$pbVAF, 
                                     use = "everything", 
                                     method = "spearman")
  
  p4 <- ggplot(df_vaf_mcf, aes(mcf, pbVAF)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, linewidth = 0.5) +
    xlim(0,1) + ylim(0,1) +
    labs(x = 'Mutant cell fraction (SCM-seq)', 
         y = 'Pseudo-bulk VAF (SCM-seq)') +
    annotate(geom = 'text',
             label = paste0('Spearman rho = ', round(cor_test_vaf_pbvaf$estimate, digits = 2), 
                            ", p = ", round(cor_test_vaf_pbvaf$p.value, digits = 3)), 
             x = 0.3, y = 0.9,
             size = 4) +
    theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.subtitle = element_text(size = 17, hjust = 0.5),
          plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))
  
  ggsave(paste0(path_results, "all_samples/plots/", t, "cor_MCF_pbVAF.png"), p4, height = 5, width = 5)
  
}

### Genotype combinations and mutation burden ###########################################

# We want to score unique genotype combinations (i.e., "clones"),
# along with completeness of each combination 
# (i.e., fraction of known genotypes across all variants for each cell)
# and mutation burden (number of mutations per cell).

df_tier <- list()
for (t in unique(genotype_all$tier)) {
  df_sample <- list()
  for (s in samples_genotype) {
    
    selected_var <- genotype_all %>%
      filter(sample == s, 
             var_keep == "yes", 
             cell_keep == "yes",
             tier == t) %>%
      dplyr::select(barcode, variant, genotype, aggregated_lineage2) %>%
      mutate(genotype = recode(genotype,
                               'mutated' = 1,
                               'wild-type' = 0))
    
    # Plot variants by number of genotyped cells
    count_var <- as.data.frame(table(selected_var$variant)) %>%
      arrange(desc(Freq))
    count_var$Var1 <- factor(count_var$Var1, levels = count_var$Var1)
    
    p1 <- count_var %>%
      ggplot(aes(x = Var1, y = Freq)) +
      geom_col(width = 0.75) +
      theme_bw() +
      ggtitle(s) +
      ylab('N genotyped cells') +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2),
            axis.title.x = element_blank(),
            panel.grid = element_blank())
    
    ggsave(paste0(path_results, s, "/plots/", t, "_n_genotyped_cells.png"), p1, width = 6, height = 5)
    
    # Total n variants/cell and completeness of genotype information
    selected_var <- selected_var %>%
      pivot_wider(names_from = variant, values_from = genotype, values_fn = mean) %>%
      column_to_rownames(var = 'barcode')
    selected_var[is.na(selected_var)] <- 2
    selected_var <- selected_var %>%
      mutate(
        percent_known = rowSums(selected_var != 2)/ncol(selected_var), # fraction of variants with valid genotype
        n_mut = apply(selected_var == 1, 1, sum)
      )
    
    # Categorize N variants per cell into discrete bins
    selected_var <- selected_var %>% 
      mutate(
        n_mut_bin = case_when(
          selected_var$n_mut == 0 ~ "0",
          selected_var$n_mut == 1 ~ "1",
          selected_var$n_mut == 2 ~ "2",
          selected_var$n_mut == 3 ~ "3",
          selected_var$n_mut > 3 ~ ">3"
        ),
        tier = t,
        sample = s)
    
    # Add frequency of "clones" (identified by unique genotype combinations, including missing data)
    x <- selected_var %>%
      dplyr::select(-c(percent_known, n_mut, n_mut_bin, tier, sample, aggregated_lineage2)) %>%
      table() %>%
      as.data.frame(stringsAsFactors = F) %>%
      filter(Freq > 0) %>%
      arrange(desc(Freq))
    x[, 1:ncol(x)] <- sapply(x[, 1:ncol(x)], as.numeric)
    x$id_clone <- paste0('clone_', 1:nrow(x))
    
    selected_var <- selected_var %>% rownames_to_column(var = "barcode")
    selected_var <- full_join(selected_var, x)
    
    df_sample[[s]] <- selected_var[,c("barcode", "percent_known", "n_mut", "n_mut_bin", "tier", "sample")]
    
    write_csv(selected_var, paste0(path_results, s, "/tables/", t, "_clones.csv"))
    
    # plot cell counts by fraction of known genotypes
    tab_known <- table(selected_var$percent_known) %>% as.data.frame() %>% mutate(sample = s)
    colnames(tab_known) <- c("percent_known", "n_cells", "sample")
    tab_known$percent_known <- as.numeric(as.character(tab_known$percent_known))
    tab_known$percent_known <- round(tab_known$percent_known, digits = 2)
    tab_known$percent_known <- factor(tab_known$percent_known, levels = unique(tab_known$percent_known))
    
    p2 <- ggplot(tab_known, aes(x = percent_known, y = n_cells, label = n_cells)) +
      geom_col(width = 0.75) +
      theme_bw() +
      ggtitle(s) +
      xlab('Fraction of known genotypes by cell') + ylab('N cells') +
      geom_text(data = tab_known, nudge_y = 15) +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2),
            panel.grid = element_blank())
    ggsave(paste0(path_results, s, "/plots/", t, "_fraction_genotyped_cells.png"), p2, width = 6, height = 5)
  }
  df_tier[[t]] <- do.call(rbind, df_sample)
}

# Transfer genotype completeness information to genotype table
df_counts <- do.call(rbind, df_tier)

genotype_all <- full_join(genotype_all, df_counts, by = c("barcode", "sample", "tier"))
write_csv(genotype_all, paste0(path_results, "all_samples/tables/genotype_all.csv"))

### Summary statistics ###########################################

# Cells from 10x dataset after QC, n
n_cells_sample <- 
  table(final_seurat$sample[which(final_seurat$sample %in% samples_genotype)]) %>%
  as.data.frame(responseName = "n_cells_sample") %>%
  rename("sample" = "Var1")

for (t in unique(genotype_all$tier)) {
  
  # Genotyped cells, n 
  a <- genotype_all %>%
    filter(cell_keep == "yes" & 
             var_keep == "yes" &
             tier == t) %>%
    distinct(barcode, .keep_all = T) %>%
    group_by(sample) %>%
    summarise(n_genotyped_cells = n())
  
  # Total mutant cells (over genotyped cells)
  b <- genotype_all %>%
    filter(cell_keep == "yes" & 
             var_keep == "yes" & 
             tier == t) %>%
    distinct(barcode, .keep_all = T) %>%
    group_by(sample, genotype) %>%
    summarise(tot_mutant_cells = n()) %>%
    pivot_wider(names_from = genotype, values_from = tot_mutant_cells)
  
  # Mutations per cell, median (range)
  c <- genotype_all %>%
    filter(cell_keep == "yes" & 
             var_keep == "yes" & 
             tier == t) %>%
    distinct(barcode, .keep_all = T) %>%
    group_by(sample) %>%
    summarise(median_n_mut = median(n_mut), 
              min_n_mut = min(n_mut),
              max_n_mut = max(n_mut)) 
  
  # Mutant and wild-type cells per variant, median (range)
  d <- genotype_all %>%
    filter(cell_keep == "yes" & 
             var_keep == "yes" &
             tier == t) %>%
    group_by(sample, variant, genotype) %>%
    summarise(n_cells = n())
  d <- d %>%
    group_by(sample, genotype) %>%
    summarise(median_variant = median(n_cells), 
              min_n_cells_variant = min(n_cells),
              max_n_cells_variant = max(n_cells)) %>%
    pivot_wider(names_from = genotype, values_from = c(median_variant, min_n_cells_variant, max_n_cells_variant))
  
  final_tab <- Reduce(full_join, list(n_cells_sample, a,b,c,d))
  
  write_csv(final_tab, paste0(path_results, "all_samples/tables/", t, "_summary.csv"))
  
  ## pie chart
  genotype_all <- read_csv(paste0(path_results, "all_samples/tables/genotype_all.csv"))
  
  dat_ <- genotype_all %>%
    filter(cell_keep == "yes", 
           n_mut_bin != 0,
           tier == t) %>%
    dplyr::select(barcode, sample, n_mut_bin) %>%
    distinct() %>%
    group_by(n_mut_bin) %>%
    summarise(Freq = n())
  
  dat_$n_mut_bin <- factor(dat_$n_mut_bin, levels = c('1', '2', '3', '>3'))
  
  # set colors
  cols <- c("#FFEDA0", "#FD8D3C", "#E31A1C", "#800026")
  
  # pie
  pdf(paste0(path_results, "all_samples/plots/", t, "_pie_n_mutations.pdf"), width = 8, height = 10)
  
  pie(dat_$Freq, labels = c('1', '2', '3', '>3'), radius = 0.8, 
      clockwise = T, col = cols, border = T, cex = 2) 
  legend(x = "bottom", legend = c('1', '2', '3', '>3'), fill = cols,
         title = "N mutations/cell", cex = 1.5, horiz = T)
  dev.off()
  
}

# Transfer genotype completeness information to seurat object (Tier2)
geno_data_all <- list()
for (s in samples_genotype) {
  geno_data <- genotype_all %>% 
    filter(sample == s,
           tier == "Tier2",
           var_keep == "yes",
           cell_keep == "yes") %>% 
    mutate(barcode = paste0(barcode, '-', s)) %>%
    dplyr::select(c(barcode, sample, percent_known, n_mut, n_mut_bin)) %>%
    distinct(barcode, .keep_all = T)
  
  geno_data_all[[s]] <- geno_data
}

geno_data_all <- Reduce(full_join, geno_data_all)

final_seurat@meta.data <- left_join(final_seurat@meta.data, geno_data_all, by = c('barcode', 'sample'))
rownames(final_seurat@meta.data) <- final_seurat$barcode

final_seurat$n_mut_bin <- factor(final_seurat$n_mut_bin, levels = c("0", "1", "2", "3", ">3"))

saveRDS(final_seurat, paste0(path_data, "final_seurat.rds"))

# Plot UMAP of mutation burden
p = DimPlot(final_seurat[,which(!is.na(final_seurat$n_mut_bin))], 
        reduction = 'umap', 
        group.by = "n_mut_bin", 
        na.value = "NA", 
        cols = c("#99D8C9", "#FFEDA0", "#FD8D3C", "#E31A1C", "#800026"),
        order = rev(c("0", "1", "2", "3", ">3")),
        pt.size = 0.1,) +
  theme_bw() +
  labs(color = "N mutations/cell") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(1,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, "all_samples/plots/umap_nmut_bin.png"), p, width = 6, height = 5)

