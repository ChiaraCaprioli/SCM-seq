
## Differently from Milo, DCATs doesn't provide any heuristics to choose the "right" number of k

## Set up
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(DCATS)
  library(ggpubr)
  })

# set paths
path_main <- '/Users/ieo4874/Desktop/working/SCMseq/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/analysis/lineage/')

load(paste0(path_data, "settings.RData"))

# load Seurat object 
final_seurat <- readRDS(paste0(path_data, 'final_seurat.rds'))
final_seurat$sample <- factor(final_seurat$sample, levels = unique(final_seurat$sample))

## Test DA
### single lineages
contrast_list <- list(
  "SRSF2mut vs SRSF2wt" = list("SRSF2mut", "SRSF2wt"),
  "SRSF2mut vs healthy BM" = list("SRSF2mut", "hBM"),
  "SRSF2wt vs healthy BM" = list("SRSF2wt", "hBM")
)

L <- list()
for (c in names(contrast_list)) {
  seurat <- final_seurat[,which(final_seurat$cohort %in% c(contrast_list[[c]][[1]], contrast_list[[c]][[2]]))]
  
  seurat$sample <- 
    factor(seurat$sample,
           levels = levels(final_seurat$sample)[which(levels(final_seurat$sample) %in% unique(seurat$sample))])
  
  seurat@graphs$RNA_snn <- final_seurat@graphs$RNA_snn[colnames(seurat), colnames(seurat)]
  
  ct <- names(which(table(seurat$lineage) >= 5))
  
  knn_mat = knn_simMat(seurat@graphs$RNA_snn, seurat$lineage[which(seurat$lineage %in% ct)])
  
  count_mat <- table(as.character(seurat$sample[which(seurat$lineage %in% ct)]),
                     as.character(seurat$lineage[which(seurat$lineage %in% ct)]))
  count_mat <- count_mat[match(levels(seurat$sample), rownames(count_mat)),]

  design_mat = seurat@meta.data %>%
    dplyr::select(sample, cohort, age, sex) %>%
    distinct() 
  rownames(design_mat) <- NULL
  design_mat <- design_mat %>% column_to_rownames(var = "sample")
  
  res <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
  
  stat.test <- data.frame(
    group1 = contrast_list[[c]][[1]],
    group2 = contrast_list[[c]][[2]],
    p = res$LRT_pvals,
    fdr = res$fdr
  ) %>% 
    rownames_to_column("lineage")
  
  L[[c]] <- stat.test
}

test_all <- do.call(rbind, L)

write_csv(test_all, paste0(path_results, "tables/dcats_test_all.csv"))

### aggregated lineages
L <- list()
for (c in names(contrast_list)) {
  seurat <- final_seurat[,which(final_seurat$cohort %in% c(contrast_list[[c]][[1]], contrast_list[[c]][[2]]))]
  
  seurat$sample <- 
    factor(seurat$sample,
           levels = levels(final_seurat$sample)[which(levels(final_seurat$sample) %in% unique(seurat$sample))])
  
  seurat@graphs$RNA_snn <- final_seurat@graphs$RNA_snn[colnames(seurat), colnames(seurat)]
  
  ct <- names(which(table(seurat$aggregated_lineage2) >= 5))
  knn_mat = knn_simMat(seurat@graphs$RNA_snn, seurat$aggregated_lineage2[which(seurat$aggregated_lineage2 %in% ct)])
  
  count_mat <- table(as.character(seurat$sample[which(seurat$aggregated_lineage2 %in% ct)]),
                     as.character(seurat$aggregated_lineage2[which(seurat$aggregated_lineage2 %in% ct)]))
  count_mat <- count_mat[match(levels(seurat$sample), rownames(count_mat)),]
  
  design_mat = seurat@meta.data %>%
    dplyr::select(sample, cohort, age, sex) %>%
    distinct() 
  rownames(design_mat) <- NULL
  design_mat <- design_mat %>% column_to_rownames(var = "sample")
  
  res <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
  
  stat.test <- data.frame(
    group1 = contrast_list[[c]][[1]],
    group2 = contrast_list[[c]][[2]],
    p = res$LRT_pvals,
    fdr = res$fdr
  ) %>% 
    rownames_to_column("aggregated_lineage2")
  
  L[[c]] <- stat.test
}

test_all <- do.call(rbind, L)

write_csv(test_all, paste0(path_results, "tables/dcats_test_all_aggr_lineage.csv"))

## plot significance of change in proportions across cohorts
test_all <- read_csv(paste0(path_results, "tables/dcats_test_all_aggr_lineage.csv"))
test_all_cohort <- test_all[,1:7]
colnames(test_all_cohort)[7] <- "fdr"
test_all_cohort$fdr.signif <- case_when(
  test_all_cohort$fdr >= 0.1 ~ "ns",
  0.05 < test_all_cohort$fdr & test_all_cohort$fdr < 0.1 ~ "*",
  0.01 < test_all_cohort$fdr & test_all_cohort$fdr < 0.05 ~ "**",
  0.001 < test_all_cohort$fdr & test_all_cohort$fdr < 0.01 ~ "***",
  0.0001 < test_all_cohort$fdr & test_all_cohort$fdr < 0.001 ~ "****",
  0 < test_all_cohort$fdr & test_all_cohort$fdr < 0.0001 ~ "*****",
)
test_all_cohort$aggregated_lineage2 <- str_replace_all(test_all_cohort$aggregated_lineage2, "_", " ")

table_samples_by_lineage <- read_csv(paste0(path_results, 'tables/table_sample_by_aggr_lineage.csv')) 

tab_prop <- left_join(
  table_samples_by_lineage, 
  distinct(final_seurat@meta.data[,c("sample", "cohort")]), 
  by = "sample"
  ) %>%
  relocate("cohort", .after = "sample") %>%
  pivot_longer(cols = 4:ncol(.), names_to = "lineage", values_to = "count") %>%
  mutate(proportion = 100*(count / total_cell_count),
         lineage = str_replace_all(lineage, "_", " "))

L_plots <- list()
for (l in unique(tab_prop$lineage[(which(tab_prop$count >= 20))])) {
  p = ggbarplot(
    tab_prop[which(tab_prop$lineage == l),], 
    x = "cohort", 
    y = "proportion",
    fill = "cohort", 
    add = "point", 
    add.params = list(group = "cohort"),
    position = position_dodge(0.8)) + 
    theme_bw() +
    ylim(0,100) +
    ggtitle(l) +
    ylab("% cells") +
    scale_fill_manual(values = colors$cohort) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)),
                       breaks = seq(0, 100, by = 25)) +
    theme(aspect.ratio = 1.5,
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2, size = 18),
          axis.title.y = element_text(size = 17),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 17)) +
    stat_pvalue_manual(
      test_all_cohort[which(test_all_cohort$aggregated_lineage2 == l),], 
      label = "fdr.signif", 
      tip.length = 0.01,
      size = 5,
      bracket.shorten = 0.1,
      bracket.nudge.y = 0.1,
      step.increase = 0.05,
      step.group.by = "aggregated_lineage2",
      y.position = c(80, 85, 80)
    )
  L_plots[[l]] <- p
}

ggsave(paste0(path_results, "plots/dcat_prop_aggr_lineage_hsc_my.png"),
       L_plots[[1]] + L_plots[[2]] + L_plots[[7]] + 
         patchwork::plot_layout(ncol = 3, widths = 4, heights = 4, guides = "collect"),
       width = 10, height = 4
)

ggsave(paste0(path_results, "plots/dcat_prop_aggr_lineage_ly.png"),
       L_plots$`T cells` + L_plots$`NK cells` + L_plots$`B mature` + L_plots$`B early` +
         patchwork::plot_layout(ncol = 2, widths = 4, heights = 4, guides = "collect"),
       width = 12, height = 8
)

# plot higher granularity

table_cohort_by_lineage <- final_seurat@meta.data %>%
  filter(!aggregated_lineage2 %in% c("T_cells", "NK_cells", "B_early", "B_mature", "Stromal")) %>%
  group_by(cohort, lineage) %>%
  summarize(count = n()) %>%
  spread(lineage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cohort', 'total_cell_count', everything())) %>%
  arrange(factor(cohort, levels = levels(final_seurat@meta.data$cohort)))

table_cohort_by_lineage %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cohort') %>%
  filter() %>%
  mutate(cohort = factor(cohort, levels = levels(final_seurat@meta.data$cohort))) %>%
  ggplot(aes(cohort, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill', 
           stat = 'identity', 
           width = 0.8,
           color = "black") +
  scale_fill_manual(name = 'Lineage', values = colors$lineage[1:17]) +
  scale_y_continuous(name = "% HSCs, progenitors and myeloid cells", 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'on') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 30, vjust = 1, hjust = 1), 
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(paste0(path_results, 'plots/composition_cohort_lineage_hsc_my.png'), width = 8, height = 8)

table_cohort_by_lineage <- final_seurat@meta.data %>%
  filter(aggregated_lineage2 %in% c("T_cells", "NK_cells", "B_early", "B_mature")) %>%
  group_by(cohort, lineage) %>%
  summarize(count = n()) %>%
  spread(lineage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cohort', 'total_cell_count', everything())) %>%
  arrange(factor(cohort, levels = levels(final_seurat@meta.data$cohort)))

table_cohort_by_lineage %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cohort') %>%
  filter() %>%
  mutate(cohort = factor(cohort, levels = levels(final_seurat@meta.data$cohort))) %>%
  ggplot(aes(cohort, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill', 
           stat = 'identity', 
           width = 0.8,
           color = "black") +
  scale_fill_manual(name = 'Lineage', values = colors$lineage[18:41]) +
  scale_y_continuous(name = "% T, NK and B cells", 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 30, vjust = 1, hjust = 1), 
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(paste0(path_results, 'plots/composition_cohort_lineage_ly.png'), width = 10, height = 8)



