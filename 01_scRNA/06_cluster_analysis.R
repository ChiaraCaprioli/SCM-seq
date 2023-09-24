---
title: "05_cluster_analysis"
author: "Chiara Caprioli"
date: "July 6th 2023"
---

**Aims:**
1. pick chosen clustering k/resolution and run new kNN/UMAP according to chosen k
2. QC by cluster and remove residual bad cells
3. save high-quality list of CBs with metadata
4. new embedding without bad cells 
5. composition by clusters and cluster composition according to different variables (sample, cohort, lineage, cell cycle)
6. compute marker genes by cluster and perform functional annotation

## Set up
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(viridis)
  library(SeuratObject)
  library(enrichR)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  })

# set paths
path_main <- '/Users/ieo4874/Desktop/working/SCMseq/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/analysis/cluster/')

load(paste0(path_data, "settings.RData"))

## Load Seurat object with clustering
seurat_k <- readRDS(paste0(path_data, 'seurat_k.rds'))

## Load Seurat object and attach chosen clustering resolution
final_seurat <- readRDS(paste0(path_data, "final_seurat.rds"))
final_seurat$RNA_snn_k_100_res_0.5 <- seurat_k$RNA_snn_k_100_res_0.5

## Run UMAP using same k as that used for preferred clustering solution (k = 50, res = 0.5).
# This was chosen based on mean silhouette and ARI between clusters and aggregated lineage.
# At this stage we are mostly interested in finding a suitable SNN, and not in the fine granularity 
# of the clusters (which instead will be investigated by subclustering)

final_seurat <- FindNeighbors(
  final_seurat,
  reduction = "scanorama",
  dims = 1:30,
  k.param = 100,
  do.plot = F,
  compute.SNN = TRUE
)

final_seurat <- RunUMAP(
  final_seurat,
  reduction = "scanorama",
  dims = 1:30,
  n.components = 2,
  seed.use = 100
)

## Plot UMAP of clusters
p1 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'RNA_snn_k_100_res_0.5', 
              cols = colors$discrete, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Cluster") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_clusters.png'), p1,
       width = 6, height = 5)

DimPlot(final_seurat, reduction = 'umap', 
        group.by = 'aggregated_lineage2', 
        cols = colors$aggregated_lineage2, pt.size = 0.1, 
        shuffle = T, seed = 123) +
  theme_bw() +
  #labs(color = "Cluster") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

## Set identities based on chosen clustering solution
Idents(final_seurat) <- final_seurat$RNA_snn_k_100_res_0.5
final_seurat$RNA_snn_k_100_res_0.5 <- factor(final_seurat$RNA_snn_k_100_res_0.5,
                                             levels = min(final_seurat$RNA_snn_k_100_res_0.5):max(final_seurat$RNA_snn_k_100_res_0.5))

### QC metrics
temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally()

p1 <- ggplot() +
  gghalves::geom_half_violin(
    data = final_seurat@meta.data, aes(RNA_snn_k_100_res_0.5, 
                                       nCount_RNA, 
                                       fill = RNA_snn_k_100_res_0.5),
    side = "l", show.legend = FALSE, trim = FALSE) +
  gghalves::geom_half_boxplot(
    data = final_seurat@meta.data, aes(RNA_snn_k_100_res_0.5, nCount_RNA, fill = RNA_snn_k_100_res_0.5),
    side = "r", outlier.color = NA, width = 0.4, show.legend = FALSE) +
  geom_text(
    data = temp_labels,
    aes(x = RNA_snn_k_100_res_0.5, y = -Inf, label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)), vjust = -1),
    color = "black", size = 2) +
  scale_fill_manual(name = 'Cluster', values = colors$discrete) +
  scale_y_continuous(name = "Number of transcripts", labels = scales::comma,
                     expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = 'black'))

p2 <- ggplot() +
  gghalves::geom_half_violin(
    data = final_seurat@meta.data, aes(RNA_snn_k_100_res_0.5, nFeature_RNA, fill = RNA_snn_k_100_res_0.5),
    side = "l", show.legend = FALSE, trim = FALSE) +
  gghalves::geom_half_boxplot(
    data = final_seurat@meta.data, aes(RNA_snn_k_100_res_0.5, nFeature_RNA, fill = RNA_snn_k_100_res_0.5),
    side = "r", outlier.color = NA, width = 0.4, show.legend = FALSE) +
  geom_text(
    data = temp_labels,
    aes(x = RNA_snn_k_100_res_0.5, y = -Inf, label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)), vjust = -1),
    color = "black", size = 2) +
  scale_y_continuous(name = "Number of expressed genes", labels = scales::comma,
                     expand = c(0.08, 0)) +
  scale_fill_manual(name = 'Cluster', values = colors$discrete) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = 'black'))

ggsave(paste0(path_results, 'plots/n_transcripts_ex_genes_by_cluster.png'), p1 + p2 + 
         plot_layout(ncol = 1, nrow = 2, widths = 10, heights = 3, guides = 'collect'),
       width = 10, height = 6)

## Remove bad cells
# cluster 12: few transcripts; stromal or inconsistent lineage; mostly specific to sAML1.

# save bad cells
cells_cluster_12 <- colnames(final_seurat)[which(final_seurat$RNA_snn_k_100_res_0.5 == '12')]  
saveRDS(cells_cluster_12, paste0(path_data, 'bad_cells_cluster_12.rds'))

# discard cells from seurat object
final_seurat <- final_seurat[,which(!colnames(final_seurat) %in% cells_cluster_12)]

# save polished seurat object
saveRDS(final_seurat, paste0(path_data, 'final_seurat.rds'))

# save high-quality list of CBs with metadata
df_meta <- final_seurat@meta.data
write_csv(df_meta, paste0(path_data, 'high_quality_cb_meta.csv')) 

# new embedding without bad cells
final_seurat <- FindNeighbors(
  final_seurat,
  reduction = "scanorama",
  dims = 1:30,
  k.param = 50,
  do.plot = F,
  compute.SNN = TRUE
)

final_seurat <- RunUMAP(
  final_seurat,
  reduction = "scanorama",
  dims = 1:30,
  n.components = 2,
  seed.use = 100
)

saveRDS(final_seurat, paste0(path_data, "final_seurat.rds"))

## Plot wild UMAPs
p1 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'sample', 
              cols = colors$sample, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Sample") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_sample.png'), p1,
       width = 6, height = 5)

p2 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'seq_run', 
              cols = rev(colors$discrete), pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "seq_run") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_seq_run.png'), p2,
       width = 6, height = 5)

p3 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'cohort', 
              cols = colors$cohort, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Cohort") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_cohort.png'), p3,
       width = 6, height = 5)

p4 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'aggregated_lineage2', 
              cols = colors$aggregated_lineage2, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Aggregated lineage") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_aggr_lineage.png'), p4,
       width = 6, height = 5)

p5 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'lineage', 
              cols = colors$lineage, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Lineage") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_lineage.png'), p5,
       width = 10, height = 5)

p6 <- FeaturePlot(final_seurat, 
                  reduction = 'umap', 
                  features = 'pt.Myelocytes', 
                  pt.size = 0.1) +
  theme_bw() +
  labs(color = "Myeloid pseudotime") +
  scale_color_viridis_c(option = "magma") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_pt_myeloid.png'), p6,
       width = 6, height = 5)

p7 <- DimPlot(final_seurat, reduction = 'umap', 
              group.by = 'RNA_snn_k_100_res_0.5', 
              cols = colors$discrete, pt.size = 0.1, 
              shuffle = T, seed = 123) +
  theme_bw() +
  labs(color = "Cluster") +
  coord_fixed(ratio = 1, clip = "on") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        legend.position="right",
        legend.justification="left", 
        legend.box.spacing = unit(5, "pt"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

ggsave(paste0(path_results, 'plots/umap_clusters_clean.png'), p7,
       width = 6, height = 5)

## Abundance
### Sample/cohort composition by clusters

# percentage sample
table_samples_by_clusters <- final_seurat@meta.data %>%
  group_by(sample, RNA_snn_k_100_res_0.5) %>%
  summarize(count = n()) %>%
  spread(RNA_snn_k_100_res_0.5, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('sample', 'total_cell_count', everything())) %>%
  arrange(factor(sample, levels = levels(final_seurat@meta.data$sample)))

write_csv(table_samples_by_clusters, 
          paste0(path_results, 'tables/table_samples_by_clusters.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(sample) %>%
  tally()

p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(final_seurat$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill', 
           color = "black",
           linewidth = 0.2,
           stat = 'identity', 
           width = 0.75) +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 3
  ) +
  scale_fill_manual(name = 'Cluster', values = colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 30, vjust = 1, hjust = 1), 
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_samples_clusters.png'),
       p1, width = 8, height = 4)

# percentage cohort
table_conditions_by_clusters <- final_seurat@meta.data %>%
  group_by(cohort, RNA_snn_k_100_res_0.5) %>%
  summarize(count = n()) %>%
  spread(RNA_snn_k_100_res_0.5, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cohort', 'total_cell_count', everything())) %>%
  arrange(factor(cohort, levels = levels(final_seurat@meta.data$cohort)))

write_csv(table_conditions_by_clusters, 
          paste0(path_results, 'tables/table_conditions_by_clusters.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(cohort) %>%
  tally()

p2 <- table_conditions_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cohort') %>%
  mutate(cohort = factor(cohort, levels = levels(final_seurat@meta.data$cohort))) %>%
  ggplot(aes(cohort, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill', 
           stat = 'identity', 
           width = 0.75,
           color = "black",
           linewidth = 0.2) +
  geom_text(
    data = temp_labels,
    aes(x = cohort, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 3
  ) +
  scale_fill_manual(name = 'Cluster', values = colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 30, vjust = 1, hjust = 1), 
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_condition_clusters.png'),
       p2, width = 5, height = 4)


### Cluster composition by sample and cohort

# percentage sample
table_clusters_by_samples <- final_seurat@meta.data %>%
  dplyr::rename('cluster' = 'RNA_snn_k_100_res_0.5') %>%
  group_by(cluster, sample) %>%
  summarize(count = n()) %>%
  spread(sample, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5)))

write_csv(table_clusters_by_samples, paste0(path_results, 'tables/table_clusters_by_samples.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally() %>%
  dplyr::rename('cluster' = RNA_snn_k_100_res_0.5)

p1 <- table_clusters_by_samples %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity',
           color = "black", linewidth = 0.2) +
  geom_text(
    data = temp_labels, 
    aes(x = cluster, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    color = 'black', size = 2.5, angle = 90, hjust = -0.1
  ) +
  scale_fill_manual(name = 'Sample', values = colors$sample) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 40, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_clusters_sample.png'),
       p1, width = 8, height = 4)

# percentage cohort
table_clusters_by_conditions <- final_seurat@meta.data %>%
  dplyr::rename('cluster' = 'RNA_snn_k_100_res_0.5') %>%
  group_by(cluster, cohort) %>%
  summarize(count = n()) %>%
  spread(cohort, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5)))

write_csv(table_clusters_by_conditions, paste0(path_results, 'tables/table_clusters_by_conditions.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally() %>%
  dplyr::rename('cluster' = RNA_snn_k_100_res_0.5)

p2 <- table_clusters_by_conditions %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity',
           color = "black", linewidth = 0.2) +
  geom_text(
    data = temp_labels, 
    aes(x = cluster, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    color = 'black', size = 2.5, angle = 90, hjust = -0.1
  ) +
  scale_fill_manual(name = 'Condition', values = colors$cohort) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 40, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_clusters_cohort.png'),
       p2, width = 8, height = 4)


### Cluster composition by cell lineage and cell cycle

# percentage lineage
table_clusters_by_lineage <- final_seurat@meta.data %>%
  dplyr::rename('cluster' = 'RNA_snn_k_100_res_0.5') %>%
  group_by(cluster, lineage) %>%
  summarize(count = n()) %>%
  spread(lineage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5)))

write_csv(table_clusters_by_lineage, paste0(path_results, 'tables/table_clusters_by_lineage.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally() %>%
  dplyr::rename('cluster' = RNA_snn_k_100_res_0.5)

p3 <- table_clusters_by_lineage %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity',
           color = "black", linewidth = 0.2, width = 0.75) +
  geom_text(
    data = temp_labels, 
    aes(x = cluster, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    color = 'black', size = 3, angle = 90, hjust = -0.1
  ) +
  scale_fill_manual(name = 'Lineage', values = colors$lineage) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 40, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_clusters_lineage.png'),
       p3, width = 16, height = 6)

# percentage aggregated lineage
table_clusters_by_aggregated_lineage <- final_seurat@meta.data %>%
  dplyr::rename('cluster' = 'RNA_snn_k_100_res_0.5') %>%
  group_by(cluster, aggregated_lineage2) %>%
  summarize(count = n()) %>%
  spread(aggregated_lineage2, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5)))

write_csv(table_clusters_by_aggregated_lineage, paste0(path_results, 'tables/table_clusters_by_aggregated_lineage.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally() %>%
  dplyr::rename('cluster' = RNA_snn_k_100_res_0.5)

p3 <- table_clusters_by_aggregated_lineage %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity',
           color = "black", linewidth = 0.2, width = 0.75) +
  geom_text(
    data = temp_labels, 
    aes(x = cluster, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    color = 'black', size = 3, angle = 90, hjust = -0.1
  ) +
  scale_fill_manual(name = 'Aggregated lineage', values = colors$aggregated_lineage2) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 40, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_clusters_aggregated_lineage.png'),
       p3, width = 8, height = 4)

# percentage cell cycle
table_clusters_by_cell_cycle <- final_seurat@meta.data %>%
  dplyr::rename('cluster' = 'RNA_snn_k_100_res_0.5') %>%
  group_by(cluster, cell_cycle_seurat) %>%
  summarize(count = n()) %>%
  spread(cell_cycle_seurat, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5)))

write_csv(table_clusters_by_cell_cycle, paste0(path_results, 'tables/table_clusters_by_cell_cycle.csv'))

temp_labels <- final_seurat@meta.data %>%
  group_by(RNA_snn_k_100_res_0.5) %>%
  tally() %>%
  dplyr::rename('cluster' = RNA_snn_k_100_res_0.5)

p5 <- table_clusters_by_cell_cycle %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(final_seurat@meta.data$RNA_snn_k_100_res_0.5))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity',
           color = "black", linewidth = 0.2) +
  geom_text(
    data = temp_labels, 
    aes(x = cluster, y = Inf, 
        label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    color = 'black', size = 3, angle = 90, hjust = -0.1
  ) +
  scale_fill_manual(name = 'Cell cycle phase', values = colors$cell_cycle) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 13, color = 'black'),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 40, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(path_results, 'plots/composition_clusters_cell_cycle.png'),
       p5, width = 8, height = 4)

## Functional annotation
### Marker genes by chosen resolution
marker_cluster <- FindAllMarkers(
  final_seurat, assay = "RNA",
  test.use = "wilcox", slot = "data",
  only.pos = F, random.seed = 123,
  min.pct = 0.25, min.diff.pct = -Inf,
  logfc.threshold = 0.5
)

write_csv(marker_cluster, paste0(path_results, 'tables/marker_cluster.csv'))

# check number of significant markers per cluster
marker_list <- list()
for (i in unique(marker_cluster$cluster)) {
  x <- marker_cluster %>%
    filter(cluster == i & p_val_adj < 0.05) %>%
    nrow()
  df_ <- data.frame(
    cluster = i,
    n_sig_markers = x
  )
  marker_list[[i]] <- df_
}

do.call(rbind, marker_list)

# select top 200 positive markers by cluster for functional analysis
selected_markers_top <- marker_cluster %>% 
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  dplyr::slice(1:200)

selected_markers_top$cluster <- factor(selected_markers_top$cluster, levels = unique(selected_markers_top$cluster))

write_csv(selected_markers_top, paste0(path_results, 'tables/top200_marker_cluster.csv'))

### ORA
# select databases of interest
dbs_selected <- c("GO_Molecular_Function_2021",
                  "GO_Cellular_Component_2021",
                  "GO_Biological_Process_2021",
                  "KEGG_2021_Human",
                  "MSigDB_Hallmark_2020",
                  "Reactome_2022",
                  "WikiPathway_2021_Human"
)

# query
for (c in levels(selected_markers_top$cluster)) {
  
  dat_ <- selected_markers_top %>% filter(cluster == c)
    
    if (length(unique(dat_$gene)) > 2) {
      enriched <- enrichr(dat_$gene, dbs_selected)
      enriched <- do.call(rbind, enriched) %>% 
        rownames_to_column(var = "Database") %>%
        mutate(cluster = c) %>%
        filter(Adjusted.P.value < 0.05) %>%
        arrange(Adjusted.P.value)
      write_csv(enriched, paste0(path_results, "tables/enrichr_up_", c, ".csv")) 
    } else {
      next
    }
  }

# plot top50 up
for (i in x) {
  for (c in unique(selected_markers_top$cluster)) {
    top_up <- read_csv(paste0(path_results, "tables/enrichr_up_", c, ".csv")) %>%
      filter(Adjusted.P.value < 0.05) %>%
      arrange(Adjusted.P.value) %>%
      dplyr::slice(1:50) %>%
      mutate(cluster = c,
             Adjusted.P.value = -log10(Adjusted.P.value)) %>%
      distinct(Term, .keep_all = T)
    
    top_up$Term <- factor(top_up$Term, levels = top_up$Term)
    
    p = ggplot(top_up, aes(x = Adjusted.P.value, y = Term)) +
      geom_col() +
      theme_bw() +
      scale_y_discrete(limits=rev) +
      labs(title = paste0("Cluster ", c)) +
      labs(x = "-log10(FDR)") +
      theme(axis.text = element_text(size = 8, color = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 13),
            panel.grid = element_blank(),
            plot.title = element_text(size = 13, hjust = 0.5, vjust = 2, face = "bold")
            )
    
    ggsave(paste0(path_results, "plots/cluster_", c, "_top50_path_up.png"), p, width = 10, height = 10)  
  }
} 

### GSEA
HALLMARK <- msigdbr(species = 'Homo sapiens', category = 'H') %>%
  dplyr::select(gs_name, gene_symbol)

GO <- msigdbr(species = 'Homo sapiens', category = 'C5') %>%
  dplyr::select(gs_name, gene_symbol)

ONCOGENIC <- msigdbr(species = 'Homo sapiens', category = 'C6') %>%
  dplyr::select(gs_name, gene_symbol)

IMMUNESIGDB <- msigdbr(species = 'Homo sapiens', category = 'C7') %>%
  dplyr::select(gs_name, gene_symbol)

TERM2GENE_list <- list(HALLMARK, GO, ONCOGENIC, IMMUNESIGDB)
names(TERM2GENE_list) <- c("HALLMARK", "GO", "ONCOGENIC", "IMMUNESIGDB")

gsea_all <- list()
gsea_all2 <- list()
for (c in levels(marker_cluster$cluster)) {
  
  # get ranked gene list
  x <- marker_cluster %>% filter(cluster == c)
  original_gene_list <- x$avg_log2FC
  names(original_gene_list) <- x$gene
  gene_list <- na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  for (i in names(TERM2GENE_list)) {
    
    # log message
    message(paste0(format(Sys.time(), '[%Y-%m-%d %H:%M:%S]'), 
                   " Doing GSEA by ", i, " for cluster ", c))
    
    # run GSEA
    gse <- GSEA(
      gene_list, TERM2GENE = TERM2GENE_list[[i]],
      minGSSize = 5, maxGSSize = 800,
      eps = 1e-10, pvalueCutoff = 0.1,
      pAdjustMethod = "fdr",
      by = "fgsea", scoreType = "pos",
      verbose = F
    )
    
    # save
    if (nrow(gse@result != 0)) {
      
      gse@result$cluster <- c
      gsea_all[[i]] <-  gse@result
      
    } else {
      next
    }
  }
  
  gsea_all_df <- do.call(rbind, gsea_all)
  gsea_all2[[c]] <-  gsea_all_df
}

gsea_complete <- do.call(rbind, gsea_all2)

write_csv(gsea_complete, paste0(path_results, 'tables/gsea_all.csv'))

sessionInfo()
