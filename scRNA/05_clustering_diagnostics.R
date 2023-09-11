
## Set up
library(tidyverse)
library(Seurat)
library(patchwork)
library(viridis)
library(SeuratObject)
library(clustree)
library(scran)
library(aricode)
library(cluster)
library(igraph)
library(bluster)
library(ComplexHeatmap)
library(log4r)

## set paths and settings
path_main <- '/hpcnfs/scratch/temporary/ccaprioli/'
path_results <- '/hpcnfs/scratch/temporary/ccaprioli/clustering_diagnostics/'
load(paste0(path_main, "settings.RData"))

## set log files
logfile = "log.txt"
file_appender = file_appender(logfile, 
                                 append = TRUE, 
                                 layout = default_log_layout())

my_logger <- logger(threshold = "INFO", appenders = file_appender)

## Load Seurat object 
final_seurat <- readRDS(paste0(path_main, 'seurat_k.rds'))

## clustering to evaluate
k_list <- c(5,10,15,30,50,100)
res_list <- seq(0.5, 2, by = 0.5)

## PART1: within a given k/res combination: evaluate clusters 
start_time <- Sys.time()
info(my_logger, paste0("Evaluating part 1"))

### silhouette
set.seed(123)

L <- list()
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    sil.approx <- bluster::approxSilhouette(
      Embeddings(final_seurat[['scanorama']]), 
      clusters = final_seurat@meta.data[[res]]
    )
    sil.data <- as.data.frame(sil.approx)
    sil.data$closest <- factor(ifelse(sil.data$width > 0, 
                                      final_seurat@meta.data[[res]], sil.data$other))
    sil.data$cluster <- as.factor(final_seurat@meta.data[[res]])
    
    p = ggplot(sil.data, aes(x = cluster, y=width, group = cluster)) +
      geom_boxplot() +
      geom_hline(yintercept = mean(sil.data$width), color = "red", linetype = 2) +
      theme_bw() +
      ggtitle(paste0("k_", k, "_res_", r)) +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_results, "k", k, "_res_", r, "_silhouette.png"),
           p, height = 8, width = 8)
  }
}

### purity
L <- list()
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    pur_data <- neighborPurity(  Embeddings(final_seurat[['scanorama']]), 
                                 clusters = final_seurat@meta.data[[res]])
    
    pur_data <- as.data.frame(pur_data)
    pur_data$maximum <- factor(pur_data$maximum)
    pur_data$cluster <- as.factor(final_seurat@meta.data[[res]])
    
    p = ggplot(pur_data, aes(x = cluster, y=purity, group = cluster)) +
      geom_boxplot() +
      geom_hline(yintercept = mean(pur_data$purity), color = "red", linetype = 2) +
      theme_bw() +
      ggtitle(paste0("k_", k, "_res_", r)) +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_results, "k", k, "_res_", r, "_purity.png"),
           p, height = 8, width = 8)
  }
}

### similarity among clusters 
sce <- as.SingleCellExperiment(final_seurat)
reducedDim(sce, 'scanorama') <- Embeddings(final_seurat[['scanorama']])[, 1:30]
g <- scran::buildSNNGraph(sce, use.dimred = 'scanorama')

for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    ratio <- bluster::pairwiseModularity(g, final_seurat@meta.data[[res]], as.ratio = TRUE)
    ratio_to_plot <- log10(ratio+1)
    x = ratio_to_plot %>%
      as_tibble() %>%
      rownames_to_column(var = 'cluster_1') %>%
      pivot_longer(
        cols = 2:ncol(.),
        names_to = 'cluster_2',
        values_to = 'probability') %>%
      mutate(
        cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
        cluster_2 = factor(cluster_2, levels = unique(cluster_2))) 
    p = ggplot(x, aes(cluster_2, cluster_1, fill = probability)) +
      geom_tile(color = 'white') +
      geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
      scale_x_discrete(name = 'Cluster', position = 'top') +
      scale_y_discrete(name = 'Cluster') +
      ggtitle(paste0("k_", k, "_res_", r)) +
      scale_fill_gradient(
        name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
        guide = guide_colorbar(
          frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
          title.theme = element_text(hjust = 1, angle = 90),
          barwidth = 0.75, barheight = 10)) +
      coord_fixed() +
      theme_bw() +
      theme(
        legend.position = 'right',
        panel.grid.major = element_blank())
    
    ggsave(paste0(path_results, "k", k, "_res_", r, "_cluster_similarity.png"), p, height = 6, width = 7)
  }
}

end_time <- Sys.time()

info(my_logger, 
     paste0("Part 1 successfully done in ", 
            floor(as.numeric(end_time - start_time, units = "mins")), " mins")
)

## PART2: within a given k: cluster trees across different res
start_time <- Sys.time()
info(my_logger, paste0("Evaluating part 2"))

for (k in k_list) {
  p1 <- clustree(final_seurat, 
                 prefix = paste0("RNA_snn_k_", k, "_res_"),
                 node_text_size = 0) +
    scale_color_manual(values = colors$discrete) +
    theme(legend.position = "right", 
          legend.text = element_text(size = 20))
  ggsave(paste0(path_results, "k", k, "_clustree_res.png"), p1, width = 20, height = 20)
}

end_time <- Sys.time()

info(my_logger, 
     paste0("Part 2 successfully done in ", 
            floor(as.numeric(end_time - start_time, units = "mins")), " mins")
)

## PART3: all k/res
start_time <- Sys.time()
info(my_logger, paste0("Evaluating part 3"))

### mean silhouette vs n clusters
cluster_all <- final_seurat@meta.data %>% dplyr::select(starts_with("RNA_snn_k"))

n_clusters <- apply(cluster_all,2,max) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "k_res") %>%
  dplyr::rename("n_clusters" = 2)

distance_matrix <- Rfast::Dist(Embeddings(final_seurat[['scanorama']])[, 1:30])

L <- list()
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    clusters <- final_seurat@meta.data[res]
    silhouette <- silhouette(as.numeric(clusters[,1]), dist = distance_matrix)
    final_seurat@meta.data$silhouette_score <- silhouette[,3]
    mean_silhouette_score <- mean(final_seurat@meta.data$silhouette_score)
    L[[res]] <- mean_silhouette_score
  }
}

mean_sil <- do.call(rbind, L) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "k_res")

df_k <- full_join(n_clusters, mean_sil, by = "k_res")
colnames(df_k) <- c("k_res", "n_clusters", "mean_sil")
df_k <- df_k %>% arrange(desc(mean_sil))
df_k$k_res <- factor(df_k$k_res)

p1 <- ggplot(df_k, aes(reorder(k_res, rev(mean_sil)), mean_sil)) +
  geom_point() +
  geom_line(group = "resolution") +
  ylab('Mean silhouette score') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, color = 'black', angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))

p2 <- ggplot(df_k, aes(reorder(k_res, rev(mean_sil)), n_clusters)) +
  geom_point() +
  geom_line(group = "resolution") +
  ylab('N clusters') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, color = 'black', angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))

ggsave(paste0(path_results, "mean_sil_n_clusters.png"),
       p1 + p2 +
         plot_layout(ncol = 1, heights = 4, widths = 4),
       height = 8, width = 8)

### ARI
my_clusters <- final_seurat@meta.data %>% dplyr::select(starts_with("RNA_snn_k_"))

res <- matrix(nrow = ncol(my_clusters), ncol = ncol(my_clusters))

for (i in (1:ncol(my_clusters))) {
  
  for (j in (i:ncol(my_clusters))){
    
    res[i,j] <- res[j,i] <- pairwiseRand(my_clusters[,i], my_clusters[,j], mode="index")
    
  }
}

colnames(res) <- colnames(my_clusters)
rownames(res) <- colnames(my_clusters)

pdf(paste0(path_results, "ARI_heat.pdf"), width = 20, height = 20)
Heatmap(res, 
        name = "ARI",
        cluster_rows = T,
        cluster_columns = T,
        show_row_dend = F,
        show_column_dend = F,
        show_row_names = T,
        row_names_side = "left",
        show_column_names = T,
        column_names_side = "top",
        col = viridisLite::mako(n = 50, direction = 1),
        #rect_gp = gpar(col = "white", lwd = 0.5),
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 4),
        width = unit(15, "cm"), height = unit(15, "cm")
)
dev.off()

end_time <- Sys.time()

info(my_logger, 
     paste0("Part 3 successfully done in ", 
            floor(as.numeric(end_time - start_time, units = "mins")), " mins")
)

## PART4: relationship between clustering and independent lineage assignment
start_time <- Sys.time()
info(my_logger, paste0("Evaluating part 4"))

### fine lineage
#### heatmap
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    
    mat_cl_lin <- final_seurat@meta.data %>% dplyr::select(c(res, lineage))
    mat_cl_lin <- table(as.character(mat_cl_lin[[res]]), as.character(mat_cl_lin[["lineage"]]))
    mat_cl_lin <- t(mat_cl_lin)/apply(mat_cl_lin,2,sum) # fraction of cells in the same lineage across clusters
    
    pdf(paste0(path_results, "k", k, "_res_", r, "_cluster_lineage.pdf"), width = 12, height = 12)
    h = Heatmap(mat_cl_lin, 
                name = "% cells",
                column_title = res,
                cluster_rows = T,
                cluster_columns = T,
                show_row_dend = F,
                show_column_dend = F,
                show_row_names = T,
                row_names_side = "left",
                show_column_names = T,
                column_names_side = "bottom",
                col = viridisLite::mako(n = 50, direction = -1),
                #rect_gp = gpar(col = "white", lwd = 0.5),
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                column_names_rot = 90,
                width = unit(10, "cm"), height = unit(10, "cm")
    )
    print(h)
    dev.off()
  }
}

#### ARI
x <- list()
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    
    df_ARI <- data.frame(
      cluster = unlist(final_seurat[[res]]),
      lineage = final_seurat$lineage
    )
    ari <- ARI(df_ARI$cluster, df_ARI$lineage)
    x[[res]] <- ari
  }
}

x <- do.call(rbind, x) %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster") %>%
  dplyr::rename("ARI" = "V1")

ggplot(x, aes(x = reorder(cluster, ARI), y = ARI)) +
  geom_point() +
  geom_line(group = "ARI") +
  theme_bw() +
  ggtitle("Agreement between clustering and cell lineage") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.7,
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2))

ggsave(paste0(path_results, "ari_cluster_lineage.png"),
       width = 6, height = 4)

### aggregated lineage
#### heatmap
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    
    mat_cl_lin <- final_seurat@meta.data %>% dplyr::select(c(res, aggregated_lineage2))
    mat_cl_lin <- table(as.character(mat_cl_lin[[res]]), as.character(mat_cl_lin[["aggregated_lineage2"]]))
    mat_cl_lin <- t(mat_cl_lin)/apply(mat_cl_lin,2,sum) # fraction of cells in the same lineage across clusters
    
    pdf(paste0(path_results, "k", k, "_res_", r, "_cluster_aggr_lineage.pdf"), width = 12, height = 12)
    h = Heatmap(mat_cl_lin, 
                name = "% cells",
                column_title = res,
                cluster_rows = T,
                cluster_columns = T,
                show_row_dend = F,
                show_column_dend = F,
                show_row_names = T,
                row_names_side = "left",
                show_column_names = T,
                column_names_side = "bottom",
                col = viridisLite::mako(n = 50, direction = -1),
                #rect_gp = gpar(col = "white", lwd = 0.5),
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                column_names_rot = 90,
                width = unit(10, "cm"), height = unit(10, "cm")
    )
    print(h)
    dev.off()
  }
}

#### ARI
x <- list()
for (k in k_list) {
  for (r in res_list) {
    res = paste0("RNA_snn_k_", k, "_res_", r)
    
    df_ARI <- data.frame(
      cluster = unlist(final_seurat[[res]]),
      lineage = final_seurat$aggregated_lineage2
    )
    ari <- ARI(df_ARI$cluster, df_ARI$lineage)
    x[[res]] <- ari
  }
}

x <- do.call(rbind, x) %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster") %>%
  dplyr::rename("ARI" = "V1")

ggplot(x, aes(x = reorder(cluster, ARI), y = ARI)) +
  geom_point() +
  geom_line(group = "ARI") +
  theme_bw() +
  ggtitle("Agreement between clustering and cell lineage") +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.7,
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2))

ggsave(paste0(path_results, "ari_cluster_aggr_lineage.png"),
       width = 6, height = 4)


end_time <- Sys.time()

info(my_logger, 
     paste0("Part 4 successfully done in ", 
            floor(as.numeric(end_time - start_time, units = "mins")), " mins")
)

# TO DO: stability by bootstrapping

write.table(readLines(logfile), 
            file = paste0(path_main, logfile), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)
