
library(ComplexHeatmap)

# get top n genes
k = "B_mature"
c = "SRSF2wt_vs_hBM"

top_up <- final_res %>%
  filter(pass == "pass" & 
           logFC > 0 &
           cluster_id == k & 
           contrast == c) %>%
  arrange(desc(logFC)) %>%
  slice(1:30)

top_down <- final_res %>%
  filter(pass == "pass" & 
           logFC < 0 &
           cluster_id == k & 
           contrast == c) %>%
  arrange(logFC) %>%
  slice(1:30)

top <- rbind(top_up, top_down) %>%
  pull(gene)

sub_sce <- sce[,which(sce$cluster_id == k & sce$group_id %in% c("SRSF2mut", "SRSF2wt"))]
exp_matrix <- as.matrix(logcounts(sub_sce)[top,])

L <- list()
for (s in unique(meta_col$sample)) {
  x <- rowMeans(exp_matrix[,sub_sce$sample_id == s])
  L[[s]] <- x              
}
exp_matrix <- do.call(rbind, L)
exp_matrix <- scale(exp_matrix)

meta_col <- data.frame(
  group = sub_sce$group_id,
  sample = sub_sce$sample_id
) %>%
  distinct()
ta <- columnAnnotation(
  Cohort = meta_col$group,
  Sample = meta_col$sample,
  show_annotation_name = F,
  show_legend = T
  )




Heatmap(
  t(exp_matrix),
  cluster_rows = T,
  cluster_columns = F,
  show_column_names = F,
  top_annotation = ta
)




