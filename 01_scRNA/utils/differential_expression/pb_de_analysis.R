
# Differential state analysis 
# Chiara Caprioli
# Nov 2nd 2023

# Load libraries
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(ComplexUpset)
library(cowplot)

# run_pb_dea
## Given some cluster identities shared across groups, performs DSA for each cluster across groups.
## Adapted from Cromwell et al., https://doi.org/10.1038/s41467-020-19894-4 (muscat r package).

run_pb_dea <- function(
  sce = sce,
  cluster_id = cluster_id,
  group_id = group_id,
  sample_id = sample_id,
  plotMDS = TRUE,
  path_save = path_save,
  method = method,
  min_cells = min_cells, # minimum n cells in a given cluster-sample to consider the sample for testing
  p.adj = p.adj, # thr for pass
  logFC = logFC, # thr for pass 
  min_freq = min_freq, # thr for pass 
  make_plots = TRUE
) {
  
  # remove lowly expressed genes
  sce <- sce[rowSums(counts(sce) >= 1) >= min_cells, ]
  
  # compute sum-factors & normalize
  sce <- computeLibraryFactors(sce)
  sce <- logNormCounts(sce)
  
  # set variables
  sce$cluster_id <- as.character(sce[[cluster_id]])
  sce$group_id <- as.character(sce[[group_id]])
  sce$sample_id <- as.character(sce[[sample_id]])
  
  # prep
  sce <- prepSCE(
    sce, 
    kid = "cluster_id", 
    gid = "group_id",  
    sid = "sample_id",   
    drop = FALSE
  )  
  
  nk <- length(kids <- levels(sce$cluster_id))
  ns <- length(sids <- levels(sce$sample_id))
  names(kids) <- kids; names(sids) <- sids
  
  sce$group_id <- factor(sce$group_id, levels = levels(sce[[group_id]]))
  sce$cluster_id <- factor(
    sce$cluster_id, 
    levels = levels(sce[[cluster_id]])[which(levels(sce[[cluster_id]]) %in% unique(sce$cluster_id))]
  )
  sce$sample_id <- factor(sce$sample_id, levels = levels(sce[[sample_id]]))
  
  # aggregate counts for each sample in each cluster
  pb <- aggregateData(
    sce,
    assay = "counts", fun = "sum",
    by = c("cluster_id", "sample_id")
  )
  
  # pseudobulk-level MDS plot
  if (plotMDS == TRUE) {
    pb_mds <- pbMDS(pb) +
      scale_color_manual(values = colors[[cluster_id]]) +
      theme_bw() +
      theme(
        panel.grid = element_blank()
      )
    ggsave(paste0(path_save, "plots/pb_mds.png"), pb_mds, width = 6, height = 5)
  }
  
  # sample-level analysis
  ## construct design matrix
  ei <- metadata(sce)$experiment_info
  ei$group_id <- factor(ei$group_id, levels = levels(sce[[group_id]]) )
  mm <- model.matrix(~ 0 + ei$group_id)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
  
  # to do: add other covariates, such as age
  
  ## construct contrast matrix
  source(paste0("~/Desktop/working/SCMseq/scripts/utils/differential_expression/make_all_contrasts.R"))
  contrast <- make_all_contrasts(colnames(mm), levels(sce[[group_id]]))
  
  ## run DS analysis
  res <- pbDS(
    pb, 
    design = mm, 
    contrast = contrast,
    min_cells = min_cells,
    filter = "none",
    method = method
  ) 
  
  ## tidy format 
  final_res <- resDS(sce, res, bind = "row")
  final_res$cluster_id <- factor(final_res$cluster_id, levels = levels(sce$cluster_id))
  final_res$contrast <- factor(final_res$contrast, levels = unique(final_res$contrast))
  
  ## compute frequencies 
  deg_list <- row.names(sce) %in% final_res$gene
  
  ### by group 
  L <- list()
  for (i in levels(sce$group_id)) {
    frq <- rowSums(counts(sce)[deg_list, which(sce$group_id == i)] > 0) / length(which(sce$group_id == i))
    L[[i]] <- frq
  }
  frq_tbl_group <- t(do.call(rbind, L)) %>%
    as.data.frame() 
  colnames(frq_tbl_group) <- paste0(colnames(frq_tbl_group), "_freq")
  frq_tbl_group$gene <- rownames(frq_tbl_group)
  
  ### by sample
  L <- list()
  for (i in levels(sce$sample_id)) {
    frq <- rowSums(counts(sce)[deg_list, which(sce$sample_id == i)] > 0) / length(which(sce$sample_id == i))
    L[[i]] <- frq
  }
  frq_tbl_sample <- t(do.call(rbind, L)) %>%
    as.data.frame() 
  colnames(frq_tbl_sample) <- paste0(colnames(frq_tbl_sample), "_freq")
  frq_tbl_sample$gene <- rownames(frq_tbl_sample)
  
  final_res <- Reduce(full_join, list(final_res, frq_tbl_group, frq_tbl_sample))
  
  ## mark passing DEGs
  # keep genes with interesting p.adj and logFC, and expressed in
  # at least 10% in at least 2 samples of either group within the appropriate contrast
  L_c_k <- list()
  for (k in levels(final_res$cluster_id)) {
    
    L_c <- list()
    for (c in levels(final_res$contrast)) {
      
      data_c_k <- final_res[which(
        final_res$p_adj.loc < p.adj & 
          abs(final_res$logFC) > logFC &
          final_res$cluster_id == k &
          final_res$contrast == c
      ),]
      
      group1 <- str_split_i(c, pattern = "_vs_", 1)
      group2 <- str_split_i(c, pattern = "_vs_", 2)
      
      samples_g1 <- as.character(unique(sce$sample_id[which(sce$group_id == group1)]))
      samples_g2 <- as.character(unique(sce$sample_id[which(sce$group_id == group2)]))
      
      data_c_k <- data_c_k %>%
        dplyr::select(c(gene, contains(c(samples_g1, samples_g2))))
      
      data_c_k$pass_g1 <- rowSums(data_c_k[,paste0(samples_g1, "_freq")] > min_freq) 
      data_c_k$pass_g2 <- rowSums(data_c_k[,paste0(samples_g2, "_freq")] > min_freq) 
      
      data_c_k$pass <- ifelse(
        data_c_k$pass_g1 >= 2 | data_c_k$pass_g2 >= 2,
        "pass", "no_pass"
      )
      
      data_c_k <- data_c_k %>%
        dplyr::select(c(gene, pass)) %>%
        mutate(contrast = c)
      
      L_c[[c]] <- data_c_k
    }
    
    res <- do.call(rbind, L_c) %>%
      as.data.frame() %>%
      mutate(cluster_id = k)
    
    L_c_k[[k]] <- res
    
  }
  
  mark <- do.call(rbind, L_c_k) %>% distinct()

  final_res <- full_join(final_res, mark)
  final_res$pass <- ifelse(is.na(final_res$pass), "no_pass", final_res$pass)
  final_res$cluster_id <- factor(final_res$cluster_id, levels = levels(sce$cluster_id))
  final_res$contrast <- factor(final_res$contrast, levels = unique(final_res$contrast))
  
  ## count DEGs
  L <- list()
  for (c in levels(final_res$contrast)) {
    tbl <- final_res[which(final_res$contrast == c & final_res$pass == "pass"),]
    n_de_cluster <- table(tbl$cluster_id)
    L[[c]] <- n_de_cluster
  }
  
  summary_all <- as.data.frame(do.call(rbind, L)) %>%
    mutate(method = method)
  
  # save DEGs and summary
  write_csv(final_res, paste0(path_save, "tables/", method, "_pb_", cluster_id, group_id, ".csv"))
  write_csv(summary_all, paste0(path_save, "tables/", method, "_n_degs_", cluster_id, group_id, ".csv"))
  
  # plots
  if (make_plots == TRUE) {
   
  colors$pass <- setNames(
    c("darkred", "steelblue", "grey"),
    c("up_pass", "down_pass", "no_pass")
  )
  
  ## Volcano
  for (k in levels(final_res$cluster_id)) {
    plot_list <- list()
    for (c in levels(final_res$contrast)) {
      data <- final_res %>%
        filter(cluster_id == k & contrast == c) %>%
        mutate(reg = factor(ifelse(logFC > 0, "up", "down"), levels = c("up", "down"))) %>%
        unite(col = "reg_pass", c("reg", "pass"), sep = "_") 
      data$reg_pass <- ifelse(!data$reg_pass %in% c("up_pass", "down_pass"), "no_pass", data$reg_pass)
      data$reg_pass <- factor(data$reg_pass, levels = c("up_pass", "down_pass", "no_pass"))
      
      p = ggplot(data, aes(logFC, -log10(p_adj.loc), color = reg_pass)) + 
        geom_point(size = 0.15, alpha = 0.8) + 
        geom_hline(yintercept = -log10(p.adj), linetype = 2, lwd = 0.2) +
        geom_vline(xintercept = 0, linetype = 2, lwd = 0.2) +
        xlim(c(-max(abs(data$logFC)), max(abs(data$logFC)))) +
        theme_bw() +
        xlab("logFC") +
        ylab("-log10(adjusted p-value)") + 
        scale_color_manual(values = colors$pass) +
        ggrepel::geom_text_repel(
          aes(label = ifelse(
            p_adj.loc < 0.05 & logFC > 2 & reg_pass != "no_pass" | 
              p_adj.loc < 0.05 & logFC < -2 & reg_pass != "no_pass", 
            gene, "")),
          color = "black", size = 2.8) +
        labs(title = str_replace_all(c, "_", " ")) +
        coord_cartesian(clip = "off") +
        theme(panel.grid = element_blank(),
              aspect.ratio = 1,
              axis.text = element_text(color = 'black'),
              plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
              axis.title = element_text(size = 12),
              legend.position = 'none') 
      plot_list[[c]] <- p
      
    }
    
    volcano_p <- plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "vh")
    ggsave(
      paste0(path_save, "plots/", method, "_volcano_", k, ".png"), 
      volcano_p +
        patchwork::plot_annotation(
          title = str_replace_all(k, "_", " "),
          theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 2.8))
        ),
      height = 4, width = 13)
  }
  
    ## Between-cluster concordance by group
    data <- final_res[which(final_res$pass == "pass"),]
    
    L_degs <- list(
      "up" = data$gene[which(data$logFC > 0)],
      "down" = data$gene[which(data$logFC < 0)]
    )
    
    for (i in names(L_degs)) {
      data_sub <- data[data$gene %in% L_degs[[i]],c("gene", "cluster_id", "contrast")] %>%
        pivot_wider(
          names_from = cluster_id, values_from = cluster_id
        )
      cluster <- setdiff(colnames(data_sub), c("gene", "contrast"))
      data_sub <- data_sub %>%
        mutate(across(all_of(cluster), ~ ifelse(is.na(.), FALSE, TRUE)))
      
      plot_list <- list()
      for (c in levels(final_res$contrast)) {
        p <- upset(
          data_sub[which(data_sub$contrast == c),], 
          cluster, 
          name = c,
          height_ratio=1,
          width_ratio=0.3,
          keep_empty_groups = TRUE,
          wrap = F,
          base_annotations = list(
            'Intersection size' = intersection_size(counts=FALSE)
          ),
          sort_sets = FALSE
        ) 
        plot_list[[c]] <- p
      }
      upset_p <- plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "vh")
      ggsave(paste0(path_save, "plots/", method, "_upset_", i, ".png"), upset_p, height = 4, width = 18)
      
    }
  }
    
  ## to do: heatmap
    
}
