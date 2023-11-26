# DA by milo

# we run milo on a batch-corrected latent space 
## (as the integration makes sure that the KNN graph captures cell states and not batch effects),
## then we add the batch covariate to make sure we don't call false positives 
## because of local incomplete batch correction.

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(miloR)
  library(Seurat)
  library(SingleCellExperiment)
})

# set paths
path_main <- '/Users/ieo4874/Desktop/working/SCMseq/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/all_samples/06_lineage_analysis/differential_abundance/plots/')

# settings
load(paste0(path_data, "settings.RData"))

# define function
run_milo <- function(
  seurat = seurat, 
  contrast = contrast, 
  reduction = reduction, 
  ndims = ndims, 
  alpha = alpha,
  lineage = lineage,
  only_frequent = TRUE,
  age = FALSE,
  seed = 01010111
) {
  
  # path to save results
  if (age == TRUE) {
      if(!dir.exists(paste0(path_results, "age_correction"))){
        dir.create(paste0(path_results, "age_correction"))
      }
    path_milo <- paste0(path_results, "age_correction/")
  } else {
    if(!dir.exists(paste0(path_results, "no_age_correction"))){
      dir.create(paste0(path_results, "no_age_correction"))
    }
    path_milo <- paste0(path_results, "no_age_correction/")
  }
  
  if(!dir.exists(paste0(path_milo, contrast))) {
    dir.create(path = paste0(path_milo, contrast))
  }
  path_save <- paste0(path_milo, contrast)
  
  # set seed
  set.seed(seed)
  
  # subset seurat to contrast of interest
  obj <- seurat[,which(
    seurat$cohort %in% 
      c(str_split_1(contrast, pattern = "_")[1], str_split_1(contrast, pattern = "_")[2])
  )]
  obj$cohort <- factor(obj$cohort, levels = unique(obj$cohort))
  obj$sample <- factor(obj$sample, levels = unique(obj$sample))
  
  # convert seurat to sce
  sce <- as.SingleCellExperiment(obj, assay = 'RNA')
  
  # create Milo object
  data_milo <- Milo(sce)
  
  # construct de novo KNN graph
  ## From Milo developers, it is recommended to choose k such that k≥Sx5, 
  ## where S is the number of experimental samples, or to have a distribution 
  ## peaking between 50 and 100, or mean average > 5 x n samples
  
  k = 5*length(unique( data_milo@colData$sample ))
  data_milo <- buildGraph(data_milo, k = k, d = ndims, reduced.dim = reduction)
  
  # Define representative neighbourhoods on the KNN graph
  data_milo <- makeNhoods(data_milo, 
                          prop = 0.1, 
                          k = k, d = ndims, 
                          refined = TRUE, 
                          refinement_scheme = "graph",
                          reduced_dims = reduction)
  
  # Check neighbourhood sizes 
  ## This heuristics serves to evaluate whether the value of k used for graph building is appropriate, 
  ## as neighbourhood sizes affect the power of DA testing.
  data_milo <- countCells(data_milo, 
                          meta.data = as.data.frame(colData(data_milo)), 
                          sample = "sample")
  
  test <- as.data.frame(data_milo@nhoodCounts) 
  test$size <- rowSums(test)
  
  p = ggplot(test, aes(size)) +
    geom_histogram(bins = 50) +
    theme_bw() +
    xlab("n cells in neighbourhoods") +
    ggtitle(paste0("k = ", k)) +
    geom_vline(xintercept = mean(test$size), linetype = 2, color = 'red') +
    geom_vline(xintercept = 5*length(unique( data_milo@colData$sample )), linetype = 2, color = 'orange') +
    scale_x_continuous(breaks = c(0,50,100,max(test$size))) +
    annotate('rect', xmin = 50, xmax = 100, 
             ymin = 0, ymax = Inf, 
             alpha = 0.2, fill = 'lightblue') +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(paste0(path_save, "/milo_NhoodSizeHist.png"), p, width = 5, height = 5)
  
  # prepare design matrix 
  data_design <- data.frame(colData(data_milo))[,c("sample", "cohort", "age", "sex")] 
  data_design <- distinct(data_design)
  rownames(data_design) <- data_design$sample
  data_design <- data_design[match(colnames(nhoodCounts(data_milo)), rownames(data_design)),]
  data_design$cohort <- factor(
    data_design$cohort, 
    levels = c(str_split_1(contrast, pattern = "_")[1],str_split_1(contrast, pattern = "_")[2])
  )
  
  # set contrast
  contr <- paste0(
    colnames(data_design)[2],str_split_1(contrast, pattern = "_")[1],
    "-",
    colnames(data_design)[2],str_split_1(contrast, pattern = "_")[2]
  )
  
  # test DA
  if (age == TRUE) {
    da_results <- testNhoods(data_milo,
                             design = ~ 0 + cohort + age, 
                             design.df = data_design, 
                             model.contrasts = contr,
                             fdr.weighting = 'graph-overlap')
  } else {
    da_results <- testNhoods(data_milo,
                             design = ~ 0 + cohort, 
                             design.df = data_design, 
                             model.contrasts = contr,
                             fdr.weighting = 'graph-overlap')
  }

  if (length(which(da_results$SpatialFDR < alpha)) > 1) {
    
    # Plot neighbourhood graph
    data_milo <- buildNhoodGraph(data_milo)
    
    p <- 
      plotNhoodGraphDA(
        data_milo, 
        da_results, 
        layout = "UMAP",
        alpha = alpha
      ) +
      labs(
        title = paste0(
          str_split_1(contrast, pattern = "_")[1]," Vs ",str_split_1(contrast, pattern = "_")[2]
        ), 
        subtitle = paste0("FDR < ", alpha)) +
      scale_fill_gradient2(name = "logFC",
                           low = "steelblue",
                           mid = "white",
                           high = "darkred",
                           midpoint = 0,
                           space = "Lab",
                           na.value = "lightgrey") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2),
            plot.subtitle = element_text(hjust = 0.5, vjust = 1))
    
    ggsave(paste0(path_save, "/da_milo_nhood.png"), p,
           width = 7, height = 7)
    
    # check p values
    p1 <- ggplot(da_results, aes(PValue)) + 
      geom_histogram(bins = floor(nrow(da_results)/10)) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    p2 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
      geom_point(alpha = 0.5, size = 0.5) +
      geom_hline(yintercept = -log10(alpha), color = "red", linetype = 2) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_save, "/da_milo_sig.png"), 
           p1 + p2 +
             patchwork::plot_layout(ncol = 2, widths = 5, heights = 4),
           width = 10, height = 4)
    
    # Find the most abundant cell type within cells in each neighbourhood
    da_results <- annotateNhoods(data_milo, da_results, coldata_col = lineage)
    
    p = ggplot(da_results, aes(x = da_results[[9]])) + 
      geom_histogram(bins = floor(nrow(da_results)/10)) +
      geom_vline(xintercept = 0.7, linetype = 2, color = 'red') +
      theme_bw() +
      xlab("lineage") +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_save, "/da_milo_lineage_fraction.png"), p,
           width = 7, height = 7)
    
    da_results$lineage <- ifelse(da_results[[9]] < 0.7, "Mixed", da_results[[8]])
    da_results$lineage <- factor(
      da_results$lineage,
      levels = c(levels(obj@meta.data[[lineage]])[which(levels(obj@meta.data[[lineage]]) %in% unique(da_results$lineage))],"Mixed")
    )
    
    if (lineage == "aggregated_lineage2") {
      freq <- "frequent_aggr_lineage"
    } else {
      freq <- "frequent_lineage"
    }
    
    da_results$keep <- ifelse(
      da_results$lineage %in% unique(obj@meta.data[[lineage]][which(obj@meta.data[[freq]] == TRUE)]),
      "TRUE", "FALSE"
    )
    
    if (only_frequent) {
      plot_data <- da_results %>%
        filter(lineage != "Mixed" & keep == TRUE)
    } else {
      plot_data <- da_results %>%
        filter(lineage != "Mixed")
    }
    
    p = plot_data %>%
      plotDAbeeswarm(group.by = "lineage", alpha = alpha) +
      ggtitle("Abundance of cell-type neighbourhood") +
      scale_colour_gradient2(low = "steelblue",
                             mid = "white",
                             high = "darkred",
                             midpoint = 0,
                             space = "Lab",
                             na.value = "lightgrey") +
      scale_x_discrete(limits = rev) +
      ylim(-10,10) +
      coord_cartesian(clip = "off") +
      coord_flip(clip = "off") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 25, size = 13),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        plot.margin = unit(c(3.5,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(linewidth = 0.5)
      ) +
      annotation_custom(grob = grid::linesGrob(
        arrow = arrow(type="open", ends="both", length = unit(3,"mm")), 
        gp = grid::gpar(col = "black", lwd = 1)
      ), 
      xmin = length(unique(plot_data$lineage))+1, 
      xmax = length(unique(plot_data$lineage))+1, 
      ymin = -10, ymax = 10) +
      annotation_custom(
        grob = grid::textGrob(
          label = str_split_1(contrast, pattern = "_")[1], 
          hjust = 0.5, vjust = -0.5, 
          rot = 0, gp = grid::gpar(col = "black")
        ), 
        xmin = length(unique(plot_data$lineage))+1.5, 
        xmax = length(unique(plot_data$lineage))+1.5, 
        ymin = 6, ymax = 10) +
      annotation_custom(
        grob = grid::textGrob(
          label = str_split_1(contrast, pattern = "_")[2],
          hjust = 0.5, vjust = -0.5, 
          rot = 0, gp = grid::gpar(col = "black")
        ), 
        xmin = length(unique(plot_data$lineage))+1.5, 
        xmax = length(unique(plot_data$lineage))+1.5, 
        ymin = -10, ymax = -6) 
    
    ggsave(paste0(path_save, "/da_results.png"), p, width = 5, height = 6)
    
  } else {
    message(paste0("No nhoods with significant FDR for contrast: ", contrast))
  }
}

