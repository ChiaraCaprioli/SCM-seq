---
title: "Annotate and select variants variants for target enrichment"
author: "Chiara Caprioli"
date: "Nov 24th 2023"
---
  
#To inform the selection of targets for mutation-specific enrichment, 
#we need to further annotate candidate variants according to:
#- expression in scRNAseq data
#- clonal structure
#
#We need to feed a sample sheet with the following structure:
#|sample_id|path_annotation|path_tenx_matrix|path_pyclone|
#|:---:|:---:|:---:|:---:|
#|sample1|/path2/annotation/somatic/sample1/|/path2/sample1/filtered_feature_bc_matrix/|/path2/sample1.pyclone.output.tsv|
  
# Load libraries
library(tidyverse)
library(maftools)
library(Seurat)
library(ggpubr)
library(gridExtra)

# Load sample sheet
path_data <- "/Users/ieo4874/Desktop/working/SCMseq/data/wes/"
sample_sheet <- read_tsv(paste0(path_data, "sample_sheet.tsv"), show_col_types = FALSE)
   
# Annotation
L_selected <- list()
L_plot <- list()

for (i in 1:length(rownames(sample_sheet))) {
  
  message(
    paste0("Processing...", as.character(sample_sheet[i,1]))
  )
  
  ## Load ALL variants
  maf <- read.maf(
    paste0(
      as.character(sample_sheet[i,2]),
      list.files(
        as.character(sample_sheet[i,2]),
        pattern = paste0(".*maf", "$")
      )
    ),
    verbose = F
  )
  var_data <- rbind(
    as.data.frame(maf@data),
    as.data.frame(maf@maf.silent)
  )
  
  ## Label candidate targets (previously selected by custom filters)
  targets <- read_tsv(
    paste0(
      as.character(sample_sheet[i,2]),
      list.files(
        as.character(sample_sheet[i,2]),
        pattern = "*final" # filtered .tsv from variantalker
      )
    ),
    show_col_types = FALSE
  ) 
  var_data <- var_data[,colnames(targets)]
  
  targets <- targets %>%
    dplyr::select(c("Hugo_Symbol","Chromosome","Start_Position","End_Position")) %>%
    mutate("candidate_target" = "PASS")
  
  var_data <- full_join(var_data, targets) 
  var_data$candidate_target <- ifelse(!is.na(var_data$candidate_target), "PASS", "NO_PASS")
  
  ## Annotate expression
  if (nrow(var_data) > 0) {
    
    # Initialize the Seurat object with raw counts
    sc_data <- Read10X(data.dir = as.character(sample_sheet[i,3]))
    seurat <- CreateSeuratObject(counts = sc_data, min.cells = 0, min.features = 0) 
    
    # Subset Seurat to variants
    seurat <- seurat[unique(var_data$Hugo_Symbol),]
    
    # Expression stats
    exp_mat <- as.matrix(GetAssayData(seurat))
    
    #### add non expressed genes
    if (length(unique(var_data$Hugo_Symbol)) != nrow(exp_mat)) {
      missing_matrix <- matrix(
        data = 0,
        nrow = length(which(!unique(var_data$Hugo_Symbol) %in% rownames(exp_mat))),
        ncol = ncol(exp_mat)
      )
      colnames(missing_matrix) <- colnames(exp_mat)
      rownames(missing_matrix) <- unique(var_data$Hugo_Symbol)[
        which(!unique(var_data$Hugo_Symbol) %in% rownames(exp_mat))
      ] 
      exp_mat <- rbind(exp_mat, missing_matrix)
    }
    
    #### summary stats
    exp_tab <- t(apply(exp_mat,1,summary))
    
    exp_tab <- as.data.frame(exp_tab) %>%
      mutate(
        fraction_non_zero_cells = rowSums(exp_mat > 0)/ncol(exp_mat),
        mean_expression_non_zero_cells = rowSums(exp_mat) / rowSums(exp_mat > 0)
      ) %>%
      rownames_to_column(var = "Hugo_Symbol") 
    exp_tab$mean_expression_non_zero_cells[is.nan(exp_tab$mean_expression_non_zero_cells)] <- 0
    
    #### normalize mean expression to enable cross-sample comparison
    exp_tab$norm_mean <- scales::rescale(exp_tab$Mean, 0:1)
    exp_tab$norm_mean_exp_non_zero <- round(scales::rescale(exp_tab$mean_expression_non_zero_cells, 0:1), digits = 3)
    
    #### add label for candidate targets
    add_target <- var_data[,c("Hugo_Symbol", "candidate_target")] %>% 
      filter(candidate_target == "PASS") %>%
      distinct()
    exp_tab <- full_join(exp_tab, add_target) 
    exp_tab <- full_join(exp_tab, var_data[,c("Hugo_Symbol", "tumor_f")]) 
    exp_tab$candidate_target <- ifelse(!is.na(exp_tab$candidate_target), "PASS", "NO PASS")
    exp_tab$candidate_target <- factor(exp_tab$candidate_target, levels = c("PASS", "NO PASS"))
    exp_tab <- exp_tab %>% arrange(desc(tumor_f))
    exp_tab$Hugo_Symbol <- factor(exp_tab$Hugo_Symbol, levels = unique(exp_tab$Hugo_Symbol))
    
    #### save summary
    write.csv(
      exp_tab,  
      paste0(
        as.character(sample_sheet[i,2]),
        as.character(sample_sheet[i,1]),
        "_expression_summary.csv"
      ),
      row.names = T
    )
    
    #### plot relationship between VAF, fraction of expressing cells and candidate targets
    # barplot (not very clear especially when including a lot of mutations)
    #p1 = ggplot(exp_tab, aes(x = Hugo_Symbol)) +
    #  geom_col(aes(y = fraction_non_zero_cells, fill = candidate_target), width = 0.5) +
    #  geom_point(aes(y = tumor_f)) +
    #  geom_line(aes(y = tumor_f), group = "tumor_f", color = "black") +
    #  scale_y_continuous(
    #    name = "Fraction of cells expressing gene",
    #    sec.axis = dup_axis(name="VAF"),
    #    limits = c(0,1)
    #  ) +
    #  theme_bw() +
    #  scale_fill_manual(name = "candidate target", values = c("red", "grey")) +
    #  ggtitle(as.character(sample_sheet[i,1])) +
    #  theme(panel.grid = element_blank(),
    #        axis.title.x = element_blank(),
    #        aspect.ratio = 1,
    #        axis.text.x = element_text(angle = 45, hjust = 1))
    #
    #L_plot1[[i]] <- p1
    
    # spearman correlation
    p <- ggscatter(
      data = exp_tab, 
      x = "tumor_f",
      y = "fraction_non_zero_cells",
      color = "candidate_target",
      palette = c("red", "grey"),
      xlab = "VAF",
      ylab = "Fraction non zero cells",
      cor.coef = TRUE,
      cor.method = "spearman",
      cor.coef.coord = c(0.05, 1.08),
      cor.coef.size = 5,
      add = "reg.line",
      add.params = list(color = "blue", linetype = 2)
    ) +
      theme_bw() +
      scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
      scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
      ggtitle(as.character(sample_sheet[i,1])) +
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)
      )
    L_plot[[i]] <- p
    
    # Label variants based on expression
    L <- list()
    for (v in unique(var_data$Hugo_Symbol)) {
      tab <- data.frame(
        Hugo_Symbol = v,
        Expression_fraction_non_zero_cells = exp_tab[which(exp_tab$Hugo_Symbol == v), "fraction_non_zero_cells"],
        Expression_mean_non_zero_cells_norm = exp_tab[which(exp_tab$Hugo_Symbol == v), "norm_mean_exp_non_zero"],
        Expression_pass = ifelse(
          exp_tab[which(exp_tab$Hugo_Symbol == v), "fraction_non_zero_cells"] > 0, 
          "PASS", "NO_PASS"
        )
      )
      L[[v]] <- tab
    }
    exp_res <- do.call(rbind, L) %>% distinct()
    var_data <- full_join(var_data, exp_res) %>% 
      arrange(desc(tumor_f))
    
    # Combine VAF (cellular prevalence) and expression
    ## here we compute cellular prevalence assuming no CNAs exist; 
    ## when available, we can also use pyclone output (challenging for indels)
    var_data <- var_data %>%
      mutate(
        cell_prevalence = ifelse((tumor_f*2) < 1, (tumor_f*2), 1), # max cellular prevalence = 1
        combined_vaf_expr = cell_prevalence*Expression_fraction_non_zero_cells
      ) %>%
      dplyr::select(-cell_prevalence)

    # Label final targets
    ## Old criteria:
    # - variant type
    # - VAF ≥ 0.1
    # - fraction expressing cells > 0
    var_data$selected_old <- ifelse(
      var_data$candidate_target == "PASS" & 
        var_data$tumor_f >= 0.1 &   
        var_data$Expression_pass == "PASS",
      "PASS", "NO_PASS"
    )
    
    ## New criteria:
    # - variant type
    # - combined VAF and expression, arbitrary threshold
    var_data$selected_new <- ifelse(
      var_data$candidate_target == "PASS" & 
        var_data$combined_vaf_expr >= 0.001, # i.e., 1/1000 cells express the target, or all cells express the target with VAF 0.1 
      "PASS", "NO_PASS"
    )
    
    # Counts of selected targets by old and new thresholds
    x <- cbind(table(var_data$selected_old), table(var_data$selected_new), as.character(sample_sheet[i,1]))
    colnames(x) <- c("selected_old", "selected_new", "sample")
    x <- as.data.frame(x) %>% rownames_to_column(var = "pass")
    L_selected[[i]] <- x
    
    ### Annotate clonal structure (if available)
    if (!is.na(sample_sheet[i,4])) {
      # Load and format pyclone output
      pyclone_out <- read_tsv(
        as.character(sample_sheet[i,4])
      ) %>%
        dplyr::select(c(mutation_id, cluster_id, cellular_prevalence)) %>%
        mutate(mutation_id = gsub("[^[:alnum:]]", "_", pyclone_out$mutation_id))
      
      pyclone_out$Chromosome <- str_split_i(pyclone_out$mutation_id, pattern = "_", 1)
      pyclone_out$Start_Position <- as.numeric(str_split_i(pyclone_out$mutation_id, pattern = "_", 2))
      pyclone_out$End_Position <- as.numeric(str_split_i(pyclone_out$mutation_id, pattern = "_", 3))
      pyclone_out$Hugo_Symbol <- str_split_i(pyclone_out$mutation_id, pattern = "_", -1)
      
      # Label variants based on clonality
      var_data <- left_join(var_data, pyclone_out[,2:ncol(pyclone_out)])
    }
    
    # Save final annotation 
    write_tsv(var_data, paste0(path_data, "annotation_filtered/", as.character(sample_sheet[i,1]), ".tsv" ))
    
  } else {
    next
  }
}

# Save counts of selected targets
count_targets <- do.call(rbind, L_selected)
write_tsv(count_targets, paste0(path_data, "annotation_filtered/count_targets.tsv"))

# Save plots
pdf(paste0(path_data, "annotation_filtered/vaf_expression_spearman.pdf"), width = 19, height = 15)
grid.arrange(grobs = L_plot, ncol = 3)
dev.off()
