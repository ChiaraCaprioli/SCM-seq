---
title: "Annotate target variants"
author: "Chiara Caprioli"
date: "Oct 12th 2023"
---

/*
To inform the selection of targets for mutation-specific enrichment, 
we need to further annotate candidate variants according to:
- expression in scRNAseq data
- clonal structure

We need to feed a sample sheet with the following structure:
|sample_id|path_annotation|path_tenx_matrix|path_pyclone|
|:---:|:---:|:---:|:---:|
|sample1|/path2/annotation/somatic/sample1/|/path2/sample1/filtered_feature_bc_matrix/|/path2/sample1.pyclone.output.tsv|
*/
  
# Load libraries
library(tidyverse)
library(maftools)
library(Seurat)

# Load sample sheet
sample_sheet <- read_tsv("~/Desktop/sample_sheet.tsv")

for (i in 1:length(rownames(sample_sheet))) {
  
  ## Load variants
  maf <- read.maf(
    paste0(
      as.character(sample_sheet[i,2]),
      list.files(
        as.character(sample_sheet[i,2]),
        pattern = paste0(".*maf", "$")
        )
      )
    )
  var_data <- rbind(
    as.data.frame(maf@data),
    as.data.frame(maf@maf.silent)
  )
  
  ## Label candidate targets based on custom filters
  targets <- read_tsv(
    paste0(
      as.character(sample_sheet[i,2]),
      list.files(
        as.character(sample_sheet[i,2]),
        pattern = "*final" # filtered .tsv from variantalker
        )
      )
    ) %>%
    dplyr::select(c("Hugo_Symbol","Chromosome","Start_Position","End_Position")) %>%
    mutate("candidate_target" = "PASS")
  
  var_data <- full_join(var_data, targets)

  ## Annotate expression
  if (nrow(var_data) > 0) {
    
    # Initialize the Seurat object with raw counts
    sc_data <- Read10X(data.dir = as.character(sample_sheet[i,3]))
    seurat <- CreateSeuratObject(counts = sc_data, min.cells = 0, min.features = 0) # QC applied later
   
    # Subset Seurat to variants
    seurat <- seurat[unique(var_data$Hugo_Symbol),]
    
    # Expression stats
    exp_mat <- GetAssayData(seurat)
    exp_tab <- t(apply(exp_mat,1,summary)) 
    exp_tab <- as.data.frame(exp_tab) %>%
      mutate(perc_non_zero_cells = round(rowSums(exp_mat > 0)/ncol(exp_mat), digits = 2)) %>%
      rownames_to_column(var = "Hugo_Symbol") %>%
      arrange(desc(perc_non_zero_cells))
    exp_tab <- left_join(exp_tab, var_data[,c("Hugo_Symbol", "candidate_target")])
    
    exp_tab <- exp_tab[!duplicated(exp_tab$Hugo_Symbol) | 
                     !is.na((duplicated(exp_tab$Hugo_Symbol) & exp_tab$candidate_target == "PASS")),]
    exp_tab$candidate_target <- ifelse(is.na(exp_tab$candidate_target), "NO PASS", "PASS")
    exp_tab$candidate_target <- factor(exp_tab$candidate_target, levels = c("PASS", "NO PASS"))
    exp_tab$Hugo_Symbol <- factor(exp_tab$Hugo_Symbol, levels = unique(exp_tab$Hugo_Symbol))
    
    write.csv(
      exp_tab,  
      paste0(
        as.character(sample_sheet[i,2]),
        as.character(sample_sheet[i,1]),
        "_expression_summary.csv"
        ),
      row.names = T
      )
    
    p = ggplot(exp_tab, aes(Hugo_Symbol, perc_non_zero_cells, fill = candidate_target)) +
      geom_col(width = 0.75) +
      theme_bw() +
      scale_fill_manual(name = "candidate target", values = c("red", "grey")) +
      ggtitle(as.character(sample_sheet[i,1])) +
      ylab("% cells expressing gene") +
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            aspect.ratio = 1,
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      paste0(
        as.character(sample_sheet[i,2]),
        as.character(sample_sheet[i,1]),
        "_expression_target.png"
        ),
      p)
    
    # Label variants based on expression
    L <- list()
    for (v in var_data$Hugo_Symbol) {
      tab <- data.frame(
        Hugo_Symbol = v,
        Expression = ifelse(exp_tab$perc_non_zero_cells[v] > 0, "PASS", "NO_PASS")
      )
      L[[v]] <- tab
    }
    exp_res <- do.call(rbind, L)
    var_data <- full_join(var_data, exp_res)
    
    ## Annotate clonal structure  
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
    
    # Save final annotation 
    write_tsv(
      var_data,
      paste0(
        as.character(sample_sheet[i,2]),
        list.files(
          as.character(sample_sheet[i,2]),
          pattern = "*final"
        )
      )
    )
    
  } else {
    next
  }
}
