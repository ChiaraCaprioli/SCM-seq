# Gene set enrichment analysis 
# Chiara Caprioli
# Nov 2nd 2023

# Set up
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

HALLMARK <- msigdbr(species = 'Homo sapiens', category = 'H') %>%
  dplyr::select(gs_name, gene_symbol)

GO <- msigdbr(species = 'Homo sapiens', category = 'C5') %>%
  dplyr::select(gs_name, gene_symbol)

IMMUNESIGDB <- msigdbr(species = 'Homo sapiens', category = 'C7') %>%
  dplyr::select(gs_name, gene_symbol)

ONCOGENIC <- msigdbr(species = 'Homo sapiens', category = 'C6') %>%
  dplyr::select(gs_name, gene_symbol)

gseaSets <- list(HALLMARK, GO, IMMUNESIGDB, ONCOGENIC)
names(gseaSets) <- c("HALLMARK", "GO", "IMMUNESIGDB", "ONCOGENIC")

# Run over-representation analysis
run_gsea <- function(
  path_degs = path_degs,
  path_save = path_save
) {
  
  # create directory
  if (!dir.exists(paste0(path_save, "GSEA/"))) {
    dir.create(paste0(path_save, "GSEA/"))
  }
  
  # get DE
  res_de <- read.csv(path_degs)
  data <- res_de[which(res_de$pass == "pass"),]
  
  # run GSEA
  
  for (g in names(gseaSets)) {
    
    for (k in unique(data$cluster_id)) {

      for (c in unique(data$contrast)) {

        data_c_k <- data %>%
          filter(
            cluster_id == k &
              contrast == c &
              p_adj.loc < 0.1
            )
        
        if (length(which(data_c_k$gene %in% gseaSets[[g]]$gene_symbol)) > 2) {
          
          # get ranked gene list
          original_gene_list <- data_c_k$logFC
          names(original_gene_list) <- data_c_k$gene
          gene_list <- na.omit(original_gene_list)
          gene_list = sort(gene_list, decreasing = TRUE)
          
          # run GSEA
          gse <- GSEA(
            gene_list, TERM2GENE = gseaSets[[g]],
            minGSSize = 5, maxGSSize = 800,
            eps = 0, pvalueCutoff = 0.1,
            pAdjustMethod = "BH",
            by = "fgsea", scoreType = "pos"
          )
          
          # save
          if (nrow(gse@result != 0)) {
            write_csv(gse@result, paste0(path_save, "GSEA/GSEA_", g, "_", k,'_', c, ".csv"))
          }
        } else {
          next
        }
      }
    }
  }

  
}




