# Over-representation analysis 
# Chiara Caprioli
# Nov 2nd 2023

# Set up
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(enrichR)
dbs_selected <- c(
  "GO_Molecular_Function_2021",
  "GO_Cellular_Component_2021",
  "GO_Biological_Process_2021",
  "KEGG_2021_Human",
  "MSigDB_Hallmark_2020",
  "MSigDB_Oncogenic_Signatures",
  "Reactome_2022",
  "WikiPathway_2021_Human",
  "GTEx_Aging_Signatures_2021"
  )

# Run over-representation analysis
run_ora <- function(
  path_degs = path_degs,
  path_save = path_save,
  make_plots = TRUE
) {
  
  # create directory
  if (!dir.exists(paste0(path_save, "ORA/"))) {
    dir.create(paste0(path_save, "ORA/"))
  }
  
  # get DE
  res_de <- read.csv(path_degs)
  res_de <- res_de[which(res_de$pass == "pass"),]
  
  # run ORA
  ## query up
  ### filter upregulated DEGs
  data <-  res_de[which(res_de$logFC > 0),]
  
  ### query
  for (k in unique(data$cluster_id)) {
    L_c <- list()
    for (c in unique(data$contrast)) {
      data_c_k <- data[which(data$cluster_id == k & data$contrast == c),]
      
      if (length(unique(data_c_k$gene)) > 5) {
        enriched <- enrichr(data_c_k$gene, dbs_selected)
        enriched <- do.call(rbind, enriched) %>% 
          rownames_to_column(var = "Database") %>%
          mutate(cluster_id = k, contrast = c, regulation = "up")
        L_c[[c]] <-  enriched 
      } else {
        next
      }
      res <- do.call(rbind, L_c)
      write_csv(res, paste0(path_save, "ORA/enrichr_up_", k, ".csv"))
    }
  }
  
  ## query down
  ### filter downregulated DEGs
  data <-  res_de[which(res_de$logFC < 0),]
  
  ### query
  for (k in unique(data$cluster_id)) {
    L_c <- list()
    for (c in unique(data$contrast)) {
      data_c_k <- data[which(data$cluster_id == k & data$contrast == c),]
      
      if (length(unique(data_c_k$gene)) > 5) {
        enriched <- enrichr(data_c_k$gene, dbs_selected)
        enriched <- do.call(rbind, enriched) %>% 
          rownames_to_column(var = "Database") %>%
          mutate(cluster_id = k, contrast = c, regulation = "down")
        L_c[[c]] <-  enriched 
      } else {
        next
      }
      res <- do.call(rbind, L_c)
      write_csv(res, paste0(path_save, "ORA/enrichr_down_", k, ".csv"))
    }
  }
  
  # plot
  if (make_plots == TRUE) {
    
    for (k in unique(data$cluster_id)) {
      
      if (file.exists(paste0(path_save, "ORA/enrichr_up_", k, ".csv"))) {
        
        data_k_up <- read_csv(paste0(path_save, "ORA/enrichr_up_", k, ".csv"))
        data_k_down <- read_csv(paste0(path_save, "ORA/enrichr_down_", k, ".csv"))
        
        for (c in unique(data$contrast)) {
          
          if (
            any(unique(data_k_up$contrast) %in% c) & 
            any(unique(data_k_down$contrast) %in% c)
          ) {
            
            top_up <- data_k_up[which(data_k_up$cluster_id == k & data_k_up$contrast == c),] %>%
              filter(Adjusted.P.value < 0.1) %>%
              arrange(Adjusted.P.value) %>%
              dplyr::slice(1:25) %>%
              mutate(Adjusted.P.value = -log10(Adjusted.P.value))
            
            top_down <- data_k_down[which(data_k_down$cluster_id == k & data_k_down$contrast == c),] %>%
              filter(Adjusted.P.value < 0.1) %>%
              arrange(Adjusted.P.value) %>%
              dplyr::slice(1:25) %>%
              mutate(Adjusted.P.value = -log10(Adjusted.P.value)*-1) %>%
              arrange(desc(Adjusted.P.value)) 
            
            top <- rbind(top_up, top_down) %>%
              distinct(Term, .keep_all = T)
            top$Term <- factor(top$Term, levels = top$Term)
            top$regulation <- factor(top$regulation, levels = c("up", "down"))
            
            p = ggplot(top, aes(x = Adjusted.P.value, y = Term, fill = regulation)) +
              geom_col() +
              theme_bw() +
              scale_y_discrete(limits=rev) +
              scale_fill_manual(name = "Enrichment",
                                values = c("darkred", "steelblue"),
                                labels = c("upregulated", "downregulated")) +
              labs(title = k, subtitle = c) +
              labs(x = "-log10(FDR)") +
              theme(axis.text = element_text(size = 8, color = "black"),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 13),
                    panel.grid = element_blank(),
                    plot.title = element_text(size = 13, hjust = 0.5, vjust = 2, face = "bold"),
                    plot.subtitle = element_text(size = 11, hjust = 0.5, vjust = 2, face = "plain"))
            
            ggsave(paste0(path_save, "ORA/", k, "_", c, "_top50_path.png"), p, width = 10, height = 10)
          } else {
            next
          }
        }
      } else {
        next
      }
    }
  }
  
}




