---
title: "04_clustering"
author: "Chiara Caprioli"
date: "July 6th 2023"
---

**Aims:**
compute Leiven clusters spanning different k and resolution

***Computationally intensive***

## Set up
library(tidyverse)
library(Seurat)
library(reticulate)
library(log4r)
library(SeuratObject)

## set paths and settings
path_main <- '/hpcnfs/scratch/temporary/ccaprioli/'
load(paste0(path_main, "settings.RData"))

## set log files
logfile = "log.txt"
file_appender = file_appender(logfile, 
                                 append = TRUE, 
                                 layout = default_log_layout())

my_logger <- logger(threshold = "INFO", appenders = file_appender)

## Load Seurat object 
final_seurat <- readRDS(paste0(path_main, 'final_seurat.rds'))

## Compute clusters by Leiden algorithm
use_condaenv("r-reticulate")

set.seed(123)

k_list <- c(5,10,15,30,50,100)
res_list <- seq(0.5, 2, by = 0.5)

DefaultAssay(final_seurat) <- "RNA"

res_L <- list() 
for (k in k_list) {
  
  k_start_time <- Sys.time()
  
  for (r in res_list) {
    
    info(my_logger, paste0("Evaluating k ", k, " and res ", r))
    
    final_seurat <- FindNeighbors(
      final_seurat,
      reduction = 'scanorama', 
      k.param = k,
      dims = 1:30,
      do.plot = T)
    
    final_seurat <- final_seurat %>%
      FindClusters(resolution = r, 
                   algorithm = 4, # 4 = leiden
                   method = "igraph", 
                   random.seed = 123) 
    
    res_L[[paste0("RNA_snn_k_", k, "_res_", r)]] <- final_seurat@meta.data[,paste0("RNA_snn_res.", r)]
    
  }
  
  k_end_time <- Sys.time()
  
  info(my_logger, 
       paste0("k ", k, " successfully done in ", 
              floor(as.numeric(k_end_time - k_start_time, units = "mins")), " mins")
       )
  
  write.table(readLines(logfile), 
              file = paste0(path_main, logfile), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE)
  
}

cluster_all <- do.call(rbind, res_L) %>% t()
final_seurat@meta.data <- cbind(final_seurat@meta.data, cluster_all)

saveRDS(final_seurat, paste0(path_main, 'seurat_k.rds'))
