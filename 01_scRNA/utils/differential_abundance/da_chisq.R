suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(corrplot)
  library(ggpubr)
})

# chi square on most abundant cell types
forcluster <- table(as.character(final_seurat$lineage), 
                    as.character(final_seurat$cohort))

forcluster <- forcluster[rownames(forcluster) %in% ct_keep,]

forcluster <- forcluster[
  match(levels(final_seurat$lineage)[which(levels(final_seurat$lineage) %in% ct_keep)], 
        rownames(forcluster)),
  match(levels(final_seurat$cohort), colnames(forcluster))
]

x <- chisq.test(forcluster) 
# X-squared = 21265, df = 12, p-value < 2.2e-16

pdf(paste0(path_results, "plots/chisq_lineage_cohort.pdf"), width = 6, height = 7)
corrplot(
  x$residuals,
  title = paste0(x$method, ": p-value < 2.2e-16"),
  is.cor = FALSE,
  mar = c(0, 0, 3, 0),
  addgrid.col = NULL,
  addCoef.col = NULL,
  addCoefasPercent = FALSE,
  order = "original", 
  tl.cex = 1.3,
  tl.col = "black",
  tl.offset = 0.4,
  tl.srt = 90,
  cl.cex = 1,
  cl.ratio = 0.6,
  cl.align.text = "r",
  cl.offset = -1
)
dev.off()

# chi square on aggregated cell types
forcluster <- table(as.character(final_seurat$aggregated_lineage), 
                    as.character(final_seurat$cohort))

forcluster <- forcluster[
  match(levels(final_seurat$aggregated_lineage)[which(levels(final_seurat$aggregated_lineage)!= "other")], rownames(forcluster)),
  match(levels(final_seurat$cohort), colnames(forcluster))
]

x <- chisq.test(forcluster) 
#X-squared = 24430, df = 10, p-value < 2.2e-16

pdf(paste0(path_results, "plots/chisq_aggr_lineage_cohort.pdf"), width = 9, height = 10)
corrplot(
  x$residuals,
  is.cor = FALSE,
  mar = c(0, 0, 0, 0),
  addgrid.col = NULL,
  addCoef.col = NULL,
  addCoefasPercent = FALSE,
  order = "original", 
  tl.cex = 1.3,
  tl.col = "black",
  tl.offset = 0.4,
  tl.srt = 90,
  cl.cex = 1,
  cl.ratio = 0.6,
  cl.align.text = "r",
  cl.offset = 2
)
dev.off()