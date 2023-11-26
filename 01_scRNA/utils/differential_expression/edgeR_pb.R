
library(scater)
library(scran)
library(batchelor)
library(edgeR)
library(tidyverse)
library(patchwork)
library(DT)


#Larger counts are more amenable to standard DE analysis pipelines designed 
#for bulk RNA-seq data. Normalization is more straightforward and certain statistical 
#approximations are more accurate e.g., the saddlepoint approximation for 
#quasi-likelihood methods or normality for linear models. 
#Collapsing cells into samples reflects the fact that our biological 
#replication occurs at the sample level (Lun and Marioni 2017). 
#Each sample is represented no more than once for each condition, 
#avoiding problems from unmodelled correlations between samples. 
#Supplying the per-cell counts directly to a DE analysis pipeline would imply 
#that each cell is an independent biological replicate, which is not true 
#from an experimental perspective. 
#(A mixed effects model can handle this variance structure but involves 
#  extra statistical and computational complexity for little benefit, see Crowell et al. (2019).) 
#Variance between cells within each sample is masked, provided it 
#does not affect variance across (replicate) samples. This avoids penalizing 
#DEGs that are not uniformly up- or down-regulated for all cells in all 
#samples of one condition. Masking is generally desirable as DEGs - unlike marker genes - 
#  do not need to have low within-sample variance to be interesting, 
#e.g., if the treatment effect is consistent across replicate populations 
#but heterogeneous on a per-cell basis. 
#(Of course, high per-cell variability will still result in weaker DE 
#if it affects the variability across populations, while homogeneous 
#per-cell responses will result in stronger DE due to a larger population-level 
#log-fold change. These effects are also largely desirable.)

sce <- as.SingleCellExperiment(final_seurat[,which(final_seurat$cohort != "SRSF2wt")])
# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.

sce$batch <- sce$sample
sce$SampleName <- sce$sample
sce$SampleGroup <- sce$cohort

columnsToUse <- c("batch", "SampleName", "SampleGroup", "aggregated_lineage", "age")
colData(sce) <- colData(sce) %>% data.frame() %>% dplyr::select(all_of(columnsToUse)) %>% DataFrame
summed <- aggregateAcrossCells(sce, 
                               id = DataFrame(
                               label = sce$aggregated_lineage,
                               sample = sce$SampleName
                               )
)
colData(summed) 

labelToGet <- "HSC_progenitors"
current <- summed[,summed$label==labelToGet]
colData(current)

countsToUse <- counts(current)
colnames(countsToUse) <- colData(current)$SampleName
y <- DGEList(countsToUse, samples=colData(current))
y
discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)

#Genes are discarded if they are not expressed above a log-CPM threshold in a minimum number of samples 
#(determined from the size of the smallest treatment group in the experimental design).
keep <- filterByExpr(
  y, 
  group=current$SampleGroup,
  min.count = 5, min.total.count = 10, min.prop = 0.5)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

par(mfrow=c(3,3))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}

plotMDS(cpm(y, log=TRUE), 
        col = as.numeric(y$samples$SampleGroup)
)

y$samples$SampleGroup <- factor(y$samples$SampleGroup)
y$samples$SampleGroup <- relevel(y$samples$SampleGroup, ref = "hBM")

design <- model.matrix(~SampleGroup+age, y$samples)

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

res <- glmQLFTest(fit, coef=ncol(design))
res_fdr <- topTags(res, n = nrow(res$table)) 
View(res_fdr$table %>% filter(FDR < 0.1))


