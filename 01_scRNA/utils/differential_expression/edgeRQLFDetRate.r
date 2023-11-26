#do_DE(): given some meta, some interesting var (with a ref) level
#and optionally some batch variable, perform DE with the cdr edgeR framework, and save results
do_DE <- function(M_test, m_test, batch_cov, var, ref) {

  #After gene filtering, cells with very low counts can have library NULL sizes (e.g.)
  #Filter our those cells 
  cells_to_filter <- colSums(M_test) >= 1
  M_test <- M_test[, cells_to_filter]
  m_test <- m_test[cells_to_filter, ]

  #Perpare DGEList object with size factors
  y <- DGEList(counts = M_test, samples = m_test) 
  y <- calcNormFactors(y) 
  
  #Calculate gene detection rate
  cdr <- scale(colMeans(M_test > 0))
  
  #If necessary prepare the batch variable
  if (batch_cov == TRUE) {
      y$samples$sequencing_run <- as.factor(y$samples$sequencing_run)
      batch <- y$samples$sequencing_run
  }
  
  ##Produce my interesting condition column, and set the reference
  y$samples[[var]] <- relevel(y$samples[[var]], ref = ref)
  BEH <- y$samples[[var]]
  
  ##Prepare the design matrix
  if (batch_cov == TRUE) {
      design <- model.matrix(~ cdr + batch + BEH) } else {
      design <- model.matrix(~ cdr + BEH)
  }

  #Estimate dispersions
  y <- estimateDisp(y, design = design)
  
  #Fit a per-gene NB model with the obtained dispersions 
  fit <- glmQLFit(y, design = design)
  
  #Use the quasi-likelihood F-test to compute pvalues for each contrast
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  
  #Filter up 
  de <- tt$table
  up <- de %>% filter(logFC > 0, FDR <= 0.1) %>% na.omit()
  
  if (dim(up)[1] > 300) {
    up <- up[1:300, ] %>% mutate(gene = row.names(.))
  } else {
    up <- up %>% mutate(gene = row.names(.))
  }

  #Reformat: pct.1, pct.2.
  if (dim(up)[1] >= 1) {

    #Create DEGs and cells masks
    DEGs_mask <- row.names(M_test) %in% up$gene
    int_cells_mask <- BEH != ref

    #Create pct.1 pct.2
    up$pct.1 <- rowSums(M_test[DEGs_mask, int_cells_mask] > 0) / sum(int_cells_mask) 
    up$pct.2 <- rowSums(M_test[DEGs_mask, !int_cells_mask] > 0) / sum(!int_cells_mask) 

    #Return de and up
    l <- list(de, up)
    names(l) <- c('all_genes', 'top_300_up')
    return(l) 

  } else {
    print('No markers found!')
  }
}

