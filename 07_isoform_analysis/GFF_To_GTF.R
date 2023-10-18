## GFF to GTF

#####GFF to GTF for new samples
samplelist=c("AML2","AML3","hBM1","hBM2","hBM3") ##add other samples here
for (ii in 1:length(samplelist)){
  print(samplelist[ii])
  gff <- rtracklayer::import(paste0("/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/",samplelist[ii],"/FLAMES_out/isoform_annotated.filtered.gff3"))
  head(gff)
  gff <- as.data.frame(gff) %>%
    separate(transcript_id, into = c('transcript_id', 'stuff'), sep = 'transcript:') %>%
    dplyr::select(-stuff) %>%
    separate(exon_id, into = c('stuff', 'exon_id'), sep = 'exon:') %>%
    dplyr::select(-stuff) %>% 
    dplyr::rename('isoform_id' = 'transcript_id') %>%
    dplyr::select(c('seqnames', 'start', 'end', 'width', 'strand', 
                    'source', 'type', 'ID', 
                    'gene_id', 'isoform_id', 'exon_id', 'rank'))
  
  gr <- makeGRangesFromDataFrame(gff,
                                 seqnames.field = 'seqnames',
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand",
                                 keep.extra.columns=TRUE)
  
  ### Map regions to gene
  geneRanges <- function(db, column="gene_id") {
    g <- GenomicFeatures::genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
  
  gns <- geneRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, column="gene_id")
  
  splitByOverlap <- function(query, subject, column="gene_id", ...) {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }
  
  gene_region <- splitByOverlap(gns, gr, "gene_id") 
  gene_region <- paste(gene_region, collapse=", ")
  
  gff <- cbind(gff, gene_region)
  
  gff <- gff %>% separate(gene_region, into = c('gene_region', 'alternative'), sep = ',') 
  gff <- gff %>% rownames_to_column(var = 'index')
  
  ### Handle unmapped regions 
  unmapped <- gff %>% dplyr::filter(gene_region == '') #could not associate some ranges to gene regions
  gff <- gff %>% dplyr::filter(gene_region != '') 
  gff$gene_region <- factor(gff$gene_region, levels = unique(gff$gene_region))
  
  gff_new <- data.frame()
  for (i in unique(gff$gene_region)) {
    x = gff %>% dplyr::filter(gene_region == i)
    x = x %>% mutate(gene_id = rep(gene_id[1], nrow(x)))
    gff_new = bind_rows(x, gff_new)
  }
  
  unmapped2 <- gff_new %>% dplyr::filter(is.na(gene_id))
  unmapped <- rbind(unmapped, unmapped2)
  
  gff_new <- gff_new %>%
    dplyr::filter(!is.na(gene_id)) %>%
    separate(gene_id, into = c('gene_id', 'stuff'), sep = 15, extra = 'merge')
  
  gff_new$gene_name = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                            keys = gff_new$gene_id, #Column containing query gene ids
                                            column = "SYMBOL",
                                            keytype = "ENSEMBL",
                                            multiVals = 'first')
  
  gff_new <- gff_new %>% 
    unite(gene_id, c('gene_id', 'stuff'), sep = '', remove = T)
  
  unmapped3 <- gff_new[which(is.na(gff_new$gene_name)),]
  unmapped <- unmapped %>% mutate(gene_name = NA)
  unmapped <- rbind(unmapped, unmapped3)
  
  gff <- gff_new %>% 
    dplyr::filter(!is.na(gene_name))
  
  ### Collapse
  gff$gene_name <- factor(gff$gene_name, levels = unique(gff$gene_name))
  
  gff_2 <- data.frame()
  for (i in levels(gff$gene_name)) {
    x = gff %>% dplyr::filter(gene_name == i) %>%
      fill(isoform_id, .direction = "down")
    gff_2 = bind_rows(x, gff_2)
  }
  
  gff_ready <- gff_2
  
  #########################################
  
  ###Troubleshoot and add unmapped annotations
  # map again ranges to genes, using a different method (biomaRt::getBM())
  to_map <- unmapped[which(unmapped$type == 'gene'),]
  
  regions <- paste(to_map$seqnames, to_map$start, to_map$end, paste0(to_map$strand, 1), sep=":") %>%
    as.data.frame()
  regions <- regions %>% separate('.', into = c('discard', 'keep'), sep = 'chr')
  regions <- as.list(regions$keep)
  
  hsapiens.anno <- useEnsembl(biomart = "genes", 
                              dataset = "hsapiens_gene_ensembl")
  
  df1 <- data.frame()
  
  for (i in regions) {
    ranges = print(i)
    x = getBM(
      attributes = c("hgnc_symbol", "chromosome_name", 
                     "start_position", "end_position", "strand"),
      filters = c("chromosomal_region"),
      values=ranges,
      mart=hsapiens.anno,
      uniqueRows = T,
      useCache = T
    )
    x <- x %>% mutate(chromosome_name = as.character(chromosome_name))
    if (nrow(x) == 0) {
      
      x = data.frame(
        hgnc_symbol = NA,
        chromosome_name = NA,
        start_position = NA,
        end_position = NA,
        strand = NA
      )
      
    } else {
      x = x
    }
    x = cbind(ranges, x) 
    df1 = bind_rows(x, df1)
  }
  
  unmapped_remapped <- df1
  
  unmapped <- unmapped %>%
    mutate(ranges = paste(unmapped$seqnames, unmapped$start, unmapped$end, paste0(unmapped$strand, 1), sep=":") ) %>%
    separate(ranges, into = c('discard', 'ranges'), sep = 'chr') %>%
    dplyr::select(-'discard')
  
  new <- full_join(unmapped, unmapped_remapped, by = c('ranges', 'start' = 'start_position', 'end' = 'end_position'))
  new <- new %>% dplyr::filter(!is.na(new$index))
  new <- new %>% distinct(index, .keep_all = T)
  
  new <- new %>% mutate(gene_name = hgnc_symbol) %>%
    dplyr::select(-c('hgnc_symbol', 'chromosome_name', 'strand.y', 'ranges'))
  
  colnames(new)[6] <- 'strand'
  
  identical(colnames(new), colnames(gff_ready)) # needs to be TRUE
  
  # fill gene_id and isoform_id, down
  
  new <- new %>% 
    fill(gene_id, .direction = "down") 
  
  new$gene_id <- factor(new$gene_id, levels = unique(new$gene_id))
  
  new_2 <- data.frame()
  for (i in levels(new$gene_id)) {
    x = new %>% dplyr::filter(gene_id == i) %>%
      fill(isoform_id, .direction = "down")
    new_2 = rbind(x, new_2)
  }
  
  gff_ready <- rbind(gff_ready, new_2)
  
  gff_ready <- gff_ready[,c('seqnames', 'start', 'end', 'width', 'strand', 
                            'source', 'type', 'gene_id', 'isoform_id',  'gene_name')]
  
  colnames(gff_ready)[9] <- 'transcript_id'
  head(gff_ready)
  rtracklayer::export(gff_ready, paste0("/hpcnfs/scratch/PGP/SCMseq/isoform_analysis/GTF_FASTA/",samplelist[ii],"ErlyLateHSC_.gtf"), format = 'gtf')
}

