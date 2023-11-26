```{r}
## reading the dataframe
adf<-read.csv("/hpcnfs/scratch/PGP/SCMseq/downstream/integration/full_dataset.csv",stringsAsFactors = F) 
samples_ne=c("AML2","AML3") 
for (i in samples_ne){
  ndf<-adf %>% dplyr::filter(sample==i) ndf$barcode<- str_split(ndf$barcode,pattern = "-",n = 2,simplify = T)[,1] write_delim(as.data.frame(ndf$barcode),paste0("/hpcnfs/scratch/PGP/SCMseq/AML_3wtSample/",i,"/short_bardcodes.txt"),col_names = F)
}
```
