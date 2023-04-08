# Set working directory at containing respective results of previous processes("~stall_no_rmChimera.rds")

output.path=".."
## ======================================== ##

# -- If you want sum up sequence read, you select "sum".
type="sum"

# reference = ** # write the absolute path where "reference" is
### (note) we use "silva NR99 v138.1 database" as the reference
##############################################

# -- Loading packages
library(dada2)
library(seqinr)

## -- Loading data table
files <- list.files()
tables=lapply(files, readRDS)

## -- Combine data table
stall <- mergeSequenceTables(tables=tables, repeats=type)
seqtab.nochim <- removeBimeraDenovo(stall, method="consensus", multithread=TRUE, verbose=TRUE)

taxa <- assignTaxonomy(seqtab.nochim, reference, multithread=TRUE)
taxa[is.na(taxa)] <- "Unidentified"

if(all(rownames(taxa)==colnames(seqtab.nochim))){
  
  seq.mat <- cbind(colnames(seqtab.nochim),sprintf('X_%s', formatC(1:ncol(seqtab.nochim), width = nchar(ncol(seqtab.nochim)), flag = "0")))
  write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/merge_clus_seq.fasta", output.path) )

  colnames(seqtab.nochim) <- rownames(taxa) <- seq.mat[,2]
  ## ======================================== ##
  write.table(cbind(sample=rownames(seqtab.nochim), seqtab.nochim), sprintf("%s/merge_seqtab.txt", output.path), sep="\t", quote=F, row.names=F) # CHANGE ME to where you want sequence table saved
  saveRDS(seqtab.nochim,  sprintf("%s/merge_seqtab.rds", output.path))
  saveRDS(taxa,  sprintf("%s/merge_taxonomylist.rds", "../table"))
  
  
}else{
  stop("rownames(taxa) and colnames(seqtab.nochim) are not match")
}


##############################################
