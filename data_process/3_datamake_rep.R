library(vegan); library(stringr)
source("functions/functions.R")

#-- read data from previous dir
Aseqtab <- readRDS("table/ASV_seqtab_5000.rds")
Rseqtab <- readRDS("table/List_seqtab_5000.rds")
taxtable <- readRDS("table/taxonomylist_5000.rds")
exp <- read.csv("table/expdata.csv")

ti <- c(2,4,6,8,10,12)

#- merge_rep_data
merge_rep <- function(seqtab,exp=exp,ti=ti){
  st <- cbind(Sample_name = rownames(seqtab),as.data.frame(seqtab))
  data <- merge(exp,st,by="Sample_name"); data <- data[(data$Day %in% ti),]
  rownames(data) <- as.character(data$Sample_name)
  ds <- list()
  for(i in 1:length(ti)){
    d <- data[(data$Day == ti[i]),]
    ds[[i]] <- as.numeric(str_sub(d$Sample_name,-4,-1))%%384
  }
  key <- intersect(ds[[1]],ds[[2]])
  for(i in 3:length(ti)){
    key <- intersect(key,ds[[i]])
  }
  data <- data[((as.numeric(str_sub(data$Sample_name,-4,-1))%%384) %in% key),]
  ans <- list(data,key,ti); return(ans)
}
merge_d <- function(seqtab,exp=exp,ti=ti){
  st <- cbind(Sample_name = rownames(seqtab),as.data.frame(seqtab))
  data <- merge(exp,st,by="Sample_name"); data <- data[(data$Day %in% ti),]
  rownames(data) <- as.character(data$Sample_name)
  
  ans <- list(data,ti); return(ans)
}

###### return genus, family table#####
out1 <- list(); out2 <- list()
for(i in 1:length(Rseqtab)){
  hoge1 <- merge_rep(seqtab=Rseqtab[[i]],exp=exp,ti=ti)
  out1[[i]] <- hoge1[[1]]
  hoge2 <- merge_d(seqtab=Rseqtab[[i]],exp=exp,ti=ti)
  out2[[i]] <- hoge2[[1]]
}

a_out <- merge_rep(seqtab=Aseqtab,exp=exp,ti=ti)
out1[[(length(Rseqtab)+1)]] <- list(a_out[[2]],a_out[[3]])
a_out2 <- merge_d(seqtab=Aseqtab,exp=exp,ti=ti)
out2[[(length(Rseqtab)+1)]] <- a_out2[[2]]

saveRDS(out1,sprintf("table/%s_data_%s.rds","List","rep"))
saveRDS(out2,sprintf("table/%s_data_%s.rds","List","5000"))
