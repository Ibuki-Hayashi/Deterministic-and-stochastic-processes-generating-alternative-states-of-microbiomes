library(vegan); library(stringr)
source("functions/functions.R")

#-- read data from previous dir
exp <- read.csv("table/expdata.csv")
seqtab <- readRDS("table/merge_seqtab.rds")
taxtable <- readRDS("table/merge_taxonomylist.rds")

#-- extract only inoculum rows
seq <- cbind(Sample_name=rownames(seqtab),as.data.frame(seqtab))
k <- seq[(str_sub(seq$Sample_name,1,4) != "HS2_"),]
data1 <- merge(exp,k,by="Sample_name")
data <- data1[(rowSums(data1[,(6:ncol(data1))])>9999),6:ncol(data1)]
rownames(data) <- data1$Sample_name[(rowSums(data1[,(6:ncol(data1))])>9999)]

###### seqtab ###########
#- rarefuction by n reads(without stdDNA)
raref_i <- function(seqtab,taxatable,n,seed=1111){
  std_out <- seqtab[,(colnames(seqtab) %in% rownames(taxtable[(str_sub(taxtable[,2],1,3) != "STD"),]))]
  use <- std_out[(rowSums(std_out)>(n-1)),]
  set.seed(seed)
  ans1 <- rrarefy(use,n); rownames(ans1) <- rownames(use)
  ans2 <- sprintf("table/Ino_ASV_seqtab_%s.rds",as.character(n))
  ans3 <- taxatable[(rownames(taxatable) %in% colnames(ans1)),]
  ans4 <- sprintf("table/Ino_taxonomylist_%s.rds",as.character(n))
  return(list(ans1,ans2,ans3,ans4))
}

number <- 10000 ## important!!!!!!!!!!!!!!!!!!!!!! the number of rarefaction

seq_ans <- raref_i(seqtab=data,taxatable=taxtable,n=number)

###### return ans#####
ans <- taxaratio(seq_ans[[1]],seq_ans[[3]])
ans2 <- list()
for(i in 1:length(ans)){
  ans2[[i]] <- ans[[i]]
}
ans2[[(length(ans)+1)]] <- seq_ans[[1]]

saveRDS(ans2,sprintf("table/%s_seqtab_%s.rds","ListIno",as.character(number)))
saveRDS(seq_ans[[1]],seq_ans[[2]]); saveRDS(seq_ans[[3]],seq_ans[[4]])
