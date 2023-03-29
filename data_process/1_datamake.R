library(vegan); library(stringr)
source("functions/functions.R")

#-- read data from previous dir
seqtab <- readRDS("table/merge_seqtab.rds")
taxtable <- readRDS("table/merge_taxonomylist.rds")

###### seqtab ###########
#- rarefuction by n reads(without stdDNA)
raref <- function(seqtab,taxatable,n,seed=1111){
  std_out <- seqtab[,(colnames(seqtab) %in% rownames(taxtable[(str_sub(taxtable[,2],1,3) != "STD"),]))]
  use <- std_out[(rowSums(std_out)>(n-1)),]
  set.seed(seed)
  ans1 <- rrarefy(use,n); rownames(ans1) <- rownames(use)
  ans2 <- sprintf("table/ASV_seqtab_%s.rds",as.character(n))
  ans3 <- taxatable[(rownames(taxatable) %in% colnames(ans1)),]
  ans4 <- sprintf("table/taxonomylist_%s.rds",as.character(n))
  return(list(ans1,ans2,ans3,ans4))
}

number <- 5000 ## important!!!!!!!!!!!!!!!!!!!!!! the number of rarefaction

seq_ans <- raref(seqtab=seqtab,taxatable=taxtable,n=number)

###### return genus, family table#####
ans <- taxaratio(seq_ans[[1]],seq_ans[[3]])
ans2 <- list()
for(i in 1:length(ans)){
  ans2[[i]] <- ans[[i]]
}
ans2[[(length(ans)+1)]] <- seq_ans[[1]]

saveRDS(ans2,sprintf("table/%s_seqtab_%s.rds","List",as.character(number)))
saveRDS(seq_ans[[1]],seq_ans[[2]]); saveRDS(seq_ans[[3]],seq_ans[[4]])
