#!/bin/bash

########################################################################
## -- Make PICRUST file by R

cat <<RRR > data_process/script_pic.R

library(seqinr); library(stringr); library(vegan)
#- 
writefastafile=function(seqlist,namelist=seqlist,outfn="out"){
  isbadarg=c(FALSE, "")
  if(!is.vector(namelist))isbadarg=c(TRUE, "argument ‘namelist’ is wrong type.")
  else if(!is.vector(seqlist))isbadarg=c(TRUE, "extra argument ‘seqlist’ is wrong type.")
  else if(!is.character(outfn))isbadarg=c(TRUE, "extra argument ‘outfn’ is wrong type.")
  else if(length(namelist)!=length(seqlist))isbadarg=c(TRUE, "Difference in vector length between seqlist and namelist")
  if(isbadarg[1]){
    warning(isbadarg[2])
    return()
  }
  library(seqinr)
  for(i in 1:length(namelist)){
    if(i == 1) mode = 'w'
    else mode = 'a'
    write.fasta(sequences = as.vector(seqlist[i]),names=as.vector(namelist[i]),file.out = paste(outfn, ".fa",sep = ''), open=mode)
  }
}
#- 

seqtab <- readRDS("table/merge_seqtab.rds")
fasta <- read.fasta("table/merge_clus_seq.fasta",as.string=T)
taxa <- readRDS("table/merge_taxonomylist.rds")

taxa_ws <- taxa[(str_sub(taxa[,2],1,3) != "STD"),]

hoge <- rownames(taxa) %in% rownames(taxa_ws)

atgc <- c(); name <- rownames(taxa)
for(i in 1:length(fasta)){
  atgc[i] <- fasta[[i]][[1]]
}
ATGC <- toupper(atgc)

writefastafile(ATGC[hoge],namelist=name[hoge],outfn="table/forPICRUST")

seq2 <- seqtab[,hoge]
seq2 <- seq2[(rowSums(seq2)>4999),]
set.seed(1234)
seq2 <- rrarefy(seq2,5000)
tseq <- as.data.frame(t(seq2))

ans <- cbind(sample=rownames(tseq),tseq)

write.table(ans ,file="table/forPICRUST.tsv",sep="\t",quote=T,row.names = F)
RRR

Rscript data_process/script_pic.R

###########################
## -- Options

fasta=table/forPICRUST.fa
abundance=table/forPICRUST.tsv

outputName=table/picrusts2
mkdir ${outputName}

###########################
conda activate picrust2 #- Please install conda & picrust if you need

picrust2_pipeline.py -v > $outputName/log.txt
picrust2_pipeline.py -s $fasta -i $abundance -o $outputName/result -p 32 >> $outputName/log.txt

