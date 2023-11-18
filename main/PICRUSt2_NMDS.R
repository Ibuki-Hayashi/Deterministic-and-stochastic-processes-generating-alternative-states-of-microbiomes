library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggnewscale); library(reshape2);library(ggstar)
library(ggtext);library(ggforce)
source("functions/functions.R")

pic <- read.table("table/path_abun_unstrat.tsv",sep="\t",header=T)
data <- readRDS("table/List_data_5000.rds")
data_pic <- read.csv("table/SupplementaryDataS3.csv")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- data[[length(tax)+1]]

##### merge #####
ptab <- as.data.frame(t(pic[,2:ncol(pic)]))
p_data <- cbind(rownames(ptab),ptab); colnames(p_data) <- c("Sample_name",as.character(pic$pathway)) 
use <- merge(exp,p_data,by="Sample_name")
ds <- list()
for(i in 1:length(ti)){
  d <- use[(use$Day == ti[i]),]
  ds[[i]] <- as.numeric(str_sub(d$Sample_name,-4,-1))%%384
}
key <- intersect(ds[[1]],ds[[2]])
for(i in 3:length(ti)){
  key <- intersect(key,ds[[i]])
}
use <- use[((as.numeric(str_sub(as.character(use$Sample_name),-4,-1))%%384) %in% key),]
use <- use[use$Day %in% ti,]
################
#-- picNMDS
# -- do NMDs
NMDS_1 <- function(data,exp,Medium,ti){ #data include experiment cols
  #- separate exp&table
  expnum <- ncol(exp)
  exp <- data[,(1:expnum)]; table <- data[,((expnum+1):ncol(data))]
  #- NMDS
  set.seed(1234)
  nmds <- metaMDS(table,trymax=20)
  d1 <- cbind(exp,nmds$point)
  colnames(d1) <- c(colnames(exp),c("NMDS1","NMDS2"))
  #- NMDS points remove outliers
  n1 <- remove.outliers(d1$NMDS1); n2 <- remove.outliers(d1$NMDS2)
  data <- d1[((d1$NMDS1 %in% n1[[1]])&(d1$NMDS2 %in% n2[[1]])),]
  #- align rep & separate Day & Medium
  ds <- list()
  for(i in 1:length(ti)){
    d <- data[(data$Day == ti[i]),]
    ds[[i]] <- as.numeric(str_sub(d$Sample_name,-4,-1))%%384
  }
  key <- intersect(ds[[1]],ds[[2]])
  for(i in 3:length(ti)){
    key <- intersect(key,ds[[i]])
  }
  data2 <- data[((as.numeric(str_sub(data$Sample_name,-4,-1))%%384) %in% key),]
  
  d_list <- list()
  for(i in 1:length(Medium)){
    d_list[[i]] <- data2[(data2$Medium == Medium[i]),]
  }
  #- plot
  N_all <- ggplot(data)+
    geom_point(aes(x=NMDS1, y=NMDS2, fill=Medium, alpha=(log(data$Day,512)+0.55), stroke=0.08), shape=21, size=2)+
    geom_point(aes(x=NMDS1, y=NMDS2, fill=Medium, stroke=0.14), color="black", shape=1, size=2)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7),legend.position="none")
  N_leg <- ggplot(data)+
    geom_point(aes(x=NMDS1, y=NMDS2, fill=Medium, alpha=Day, stroke=0.15), color="black", size=2)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7))
  N_leg <- g_legend(N_leg)
    
  return(list(N_all,N_leg))
}
N_1 <- NMDS_1(data=use,exp=exp,Medium=Medium,ti=ti)

##################
###-- output --###
##################

saveRDS(N_1, "pic_nmds.rds")
