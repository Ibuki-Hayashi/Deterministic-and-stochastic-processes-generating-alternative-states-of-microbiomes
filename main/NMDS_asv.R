library(ggplot2); library(vegan);library(gridExtra);library(ggstar);library(lemon)library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
source("functions/functions.R")

# -- read data from previous dir
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/List_data_5000.rds")
data2 <- read.csv("table/SupplementaryDataS3.csv")

#########################################
########### --    NMDS    -- ############
#########################################
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC"); ti <- data[[length(tax)+1]]
sh <- c(1, 13, 15, 11, 12, 14, 29, 2, 27)
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
  data$Medium <- factor(data$Medium,levels=Medium)
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
  ans <- list(data,d_list); return(ans)
}
#- plot
N_plots <- function(data,d_list,Medium,ti,sh=sh){
  N_all <- ggplot(data)+
    geom_point(aes(x=NMDS1, y=NMDS2, fill=Medium, alpha=(log(data$Day,512)+0.55),stroke=0.08), color="black", shape=21,size=2)+
    geom_point(aes(x=NMDS1, y=NMDS2, stroke=0.14), color="black", shape=1,size=2)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7),legend.position="none")
  return(N_all)
}
N_legends <- function(data,d_list,Medium,ti,sh=sh){
  N_all <- ggplot(data)+
    geom_star(aes(x=NMDS1, y=NMDS2,fill=Medium,alpha=Day,starstroke=0.25), color="black", starshape=15, size=1.8)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7))
  ans <- g_legend(N_all); return(ans)
}

# --output
alls1 <- NMDS_1(data=data[[5]],exp=exp,Medium=Medium,ti=ti)
asvNMDS <- N_plots(data=alls1[[1]],d_list=alls1[[2]],Medium=Medium,ti=ti,sh=sh)
asvleg <- N_legends(data=alls1[[1]],d_list=alls1[[2]],Medium=Medium,ti=ti,sh=sh)
