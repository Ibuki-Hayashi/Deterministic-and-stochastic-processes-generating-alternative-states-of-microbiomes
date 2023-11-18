library(ggplot2); library(vegan);library(gridExtra); library(cowplot); library(lemon)
library(tidyr); library(stringr);library(ggnewscale); library(reshape2);library(ggtext)
source("functions/functions.R")

#-- read data from previous dir
exp <- read.csv("table/expdata.csv")
pic <- read.table("table/path_abun_unstrat.tsv",sep="\t",header=T)
data <- readRDS("table/List_data_5000.rds")

##### merge #####
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- data[[length(tax)+1]]
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

# -- listにDay(ti)とMedium(med)に分けて格納
sp <- function(data,ti=ti,med=med){
  ans <- list()
  for(i in 1:length(ti)){
    ans[[i]] <- list()
    for(j in 1:length(med)){
      ans[[i]][[j]] <- data[(data$Day == ti[i])&(data$Medium == med[j]),]
    }
  }
  return(ans)
}
pic_d <- sp(data=use,ti=ti,med=Medium)

# -- plot function(beta-diversityを計算してDay数だけのgg_boxplotを生成)
pl <- function(list,med=med,ti=ti,method="bray"){
  sample_beta <- function(seqtab,exp=exp,method=method){
    times <- (nrow(seqtab)*(nrow(seqtab)-1))/2
    exp_length <- ncol(exp); ans <- rep(NA,times); count <- 0
    for(i in 1:nrow(seqtab)){
      for(j in 1:nrow(seqtab)){
        if(i > j){
          count <- count+1
          ans[count] <- vegdist(rbind(seqtab[i,(exp_length+1):ncol(seqtab)],seqtab[j,(exp_length+1):ncol(seqtab)]),method=method)
        }
      }
    }
    return(ans)
  }
  aa <- list(); for(i in 1:length(list)){
    ml <- list()
    ms <- c(); bs <- c(); for(j in 1:length(list[[i]])){
      bs <- sample_beta(seqtab=list[[i]][[j]],exp=exp,method=method)
      ms <- rep(med[j],nrow(list[[i]][[j]])*(nrow(list[[i]][[j]])-1)/2)
      ml[[j]] <- data.frame(Medium=ms,beta=bs)
    }
    bms <- rbind(ml[[1]],ml[[2]])
    for(j in 3:length(list[[i]])){
      bms <- rbind(bms,ml[[j]])
    }
    aa[[i]] <- data.frame(bms,Day=rep(ti[i],nrow(bms)))
  }
  ans <- rbind(aa[[1]],aa[[2]])
  for(i in 3:length(list)){
    ans <- rbind(ans,aa[[i]])
  }
  return(ans)
}

plpl <- function(mat){
  ans <- ggplot(mat,aes(x=beta))+
    geom_histogram(binwidth=0.005)+theme_light()+facet_grid(Medium~Day)+
    xlab("Difference in community structure (<i>&beta;</i>-diversity)<br>in each pair of replicate samples")+
    ylab("Frequency")+xlim(c(0,1.0))+ylim(c(0,150))+
    #scale_x_continuous(expand=c(0,0,1,0))+
    theme(axis.title=element_text(size=7),axis.text=element_text(size=5),
          axis.title.x=element_markdown(),
          panel.grid.minor=element_blank(),legend.position="none",panel.spacing=unit(3,"pt"),
          strip.text=element_text(size=7),strip.text.x=element_text(color="black",margin=margin(0.5,0,0.5,0,unit="pt")),
          panel.background=element_blank(),strip.text.y=element_text(color="black",angle=0,margin=margin(0.5,1,0.5,1,unit="pt")))
  return(ans)
}

pl_p <- pl(list=pic_d,med=Medium,ti=ti,method="bray")
plpl_p <- plpl(mat=pl_p)

#-- output
pdf("FigS17.pdf",height=10.8,width=7.2)
plpl_p
dev.off()
