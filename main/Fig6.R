library(ggplot2); library(vegan);library(gridExtra);library(ggstar)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon)
source("functions/functions.R")

#-- read data from previous dir
exp <- read.csv("table/expdata.csv")
L_tab <- readRDS("table/List_data_5000.rds")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- L_tab[[length(tax)+1]]
taxtable <- readRDS("table/taxonomylist_5000.rds")

# -- make data
sn <- function(seq){return(cbind(Sample_name=rownames(seq),as.data.frame(seq)))}
L_mer <- L_tab
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

L_dlist <- list()
for(i in 1:length(tax)){
  L_dlist[[i]] <- sp(data=L_mer[[i]],ti=ti,med=Medium)
}
# -- plot function(calculate beta-diversity & geom_boxplot)
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
    geom_histogram(binwidth=0.02)+theme_light()+facet_grid(Medium~Day)+
    xlab("Difference in community structure (<i>&beta;</i>-diversity)<br>in each pair of replicate samples")+
    ylab("Frequency")+
    theme(axis.title=element_text(size=7),axis.text=element_text(size=5),
          axis.title.x=element_markdown(),
          panel.grid.minor=element_blank(),legend.position="none",panel.spacing=unit(3,"pt"),
          strip.text=element_text(size=7),strip.text.x=element_text(color="black",margin=margin(0.5,0,0.5,0,unit="pt")),
          panel.background=element_blank(),strip.text.y=element_text(color="black",angle=0,margin=margin(0.5,1,0.5,1,unit="pt")))
  return(ans)
}
pl_list <- pl(list=L_dlist[[5]],med=Medium,ti=ti,method="bray")
plpl_list <- plpl(mat=pl_list)
#-- output
pdf("Fig6.pdf",height=10.5,width=7.2)
plpl_list
dev.off()
