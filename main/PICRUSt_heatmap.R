library(ggplot2); library(vegan);library(gridExtra);library(ggstar);library(lemon)
library(tidyr); library(stringr);library(ggforce);library(cowplot); library(ggtext)
library(cluster); library(clusterSim); library(reshape2); library(RColorBrewer)
library(philentropy); library(NbClust); library(factoextra); library(Rmisc); library(concaveman)
source("functions/functions.R")

#-- information of directory
pic <- read.table("table/path_abun_unstrat_descrip.tsv",sep="\t",header=T)
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/List_data_rep.rds")
ptab <- as.data.frame(t(pic[,3:ncol(pic)]))
p_data <- cbind(rownames(ptab), ptab); colnames(p_data) <- c("Sample_name",as.character(pic$description)) 
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- c(2,4,6,8,10,12)

################################################################################
#-- k-means
data2 <- list()
for(i in 1:5){
  data2[[i]] <- list()
  for(j in 1:length(ti)){
    data2[[i]][[j]] <- list()
    for(k in 1:length(Medium)){
      data2[[i]][[j]][[k]] <- data[[i]][(data[[i]]$Day==ti[j])&(data[[i]]$Medium==Medium[k]),]
      data2[[i]][[j]][[k]] <- data2[[i]][[j]][[k]][,6:ncol(data2[[i]][[j]][[k]])]
    }
  }
}

kMed2gg <- function(data){
  BCd <- vegdist(data, method="bray")
  Clu <- list()
  max_Si <- -2000000
  for(i in 2:(nrow(data)-1)){
    Clu[[i-1]] <- pam(BCd, k=i, diss=T)$clustering
    BC_Si <- pam(BCd, k=i, diss=T)$silinfo$avg.width
    if(BC_Si>max_Si){
      sii <- i; max_Si <- BC_Si
    }
  }
  ans <- data.frame(Sample_name=rownames(data),
                    Si_Cluster=Clu[[sii-1]])
  return(ans)
}

ggans <- list()
for(i in 3:5){
  ggans[[i]] <- list()
  for(j in 1:length(ti)){
    ggans[[i]][[j]] <- list()
    for(k in 1:length(Medium)){
      ggans[[i]][[j]][[k]] <- kMed2gg(data=data2[[i]][[j]][[k]])
      message(sprintf("(i,j,k)=(%s,%s,%s)/(5,6,8) was finished",i,j,k))
    }
  }
}
uuu <- ggans

##### merge #####
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
pic_list <- list(); PIC_new <- list()
for(i in 3:length(tax)){
  pic_list[[i]] <- list(); PIC_new[[i]] <- list()
  for(j in 1:length(ti)){
    pic_list[[i]][[j]] <- list()
    for(k in 1:length(Medium)){
      hoge <- use[(use$Day == ti[j])&(use$Medium == Medium[k]),]
      pic_list[[i]][[j]][[k]] <- merge(uuu[[i]][[j]][[k]], hoge, by="Sample_name")
      pic_list[[i]][[j]][[k]]$Si_Cluster <- sprintf("%s-%s", pic_list[[i]][[j]][[k]]$Medium, pic_list[[i]][[j]][[k]]$Si_Cluster)
    }
    PIC_new[[i]][[j]] <- rbind(pic_list[[i]][[j]][[1]], pic_list[[i]][[j]][[2]])
    for(k in 3:length(Medium)){
      PIC_new[[i]][[j]] <- rbind(PIC_new[[i]][[j]], pic_list[[i]][[j]][[k]])
    }
  }
}

#-- logFC
logfc <- function(data){
  data <- data[rowSums(data[,(7:ncol(data))])>0,] #data <- PIC_new[[i]][[6]]
  key <- unique(data$Si_Cluster)
  ls <- list()
  for(i in 1:length(key)){
    ls[[i]] <- (colSums(data[(data$Si_Cluster == key[i]), (7:ncol(data))])+0.0001)/nrow(data[(data$Si_Cluster == key[i]), (7:ncol(data))])
  }
  ans <- as.data.frame(matrix(NA, ncol=length(key), nrow=length(ls[[1]]))); colnames(ans) <- key
  for(i in 1:length(key)){
    ans[,i] <- ls[[i]]
  }
  la <- log(ans/(rowSums(ans)/length(key)), 10)
  rownames(la) <- colnames(data[, (7:ncol(data))])
  
  return(list(ans, la))
}

picPC <- list(); picFC <- list()
for(i in 3:length(tax)){
  picPC[[i]] <- logfc(PIC_new[[i]][[6]])[[2]]
  a <- diversity(logfc(PIC_new[[i]][[6]])[[1]], index="shannon")
  
  thr <- 30
  picFC[[i]] <- picPC[[i]][(a < a[order(a, decreasing = F)][thr+1]),]
  picFC[[i]] <- cbind(Pathway=rownames(picFC[[i]]), picFC[[i]])
  picFC[[i]] <- melt(picFC[[i]]); colnames(picFC[[i]]) <- c("Pathway", "Cluster", "logFC")
}

# --plot
ppp <- function(data, title){
  ans <- ggplot(data, aes(x=Cluster, y=Pathway, fill=logFC))+
    geom_tile()+theme_bw()+ggtitle(title)+
    theme(plot.background = element_blank(),
          plot.title = element_text(size=7),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.position="none",
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = "white", colour = "white"),
          axis.text.x = element_text(size=6, angle=60, hjust=1, vjust=1),
          axis.text.y = element_text(size=6))+
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")[4:9]), na.value = "white")
  return(ans)
}

ppp_l <- function(data){
  ans <- ggplot(data, aes(x=Cluster, y=Pathway, fill=logFC))+
    geom_tile()+theme_bw()+
    theme(plot.background = element_blank(),
          plot.title = element_text(size=10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = "white", colour = "white"),
          axis.text.x = element_text(size=6, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size=6, vjust = 0.5, hjust = 1))+
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")[4:9]), na.value = "white")
  return(ans)
}

ff <- c("Family", "Genus", "ASV")

#-- ggplot2 object of PICRUSt heatmap
plots <- list()
for(i in 3:length(tax)){
  plots[[i]] <- ppp(data=picFC[[i]], title=ff[i-2])
}
