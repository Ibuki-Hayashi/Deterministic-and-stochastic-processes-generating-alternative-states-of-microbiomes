library(ggplot2); library(vegan);library(gridExtra);library(ggstar)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon)
source("functions/functions.R")

#-- read data from previous dir
data <- readRDS("table/List_data_rep.rds")
exp <- read.csv("table/expdata.csv")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC"); ti <- data[[length(tax)+1]][[2]]
sh <- c(1, 13, 15, 11, 12, 14, 29, 2, 27)

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
  data <- d1[((d1$NMDS1 %in% n1[[1]])&(d1$NMDS2 %in% n2[[1]])),]; data$Day <- as.factor(data$Day)
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
N_plots1 <- function(data,d_list,Medium,ti,sh=sh){
  N_sepmedarrow1 <- list()
  for(i in 1:length(Medium)){
    pd <- d_list[[i]]
    N_sepmedarrow1[[i]] <- ggplot(pd)+
      geom_star(aes(x=NMDS1, y=NMDS2, starshape=Day, starstroke=0.15),fill=palettes(Medium)[i],size=1.5)+
      ggtitle(as.character(Medium[i])) + xlab("NMDS1") + ylab("NMDS2")+theme_light()
    for(j in 1:(length(unique(pd$Day))-1)){
      lines <- data.frame(
        x = pd[(pd$Day==unique(pd$Day)[j]),]$NMDS1,
        y = pd[(pd$Day==unique(pd$Day)[j]),]$NMDS2,
        xend = pd[(pd$Day==unique(pd$Day)[j+1]),]$NMDS1,
        yend = pd[(pd$Day==unique(pd$Day)[j+1]),]$NMDS2
      )
      N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + geom_link(data=lines, color="grey70",
                                                             aes(x=x, y=y, xend=xend, yend=yend),size=0.00001,
                                                             arrow=arrow(angle=7.5,length=unit(0.2,"cm"),type="closed"))
    }
    N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + theme(legend.position = "none",text=element_text(size=7))
  }
  return(N_sepmedarrow1)
}
Leg_plots <- function(data,d_list,Medium,ti,sh=sh){
  N_sepmedarrow1 <- list()
  for(i in 1:length(Medium)){
    pd <- d_list[[i]]
    N_sepmedarrow1[[i]] <- ggplot(pd)+
      geom_star(aes(x=NMDS1, y=NMDS2, starshape=Day, starstroke=0.25), fill=palettes(Medium)[i],size=1.5)+
      ggtitle(as.character(Medium[i])) + xlab("NMDS1") + ylab("NMDS2")+theme_light()
    N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + theme(text=element_text(size=7))
    N_sepmedarrow1[[i]] <- g_legend(N_sepmedarrow1[[i]])
  }
  return(N_sepmedarrow1)
}
# --output
alls1 <- list(); alls2 <- list(); alls3 <- list()
for(i in 3:length(tax)){
  alls1[[i]] <- NMDS_1(data=data[[i]],exp=exp,Medium=Medium,ti=ti)
  alls2[[i]] <- N_plots1(data=alls1[[i]][[1]],d_list=alls1[[i]][[2]],Medium=Medium,ti=ti,sh=sh)
  alls3[[i]] <- Leg_plots(data=alls1[[i]][[1]],d_list=alls1[[i]][[2]],Medium=Medium,ti=ti,sh=sh)
}

ASN <- plot_grid(alls2[[5]][[1]],alls2[[5]][[2]],alls2[[5]][[3]],alls2[[5]][[4]],
                 alls2[[5]][[5]],alls2[[5]][[6]],alls2[[5]][[7]],alls2[[5]][[8]],ncol=2,nrow=4)

pdf("Fig5.pdf",height=10.5,width=7.2)
plot_grid(ASN,plot_grid(alls3[[5]][[1]],NA,NA,ncol=1,rel_heights=c(0.7,1,1)),ncol=2,rel_widths=c(1,0.3))
dev.off()
