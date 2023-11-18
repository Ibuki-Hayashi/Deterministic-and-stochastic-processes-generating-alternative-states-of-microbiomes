library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale)
source("functions/functions.R")

#-- read data from previous dir
data <- readRDS("table/List_data_5000.rds")
exp <- read.csv("table/expdata.csv")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- data[[length(tax)+1]]

# -- taxa num
seqk <- function(data=data,exp=exp){expnum<-ncol(exp);return(data[,((expnum+1):ncol(data))])}
tanum <- function(data=data,exp=exp){
  taxa01 <- function(data){
    ans1 <- c()
    for(i in 1:nrow(data)){
      ans1[i] <- 0
      for(j in 1:ncol(data)){
        if(data[i,j] != 0){(ans1[i] <- ans1[i]+1)}
      }
    }
    return(ans1)
  }
  return(taxa01(data=seqk(data=data,exp=exp)))
}

# -- all matrix make
t_list <- list()
for(i in 1:length(tax)){
  mat <- as.data.frame(matrix(NA,nrow=nrow(data[[i]]),ncol=6))
  mat[,1] <- data[[i]]$Medium; mat[,2] <- as.factor(data[[i]]$Day); mat[,3] <- rep(data[[(length(tax)+1)]][[1]],length(data[[(length(tax)+1)]][[2]]))
  mat[,4] <- tanum(data=data[[i]],exp=exp); mat[,5] <- diversity(seqk(data[[i]],exp=exp))
  mat[,6] <- as.factor(rep(tax[i],nrow(data[[i]]))); colnames(mat) <- c("Medium","Day","rep","number","alpha","tax")
  mat[,1] <- factor(mat[,1],levels=Medium)
  
  t_list[[i]] <- mat
}

# -- violinplot
#box_med
t_box_jit_n <- function(ddd, title){
  box <- ggplot(ddd,aes(x=Day,y=number))+
    geom_violin(aes(fill=Medium),lwd=0.1)+
    geom_line(aes(x=Day, y=number, group=1), stat="summary", fun="median")+
    facet_wrap(.~Medium, ncol=4)+
    scale_fill_manual(values = palettes(unique(ddd$Medium)))+
    stat_summary(fun="median", geom="point", aes(group=Day),
                 shape=23, size=1.5, fill="white")+
    ylab(sprintf("Number of %s", title))+
    theme_light()+theme(text=element_text(size=7),axis.title.y=element_markdown())+
    theme(legend.position="none")
  return(box)
}

t_box_jit_a <- function(ddd){
  box <- ggplot(ddd,aes(x=Day,y=alpha))+
    geom_violin(aes(fill=Medium),lwd=0.1)+
    geom_line(aes(x=Day, y=alpha, group=1), stat="summary", fun="median")+
    facet_wrap(.~Medium, ncol=4)+
    scale_fill_manual(values = palettes(unique(ddd$Medium)))+
    stat_summary(fun="median", geom="point", aes(group=Day),
                 shape=23, size=1.5, fill="white")+
    ylab("<i>&alpha;</i> diversity of community<br>(Shannon's <i>H'</i>)")+
    theme_light()+theme(text=element_text(size=7),axis.title.y=element_markdown())+
    theme(legend.position="none")
  return(box)
}

legg <- function(ddd){
  box <- ggplot(ddd,aes(x=Medium,y=alpha,starshape=Day,starstroke=0.1))+
    geom_violin(aes(fill=Medium),lwd=0.1)+scale_fill_manual(values = palettes(unique(ddd$Medium)))+
    new_scale_fill()+ylab("<i>&alpha;</i> diversity of community<br>(Shannon's <i>H'</i>)")+
    geom_star(aes(fill="black",color="black",colour="device"),fill="black",
              position=position_jitterdodge(jitter.width=0.1,dodge.width=0.85),size=0.8)+
    theme_light()+theme(text=element_text(size=7),axis.title.y=element_markdown())
  ans <- g_legend(box); return(ans)
}

npl <- list(); apl <- list(); lpl <- list()
a <- c("families", "genera", "taxa of ASVs")
for(i in 3:length(tax)){
  apl[[i]] <- t_box_jit_a(ddd=t_list[[i]])
  npl[[i]] <- t_box_jit_n(ddd=t_list[[i]], title=a[i-2])
  lpl[[i]] <- legg(ddd=t_list[[i]])
}

apap <- plot_grid(apl[[5]],apl[[4]],apl[[3]],ncol=1,labels="AUTO")
apap <- plot_grid(apap,plot_grid(NA,lpl[[4]],NA,ncol=1,rel_heights=c(0.2,1,1)),rel_widths=c(1,0.15))
npnp <- plot_grid(npl[[5]],npl[[4]],npl[[3]],ncol=1,labels="AUTO")
npnp <- plot_grid(npnp,plot_grid(NA,lpl[[4]],NA,ncol=1,rel_heights=c(0.2,1,1)),rel_widths=c(1,0.15))

# -- output
pdf("FigS3.pdf",width=7.2,height=10)
apap
dev.off()
pdf("FigS4.pdf",width=7.2,height=10)
npnp
dev.off()
