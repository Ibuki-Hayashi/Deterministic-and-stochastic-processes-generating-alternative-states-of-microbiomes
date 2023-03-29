library(ggplot2); library(vegan);library(gridExtra);library(ggstar);library(lemon)library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
source("functions/functions.R")

# -- read data from previous dir
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/List_data_5000.rds")
data2 <- read.csv("table/SupplementaryDataS3.csv")

#########################################
########### -- For Fig.3A -- ############
#########################################
data2 <- data2[(data2$Day!="all"),]
data2[,3] <- as.numeric(str_sub(data2[,3],start=4))
data2[,5] <- as.character(data2[,5])

lev <- c("family","genus","ASV","PICRUSTs2")
lev2 <- paste(c("Family","Genus","ASV","Functional compositions by PICRUSTs2"),c(rep("level communities",3),""),sep=" ")

d_lb <- list()
for(i in 1:length(lev)){
  d_lb[[i]] <- data2[(data2[,2]==lev[i])&(data2[,4]=="bray"),]
}
#-- plot
plot_notall <- function(dd){
  dd[,5][dd[,5]=="glucose"] <- "Concentration<br>of glucose"
  dd[,5][dd[,5]=="leucine"] <- "Leucine"
  dd[,5][dd[,5]=="citrate"] <- "Citrate"
  FDRq <- as.character(dd$qs<0.05); FDRq <- factor(FDRq,levels=c("TRUE","FALSE"))
  ans <- ggplot(dd,aes(x=Day,y=R2))+
    geom_point(size=1.2, aes(shape=FDRq))+geom_line(aes(group=Treatment,color=Treatment),size=0.4)+
    scale_color_manual(values = palettes(c(as.factor(c(dd$Treatment)))))+
    theme_light()+theme(text=element_text(size=7))+
    scale_x_continuous(breaks=seq(0,12,2))
  ans <- ans +
    ylab(paste("<i>R^2</i>","of PERMANOVA<br>for the communities of the genus level"))+
    theme(axis.title.y=element_markdown(),legend.text=element_markdown())
  return(ans)
}
plots_na <- list()
for(i in 1:length(lev)){
  plots_na[[i]] <- plot_notall(dd=d_lb[[i]])
}
FigA <- plots_na[[2]]

#########################################
########### -- For Fig.3B -- ############
#########################################
# -- taxa num
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- data[[length(tax)+1]]
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
t_box_jit_a <- function(ddd){
  box <- ggplot(ddd,aes(x=Medium,y=alpha,starshape=Day,starstroke=0.1))+
    geom_violin(aes(fill=Medium),lwd=0.1)+scale_fill_manual(values = palettes(unique(ddd$Medium)))+
    new_scale_fill()+ylab("<i>&alpha;</i> diversity of community<br>(Shannon's <i>H'</i>)")+
    theme_light()+theme(text=element_text(size=7),axis.title.y=element_markdown())+
    theme(legend.position="none")
  return(box)
}
FigB <- t_box_jit_a(ddd=t_list[[5]])

#########################################
########### -- For Fig.3C -- ############
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
N_plots <- function(data,d_list,Medium,ti,sh=sh){
  N_all <- ggplot(data)+
    geom_star(aes(x=NMDS1, y=NMDS2,fill=Medium,starshape=Day,starstroke=0.25),size=2.5)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7),legend.position="none")
  return(N_all)
}
N_legends <- function(data,d_list,Medium,ti,sh=sh){
  N_all <- ggplot(data)+
    geom_star(aes(x=NMDS1, y=NMDS2,fill=Medium,starshape=Day,starstroke=0.25),size=2.5)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7))
  ans <- g_legend(N_all); return(ans)
}
# --output
alls1 <- NMDS_1(data=data[[5]],exp=exp,Medium=Medium,ti=ti)
FigC <- N_plots(data=alls1[[1]],d_list=alls1[[2]],Medium=Medium,ti=ti,sh=sh)
FigC_leg <- N_legends(data=alls1[[1]],d_list=alls1[[2]],Medium=Medium,ti=ti,sh=sh)

#########################################
########### -- Output-- #################
#########################################
#-- plot_grid
AB <- plot_grid(FigA,FigB,ncol=2,labels=c("A","B"),rel_widths=c(0.6,1))
C <- plot_grid(FigC,FigC_leg,ncol=2,labels=c("C",""),rel_widths=c(1,0.12))
ABCD <- plot_grid(AB,C,nrow=2,rel_heights=c(0.4,1))

#- plot
pdf("Fig3.pdf",width=7.2,height=8.5)
ABCD
dev.off()
