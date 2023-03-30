library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggnewscale); library(reshape2);library(ggstar)
library(ggtext);library(ggforce)
source("functions/functions.R")

pic <- read.table("table/path_abun_unstrat.tsv",sep="\t",header=T)
data <- readRDS("table/List_data_5000.rds")
data_pic <- read.csv("table/SupplementaryDataS3.csv")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- data[[length(tax)+1]]

#########################################
########### -- For Fig.7A -- ############
#########################################
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
  data <- d1[((d1$NMDS1 %in% n1[[1]])&(d1$NMDS2 %in% n2[[1]])),]; data$Day <- as.factor(data$Day)
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
    geom_star(aes(x=NMDS1, y=NMDS2,fill=Medium,starshape=Day,starstroke=0.25),size=2.5)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7),legend.position="none")
  N_leg <- ggplot(data)+
    geom_star(aes(x=NMDS1, y=NMDS2,fill=Medium,starshape=Day,starstroke=0.25),size=2.5)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_fill_manual(values=palettes(unique(data$Medium)))+
    theme_light()+theme(text=element_text(size=7))
  N_leg <- g_legend(N_leg)
  
  return(list(N_all,N_leg))
}
N_1 <- NMDS_1(data=use,exp=exp,Medium=Medium,ti=ti)
FigA1 <- N_1[[1]]; FigA2 <- N_1[[2]]

#########################################
########### -- For Fig.7B -- ############
#########################################
data_pic <- data_pic[(data_pic$Day!="all"),]
data_pic[,3] <- as.numeric(str_sub(data_pic[,3],start=4))
data_pic[,5] <- as.character(data_pic[,5])

lev <- c("family","genus","ASV","PICRUSTs2")

d_lb <- list()
for(i in 1:length(lev)){
  d_lb[[i]] <- data_pic[(data_pic[,2]==lev[i])&(data_pic[,4]=="bray"),]
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
    ylab(paste("<i>R^2</i>","of PERMANOVA<br>for the functional compositions"))+
    theme(axis.title.y=element_markdown(),legend.text=element_markdown())
  return(ans)
}
plots_na <- list()
for(i in 1:length(lev)){
  plots_na[[i]] <- plot_notall(dd=d_lb[[i]])
}
FigB <- plots_na[[4]]

#########################################
########### -- For Fig.7C -- ############
#########################################
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
  data <- d1[((d1$NMDS1 %in% n1[[1]])&(d1$NMDS2 %in% n2[[1]])),]; data$Day <- as.factor(data$Day)
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
  N_sepmedarrow1 <- list()
  for(i in 1:length(Medium)){
    pd <- d_list[[i]]
    N_sepmedarrow1[[i]] <- ggplot(pd)+
      geom_star(aes(x=NMDS1, y=NMDS2,starshape=Day, fill=Medium, starstroke=0.15),fill=palettes(unique(Medium))[i],size=1.5)+
      ggtitle(as.character(Medium[i])) + xlab("NMDS1") + ylab("NMDS2")+theme_light()
    for(j in 1:(length(unique(pd$Day))-1)){
      lines <- data.frame(
        x = pd[(pd$Day==unique(pd$Day)[j]),]$NMDS1,
        y = pd[(pd$Day==unique(pd$Day)[j]),]$NMDS2,
        xend = pd[(pd$Day==unique(pd$Day)[j+1]),]$NMDS1,
        yend = pd[(pd$Day==unique(pd$Day)[j+1]),]$NMDS2
      )
      N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + geom_link(data=lines,
                                                             aes(x=x, y=y, xend=xend,yend=yend),
                                                             size=0.006,color="grey70",
                                                             arrow=arrow(angle=7.5,length=unit(0.2,"cm"),
                                                                         type="closed"))
    }
    N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + theme(legend.position = "none",text=element_text(size=7))
  }
  return(N_sepmedarrow1)
}
FigC <- plot_grid(alls[[4]],alls[[7]],ncol=2)

#########################################
########### -- For Fig.7D -- ############
#########################################
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
################

plpl <- function(mat){
  ans <- ggplot(mat,aes(x=beta))+
    geom_histogram(binwidth=0.005)+theme_light()+facet_grid(Medium~Day)+
    xlab("Difference in community structure (<i>&beta;</i>-diversity)<br>in each pair of replicate samples")+
    ylab("Frequency")+xlim(c(0,1))+ylim(c(0,150))+
    theme(axis.title=element_text(size=7),panel.grid.minor=element_blank(),legend.position="none",panel.spacing=unit(3,"pt"),
          strip.text=element_text(size=7),strip.text.x=element_text(color="black",margin=margin(0.5,0,0.5,0,unit="pt")),
          panel.background=element_blank(),strip.text.y=element_text(color="black",angle=0,margin=margin(0.5,1,0.5,1,unit="pt")))+
    theme(axis.text=element_text(size=5),axis.title.x = element_markdown())
  return(ans)
}

pl_p <- pl(list=pic_d,med=Medium,ti=ti,method="bray")
FigD <- plpl(mat=pl_p[(pl_p$Medium == "GLC")|(pl_p$Medium == "HGC"),])

##################
###-- output --###
##################

#-- plot_grid
A <- plot_grid(FigA1,FigA2,nrow=1,rel_widths=c(1,0.25),labels=c("A",""))
BC <- plot_grid(FigB,FigC,nrow=1,rel_widths=c(0.7,0.8),labels=c("B","C"))
ans <- plot_grid(A,BC,FigD,nrow=3,rel_heights=c(1.2,0.45,0.35),labels=c("","","D"))

#- plot
pdf("Fig7.pdf",height=10,width=7.2)
ans
dev.off()
