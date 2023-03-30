library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale)
library(ggplate); library(tibble)
source("functions/functions.R")

#-- read data from previous dir
data <- readRDS("table/List_data_rep.rds")
exp <- read.csv("table/expdata.csv")
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")

#########################################
############# -- For NMDS -- ############
#########################################
# -- do NMDs
tax <- c("cla","ord","fam","gen","ASV")
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC"); ti <- data[[length(tax)+1]][[2]]
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
  data2 <- data2[(data2$Day=="2")|(data2$Day=="12"),]
  
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
    lines <- data.frame(
      x = pd[(pd$Day==unique(pd$Day)[1]),]$NMDS1,
      y = pd[(pd$Day==unique(pd$Day)[1]),]$NMDS2,
      xend = pd[(pd$Day==unique(pd$Day)[2]),]$NMDS1,
      yend = pd[(pd$Day==unique(pd$Day)[2]),]$NMDS2
    )
    N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + geom_link(data=lines, color="grey70",
                                                           aes(x=x, y=y, xend=xend, yend=yend),size=0.00001,
                                                           arrow=arrow(angle=7.5,length=unit(0.2,"cm"),type="closed"))
    
    N_sepmedarrow1[[i]] <- N_sepmedarrow1[[i]] + theme(legend.position = "none",text=element_text(size=7))+
      scale_x_continuous(limits = c(-2,1.2))+scale_y_continuous(limits = c(-1.2,2.2))
  }
  return(N_sepmedarrow1)
}

# --output
alls1 <- list(); alls2 <- list(); alls3 <- list()
for(i in 5){
  alls1[[i]] <- NMDS_1(data=data[[i]],exp=exp,Medium=Medium,ti=ti)
  alls2[[i]] <- N_plots1(data=alls1[[i]][[1]],d_list=alls1[[i]][[2]],Medium=Medium,ti=ti,sh=sh)
}

FigS13_nmds <- plot_grid(alls2[[5]][[1]],alls2[[5]][[2]],alls2[[5]][[3]],alls2[[5]][[4]],ncol=1,nrow=4)
FigS14_nmds <- plot_grid(alls2[[5]][[5]],alls2[[5]][[6]],alls2[[5]][[7]],alls2[[5]][[8]],ncol=1,nrow=4)

# --outlier_well
Go <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "G")&(alls1[[5]][[1]]$NMDS2 > 1),]
GLo <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "GL")&(alls1[[5]][[1]]$NMDS2 > 1),]
GCo <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "GC"),]
GLCo <- list()
GLCo[[1]] <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "GLC")&(alls1[[5]][[1]]$NMDS2 > 1),]
GLCo[[2]] <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "GLC")&(alls1[[5]][[1]]$NMDS2 < -1),]
HGo <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "HG"),]
HGLo <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "HGL")&(alls1[[5]][[1]]$NMDS2 > 1.5),]
HGCo <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "HGC")&(alls1[[5]][[1]]$NMDS2 > 1.5),]
HGLCo <- list()
HGLCo[[1]] <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "HGLC")&(alls1[[5]][[1]]$NMDS2 > 1.5),]
HGLCo[[2]] <- alls1[[5]][[1]][((alls1[[5]][[1]]$Day == "2")|(alls1[[5]][[1]]$Day == "12"))&(alls1[[5]][[1]]$Medium == "HGLC")&(alls1[[5]][[1]]$NMDS2 < -0.5)&(alls1[[5]][[1]]$NMDS1 > 0),]

outs <- list(Go,GLo,GCo,GLCo,
            HGo,HGLo,HGCo,HGLCo)

#########################################
############# -- For Well -- ############
#########################################
#-- data2plate
data2plate <- function(data,colornum,number,Med,colorr="#FF0000"){
  if(colornum==3){
    ans <- plate_plot(data=data,position=Location,value=State,plate_size=384,
                      colour=c(palettes(unique(data$Medium))[number],colorr,"blue"),
                      title=sprintf("Medium-%s Layout",Med[number]),title_size=8,scale=0.5)+
      theme(legend.text = element_text(size=7),
            legend.title = element_text(size=7),
            axis.text = element_text(size=6),
            plot.title = element_text(size=9),)
    if(number==8){
      ans <- plate_plot(data=data,position=Location,value=State,plate_size=384,
                        colour=c("blue",palettes(unique(data$Medium))[number],colorr),
                        title=sprintf("Medium-%s Layout",Med[number]),title_size=8,scale=0.5)+
        theme(legend.text = element_text(size=7),
              legend.title = element_text(size=7),
              axis.text = element_text(size=6),
              plot.title = element_text(size=9))
    }
  }else{
    if(colornum==1){
      ans <- plate_plot(data=data,position=Location,value=State,plate_size=384,
                        colour=c(palettes(unique(data$Medium))[number],palettes(unique(data$Medium))[number]),
                        title=sprintf("Medium-%s Layout",Med[number]),title_size=8,scale=0.5)+
        theme(legend.text = element_text(size=7),
              legend.title = element_text(size=7),
              axis.text = element_text(size=6),
              plot.title = element_text(size=9))
    }else{
      ans <- plate_plot(data=data,position=Location,value=State,plate_size=384,
                        colour=c(palettes(unique(data$Medium))[number],colorr),
                        title=sprintf("Medium-%s Layout",Med[number]),title_size=8,scale=0.5)+
        theme(legend.text = element_text(size=7),
              legend.title = element_text(size=7),
              axis.text = element_text(size=6),
              plot.title = element_text(size=9))
    }
  }
  return(ans)
}

#-- all plate
fa <- exp[exp$Day==2,] #big letter
fa$Location <- sapply(as.character(fa$Location),
       function(value){
         v1 <- toupper(str_sub(value,1,1))
         return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})

fat <- as_tibble(fa)

all_plate <- plate_plot(
  data=fat,
  position=Location,
  value=Medium,
  plate_size=384,
  colour=palettes(unique(fat$Medium)),
  title_size=8,
)+theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size=9))

pdf("FigS13.pdf",width=7.2,height=8)
all_plate
dev.off()

#-- respective plates
res_plate <- list(); fat2 <- fat
#-- G plate
j<-1; key <- sapply(as.character(outs[[j]]$Location),
                    function(value){
                      v1 <- toupper(str_sub(value,1,1))
                      return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=2,number=j,Med=Medium)
#-- GL plate
j<-2; key <- sapply(as.character(outs[[j]]$Location),
                    function(value){
                      v1 <- toupper(str_sub(value,1,1))
                      return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=2,number=j,Med=Medium)
#-- GC plate
j<-3
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(i%%2 ==0){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=1,number=j,Med=Medium)
#-- GLC plate
j<-4
key1 <- sapply(as.character(outs[[j]][[1]]$Location),
                function(value){
                  v1 <- toupper(str_sub(value,1,1))
                  return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
key2 <- sapply(as.character(outs[[j]][[2]]$Location),
               function(value){
                 v1 <- toupper(str_sub(value,1,1))
                 return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key1){
      vv2[i] <- "State 2"
    }else{
      if(fat2$Location[i] %in% key2){
        vv2[i] <- "State 3"
      }
    }
  }
}
fat2$State <- vv2
res_plate[[j]] <- data2plate(data=fat2,colornum=3,number=j,Med=Medium)
#-- HG plate
j<-5
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(i%%2 ==0){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=1,number=j,Med=Medium)
#-- HGL plate
j<-6; key <- sapply(as.character(outs[[j]]$Location),
                  function(value){
                    v1 <- toupper(str_sub(value,1,1))
                    return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=2,number=j,Med=Medium)
#-- HGC plate
j<-7; key <- sapply(as.character(outs[[j]]$Location),
                    function(value){
                      v1 <- toupper(str_sub(value,1,1))
                      return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key){
      vv2[i] <- "State 2"
    }
  }
}
fat2$State <- vv2; res_plate[[j]] <- data2plate(data=fat2,colornum=2,number=j,Med=Medium)
#-- HGLC plate
j<-8
key1 <- sapply(as.character(outs[[j]][[1]]$Location),
               function(value){
                 v1 <- toupper(str_sub(value,1,1))
                 return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
key2 <- sapply(as.character(outs[[j]][[2]]$Location),
               function(value){
                 v1 <- toupper(str_sub(value,1,1))
                 return(paste(v1,str_sub(value,2,nchar(value)),sep=""))})
vv2 <- rep(NA,nrow(fat2))
for(i in 1:length(vv2)){
  if(fat2$Medium[i] == Medium[j]){
    vv2[i] <- "State 1"
    if(fat2$Location[i] %in% key1){
      vv2[i] <- "State 2"
    }else{
      if(fat2$Location[i] %in% key2){
        vv2[i] <- "State 3"
      }
    }
  }
}
fat2$State <- vv2
res_plate[[j]] <- data2plate(data=fat2,colornum=3,number=j,Med=Medium)
#-- out 
FigS13_well <- plot_grid(res_plate[[1]],res_plate[[2]],res_plate[[3]],res_plate[[4]],ncol=1)
FigS14_well <- plot_grid(res_plate[[5]],res_plate[[6]],res_plate[[7]],res_plate[[8]],ncol=1)

#########################################
########### -- Output-- #################
#########################################
#-- plot_grid
NW13 <- plot_grid(FigS13_nmds,FigS13_well,ncol=2,rel_widths=c(0.6,1))
NW14 <- plot_grid(FigS14_nmds,FigS14_well,ncol=2,rel_widths=c(0.6,1))
#- plot
pdf("FigS14.pdf",width=7.2,height=11)
NW13
dev.off()
pdf("FigS15.pdf",width=7.2,height=11)
NW14
dev.off()
