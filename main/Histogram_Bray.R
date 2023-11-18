library(ggplot2); library(vegan);library(gridExtra);library(ggstar); library(ggtext)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon);library(philentropy)

#-- information of directory
exp <- read.csv("table/expdata.csv")
source("table/functions.R")
data <- readRDS("table/Data_step10_median.rds")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
KL_d <- readRDS("table/KL_calc.rds")
ti <- c(2,4,6,8,10,12)

#################################################################################
#-- Convert List to long matrix
C2L <- function(data, t=ti, med=Medium){
  maT <- list()
  for(j in 1:length(t)){
    maTT <- list()
    for(k in 1:length(med)){
      maTT[[k]] <- data[[j]][[k]]; maTT[[k]] <- cbind(as.data.frame(maTT[[k]]), cbind(Medium=rep(med[k], nrow(maTT[[k]])), Day=rep(t[j], nrow(maTT[[k]]))))
    }
    maT[[j]] <- rbind(maTT[[1]], maTT[[2]])
    for(k in 3:length(med)){
      maT[[j]] <- rbind(maT[[j]], maTT[[k]])
    }
  }
  ans <- rbind(maT[[1]], maT[[2]])
  for(j in 3:length(t)){
    ans <- rbind(ans, maT[[j]])
  }
  return(ans)
}

matL <- list()
for(i in 1:length(tax)){
  matL[[i]] <- C2L(data=data[[i]])
  matL[[i]]$Day <- factor(matL[[i]]$Day, levels=ti)
}

#####################################

onon_pl_sho <- function(mat, thr, shogen){
  if(thr>1.5){
    thc <- "red"; thw <- "#FFC0CB"
  }else if(thr>1){
    thc <- "blue"; thw <- "#87CEFA"
  }else{
    thc <- "black"; thw <- "white"
  }
  mat$Category <- factor(mat$Category, levels=c("Simulation", "Actual"))
  ans <- ggplot(mat, aes(x=ans, fill=Category))+
    geom_histogram(binwidth=0.005, position="identity", alpha=0.35)+xlab("")+ylab("")+
    scale_fill_manual(values=c("grey30", "forestgreen"))+xlim(c(0,1))+
    theme(axis.title.x=element_markdown(),
          axis.text.y=element_text(size=5),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = thw, colour = thc),
          axis.line = element_line(color = "black", linewidth=0.6), legend.position="none", panel.spacing=unit(3,"pt"))
  ans <- ans+theme(plot.margin = margin(2.5,3,-0.5,-10)) #-- manual change
  
  if(shogen==1){
    ans <- ans + theme(axis.text.x=element_blank())
  }else if(shogen==2){
    ans <- ans + theme(axis.text.x=element_text(size=5))
  }
  return(ans)
}

###########################
#- make ggplot2 object of histogram
all_gg <- list(); Fasg1 <- list(); Fasg2 <- list()
for(i in 1:length(tax)){
  all_gg[[i]] <- list(); Fasg1[[i]] <- list()
  for(j in 1:length(ti)){
    all_gg[[i]][[j]] <- list()
    for(k in 1:7){
      KL_m <- KL_d[[i]][[2]][(KL_d[[i]][[2]]$Day==ti[j])&(KL_d[[i]][[2]]$Medium==Medium[k]),1]
      all_gg[[i]][[j]][[k]] <- onon_pl_sho(mat=matL[[i]][(matL[[i]]$Medium==Medium[k])&(matL[[i]]$Day==ti[j]),],
                                       thr=KL_m, shogen=1)
    }
    KL_m <- KL_d[[i]][[2]][(KL_d[[i]][[2]]$Day==ti[j])&(KL_d[[i]][[2]]$Medium==Medium[8]),1]
    all_gg[[i]][[j]][[8]] <- onon_pl_sho(mat=matL[[i]][(matL[[i]]$Medium==Medium[8])&(matL[[i]]$Day==ti[j]),],
                                      thr=KL_m, shogen=2)
    Fasg1[[i]][[j]] <- plot_grid(all_gg[[i]][[j]][[1]], all_gg[[i]][[j]][[2]], all_gg[[i]][[j]][[3]], all_gg[[i]][[j]][[4]],
                                 all_gg[[i]][[j]][[5]], all_gg[[i]][[j]][[6]], all_gg[[i]][[j]][[7]], all_gg[[i]][[j]][[8]],
                                 ncol=1, nrow=8, scale=0.95)
  }
  Fasg2[[i]] <- plot_grid(Fasg1[[i]][[1]], Fasg1[[i]][[2]], Fasg1[[i]][[3]],
                          Fasg1[[i]][[4]], Fasg1[[i]][[5]], Fasg1[[i]][[6]], ncol=6, nrow=1, scale=0.98)
}
