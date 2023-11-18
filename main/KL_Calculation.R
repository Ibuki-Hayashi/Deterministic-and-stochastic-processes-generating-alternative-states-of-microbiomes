library(ggplot2); library(vegan);library(gridExtra);library(ggstar); library(ggtext)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon);library(philentropy)

#-- read data from previous dir
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/Data_step10_median.rds")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
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

# -- plot function(beta-diversityを計算してDay数だけのgg_boxplotを生成)
plpl <- function(mat){
  ans <- ggplot(mat,aes(x=ans, fill=Category))+
    geom_histogram(binwidth=0.005, position="identity", alpha=0.4)+theme_light()+facet_grid(Medium~Day)+
    xlab("Difference in community structure (<i>&beta;</i>-diversity)<br>in each pair of replicate samples")+
    ylab("Frequency")+xlim(c(0,1))+
    theme(axis.title=element_text(size=6.5),axis.text=element_text(size=5),
          axis.title.x=element_markdown(),
          panel.grid.minor=element_blank(),legend.position="none",panel.spacing=unit(3,"pt"),
          strip.text=element_text(size=5),strip.text.x=element_text(color="black",margin=margin(0.5,0,0.5,0,unit="pt")),
          panel.background=element_blank(),strip.text.y=element_text(color="black",angle=0,margin=margin(0.5,1,0.5,1,unit="pt")))
  return(ans)
}

plpl_list <- list()
for(i in 1:length(tax)){
  plpl_list[[i]] <- plpl(mat=matL[[i]])
}

#################################################################################
KL_plot <- function(data, t=ti, med=Medium){
  KL_value <- as.data.frame(matrix(NA, ncol=3, nrow=length(med)*length(t)))
  colnames(KL_value) <- c("KL", "Day", "Medium")
  
  a <- c()
  tt <- c(2,4,6,8,10,12)
  for(i in 1:length(t)){ a <- c(a, rep(tt[i], length(med)))}
  a <- factor(a, levels=tt)
  
  b <- rep(med, length(t))
  b <- factor(b, levels=med)
  KL_value[,2] <- a; KL_value[,3] <- b
  
  for(i in 1:nrow(KL_value)){
    use <- data[(data$Day == KL_value[i, 2])&(data$Medium == KL_value[i, 3]),]
    act <- hist(use[(use$Category=="Actual"),1], breaks=seq(0,1,0.005))$counts/sum(hist(use[(use$Category=="Actual"),1], breaks=seq(0,1,0.005))$counts)
    sim <- hist(use[(use$Category=="Simulation"),1], breaks=seq(0,1,0.005))$counts/sum(hist(use[(use$Category=="Simulation"),1], breaks=seq(0,1,0.005))$counts)
    KL_value[i,1] <- KL(rbind(sim, act))
  }
  
  ans <- ggplot(KL_value)+
    geom_line(aes(x=Day, y=KL))+facet_grid(Medium~.)+theme_light()+
    xlab("Day")+ylab("Kullback-Leibler divergence<br>between Actual data & Expectation")+
    theme(axis.title=element_text(size=6.5),axis.text=element_text(size=5),
          axis.title.y=element_markdown(),
          panel.grid.minor=element_blank(),legend.position="none",panel.spacing=unit(3,"pt"),
          strip.text=element_text(size=5),strip.text.x=element_text(color="black",margin=margin(0.5,0,0.5,0,unit="pt")),
          panel.background=element_blank(),strip.text.y=element_text(color="black",angle=0,margin=margin(0.5,1,0.5,1,unit="pt")))
  return(list(ans, KL_value))
}

plpl_list2 <- list()
for(i in 1:length(tax)){
  plpl_list2[[i]] <- KL_plot(data=matL[[i]])
}

saveRDS(plpl_list2, "table/KL_calc.rds")
