library(ggplot2); library(vegan);library(gridExtra);library(ggstar); library(ggtext)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon);library(philentropy)
source("functions/functions.R")

#-- read data from previous dir
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/Data_step10_median.rds")
tax <- c("cla","ord","Family","Genus","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
KL_d <- readRDS("table/KL_calc.rds")
ti <- c(2,4,6,8,10,12)

#-- hist plot with Medium color
hi_KL <- function(values, title, me=Medium){
  values$Medium <- factor(values$Medium, levels=me)
  cocol <- palettes(me)
  ans <- ggplot(data=values, aes(x=KL, fill=Medium))+
    geom_histogram(alpha=0.8, binwidth=0.1)+theme_light()+
    scale_fill_manual(values=cocol)+xlim(c(0,4))+
    xlab("Kullback-Leibler divergence<br>between Actual data & Expectation")+
    ylab("Frequency")+ggtitle(title)+
    theme(axis.title=element_text(size=6.5), axis.text=element_text(size=5),
        axis.title.x=element_markdown(), legend.position="none")
  
  return(ans)
}
hi_KL_leg <- function(values, title, me=Medium){
  values$Medium <- factor(values$Medium, levels=me)
  cocol <- palettes(me)
  ans <- ggplot(data=values, aes(x=KL, fill=Medium))+
    geom_histogram(alpha=0.8, binwidth=0.1)+theme_light()+
    scale_fill_manual(values=cocol)+
    xlab("Kullback-Leibler divergence<br>between Actual data & Expectation")+
    ylab("Frequency")+ggtitle(title)+
    theme(axis.title=element_text(size=6.5), axis.text=element_text(size=5),
          axis.title.x=element_markdown())
  ans <- g_legend(ans)
  return(ans)
}

KL_ll <- list()
for(i in 1:length(tax)){
  KL_ll[[i]] <- hi_KL(values=KL_d[[i]][[2]], title=sprintf("KL divergence(%s)",tax[i]))
}

KL_leg <- hi_KL_leg(values = KLH_all, title = sprintf("KL divergence(All)"))

#-- ggplot2 object of KL_histogram
ans <- plot_grid(KL_ll[[3]], KL_ll[[4]], KL_leg,
                 KL_ll[[5]], NA, NA, ncol=3, rel_widths = c(1,1,0.4))
