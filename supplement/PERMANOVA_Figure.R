library(ggplot2); library(vegan);library(gridExtra);library(ggstar);library(ggtext)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon)

#-- information of directory
data <- read.csv("table/SupplementaryDataS3.csv")
data <- data[(data$Day!="all"),]
data[,3] <- as.numeric(str_sub(data[,3],start=4))
data[,5] <- as.character(data[,5])

lev <- c("family","genus","ASV","PICRUSTs2")
lev2 <- paste(c("Family","Genus","ASV","Functional compositions by PICRUSTs2"),c(rep("level communities",3),""),sep=" ")

d_lb <- list()
for(i in 1:length(lev)){
  d_lb[[i]] <- data[(data[,2]==lev[i])&(data[,4]=="bray"),]
}

#-- plot
plot_notall <- function(dd,ti){
  dd[,5][dd[,5]=="glucose"] <- "Concentration<br>of glucose"
  dd[,5][dd[,5]=="leucine"] <- "Leucine"
  dd[,5][dd[,5]=="citrate"] <- "Citrate"
  FDRq <- as.character(dd$qs<0.05); FDRq <- factor(FDRq,levels=c("TRUE","FALSE"))
  ans <- ggplot(dd,aes(x=Day,y=R2))+
    geom_point(size=1.2, aes(shape=FDRq))+geom_line(aes(group=Treatment,color=Treatment),size=0.4)+
    scale_color_manual(values = palettes(c(as.factor(c(dd$Treatment)))))+
    theme_light()+theme(text=element_text(size=7))+
    scale_x_continuous(breaks=seq(0,12,2))
  ans <- ans + theme(legend.position="none")+
    ylab(paste("<i>R^2</i>","of PERMANOVA<br>for the communities"))+
    theme(axis.title.y=element_markdown(),legend.text=element_markdown())+
    ggtitle(ti)
  return(ans)
}
plot_leg <- function(dd){
  dd[,5][dd[,5]=="glucose"] <- "Concentration<br>of glucose"
  dd[,5][dd[,5]=="leucine"] <- "Leucine"
  dd[,5][dd[,5]=="citrate"] <- "Citrate"
  FDRq <- as.character(dd$qs<0.05); FDRq <- factor(FDRq,levels=c("TRUE","FALSE"))
  ans <- ggplot(dd,aes(x=Day,y=R2))+
    geom_point(size=1.2, aes(shape=FDRq))+geom_line(aes(group=Treatment,color=Treatment),size=0.4)+
    scale_color_manual(values = palettes(c(as.factor(c(dd$Treatment)))))+
    theme_light()+theme(text=element_text(size=7))+
    scale_x_continuous(breaks=seq(0,12,2))
  ans <- ans+
    ylab(paste("<i>R^2</i>","of PERMANOVA<br>for the communities"))+
    theme(axis.title.y=element_markdown(),legend.text=element_markdown())
  ans <- g_legend(ans)
  return(ans)
}

plots_na <- list()
for(i in 1:length(lev)){
  plots_na[[i]] <- plot_notall(dd=d_lb[[i]],ti=lev2[i])
}
plotsL <- plot_leg(dd=d_lb[[1]])

AB <- plot_grid(plots_na[[1]],plots_na[[2]],labels=c("A","B"),nrow=1)
CD <- plot_grid(plots_na[[3]],plots_na[[4]],labels=c("C","D"),nrow=1)
ABCD <- plot_grid(AB,CD,nrow=2)
ansB <- plot_grid(ABCD,plotsL,rel_widths = c(1,0.3))

#-- output
pdf("FigS7.pdf",width=7.2,height=6)
ansB
dev.off()
