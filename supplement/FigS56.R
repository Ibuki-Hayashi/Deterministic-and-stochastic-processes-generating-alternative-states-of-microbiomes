library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggforce); library(ggstar);library(ggnewscale)
source("functions/functions.R")

exp <- read.csv("table/expdata.csv")
ino_data <- readRDS("table/ListIno_seqtab_10000.rds")
data <- readRDS("table/List_data_rep.rds"); tax <- c("cla","ord","fam","gen","ASV")
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC"); ti <- data[[length(tax)+1]][[2]]

#-- make data ordered, & long data & color palettes
med_long <- function(data,ino,Nexp,top){ ##-- data include exp, med is vector of Medium, ti is vector of Day, top means the number of used Taxonomy
  exp <- data[,(1:ncol(Nexp))]; table <- data[,((ncol(exp)+1):ncol(data))]
  table <- sum_other(mat=table,num=top); data <- cbind(exp,as.data.frame(table))
  #-
  table3 <- table[(exp$Day==2),]; rownames(table3) <- as.numeric(str_sub(rownames(table3),-4,-1))%%384
  kkk <- hclust(as.dist(1-cor(t(table3))),method="average")
  #-
  ino_other <- function(itab,dtab){
    u_itab <- itab[,(colnames(itab) %in% colnames(dtab))]
    ot <- rowSums(itab[,!(colnames(itab) %in% colnames(dtab))])
    ans <- cbind(as.data.frame(u_itab),Others=ot); return(ans)
  }
  
  ino <- ino_other(itab=ino,dtab=table); ino <- cbind(Sample_name=rownames(ino),as.data.frame(ino))
  ino <- merge(Nexp,ino,by="Sample_name")
  
  tab_a <- table[,(colnames(table) != "Unidentified")&(colnames(table) != "Others")]
  key_a <- colSums(tab_a); key_a <- colnames(tab_a)[order(key_a,decreasing=T)]; key <- c(key_a,"Unidentified","Others")
  cpal <- palettes(key_a); cpal <- rev(c(cpal,"#D3D3D3","#2a333c"))
  
  Ld_table <- gather(data,Taxonomy,abundance,-c(1:ncol(exp)))
  Ld_table <- cbind(rep=(as.numeric(str_sub(Ld_table$Sample_name,-4,-1))%%384),Ld_table)
  Ld_table$Taxonomy <- factor(Ld_table$Taxonomy,levels=rev(key))
  Ld_table$rep <- factor(Ld_table$rep,levels=unique(Ld_table$rep)[kkk$order])
  
  Li_table <- gather(ino,Taxonomy,abundance,-c(1:ncol(exp)))
  Li_table$Taxonomy <- factor(Li_table$Taxonomy,levels=key[key %in% Li_table$Taxonomy])
  
  ans <- list(Ld_table,Li_table,key,cpal)
  return(ans)
}
#-- plot with Day
each_time_plot <- function(longlist,med,ti){
  each_med <- list(); each_pal <- list()
  L_d <- longlist[[1]]; L_i <- longlist[[2]]
  key <- longlist[[3]]; c_d <- longlist[[4]]
  for(i in 1:length(med)){
    each_med[[i]] <- L_d[(L_d$Medium == med[i]),]; each_pal[[i]] <- c_d[(key %in% unique(as.character(each_med[[i]]$Taxonomy)))]
  }
  ans <- list()
  for(i in 1:length(med)){
    ans[[i]] <- ggplot(each_med[[i]])+
      geom_bar(aes(x=rep, y=abundance, fill=Taxonomy),position="fill",sta="identity",width=0.9)+
      ggtitle(med[i])+xlab("Replicate communities")+ylab("Relative abundance")+
      theme(text=element_text(size=7),panel.grid.minor=element_blank(),axis.ticks=element_blank(),legend.position="none",
            panel.spacing=unit(1.5,"pt"),plot.margin=unit(c(1,2.5,1,2.5),"pt"),
            panel.grid.major=element_blank(),panel.background=element_blank(), axis.text.x=element_blank(),strip.text.x=element_text(margin=margin(0,-0.8,0,-0.8,unit="pt")))+
      facet_wrap(~Day,scales = "free_x",nrow=3)+scale_y_continuous(expand=c(0,0))+scale_x_discrete(expand=c(0,0))+
      scale_fill_manual(values=each_pal[[i]])
  }
  return(ans)
}

#-- legend
legend_plot <- function(longlist,med,ti){
  pl <- ggplot(longlist[[1]])+
    geom_bar(aes(x=Sample_name, y=abundance, fill=Taxonomy),position="fill",color="black",sta="identity",width=0.95)+
    theme(text=element_text(size=5),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(), axis.text.x=element_blank())+
    scale_fill_manual(values=longlist[[4]])+
    theme(legend.key.size = unit(5,'pt'),legend.text=element_text(face="italic",size=5),legend.margin = margin(0,0,0,0))+guides(fill=guide_legend(ncol=1,))
  ans <- g_legend(pl); return(ans)
}

#-- output
ans <- list()
t_plist <- list(); l_plist <- list()

ts <- c(10,15,20,25,30) #(cla,ord,fam,gen,asv)のtop

for(i in 3:length(tax)){
  ans[[i]] <- med_long(data=data[[i]],ino=ino_data[[i]],Nexp=exp,top=ts[i])
  t_plist[[i]] <- each_time_plot(longlist=ans[[i]],med=Medium,ti=ti)
  l_plist[[i]] <- legend_plot(longlist=ans[[i]],med=Medium,ti=ti)
}

Fab <- grid.arrange(t_plist[[3]][[1]],t_plist[[3]][[2]],t_plist[[3]][[3]],t_plist[[3]][[4]],
                    t_plist[[3]][[5]],t_plist[[3]][[6]],t_plist[[3]][[7]],t_plist[[3]][[8]],ncol=2,nrow=4)
ASb <- grid.arrange(t_plist[[5]][[1]],t_plist[[5]][[2]],t_plist[[5]][[3]],t_plist[[5]][[4]],
                    t_plist[[5]][[5]],t_plist[[5]][[6]],t_plist[[5]][[7]],t_plist[[5]][[8]],ncol=2,nrow=4)

Fl <- plot_grid(l_plist[[3]],NA,ncol=1,rel_heights=c(0.3,1))
Al <- plot_grid(l_plist[[5]],NA,ncol=1,rel_heights=c(0.5,1))

##### output #####
pdf("FigS5.pdf",height=10.8,width=7.2)
plot_grid(ASb,Al,ncol=2,rel_widths=c(1,0.3))
dev.off()

pdf("FigS6.pdf",height=10.8,width=7.2)
plot_grid(Fab,Fl,ncol=2,rel_widths=c(1,0.2))
dev.off()
