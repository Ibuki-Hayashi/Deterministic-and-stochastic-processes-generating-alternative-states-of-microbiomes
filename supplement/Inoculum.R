library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggnewscale); library(reshape2)
source("functions/functions.R")

#-- read data from previous dir
exp <- read.csv("table/expdata.csv",header=T)
seqtab <- readRDS("table/merge_seqtab.rds")
taxa.print <- readRDS("table/merge_taxonomylist.rds")
ino_data <- readRDS("table/ListIno_seqtab_10000.rds")
data <- readRDS("table/List_data_rep.rds")

#########################################
########## -- For Fig.S1AB -- ###########
#########################################
#-- extract only inoculum rows
seq <- cbind(Sample_name=rownames(seqtab),as.data.frame(seqtab))
k <- seq[(str_sub(seq$Sample_name,1,4) != "HS2_"),]
data <- merge(exp,k,by="Sample_name")
seqtab <- data[(rowSums(data[,(6:ncol(data))])>9999),6:ncol(data)]

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Nnumber of copies of standard DNA
std.bairitu <- 500 ## important!!!!!!!!!!!!!!!!!!!!!!!!!!
std.copy.n <- c( 0.1, 0.05, 0.02, 0.01,0.005)*(6.02*10^14/1000000)/std.bairitu

############################################################################
####
#### F1. Collection of helper functions for DNA extration study
#### 2017.12.1 Ushio
####

# ggplot function 1
PlotStyle <-  function(ggobject){
  return(ggobject + theme_light() + theme(axis.text.x = element_text(angle=0),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.text = element_text(size=7),
                                          axis.title = element_text(size=7),
                                          panel.background=element_rect(colour="black", fill=NA, size=0.8)))
}

# ggplot function 2
PlotStyle2 <- function(ggobject){
  return(ggobject + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                          axis.title.x = element_blank(),
                          legend.position = "none") +
           geom_jitter(shape = 16, size = 2, alpha = 0.8, width = 0.1, height = 0))
}

# Merge standard DNA sequences
MergeSTD <- function(std.i, std.data = std.table){
  index.std <- which(match(colnames(std.table), std.i) == 1)
  if(length(index.std) > 1){
    std.tmp <- rowSums(std.table[,index.std])
  }else{
    std.tmp <- std.table[,index.std]
  }
  return(std.tmp)
}
############################################################################

# -- Extract standard sequeces
detected.std.name <- unique(taxa.print[which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro"), "Phylum"])

n.std.seq <- which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro")
std.table <- seqtab[,n.std.seq]
std.taxa <- taxa.print[n.std.seq, "Phylum"]

# --  STD reads - copy number relationship
# --  Rename colnames
colnames(std.table) <- std.taxa
# --  Merge the same colnames
new.std.table <- data.frame(std_rank1 = MergeSTD(detected.std.name[1], std.data = std.table),
                            std_rank2 = MergeSTD(detected.std.name[2], std.data = std.table),
                            std_rank3 = MergeSTD(detected.std.name[3], std.data = std.table),
                            std_rank4 = MergeSTD(detected.std.name[4], std.data = std.table),
                            std_rank5 = MergeSTD(detected.std.name[5], std.data = std.table))

# -- Interpolate 0 into missing value
new.std.table[is.na(new.std.table)] <- 0

# --  Linear regression
adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared
lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1]
r2.summary <- apply(new.std.table, 1, adj.r.fun)
coef.summary <- apply(new.std.table, 1, lm.coef.fun)
new.seqtab <- as.data.frame(seqtab[,-n.std.seq]) # Make seq table without standard DNA
new.seqtab2 <- new.seqtab

hosei.seq <- new.seqtab
for(i in 1:ncol(hosei.seq)){
  hosei.seq[,i] <- hosei.seq[,i]/coef.summary
}
saveRDS(hosei.seq,file="table/seqtab_copies.rds")

colnames(new.seqtab2) <- taxa.print[-n.std.seq, "Phylum"]

# --  Visualize regression results
# --  1. R2 value distribution
g1 <- ggplot(data.frame(values = r2.summary), aes(x = values))
g1 <- g1 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2) #+ xlim(0.9,1) + ylim(0,13)
g1 <- g1 + xlab(expression(paste("R"^{2}, " values"))) + ylab("Count") #+ scale_x_continuous(breaks=seq(0.8,1,by=0.05))
g1 <- PlotStyle(g1)
FigA <- g1

# --  3. Regression examples
max.slope <- as.numeric(c(new.std.table[which.max(coef.summary),], 0))
med.slope <- as.numeric(c(new.std.table[which.min(abs(coef.summary - median(coef.summary))),], 0))
min.slope <- as.numeric(c(new.std.table[which.min(coef.summary),], 0))
slope.summary <- melt(data.frame(copy =c(std.copy.n, 0),
                                 max_slope = max.slope,
                                 med_slope = med.slope,
                                 min_slope = min.slope), id.vars = "copy")
g3 <- ggplot(slope.summary, aes(x = copy, y = value, group = variable, colour = variable))
g3 <- g3 + geom_point(size = 1) + scale_color_manual(name = "Regression slope", values = c("red3", "darkred", "black"))
g3 <- g3 + geom_smooth(method = "lm", size = 0.5, se = F)+theme_bw(base_size=7)
FigB <- g3

#########################################
########### -- For Fig.S1C -- ###########
#########################################
# -- data
data_copy <- readRDS("table/seqtab_copies.rds")
me <- mean(rowSums(data_copy)*50/50*1*10) #copies/1ul(temp)*(powersoil out)/(powersoil in)*10ul(inoculum)
copies <- rowSums(data_copy)*50/50*1*10
t <- rep("A",length(copies))
ddd <- data.frame(t,copies)

pl <- ggplot(ddd)+
  theme_light()+geom_boxplot(aes(x=t,y=copies))+
  theme(text=element_text(size=7),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())+
  ylab("DNA copies in 10 μl initial community (copies)")+xlab("Inoculum")+theme(axis.text.x=element_blank())
FigC <- pl

#########################################
######### -- For Fig.S1DEF -- ###########
#########################################
#-- make data ordered, & long data & color palettes
tax <- c("cla","ord","fam","gen","ASV")
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC"); ti <- data[[length(tax)+1]][[2]]
med_long <- function(data,ino,Nexp,top){ ##-- data include exp, med is vector of Medium, ti is vector of time, top means the number of used taxonomy
  exp <- data[,(1:ncol(Nexp))]; table <- data[,((ncol(exp)+1):ncol(data))]
  table <- sum_other(mat=table,num=top); data <- cbind(exp,as.data.frame(table))
  
  ino_other <- function(itab,dtab){
    u_itab <- itab[,(colnames(itab) %in% colnames(dtab))]
    ot <- rowSums(itab[,!(colnames(itab) %in% colnames(dtab))])
    ans <- cbind(as.data.frame(u_itab),Others=ot); return(ans)
  }
  
  ino <- ino_other(itab=ino,dtab=table); ino <- cbind(Sample_name=rownames(ino),as.data.frame(ino))
  ino <- merge(Nexp,ino,by="Sample_name"); ino_d <- ino[1,]
  ino_d[1,((ncol(exp)+1):ncol(ino))] <- colSums(ino[,((ncol(exp)+1):ncol(ino))])
  ino <- ino_d
  
  tab_a <- table[,(colnames(table) != "Unidentified")&(colnames(table) != "Others")]
  key_a <- colSums(tab_a); key_a <- colnames(tab_a)[order(key_a,decreasing=T)]; key <- c(key_a,"Unidentified","Others")
  cpal <- palettes(key_a); cpal <- c(cpal,"#D3D3D3","#2a333c")
  
  Ld_table <- gather(data,taxonomy,abundance,-c(1:ncol(exp)))
  Ld_table <- cbind(rep=(as.numeric(str_sub(Ld_table$Sample_name,-4,-1))%%384),Ld_table)
  Ld_table$taxonomy <- factor(Ld_table$taxonomy,levels=key)
  
  Li_table <- gather(ino,taxonomy,abundance,-c(1:ncol(exp)))
  Li_table$taxonomy <- factor(Li_table$taxonomy,levels=key[key %in% Li_table$taxonomy])
  
  ans <- list(Ld_table,Li_table,key,cpal)
  return(ans)
}
######
PlotStyle <-  function(ggobject){
  return(ggobject + theme_light() + theme(axis.text.x = element_text(angle=0),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.text = element_text(size=12),
                                          axis.title = element_text(size=12)))
}
##########

#-- Inoplot
Ino_plot <- function(longlist){
  L_d <- longlist[[1]]; L_i <- longlist[[2]]
  key <- longlist[[3]]; c_d <- longlist[[4]]
  
  ipal <- c_d[(key %in% unique(as.character(L_i$taxonomy)))]
  
  ans <- ggplot(L_i)+
    geom_bar(aes(x=Sample_name, y=abundance, fill=taxonomy),position="fill",color="black",sta="identity",width=0.95)+
    xlab("sample")+ylab("relative abundance")+scale_fill_manual(values=ipal)+theme(legend.position="none")+
    theme(text=element_text(size=7),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(), axis.text.x=element_blank())+
    scale_y_continuous(expand=c(0,0))+scale_x_discrete(expand=c(0,0))+xlab("sample")
  
  return(ans)
}
#-- legend
legend_plot <- function(longlist,med,ti){
  pl <- ggplot(longlist[[1]])+
    geom_bar(aes(x=Sample_name, y=abundance, fill=taxonomy),position="fill",color="black",sta="identity",width=0.95)+
    theme(text=element_text(size=6),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(), axis.text.x=element_blank())+
    scale_fill_manual(values=longlist[[4]])+
    theme(legend.key.size = unit(6,'pt'),legend.text=element_text(face="italic",size=6),legend.margin = margin(0,0,0,0))+guides(fill=guide_legend(ncol=1,))
  ans <- g_legend(pl); return(ans)
}
#-- output
ans <- list()
l_plist <- list(); i_plist <- list()
ts <- c(10,15,20,25,30) #top of (cla,ord,fam,gen,asv)

for(i in 3:length(tax)){
  ans[[i]] <- med_long(data=data[[i]],ino=ino_data[[i]],Nexp=exp,top=ts[i])
  i_plist[[i]] <- Ino_plot(longlist=ans[[i]])
  l_plist[[i]] <- legend_plot(longlist=ans[[i]],med=Medium,ti=ti)
}
s_l <- list()
s_l[[5]] <- plot_grid(i_plist[[5]],l_plist[[5]],rel_widths=c(1,0.6))
s_l[[4]] <- plot_grid(i_plist[[4]],l_plist[[4]],rel_widths=c(1,2.2))
s_l[[3]] <- plot_grid(i_plist[[3]],l_plist[[3]],rel_widths=c(1,1.4))

FigDEF <- plot_grid(s_l[[5]],s_l[[4]],s_l[[3]],nrow=1,labels=c("D","E","F"),rel_widths=c(0.5,1,0.75))

#########################################
############ -- For output-- ############
#########################################
#-- plot_grid
ABC <- plot_grid(FigA,FigB,C,nrow=1,labels=c("A","B","C"),rel_widths=c(0.8,1,0.4))
ans <- plot_grid(ABC,FigDEF,nrow=2)

pdf("FigS1.pdf",width=7.8,height=6.8)
ans
dev.off()
