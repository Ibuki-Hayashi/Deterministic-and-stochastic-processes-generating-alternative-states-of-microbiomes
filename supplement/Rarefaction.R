library(ggplot2); library(vegan);library(gridExtra)
library(cowplot); library(markdown); library(ggtext)
library(tidyr); library(stringr);library(ggforce)
library(ggstar);library(ggnewscale);library(ggsci)
library(tidyverse);library(parallel);library(phyloseq)
source("functions/functions.R")

#-- Load data
exp <- read.csv("table/expdata.csv")
seqtab <- readRDS("table/merge_seqtab.rds")
taxtable <- readRDS("table/merge_taxonomylist.rds")

#-- delete std
use_std <- seqtab[,(colnames(seqtab) %in% rownames(taxtable[(str_sub(taxtable[,2],1,3) != "STD"),]))]
use_std <- use_std[(rowSums(use_std) > 4999),]; use_std <- use_std[(str_sub(rownames(use_std),1,4) == "HS2_"),]
seq <- cbind(Sample_name=rownames(use_std),as.data.frame(use_std))
data <- merge(exp,seq,by="Sample_name")

phy_s <- data[,(ncol(exp)+1):ncol(data)]; rownames(phy_s) <- data$Sample_name
phy_t <- taxtable[!(taxtable[,3]=="c_"),]
phy_e <- data[,2:ncol(exp)]; rownames(phy_e) <- data$Sample_name

ps_all <- phyloseq(otu_table(as.matrix(phy_s), taxa_are_rows = F),
                   sample_data(phy_e),
                   tax_table(as.matrix(phy_t)))

phys <- list()
for(i in 1:length(unique(data$Day))){
  data_d <- data[(data$Day==unique(data$Day)[i]),]
  phy_s <- data_d[,(ncol(exp)+1):ncol(data)]; rownames(phy_s) <- data_d$Sample_name
  phy_t <- taxtable[!(taxtable[,3]=="c_"),]
  phy_e <- data_d[,2:ncol(exp)]; rownames(phy_e) <- data_d$Sample_name
  
  ps_all <- phyloseq(otu_table(as.matrix(phy_s), taxa_are_rows = F),
                     sample_data(phy_e),
                     tax_table(as.matrix(phy_t)))
  phys[[i]] <- ps_all
}

#-- rarefy with phyloseq
gg_rrr1 <- function(ps_all){
  ps_rr <- rarecurve(as.matrix(otu_table(ps_all)), step = 50, label = TRUE)
  hoge <- sample_data(ps_all)
  rare <- lapply(ps_rr, function(x){
    b <- as.data.frame(x)
    b <- data.frame(ASV = b[,1], raw_read = rownames(b))
    b$raw_read <- as.numeric(gsub("N", "",  b$raw_read))
    return(b)
  })
  # Label list
  names(rare) <- rownames(sample_data(ps_all))
  # Convert to the data frame
  rare_df <- map_dfr(rare, function(x) return(data.frame(x)), .id = "sample")
  meds <- c()
  for(i in 1:length(unique(rare_df$sample))){
    m <- hoge[[4]][rownames(hoge)==unique(rare_df$sample)[i]]
    mm <- rep(as.character(m),nrow(rare_df[rare_df$sample==unique(rare_df$sample)[i],]))
    meds <- c(meds,mm)
  }
  rare_dff <- cbind(rare_df,Medium=meds)
  
  return(rare_dff)
}
rare_dfs <- lapply(phys,gg_rrr1)
rare_dfs2 <- list()
for(i in 1:length(rare_dfs)){
  rare_dfs2[[i]] <- list(rare_dfs[[i]],sprintf("Day %s",as.character(unique(data$Day)[i])))
}

gg_rrr22 <- function(rare_dff){
  set.seed(1234)
  rar <- sample(unique(rare_dff[[1]]$sample),size=50)
  rar <- rare_dff[[1]][(rare_dff[[1]]$sample %in% rar),]
  ans <- ggplot(rar, aes(x = raw_read, y = ASV, color=sample))+
    geom_line()+
    xlab("Reads")+ylab("The number of ASV")+theme_bw()+
    ggtitle(rare_dff[[2]])+scale_color_igv()+
    theme(text=element_text(size=7),legend.position="none")
  return(ans)
}
rare_ggs2 <- lapply(rare_dfs2,gg_rrr22)

ansplot2 <- plot_grid(rare_ggs2[[1]],rare_ggs2[[2]],
                     rare_ggs2[[3]],rare_ggs2[[4]],
                     rare_ggs2[[5]],rare_ggs2[[6]],nrow=3,ncol=2)

# -- output
pdf("FigS2.pdf",width=7.2,height=10.5)
ansplot2
dev.off()
