library(ggplot2); library(vegan);library(gridExtra);library(ggstar); library(ggtext);library(mclust)
library(tidyr); library(stringr);library(ggforce);library(cowplot);library(lemon);library(philentropy)
library(cluster); source("function/functions.R")

#-- information of directory
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/Data_step10_median.rds")
dataN <- readRDS("table/List_data_5000.rds")
fd <- "table"
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- c(2,4,6,8,10,12)

#-- read Cluster_number
data2 <- list()
for(i in 1:5){
  data2[[i]] <- list()
  for(j in 1:length(ti)){
    data2[[i]][[j]] <- list()
    for(k in 1:length(Medium)){
      data2[[i]][[j]][[k]] <- dataN[[i]][(dataN[[i]]$Day==ti[j])&(dataN[[i]]$Medium==Medium[k]),]
      data2[[i]][[j]][[k]] <- data2[[i]][[j]][[k]][,6:ncol(data2[[i]][[j]][[k]])]
    }
  }
}

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
  matL[[i]] <- matL[[i]][(matL[[i]]$Category == "Actual"),]
}

##############################################
ans_mat <- as.data.frame(matrix(NA, ncol=14, nrow=144))
colnames(ans_mat) <- c("LevelofCommunityCompositionalData", "Day", "Medium", "TheNumberofSamples",
                       "logLL1", "BIC1", "logLL2", "BIC2", "df", "TestStatics", "P", "FDR(q)",
                       "ClusterNumber", "AveofSilhouetteCoeffecient")

#-- setting variables
ans_mat$LevelofCommunityCompositionalData <- c(rep("Family",48), rep("Genus",48), rep("ASV",48))
ans_mat$Day <- rep(c(rep(2,8), rep(4,8), rep(6,8), rep(8,8), rep(10,8), rep(12,8)), 3)
ans_mat$Medium <- rep(Medium, 18)

co1 <- 1; co2 <- 2
for(i in 1:nrow(ans_mat)){
  #-- Select use level
  if(ans_mat$LevelofCommunityCompositionalData[i] == "Family"){
    ni=3
  }else if(ans_mat$LevelofCommunityCompositionalData[i] == "Genus"){
    ni=4
  }else{ni=5}
  
  #-- 2 GMMs
  ddd <- matL[[ni]][(matL[[ni]]$Medium == ans_mat$Medium[i])&(matL[[ni]]$Day == ans_mat$Day[i]),]
  MM1 <- summary(Mclust(ddd$ans, G=co1, modelNames = "V"), parameters=T)
  MM2 <- summary(Mclust(ddd$ans, G=co2, modelNames = "V"), parameters=T)
  ans_mat[i,c(4, 5, 6, 7, 8)] <- c(MM1$n, MM1$loglik, MM1$bic,
                                   MM2$loglik, MM2$bic)
  #-- likelihood ratio test
  ans_mat[i,9] <- df12 <- MM2$df - MM1$df
  ans_mat[i,10] <- devbic <- -2*(MM1$loglik-MM2$loglik) #- MM1 is null, MM2 is alternative
  ans_mat[i,11] <- max((1 - pchisq(devbic, df12)), 1e-17)
  
  #-- Silhouette
  nj <- c(1:length(ti))[ti==ans_mat[i,2]]; nk <- c(1:length(Medium))[Medium==ans_mat[i,3]]
  BCd <- vegdist(data2[[ni]][[nj]][[nk]], method="bray"); BC_Si <- c()
  for(kkk in 2:(nrow(data2[[ni]][[nj]][[nk]])-1)){
    BC_Si[kkk-1] <- pam(BCd, k=kkk, diss=T)$silinfo$avg.width
  }
  ans_mat[i, 13] <- c(2:(nrow(data2[[ni]][[nj]][[nk]])-1))[BC_Si==max(BC_Si)]
  ans_mat[i, 14] <- BC_Si[BC_Si==max(BC_Si)]
}

#-- qvalue by FDR
ans_mat[,12] <- p.adjust(ans_mat[,11], method="BH")
write.csv(ans_mat, sprintf("%s/SD_multimodal_12.csv", fd))
###############################################

##############################################
ans_mat <- as.data.frame(matrix(NA, ncol=14, nrow=144))
colnames(ans_mat) <- c("LevelofCommunityCompositionalData", "Day", "Medium", "TheNumberofSamples",
                       "logLL1", "BIC1", "logLL2", "BIC2", "df", "TestStatics", "P", "FDR(q)",
                       "ClusterNumber", "AveofSilhouetteCoeffecient")

#-- setting variables
ans_mat$LevelofCommunityCompositionalData <- c(rep("Family",48), rep("Genus",48), rep("ASV",48))
ans_mat$Day <- rep(c(rep(2,8), rep(4,8), rep(6,8), rep(8,8), rep(10,8), rep(12,8)), 3)
ans_mat$Medium <- rep(Medium, 18)

co1 <- 1; co2 <- 3
for(i in 1:nrow(ans_mat)){
  #-- Select use level
  if(ans_mat$LevelofCommunityCompositionalData[i] == "Family"){
    ni=3
  }else if(ans_mat$LevelofCommunityCompositionalData[i] == "Genus"){
    ni=4
  }else{ni=5}
  
  #-- 2 GMMs
  ddd <- matL[[ni]][(matL[[ni]]$Medium == ans_mat$Medium[i])&(matL[[ni]]$Day == ans_mat$Day[i]),]
  MM1 <- summary(Mclust(ddd$ans, G=co1, modelNames = "V"), parameters=T)
  MM2 <- summary(Mclust(ddd$ans, G=co2, modelNames = "V"), parameters=T)
  ans_mat[i,c(4, 5, 6, 7, 8)] <- c(MM1$n, MM1$loglik, MM1$bic,
                                   MM2$loglik, MM2$bic)
  #-- likelihood ratio test
  ans_mat[i,9] <- df12 <- MM2$df - MM1$df
  ans_mat[i,10] <- devbic <- -2*(MM1$loglik-MM2$loglik) #- MM1 is null, MM2 is alternative
  ans_mat[i,11] <- max((1 - pchisq(devbic, df12)), 1e-17)
  
  #-- Silhouette
  nj <- c(1:length(ti))[ti==ans_mat[i,2]]; nk <- c(1:length(Medium))[Medium==ans_mat[i,3]]
  BCd <- vegdist(data2[[ni]][[nj]][[nk]], method="bray"); BC_Si <- c()
  for(kkk in 2:(nrow(data2[[ni]][[nj]][[nk]])-1)){
    BC_Si[kkk-1] <- pam(BCd, k=kkk, diss=T)$silinfo$avg.width
  }
  ans_mat[i, 13] <- c(2:(nrow(data2[[ni]][[nj]][[nk]])-1))[BC_Si==max(BC_Si)]
  ans_mat[i, 14] <- BC_Si[BC_Si==max(BC_Si)]
}

#-- qvalue by FDR
ans_mat[,12] <- p.adjust(ans_mat[,11], method="BH")
write.csv(ans_mat, sprintf("%s/SD_multimodal_13.csv", fd))
###############################################

##############################################
ans_mat <- as.data.frame(matrix(NA, ncol=14, nrow=144))
colnames(ans_mat) <- c("LevelofCommunityCompositionalData", "Day", "Medium", "TheNumberofSamples",
                       "logLL1", "BIC1", "logLL2", "BIC2", "df", "TestStatics", "P", "FDR(q)",
                       "ClusterNumber", "AveofSilhouetteCoeffecient")

#-- setting variables
ans_mat$LevelofCommunityCompositionalData <- c(rep("Family",48), rep("Genus",48), rep("ASV",48))
ans_mat$Day <- rep(c(rep(2,8), rep(4,8), rep(6,8), rep(8,8), rep(10,8), rep(12,8)), 3)
ans_mat$Medium <- rep(Medium, 18)

co1 <- 1; co2 <- 4
for(i in 1:nrow(ans_mat)){
  #-- Select use level
  if(ans_mat$LevelofCommunityCompositionalData[i] == "Family"){
    ni=3
  }else if(ans_mat$LevelofCommunityCompositionalData[i] == "Genus"){
    ni=4
  }else{ni=5}
  
  #-- 2 GMMs
  ddd <- matL[[ni]][(matL[[ni]]$Medium == ans_mat$Medium[i])&(matL[[ni]]$Day == ans_mat$Day[i]),]
  MM1 <- summary(Mclust(ddd$ans, G=co1, modelNames = "V"), parameters=T)
  MM2 <- summary(Mclust(ddd$ans, G=co2, modelNames = "V"), parameters=T)
  ans_mat[i,c(4, 5, 6, 7, 8)] <- c(MM1$n, MM1$loglik, MM1$bic,
                                   MM2$loglik, MM2$bic)
  #-- likelihood ratio test
  ans_mat[i,9] <- df12 <- MM2$df - MM1$df
  ans_mat[i,10] <- devbic <- -2*(MM1$loglik-MM2$loglik) #- MM1 is null, MM2 is alternative
  ans_mat[i,11] <- max((1 - pchisq(devbic, df12)), 1e-17)
  
  #-- Silhouette
  nj <- c(1:length(ti))[ti==ans_mat[i,2]]; nk <- c(1:length(Medium))[Medium==ans_mat[i,3]]
  BCd <- vegdist(data2[[ni]][[nj]][[nk]], method="bray"); BC_Si <- c()
  for(kkk in 2:(nrow(data2[[ni]][[nj]][[nk]])-1)){
    BC_Si[kkk-1] <- pam(BCd, k=kkk, diss=T)$silinfo$avg.width
  }
  ans_mat[i, 13] <- c(2:(nrow(data2[[ni]][[nj]][[nk]])-1))[BC_Si==max(BC_Si)]
  ans_mat[i, 14] <- BC_Si[BC_Si==max(BC_Si)]
}

#-- qvalue by FDR
ans_mat[,12] <- p.adjust(ans_mat[,11], method="BH")
write.csv(ans_mat, sprintf("%s/SD_multimodal_14.csv", fd))
###############################################