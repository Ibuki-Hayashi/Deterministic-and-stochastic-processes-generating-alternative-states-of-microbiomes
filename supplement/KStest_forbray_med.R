########################################################################
###-  Experiment whether Bray-curtis follows a normal distribution  -###
########################################################################

library(ggplot2); library(vegan);library(gridExtra); library(cowplot);library(lemon)
library(tidyr); library(stringr);library(ggnewscale); library(reshape2);library(ggstar)
library(ggtext);library(ggforce)
source("functions/functions.R")

#-- read data and set directory
exp <- read.csv("table/expdata.csv")
L_tab <- readRDS("table/List_data_5000.rds")
tax <- c("cla","ord","fam","gen","ASV"); Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")
ti <- L_tab[[length(tax)+1]]

####### Functions ################################################
## -- Function of Making simulated data set("N" means the number of trial, "size" means the community size of each community)
multi_rep <- function(N, size, prob){
  com <- rmultinom(n=N, size=size, prob=prob)
  return(com)
}

## -- Random walk simulation -- ##
randW <- function(data, Step=1, Read=5000){
  RRR <- function(vec, step=Step, read=Read){
    iN <- floor(runif(1, min=1, max=read)); ouT <- floor(runif(1, min=1, max=read+step-1))
    iN2 <- 0; ouT2 <- 0
    for(i in 1:length(vec)){
      iN2 <- iN2+vec[i]
      if(iN2 > iN){
        vec[i] <- vec[i]+step; break
      }
    }
    for(i in 1:length(vec)){
      ouT2 <- ouT2+vec[i]
      if(ouT2 > ouT){
        if((vec[i]>step)|(vec[i]==step)){
          vec[i] <- vec[i]-step; break
        }else{
          vec[i] <- 0
          ran <- floor(runif(1, min=1, max=length(vec)+1))
          vec[ran] <- vec[ran]+step-vec[i]; break
        }
      }
    }
    return(vec)
  }
  ans <- apply(data, 2, RRR)
  
  return(ans)
}

## -- Function of Calculating pairwise beta diversity
Cal_pBC <- function(N, com){
  #- Calculate bray
  ans <- c(); use1 <- c(); use2 <- c()
  for(i in 1:(N-1)){
    com1 <- com[,i] #- comm1
    A <- c(); B <- c()
    for(j in 1:(N-i)){
      com2 <- com[,(i+j)] #- comm2
      A[j] <- vegdist(rbind(com1,com2), method="bray")
      B[j] <- i+j
    }
    ans <- c(ans, A)
    use1 <- c(use1, rep(i, length(B))); use2 <- c(use2, B)
  }
  AAA2 <- cbind(ans, cbind(use1, use2))
  return(AAA2)
}
#######################

####### Use data #############
## -- list; Day(ti) & Medium(med)
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
L_dlist <- list()
for(i in 1:length(tax)){
  L_dlist[[i]] <- sp(data=L_tab[[i]],ti=ti,med=Medium)
}
#######################

####### Actual calculation #########
Act <- list()
for(i in 1:length(tax)){
  Act[[i]] <- list()
  for(j in 1:length(ti)){
    Act[[i]][[j]] <- list()
    for(k in 1:length(Medium)){
      A_data <- t(L_dlist[[i]][[j]][[k]][,6:ncol(L_dlist[[i]][[j]][[k]])])
      A_BC <- Cal_pBC(N=ncol(A_data), com=A_data)
      Act[[i]][[j]][[k]] <- list(median(A_BC[,1]), A_BC)
    }
  }
}

####### Simulated calculation & KStest #########
Sit <- list(); KSt <- list()
for(i in 1:5){ #-- Taxonomy
  Sit[[i]] <- list(); KSt[[i]] <- list()
  for(j in 1:6){ #-- Time
    Sit[[i]][[j]] <- list(); KSt[[i]][[j]] <- list()
    for(k in 1:8){ #-- Medium
      start <- Sys.time()
      A_data <- t(L_dlist[[i]][[j]][[k]][,6:ncol(L_dlist[[i]][[j]][[k]])])
      ave <- rowSums(A_data)/(ncol(A_data)*5000)
      S_data <- multi_rep(N=ncol(A_data), size=5000, prob=ave)
      key <- median(Cal_pBC(N=ncol(S_data), com=S_data)[,1])
      thr <- Act[[i]][[j]][[k]][[1]]; A_BC <- Act[[i]][[j]][[k]][[2]]
      
      while(key < thr){
        S_data <- randW(data=S_data, Step=10, Read=5000)
        key <- median(Cal_pBC(N=ncol(S_data), com=S_data)[,1])
      }
      
      S_BC <- Cal_pBC(N=ncol(S_data), com=S_data)

      Sit[[i]][[j]][[k]] <- list(median(S_BC[,1]), S_BC)
      KSt[[i]][[j]][[k]] <- ks.test(A_BC[,1], S_BC[,1])
      
      #-- System message
      end <- Sys.time()
      message(sprintf("Process (i,j,k)=(%s/5,%s/6,%s/8) is done: Exe. time = %s min.",
                      as.character(i), as.character(j), as.character(k), end-start))
    }
  }
}

####### data merge ############
DDD <- list()
for(i in 1:5){
  DDD[[i]] <- list()
  for(j in 1:6){
    DDD[[i]][[j]] <- list()
    for(k in 1:8){
      Value <- rbind(Act[[i]][[j]][[k]][[2]], Sit[[i]][[j]][[k]][[2]])
      Category <- c(rep("Actual", nrow(Act[[i]][[j]][[k]][[2]])), rep("Simulation", nrow(Sit[[i]][[j]][[k]][[2]])))
      
      
      DDD[[i]][[j]][[k]] <- cbind(as.data.frame(Value), Category)
    }
  }
}

saveRDS(DDD, file="table/Data_step10_median.rds")
##############################
