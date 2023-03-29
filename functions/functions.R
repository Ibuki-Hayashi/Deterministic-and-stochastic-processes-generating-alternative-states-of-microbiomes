taxaratio <- function(seqtab,taxa){
  f <- function(seqtab,taxa,taxid){
    taxa_n <- cbind(ASV=rownames(taxa),ta=taxa[,(colnames(taxa)==taxid)])
    ls <- unique(taxa_n[,2]); ans <- matrix(0,nrow=nrow(seqtab),ncol=length(ls))
    for(i in 1:length(ls)){
      key <- taxa_n[(taxa_n[,2] == ls[i]),1]
      if(length(key)>1){
        ans[,i] <- rowSums(seqtab[,(colnames(seqtab) %in% key)])
      }else{
        ans[,i] <- seqtab[,(colnames(seqtab) %in% key)]
      }
    }
    rownames(ans) <- rownames(seqtab); colnames(ans) <- ls
    return(ans)
  }
  cla <- f(seqtab=seqtab,taxa=taxa,taxid="Class")
  ord <- f(seqtab=seqtab,taxa=taxa,taxid="Order")
  fam <- f(seqtab=seqtab,taxa=taxa,taxid="Family")
  gen <- f(seqtab=seqtab,taxa=taxa,taxid="Genus")
  
  return(list(cla,ord,fam,gen))
}

sum_other <- function(mat,num){
  a <- colSums(mat)
  b <- sort(a,T)
  c <- (a>b[num])|(a==b[num])
  out <- matrix(0, nrow=nrow(mat), ncol=(sum(c)+1))
  count <- 1
  outname <- c(rep(NA,sum(c)),"Others")
  for(i in 1:ncol(mat)){
    if(c[i]==T){
      out[,count] <- mat[,i]
      outname[count] <- colnames(mat)[i]
      count <- count+1
    }
    else{
      out[,num+1] <- out[,num+1]+mat[,i]
    }
  }
  rownames(out) <- rownames(mat)
  colnames(out) <- outname
  return(out)
}

rare_merge <- function(mat,csv,rarefuction){
  a <- rrarefy(mat, rarefuction)
  b <- cbind(rownames(a),as.data.frame(a))
  colnames(b) <- c("sample_name",colnames(b))
  c <- merge(b,csv,by="sample_name")
  
  return(c)
}

###############
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib) 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

palettes <- function(x){ col.vector[1:length(unique(x))]}


###############
remove.outliers <- function(x, conf.level = 0.95)
{
  x <- x[!is.na(x)]
  del.val <- NULL
  
  while (TRUE) {
    n <- length(x)
    if (n < 3) {
      break
    }
    
    r <- range(x)
    t <- abs(r - mean(x)) / sd(x)
    q2 <- (n - 2) / ((n - 1) ^ 2 / t ^ 2 / n - 1)
    q2[q2 < 0] <- 0
    q <- sqrt(q2)
    p <- n * pt(q, n - 2, lower.tail = FALSE)
    
    if (t[1] < t[2]) {
      if (p[2] < 1 - conf.level) {
        del.val <- c(del.val, r[2])
        x <- x[x != r[2]]
        next
      }
    } else {
      if (p[1] < 1 - conf.level) {
        del.val <- c(del.val, r[1])
        x <- x[x != r[1]]
        next
      }
    }
    break
  }
  return(list(x = x, del.val = del.val))
}

forjac <- function(mat,thr=1,dataframe=F){
  ans <- mat
  for(i in 1:nrow(ans)){
    for(j in 1:ncol(ans)){
      if(ans[i,j]>(thr-1)){
        ans[i,j] <- 1
      }else{
        ans[i,j] <- 0
      }
    }
  }
  if(dataframe==T){
    ans <- as.data.frame(ans)
  }
  return(ans)
}
