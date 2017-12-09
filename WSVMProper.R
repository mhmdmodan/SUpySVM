library(tidyverse)
library(ggplot2)
library(ggthemes)
library(fpCompare)

dot <- function(x,y) {return((x %*% y)[[1]])}
`%.%` <- function(x,y) {return((x %*% y)[[1]])}

checkSameDim <- function(n, P) {
  if(!is.list(P)) {stop('P is not a list')}
  dimension <- length(n)
  for(pt in P) {
    if(length(pt) != dimension) {return(FALSE)}
  }
  return(TRUE)
}

getMaxIndex <- function(n, P) {
  maxInd <- 1
  for(i in 1:length(P)) {
    curVal <- dot(n, P[[i]])
    maxVal <- dot(n, P[[maxInd]])
    maxInd <- ifelse(curVal > maxVal, i, maxInd)
  }
  return(maxInd)
}

#Finding WRCH vertex
findVertex <- function(P, s, n, u) {
  if(!checkSameDim(n,P)) {stop('P and n are not same dims')}
  if(length(P) != length(s)) {stop('P and s are not same length')}
  
  a <- vector(mode = 'numeric', length = length(P))
  total <- 0
  while(total %!=% 1) {
    i <- getMaxIndex(n, P)
    a[i] <- min(s[i]*u, 1-total)
    total <- total + a[i]
  }
  
  vSum <- 0
  for(i in 1:length(P)) {
    vSum <- vSum + a[i]*P[[i]]
  }
  
  return(vSum)
}

#WSVM Algorithm
WSVM <- function(P, s, y, rchFactor, ep) {
  
  classPos <- which(y > 0)
  classNeg <- which(y < 0)
  
  pos <- P[classPos]
  neg <- P[classNeg]
  sPos <- s[classPos]
  sNeg <- s[classNeg]
  
  randVec1 <- genRand()
  randVec2 <- genRand()
  pPos <- findVertex(pos, sPos, randVec1, rchFactor)
  pNeg <- findVertex(neg, sNeg, randVec2, rchFactor)
  
  while(TRUE) {
    w <- pPos - pNeg
    
    vPos <- findVertex(pos, sPos, -w, rchFactor)
    vNeg <- findVertex(neg, sNeg, w, rchFactor)
    
    if((1 - dot(w, (pPos - vNeg))/dot(w,w)) < ep & (1 - dot(w, (vPos - pNeg))/dot(w,w)) < ep) {
      break
    }
    
    if(w %.% (pPos - vPos) > w %.% (vNeg - pNeg)) {
      
      q <- ((pPos - pNeg) %.% (pPos - vPos))/dot(pPos - vPos,pPos - vPos)
      pPos <- (1-q)*pPos + q*vPos
      
    } else {
      
      q <- -((pPos - pNeg) %.% (pNeg - vNeg))/dot(pNeg - vNeg,pNeg - vNeg)
      pNeg <- (1-q)*pNeg + q*vNeg
      
    }
  }
  
  return(list(w=w, bisect=pNeg+w/2, pPos=pPos, pNeg=pNeg))
}

#####################
set.seed(1345)
pts <- tibble(x = sample(1:100, size=50, replace=TRUE), y = sample(1:100, size=50, replace=TRUE), class = 1)

theorLine <- function(x) {-25/12.5*(x-50)+37.5}

for(i in 1:nrow(pts)) {
  if(pts$y[i] < theorLine(pts$x[i])) {pts[i,3] <- -1}
}

ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() + 
  stat_function(fun=theorLine) + xlim(0,100) + ylim(0,100)

ptsList <- lapply(1:nrow(pts), function(n) as.double(c(pts[n,1],pts[n,2])))
yList <- pts$class
weights <- sapply(ptsList, function(n) {1})
###########

out <- WSVM(ptsList,weights,yList,1,10^-2)

slope = -1/((out$pPos[2]-out$pNeg[2])/(out$pPos[1]-out$pNeg[1]))
line <- function(x) {slope*(x - out$bisect[1]) + out$bisect[2]}
ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() + geom_segment(aes(x=out$pNeg[1],y=out$pNeg[2],xend=out$pPos[1],yend=out$pPos[2]))+
  stat_function(fun=line) + xlim(0,100) + ylim(0,100)