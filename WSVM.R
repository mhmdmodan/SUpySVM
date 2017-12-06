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

##################
# Test this shit #
##################
testList <- list(c(3,4,5),c(2,3,5),c(6,7,4))
set.seed(1234)
testWRCHPts <- tibble(x = sample(1:50, size=10, replace=TRUE), y = sample(1:50, size=10, replace=TRUE))
ggplot(testWRCHPts, aes(x=x, y=y)) + geom_point()

pts <- lapply(1:nrow(testWRCHPts), function(n) as.double(c(testWRCHPts[n,1],testWRCHPts[n,2])))
weights <- sapply(pts, function(n) {1})

findVertex(pts, weights, c(0,1), 1)

####################
genRand <- function() {c(runif(1,-1,1),runif(1,-1,1))}


WSK <- function(P, s, y) {
  rchFactor <- 1
  ep <- 0.0000001
  classPos <- which(y > 0)
  classNeg <- which(y < 0)
  
  pos <- P[classPos]
  neg <- P[classNeg]
  
  sPos <- s[classPos]
  sNeg <- s[classNeg]
  
  pPos <- findVertex(pos, sPos, genRand(), rchFactor)
  pNeg <- findVertex(neg, sNeg, genRand(), rchFactor)
  
  while(TRUE) {
    w <- pPos - pNeg
    
    vPos <- findVertex(pos, sPos, -w, rchFactor)
    vNeg <- findVertex(neg, sNeg, w, rchFactor)
    
    if(w %.% (pPos - vPos) > w %.% (vNeg - pNeg)) {
      if(1 - dot(w, (pPos - vNeg))/dot(w,w) < ep) {
        break
      }
      
      q <- ((pNeg - pPos) %.% (pPos - vPos))/(pPos - vPos)^2
      pPos <- (1-q)*pPos + q*vPos
      print('pos')
    } else {
      
      if(1 - dot(w, (vPos - pNeg))/dot(w,w) < ep) {
        break
      }
      
      q <- ((pPos - pNeg) %.% (pNeg - vNeg))/(pNeg - vNeg)^2
      pNeg <- (1-q)*pNeg + q*vNeg
      print('neg')
    }
  }
  
  b <- .5*(w%.%pPos + w%.%pNeg)
  return(list(pPos,pNeg,w,b))
}



set.seed(123)
pts <- tibble(x = sample(1:100, size=50, replace=TRUE), y = sample(1:100, size=50, replace=TRUE), class = 1)

for(i in 1:nrow(pts)) {
  if(pts[i,1] < 37.5) {pts[i,3] <- -1}
}

ggplot(pts,aes(x=x, y=y, color = class)) + geom_point()

ptsList <- lapply(1:nrow(pts), function(n) as.double(c(pts[n,1],pts[n,2])))
yList <- pts[,3]
weights <- sapply(ptsList, function(n) {1})

WSK(ptsList,weights,yList)
