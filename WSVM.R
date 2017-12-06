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

findVertex(pts, weights, c(1,0), 1)
