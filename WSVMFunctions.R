library(tidyverse)

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

getMaxIndex <- function(n, P, a) {
  maxInd <- -1
  for(i in 1:length(P)) {
    if(a[i] == 0) {
      maxInd <- i
      break
    }
  }
  if(maxInd == -1) {
    stop('None are 0')
  }
  for(i in 1:length(P)) {
    if(a[i] == 0) {
      curVal <- dot(n, P[[i]])
      maxVal <- dot(n, P[[maxInd]])
      maxInd <- ifelse(curVal > maxVal, i, maxInd)
    }
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
    i <- getMaxIndex(n, P, a)
    a[i] <- min(s[i]*u, 1-total)
    total <- total + a[i]
  }
  
  vSum <- 0
  for(i in 1:length(P)) {
    vSum <- vSum + a[i]*P[[i]]
  }
  
  return(vSum)
}

#Get CH points
WRCH <- function(P, s, u, increment) {
  angles <- seq(0+increment, 2*pi, increment)
  
  CH <- findVertex(P, s, c(cos(increment),sin(increment)), u)
  for(angle in angles) {
    CH <- rbind(CH,findVertex(P, s, c(cos(angle),sin(angle)), u))
  }
  
  CH <- CH %>% as.tibble() %>% unique
  CH[nrow(CH)+1,] <- CH[1,]
  return(CH)
}

#Get center of points
getCenter <- function(P, s) {
  totalWeights <- sum(s)
  newS <- sapply(s, function(x){x/totalWeights})
  
  cumMean <- vector(mode='numeric', length = length(P[[1]]))
  for(i in 1:length(P)) {
    cumMean <- cumMean + P[[i]]*newS[i]
  }
  
  return(cumMean)
}

#Find a mu
findMu <- function(weights1,weights2) {
  return(1/(0.9*min(sum(weights1), sum(weights2))))
}

genRand <- function() {c(runif(1,-1,1),runif(1,-1,1))}

#WSVM Algorithm
WSVM <- function(P, s, y, rchFactor, ep, nonSep) {
  
  classPos <- which(y > 0)
  classNeg <- which(y < 0)
  
  pos <- P[classPos]
  neg <- P[classNeg]
  sPos <- s[classPos]
  sNeg <- s[classNeg]
  
  centerPos <- getCenter(pos, sPos)
  centerNeg <- getCenter(neg, sNeg)
  
  pPos <- findVertex(pos, sPos, centerNeg - centerPos, rchFactor)
  pNeg <- findVertex(neg, sNeg, centerPos - centerNeg, rchFactor)
  
  numLoops <- 0
  
  while(TRUE) {
    numLoops <- numLoops + 1
    w <- pPos - pNeg
    
    vPos <- findVertex(pos, sPos, -w, rchFactor)
    vNeg <- findVertex(neg, sNeg, w, rchFactor)
    
    if(sqrt(dot(out$w,out$w)) < nonSep) {
      stop('Hulls are overlapping!')
    }
    
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
  print(numLoops)
  return(list(w=w, bisect=pNeg+w/2, pPos=pPos, pNeg=pNeg))
}