library(tidyverse)
library(fpCompare)

dot <- function(x,y) {return((x %*% y)[[1]])}
`%.%` <- function(x,y) {return((x %*% y)[[1]])}

summ <- function(func, end=10, ...) {
  initial <- func(1, ...)
  for(i in 2:end) {
    initial <- initial + func(i, ...)
  }
  return(initial)
}

dSumm <- function(func, iEnd = 10, jEnd = 10, ...) {
  iSum <- func(1,1, ...)*0
  for(i in 1:iEnd) {
    jSum <- func(1, 1, ...)*0
    for(j in 1:jEnd) {
      jSum <- jSum + func(i, j, ...)
    }
    iSum <- iSum + jSum
  }
  return(iSum)
}

clamp <- function(expr, minim, maxim) {
  if(expr < minim) {
    return(minim)
  } else if(expr > maxim) {
    return(maxim)
  } else {
    return(expr)
  }
}

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
      maxInd <- ifelse(curVal %>=% maxVal, i, maxInd)
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
  
  return(list(pt = vSum, wt = a))
}

#Get CH points
WRCH <- function(P, s, u, increment) {
  angles <- seq(0+increment, 2*pi, increment)
  
  CH <- findVertex(P, s, c(cos(increment),sin(increment)), u)$pt
  for(angle in angles) {
    CH <- rbind(CH,findVertex(P, s, c(cos(angle),sin(angle)), u)$pt)
  }
  
  CH <- CH %>% as.tibble() %>% unique
  CH[nrow(CH)+1,] <- CH[1,]
  return(CH)
}

getMaxIndex2 <- function(whichClass, P, a, allPts, ptsClass, ptsWts) {
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
      curVal <- 0
      maxVal <- 0
      for(j in 1:length(allPts)) {
        curVal <- curVal + whichClass*ptsWts[j]*ptsClass[j]*kern(allPts[[j]],P[[i]])
        maxVal <- maxVal + whichClass*ptsWts[j]*ptsClass[j]*kern(allPts[[j]],P[[maxInd]])
      }
      maxInd <- ifelse(curVal %>=% maxVal, i, maxInd)
    }
  }
  return(maxInd)
}

findVertex2 <- function(P, s, whichClass, u, allPts, ptsClass, ptsWts) {
  if(length(P) != length(s)) {stop('P and s are not same length')}
  
  a <- vector(mode = 'numeric', length = length(P))
  total <- 0
  while(total %!=% 1) {
    i <- getMaxIndex2(whichClass, P, a, allPts, ptsClass, ptsWts)
    a[i] <- min(s[i]*u, 1-total)
    total <- total + a[i]
  }
  
  vSum <- 0
  for(i in 1:length(P)) {
    vSum <- vSum + a[i]*P[[i]]
  }
  
  return(list(pt = vSum, wt = a))
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

#kernel
kern <- function(x, y) {return(dot(x,y)^6)}

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
  
  pPosPt <- pPos$pt
  pNegPt <- pNeg$pt
  pPosWt <- pPos$wt
  pNegWt <- pNeg$wt
  
  numLoops <- 0
  
  while(TRUE) {
    numLoops <- numLoops + 1
    w <- pPosPt - pNegPt
    ptsClass <- c(rep(-1, length(neg)), rep(1, length(pos)))
    
    vPos <- findVertex2(pos, sPos, -1, rchFactor, c(neg,pos), ptsClass, c(pNegWt,pPosWt))
    vNeg <- findVertex2(neg, sNeg, 1, rchFactor, c(neg,pos), ptsClass, c(pNegWt,pPosWt))
    
    vPosPt <- vPos$pt
    vNegPt <- vNeg$pt
    vPosWt <- vPos$wt
    vNegWt <- vNeg$wt
    
    # if(sqrt(dot(w,w)) < nonSep) {
    #   stop('Hulls are overlapping!')
    # }
    
    # if((1 - dot(w, (pPosPt - vNegPt))/dot(w,w)) < ep & (1 - dot(w, (vPosPt - pNegPt))/dot(w,w)) < ep) {
    #   break
    # }
    if(numLoops > 3) {break}
    print(numLoops)
    # print(pPosPt)
    # print(pNegPt)
    if(w %.% (pPosPt - vPosPt) > w %.% (vNegPt - pNegPt)) {
      
      #q <- ((pPosPt - pNegPt) %.% (pPosPt - vPosPt))/dot(pPosPt - vPosPt,pPosPt - vPosPt)
      
      #### New q calculation
      
      allWt <- c(pNegWt, pPosWt)
      allPts <- c(neg, pos)
      ptClass <- c(rep(-1, length(neg)), rep(1, length(pos)))
      
      numerator <- 0
      for(i in 1:length(allPts)) {
        jSum <- 0
        for(j in 1:length(pos)) {
          jSum <- jSum + allWt[i]*(pPosWt[j]-vPosWt[j])*ptClass[i]*kern(allPts[[i]], pos[[j]])
        }
        numerator <- numerator + jSum
      }
      
      denominator <- 0
      for(i in 1:length(pos)) {
        jSum <- 0
        for(j in 1:length(pos)) {
          jSum <- jSum + (pPosWt[i] - vPosWt[i])*(pPosWt[j] - vPosWt[j])*kern(pos[[i]], pos[[j]])
        }
        denominator <- denominator + jSum
      }
      #print(numerator/denominator)
      q <- clamp(numerator/denominator, 0, 1)
      # print(((pPosPt - pNegPt) %.% (pPosPt - vPosPt))/dot(pPosPt - vPosPt,pPosPt - vPosPt))
      
      #pPosPt <- (1-q)*pPosPt + q*vPosPt
      
      for(i in 1:length(pPosWt)) {
        pPosWt[i] <- (1-q)*pPosWt[i] + q*vPosWt[i]
      }
      
      pPosPt <- 0
      for(i in 1:length(pPosWt)) {
        pPosPt <- pPosPt + pPosWt[i]*pos[[i]]
      }
      
      # print('pos')
    } else {
      
      #q <- -((pPosPt - pNegPt) %.% (pNegPt - vNegPt))/dot(pNegPt - vNegPt,pNegPt - vNegPt)
      
      
      #### New q calculation
      
      allWt <- c(pNegWt, pPosWt)
      allPts <- c(neg, pos)
      ptClass <- c(rep(-1, length(neg)), rep(1, length(pos)))
      
      numerator <- 0
      for(i in 1:length(allPts)) {
        jSum <- 0
        for(j in 1:length(neg)) {
          jSum <- jSum + allWt[i]*(pNegWt[j]-vNegWt[j])*ptClass[i]*kern(allPts[[i]], neg[[j]])
        }
        numerator <- numerator + jSum
      }
      
      denominator <- 0
      for(i in 1:length(neg)) {
        jSum <- 0
        for(j in 1:length(neg)) {
          jSum <- jSum + (pNegWt[i] - vNegWt[i])*(pNegWt[j] - vNegWt[j])*kern(neg[[i]], neg[[j]])
        }
        denominator <- denominator + jSum
      }
      #print(-numerator/denominator)
      q <- clamp(-numerator/denominator, 0, 1)
      # print(-((pPosPt - pNegPt) %.% (pNegPt - vNegPt))/dot(pNegPt - vNegPt,pNegPt - vNegPt))
      #pNegPt <- (1-q)*pNegPt + q*vNegPt
      
      for(i in 1:length(pNegWt)) {
        pNegWt[i] <- (1-q)*pNegWt[i] + q*vNegWt[i]
      }
      
      pNegPt <- 0
      for(i in 1:length(pNegWt)) {
        pNegPt <- pNegPt + pNegWt[i]*neg[[i]]
      }
      # print('neg')
    }
  }
  print(numLoops)
  return(list(w=w, bisect=.5*(w%.%pPosPt + w%.%pNegPt), pts = c(neg, pos), wts=c(pNegWt,pPosWt), ptClass = c(rep(-1, length(neg)), rep(1, length(pos))), 
              pos=pos,neg=neg, pPosWt = pPosWt, pNegWt = pNegWt))
}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#Old b: pNeg+w/2