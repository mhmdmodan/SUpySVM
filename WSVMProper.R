library(tidyverse)
library(ggplot2)
library(ggthemes)
library(fpCompare)

five38Mod <- theme_fivethirtyeight() + theme(legend.position = "right", 
                                             legend.direction = 'vertical',
                                             axis.title = element_text())

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
set.seed(1232234)
pts <- tibble(x = sample(1:160, size=80, replace=TRUE), y = sample(1:90, size=80, replace=TRUE), class = 1) %>% 
  if(x=)
ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() + 
  stat_function(fun=theorLine) + xlim(0,100) + ylim(0,100)
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

out <- WSVM(ptsList,weights,yList,1/5,10^-5)

slope = -1/((out$pPos[2]-out$pNeg[2])/(out$pPos[1]-out$pNeg[1]))
line <- function(x) {slope*(x - out$bisect[1]) + out$bisect[2]}
ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() + geom_segment(aes(x=out$pNeg[1],y=out$pNeg[2],xend=out$pPos[1],yend=out$pPos[2]))+
  stat_function(fun=line) + xlim(0,100) + ylim(0,100)

#############
#Test out CH finder
chPts <- WRCH(ptsList, weights, 1/5, 0.01)
ggplot(data=pts, aes(x=x,y=y)) + 
  geom_point() + 
  geom_path(data=chPts, aes(x=V1,y=V2), color="#EE6363", size=1.3) +
  theme_fivethirtyeight() + 
  xlim(0,100) +
  ylim(0,100)

ggsave(filename='test.png',width=7,height=7,units='in',dpi=72)