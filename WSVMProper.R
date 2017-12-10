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

savePlot <- function(params, index) {
  out <- params
  
  slope = -1/((out$pPos[2]-out$pNeg[2])/(out$pPos[1]-out$pNeg[1]))
  line <- function(x) {slope*(x - out$bisect[1]) + out$bisect[2]}
  toSave <- ggplot(pts,aes(x=x, y=y, color = class)) + 
    geom_path(data=negCH, aes(x=V1,y=V2), color='palegreen3', size=1.3) +
    geom_path(data=posCH, aes(x=V1,y=V2), color="palegreen3", size=1.3) +
    geom_point() + 
    geom_segment(aes(x=out$pNeg[1],y=out$pNeg[2],xend=out$pPos[1],yend=out$pPos[2]), color = "dodgerblue", size=1.4, lineend = 'round')+
    stat_function(fun=line, color = c("#7D5BA6"), size=1.0) + 
    scale_color_gradient(low="black", high="#FC6471") +
    xlim(0,1) + 
    ylim(0,1) + 
    theme_fivethirtyeight() + 
    theme(legend.position="none") +
    labs(title = 'Calculating the Hyperplane')
  suppressWarnings(ggsave(filename=paste0('Anims/WAnim/',index,'.png'), plot = toSave,width=7,height=7,units='in',dpi=200))
}

genRand <- function() {c(runif(1,-1,1),runif(1,-1,1))}

#WSVM Algorithm
WSVM <- function(P, s, y, rchFactor, ep) {
  
  classPos <- which(y > 0)
  classNeg <- which(y < 0)
  
  pos <- P[classPos]
  neg <- P[classNeg]
  sPos <- s[classPos]
  sNeg <- s[classNeg]
  
  # randVec1 <- genRand()
  # randVec2 <- genRand()
  # pPos <- findVertex(pos, sPos, randVec1, rchFactor)
  # pNeg <- findVertex(neg, sNeg, randVec2, rchFactor)

  centerPos <- getCenter(pos, sPos)
  centerNeg <- getCenter(neg, sNeg)

  pPos <- findVertex(pos, sPos, centerNeg - centerPos, rchFactor)
  pNeg <- findVertex(neg, sNeg, centerPos - centerNeg, rchFactor)

  numLoops <- 0
  
  while(TRUE) {
    numLoops <- numLoops + 1
    w <- pPos - pNeg
    
    #savePlot(list(w=w, bisect=pNeg+w/2, pPos=pPos, pNeg=pNeg), numLoops)
    
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
  print(numLoops)
  return(list(w=w, bisect=pNeg+w/2, pPos=pPos, pNeg=pNeg))
}
#####################
# SETUP
#####################
set.seed(174)
pts <- tibble(x = sample(1:100, size=80, replace=TRUE)/100, y = sample(1:100, size=80, replace=TRUE)/100, class = 1)
theorLine <- function(x) {.28/.125*(x-.55)+.5}

for(i in 1:nrow(pts)) {
  if(pts[i,2] > theorLine(pts[i,1])) {pts[i,3] <- -1}
}

ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() +
  stat_function(fun=theorLine) + xlim(0,1) + ylim(0,1)

ptsList <- lapply(1:nrow(pts), function(n) as.double(c(pts[n,1],pts[n,2])))
yList <- pts$class
weights <- sapply(ptsList, function(n) {1})
###################################
###################################
###################################

negCH <- WRCH(ptsList[which(yList<0)], weights[which(yList<0)], 1, 0.01)
posCH <- WRCH(ptsList[which(yList>0)], weights[which(yList>0)], 1, 0.01)
out <- WSVM(ptsList,weights,yList,1,10^-5)

slope = -1/((out$pPos[2]-out$pNeg[2])/(out$pPos[1]-out$pNeg[1]))
line <- function(x) {slope*(x - out$bisect[1]) + out$bisect[2]}
ggplot(pts,aes(x=x, y=y, color = class)) + 
  geom_path(data=negCH, aes(x=V1,y=V2), color='palegreen3', size=1.3) +
  geom_path(data=posCH, aes(x=V1,y=V2), color="palegreen3", size=1.3) +
  geom_point() + 
  geom_segment(aes(x=out$pNeg[1],y=out$pNeg[2],xend=out$pPos[1],yend=out$pPos[2]), color = "dodgerblue", size=1.4, lineend = 'round')+
  stat_function(fun=line, color = c("#7D5BA6"), size=1.0) + 
  scale_color_gradient(low="black", high="#FC6471") +
  xlim(0,1) + 
  ylim(0,1) + 
  theme_fivethirtyeight() + 
  theme(legend.position="none") +
  labs(title = 'Calculating the Hyperplane')
 ggsave(filename='Static/optimal.png',width=7,height=7,units='in',dpi=200)

#############
#Test out CH finder
set.seed(127)
pts <- tibble(x = sample(1:100, size=15, replace=TRUE)/100, y = sample(1:100, size=15, replace=TRUE)/100, class = 1) %>% 
  filter(!(x>.7 & x<.8))

ptsList <- lapply(1:nrow(pts), function(n) as.double(c(pts[n,1],pts[n,2])))
yList <- pts$class
weights <- sapply(ptsList, function(n) {1})

numFrames <- 150
increment <- (1/numFrames)*.9
rFac <- 1
for(i in 1:numFrames){
  chPts <- WRCH(ptsList, weights, rFac, 0.01)
  rchToSave<-ggplot(data=pts, aes(x=x,y=y)) + 
    geom_path(data=chPts, aes(x=V1,y=V2), color="#EE6363", size=1.3) +
    geom_point() + 
    theme_fivethirtyeight() + 
    xlim(0,1) +
    ylim(0,1) +
    labs(title=paste0('Convex hull reduction - μ = ',formatC(as.numeric(rFac), format = 'f', flag='0', digits = 3)))
  suppressWarnings(ggsave(filename=paste0('Anims/RCHAnim/',i,'.png'), plot = rchToSave,width=7,height=7,units='in',dpi=200)) 
  rFac <- rFac - increment
}