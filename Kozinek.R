library(tidyverse)
library(ggplot2)
library(lubridate)
library(readr)
library(stringr)
library(ggthemes)

set.seed(123)
pts <- tibble(x = sample(1:100, size=50, replace=TRUE), y = sample(1:100, size=50, replace=TRUE), class = 1)

for(i in 1:nrow(pts)) {
  if(pts[i,1] < 37.5) {pts[i,3] <- -1}
}

ggplot(pts,aes(x=x, y=y, color = class)) + geom_point()

X1 <- lapply(1:length(which(pts$class==1)), function(n) as.double(c(pts[n,1],pts[n,2],1)))
X2 <- lapply(1:length(which(pts$class==-1)), function(n) as.double(c(-pts[n,1],-pts[n,2],-1)))

X <- c(X1,X2)

testpts <- as.tibble(data.frame(do.call(rbind, X)))
ggplot(testpts,aes(x=X1, y=X2, color = X3)) + geom_point() + stat_function(fun=function(x) -.543*x) + xlim(-100,100)

normVec <- function(x) sqrt(sum(x^2))
qMin <- function(q,w,xt) {
  normVec(w*(1-q)+q*xt)
}

w <- sample(X,1)[[1]]
xt <- sample(X,1)[[1]]

checkSatisfied <- function(w) {
  for(vec in X) {
    if(checkZero((w %*% vec)[[1]]) < 0) {
      xt <<- vec
      return(TRUE)
      }
  }
  return(FALSE)
}

checkZero <- function(num) {return(ifelse(abs(num)<10^-5, 0, num))}

qvals <- vector(mode='numeric')
while(checkSatisfied(w)) {
  q <- optim(par=.5,fn=qMin,gr=NULL, w, xt,method='Brent',lower=0, upper=1)$par
  qvals <- c(qvals,q)
  w <- w*(1-q)+q*xt
}

w

wPlane <- function(x) {
  (w[1]*x + w[3])/(-w[2])
}

ggplot(pts,aes(x=x, y=y, color = class)) + geom_point() + stat_function(fun=wPlane) + xlim(0,100) + ylim(0,100)

sapply(X, function(n) w%*%n)
