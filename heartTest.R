library(tidyverse)
library(GGally)
library(caret)

source('WSVMFunctions.R')

heart <- read_csv('Data/heart.csv')

heart <- heart %>% mutate(disease = ifelse(disease != 1,-1,1))

heartPlot <- heart
heartPlot$disease <- as.factor(heart$disease)
#ggpairs(heartPlot,columns=1:(ncol(heart)-1),mapping=ggplot2::aes(colour = disease))

heart[1:10] <- as.tibble(lapply(heart[1:10], normalize))

yList <- heart$disease
ptsList <- lapply(1:nrow(heart), function(n) as.double(c(heart[n,1],
                                                         heart[n,2],
                                                         heart[n,3],
                                                         heart[n,4],
                                                         heart[n,5],
                                                         heart[n,6],
                                                         heart[n,7],
                                                         heart[n,8],
                                                         heart[n,9],
                                                         heart[n,10])))
posClass <- which(yList>0)
negClass <- which(yList<0)
weights <- sapply(ptsList, function(n) {1})

weights[posClass] <- rep(1, each=length(posClass))
weights[negClass] <- rep(5, each=length(negClass))

mu <- findMu(weights[posClass], weights[negClass])

out <- WSVM(ptsList, weights, yList, mu, 10^-2 ,10^-5)

#####Test new B

left <- dSumm(func=function(i,j){out$wts[i]*out$ptClass[i]*out$pPosWt[j]*kern(out$pts[[i]],out$pos[[j]])},
              iEnd = length(out$pts), jEnd = length(out$pos))
right <- dSumm(func=function(i,j){out$wts[i]*out$ptClass[i]*out$pNegWt[j]*kern(out$pts[[i]],out$neg[[j]])},
               iEnd = length(out$pts), jEnd = length(out$neg))
bisect <- .5*(left+right)

wts <- out$wts
pts <- out$pts
ptClass <- out$ptClass
pred <- sapply(ptsList, function(pt) {
  wVal <- 0
  for(i in 1:length(pts)) {
    wVal <- wVal + wts[i]*ptClass[i]*kern(pts[[i]], pt)
  }
  return(ifelse(wVal - bisect > 0, 1, -1))
})

confusionMatrix(pred, yList)