library(tidyverse)
library(GGally)
library(caret)

source('WSVMFunctions.R')

wine <- read_csv('Data/Wine.csv')

names(wine) <- c('type',
                 'alcohol',
                 'malicAcid',
                 'ash',
                 'alcAsh',
                 'mg',
                 'phenols',
                 'flav',
                 'nonflav',
                 'proanth',
                 'color',
                 'hue',
                 'code',
                 'proline')

smWine <- wine %>% select(type,alcohol,flav,mg) %>% mutate(type = ifelse(type != 1,-1,1))
smWinePlot <- smWine
smWinePlot$type <- as.factor(smWinePlot$type)
#ggpairs(smWinePlot,columns=2:ncol(smWinePlot),mapping=ggplot2::aes(colour = type))

smWine[2:4] <- as.tibble(lapply(smWine[2:4], normalize))

yList <- smWine$type
ptsList <- lapply(1:nrow(smWine), function(n) as.double(c(smWine[n,2],smWine[n,3],smWine[n,4])))
posClass <- which(yList>0)
negClass <- which(yList<0)
weights <- sapply(ptsList, function(n) {1})

# weights[posClass] <- rep(1/length(posClass), each=length(posClass))
# weights[negClass] <- rep(1/length(negClass), each=length(negClass))

mu <- findMu(weights[posClass], weights[negClass])

out <- WSVM(ptsList, weights, yList, mu, 10^-5 ,10^-5)

wts <- out$wts
pts <- out$pts
ptClass <- out$ptClass
pred <- sapply(ptsList, function(pt) {
  wVal <- 0
  for(i in 1:length(pts)) {
    wVal <- wVal + wts[i]*ptClass[i]*kern(pts[[i]], pt)
  }
  return(ifelse(wVal - out$bisect > 0, 1, -1))
})

confusionMatrix(pred, yList)
