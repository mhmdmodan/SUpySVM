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
ggpairs(smWinePlot,columns=2:ncol(smWinePlot),mapping=ggplot2::aes(colour = type))

smWine[2:4] <- as.tibble(lapply(smWine[2:4], normalize))

yList <- smWine$type
ptsList <- lapply(1:nrow(smWine), function(n) as.double(c(smWine[n,2],smWine[n,3],smWine[n,4])))
weights <- sapply(ptsList, function(n) {1})
mu <- findMu(weights[which(yList>0)], weights[which(yList<0)])

out <- WSVM(ptsList, weights, yList, mu, 10^-8 ,10^-5)

pred <- sapply(ptsList, function(pt) {
  ifelse((pt-out$bisect) %.% out$w > 0, 1, -1)
  })

confusionMatrix(pred, yList)
