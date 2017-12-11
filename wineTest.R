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

set.seed(123)
smWine[2:4] <- as.tibble(lapply(smWine[2:4], normalize))
smWine <- smWine[sample(nrow(smWine)),]

smWineTrain <- smWine[1:round(nrow(smWine)*.70),]
smWineTest <- smWine[-(1:round(nrow(smWine)*.70)),]

yList <- smWineTrain$type
ptsList <- lapply(1:nrow(smWineTrain), function(n) as.double(c(smWineTrain[n,2],smWineTrain[n,3],smWineTrain[n,4])))
ptsListTest <- lapply(1:nrow(smWineTest), function(n) as.double(c(smWineTest[n,2],smWineTest[n,3],smWineTest[n,4])))
weights <- sapply(ptsList, function(n) {1})
mu <- findMu(weights[which(yList>0)], weights[which(yList<0)])

out <- WSVM(ptsList, weights, yList, mu, 10^-5 ,10^-5)

pred <- sapply(ptsListTest, function(pt) {
  ifelse((pt-out$bisect) %.% out$w > 0, 1, -1)
  })

confusionMatrix(pred, smWineTest$type)
