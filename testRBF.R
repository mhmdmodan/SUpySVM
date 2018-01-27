library(tidyverse)
library(GGally)
library(ggplot2)
library(caret)

source('WSVMFunctions.R')

set.seed(123)

poop <- as.data.frame(as.list(c(genRand() / 2, 1)))
for (i in 1:200) {
    poop <- rbind(poop, as.list(c(genRand() / 2, 1)))
}
colnames(poop) <- c('x', 'y', 'class')


topCirc <- function(x) { sqrt(-(x) ^ 2 + .25 ^ 2) }
botCirc <- function(x) {-sqrt(-(x) ^ 2 + .25 ^ 2) }
poop <- as.tibble(poop)
todel <- 1
for (i in 1:nrow(poop)) {
    margin <- .08
    if (poop$x[i] > -0.25 & poop$x[i] < .25) {
        if (poop$y[i] < topCirc(poop$x[i]) & poop$y[i] > botCirc(poop$x[i])) {
            poop$class[i] <- -1
        }
        if (abs(poop$y[i] - topCirc(poop$x[i])) < margin | abs(poop$y[i] - botCirc(poop$x[i])) < margin) {
            todel <- c(todel, i)
        }
    }
}

poop <- poop[-todel,]

ggplot(poop, aes(x = x, y = y)) +
    geom_point(aes(color = class)) +
    stat_function(fun = topCirc, n = 1000, na.rm = TRUE) +
    stat_function(fun = botCirc, n = 1000, na.rm = TRUE) +
    xlim(-.5, .5) + ylim(-.5, .5)



yList <- poop$class
ptsList <- lapply(1:nrow(poop), function(n) as.double(c(poop[n, 1], poop[n, 2])))
posClass <- which(yList > 0)
negClass <- which(yList < 0)
weights <- sapply(ptsList, function(n) { 1 })

mu <- findMu(weights[posClass], weights[negClass])

out <- WSVM(ptsList, weights, yList, 1, 10 ^ -2, 10 ^ -5)


#####Test new B

bisect <- out$bisect
wts <- out$wts
pts <- out$pts
ptClass <- out$ptsClass
pred <- sapply(ptsList, function(pt) {
    wVal <- 0
    for (i in 1:length(pts)) {
        wVal <- wVal + wts[i] * ptClass[i] * kern(pts[[i]], pt)
    }
    return(ifelse(wVal - bisect > 0, 1, -1))
})

confusionMatrix(pred, yList)

dtInit <- as.data.frame(list(x = 1, y = 1, w = 1))
xVal <- -0.5
while (xVal %<=% 0.50) {
    yVal <- -0.5
    while (yVal %<=% 0.50) {
        wVal <- 0
        for (i in 1:length(pts)) {
            wVal <- wVal + wts[i] * ptClass[i] * kern(pts[[i]], c(xVal, yVal))
        }
        dtInit <- rbind(dtInit, as.data.frame(list(x = xVal, y = yVal, w = ifelse(wVal - bisect > 0, 1, -1))))

        yVal <- yVal + 0.02
    }
    xVal <- xVal + 0.02
}

ggplot(dtInit[-1,], aes(x = x, y = y)) +
    geom_raster(aes(fill = w), interpolate = FALSE) +
    scale_fill_gradient(low = 'salmon', high = 'dodgerblue') +
    geom_point(data = poop, aes(x = x, y = y, color = class)) + labs(fill = 'w-b') +
    scale_color_gradient(high = 'blue', low = 'red')
