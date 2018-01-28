library(tidyverse)
library(fpCompare)

dot <- function(x, y) { return((x %*% y)[[1]]) }
`%.%` <- function(x, y) { return((x %*% y)[[1]]) }

summ <- function(func, end = 10, ...) {
    initial <- func(1, ...) * 0
    for (i in 1:end) {
        initial <- initial + func(i, ...)
    }
    return(initial)
}

dSumm <- function(func, iEnd = 10, jEnd = 10, ...) {
    iSum <- func(1, 1, ...) * 0
    for (i in 1:iEnd) {
        jSum <- func(1, 1, ...) * 0
        for (j in 1:jEnd) {
            jSum <- jSum + func(i, j, ...)
        }
        iSum <- iSum + jSum
    }
    return(iSum)
}

clamp <- function(expr, minim, maxim) {
    if (expr < minim) {
        return(minim)
    } else if (expr > maxim) {
        return(maxim)
    } else {
        return(expr)
    }
}

checkSameDim <- function(n, P) {
    if (!is.list(P)) { stop('P is not a list') }
    dimension <- length(n)
    for (pt in P) {
        if (length(pt) != dimension) { return(FALSE) }
        }
    return(TRUE)
}

getMaxIndex <- function(n, P, a) {
    maxInd <- -1
    for (i in 1:length(P)) {
        if (a[i] == 0) {
            maxInd <- i
            break
        }
    }
    if (maxInd == -1) {
        stop('None are 0')
    }
    for (i in 1:length(P)) {
        if (a[i] == 0) {
            curVal <- dot(n, P[[i]])
            maxVal <- dot(n, P[[maxInd]])
            maxInd <- ifelse(curVal %>=% maxVal, i, maxInd)
        }
    }
    return(maxInd)
}

#Finding WRCH vertex
findVertex <- function(P, s, n, u) {
    if (!checkSameDim(n, P)) { stop('P and n are not same dims') }
    if (length(P) != length(s)) { stop('P and s are not same length') }

    a <- vector(mode = 'numeric', length = length(P))
    total <- 0
    while (total %!=% 1) {
        i <- getMaxIndex(n, P, a)
        a[i] <- min(s[i] * u, 1 - total)
        total <- total + a[i]
    }

    vSum <- summ(function(i) { a[i] * P[[i]] }, length(P))

    return(list(pt = vSum, wt = a))
}

#Get CH points
WRCH <- function(P, s, u, increment) {
    angles <- seq(0 + increment, 2 * pi, increment)

    CH <- findVertex(P, s, c(cos(increment), sin(increment)), u)$pt
    for (angle in angles) {
        CH <- rbind(CH, findVertex(P, s, c(cos(angle), sin(angle)), u)$pt)
    }

    CH <- CH %>% as.tibble() %>% unique
    CH[nrow(CH) + 1,] <- CH[1,]
    return(CH)
}

getMaxIndex2 <- function(whichClass, P, a, allPts, ptsClass, ptsWts) {
    maxInd <- -1
    for (i in 1:length(P)) {
        if (a[i] == 0) {
            maxInd <- i
            break
        }
    }
    if (maxInd == -1) {
        stop('None are 0')
    }
    for (i in 1:length(P)) {
        if (a[i] == 0) {
            curVal <- 0
            maxVal <- 0
            for (j in 1:length(allPts)) {
                curVal <- curVal + whichClass * ptsWts[j] * ptsClass[j] * kern(allPts[[j]], P[[i]])
                maxVal <- maxVal + whichClass * ptsWts[j] * ptsClass[j] * kern(allPts[[j]], P[[maxInd]])
            }
            maxInd <- ifelse(curVal %>=% maxVal, i, maxInd)
        }
    }
    return(maxInd)
}

findVertex2 <- function(P, s, whichClass, u, allPts, ptsClass, ptsWts) {
    if (length(P) != length(s)) { stop('P and s are not same length') }

    a <- vector(mode = 'numeric', length = length(P))
    total <- 0
    while (total %!=% 1) {
        i <- getMaxIndex2(whichClass, P, a, allPts, ptsClass, ptsWts)
        a[i] <- min(s[i] * u, 1 - total)
        total <- total + a[i]
    }

    vSum <- summ(function(i) { a[i] * P[[i]] }, length(P))

    return(list(pt = vSum, wt = a))
}

#Get center of points
getCenter <- function(P, s) {
    totalWeights <- sum(s)
    newS <- sapply(s, function(x) { x / totalWeights })

    cumMean <- summ(function(i) { P[[i]] * newS[i] }, length(P[[1]]))

    return(cumMean)
}

#Find a mu
findMu <- function(weights1, weights2) {
    return(1 / (0.9 * min(sum(weights1), sum(weights2))))
}

genRand <- function() { c(runif(1, -1, 1), runif(1, -1, 1)) }

#kernel
kern <- function(x, y) { return(dot(x, y) ^ 2) }

#WSVM Algorithm
WSVM <- function(P, s, y, rchFactor, ep, nonSep) {

    classPos <- which(y > 0)
    classNeg <- which(y < 0)

    pos <- P[classPos]
    neg <- P[classNeg]
    allPts <- c(neg, pos)

    sPos <- s[classPos]
    sNeg <- s[classNeg]

    centerPos <- getCenter(pos, sPos)
    centerNeg <- getCenter(neg, sNeg)

    pPosWt <- findVertex(pos, sPos, centerNeg - centerPos, rchFactor)$wt
    pNegWt <- findVertex(neg, sNeg, centerPos - centerNeg, rchFactor)$wt

    ptsClass <- c(rep(-1, length(neg)), rep(1, length(pos)))

    numLoops <- 0

    while (TRUE) {
        numLoops <- numLoops + 1
        allWt <- c(pNegWt, pPosWt)

        vPosWt <- findVertex2(pos, sPos, -1, rchFactor, c(neg, pos), ptsClass, c(pNegWt, pPosWt))$wt
        vNegWt <- findVertex2(neg, sNeg, 1, rchFactor, c(neg, pos), ptsClass, c(pNegWt, pPosWt))$wt

        allVWt <- c(vNegWt, vPosWt)

        # if(numLoops > 4) {break}
        print(numLoops)

        w.vPos_pNeg <- dSumm(function(i, j) {
            allWt[i] * ptsClass[i] * (
            max(ptsClass[j], 0) * allVWt[j] + min(ptsClass[j], 0) * allWt[j]
            ) * kern(allPts[[i]], allPts[[j]])
        }, iEnd = length(allPts), jEnd = length(allPts))

        w.pPos_vNeg <- dSumm(function(i, j) {
            allWt[i] * ptsClass[i] * (
            max(ptsClass[j], 0) * allWt[j] + min(ptsClass[j], 0) * allVWt[j]
            ) * kern(allPts[[i]], allPts[[j]])
        }, iEnd = length(allPts), jEnd = length(allPts))

        w.w <- dSumm(function(i, j) {
            allWt[i] * allWt[j] * ptsClass[i] * ptsClass[j] * kern(allPts[[i]], allPts[[j]])
        }, iEnd = length(allPts), jEnd = length(allPts))

        if (w.pPos_vNeg > w.vPos_pNeg) {

            if (1 - w.vPos_pNeg / w.w < ep) { break }

            numerator <- dSumm(function(i, j) {
                allWt[i] * (pPosWt[j] - vPosWt[j]) * ptsClass[i] * kern(allPts[[i]], pos[[j]])
            }, iEnd = length(allPts), jEnd = length(pos))

            denominator <- dSumm(function(i, j) {
                (pPosWt[i] - vPosWt[i]) * (pPosWt[j] - vPosWt[j]) * kern(pos[[i]], pos[[j]])
            }, iEnd = length(pos), jEnd = length(pos))

            q <- clamp(numerator / denominator, 0, 1)

            for (i in 1:length(pPosWt)) {
                pPosWt[i] <- (1 - q) * pPosWt[i] + q * vPosWt[i]
            }

            # print('pos')
        } else {

            if (1 - w.pPos_vNeg / w.w < ep) { break }

            numerator <- dSumm(function(i, j) {
                allWt[i] * (pNegWt[j] - vNegWt[j]) * ptsClass[i] * kern(allPts[[i]], neg[[j]])
            }, iEnd = length(allPts), jEnd = length(neg))


            denominator <- dSumm(function(i, j) {
                (pNegWt[i] - vNegWt[i]) * (pNegWt[j] - vNegWt[j]) * kern(neg[[i]], neg[[j]])
            }, iEnd = length(neg), jEnd = length(neg))

            q <- clamp(-numerator / denominator, 0, 1)

            for (i in 1:length(pNegWt)) {
                pNegWt[i] <- (1 - q) * pNegWt[i] + q * vNegWt[i]
            }

            # print('neg')
        }
    }
    print(numLoops)

    bisect <- .5 * dSumm(
    function(i, j) {
        allWt[i] * ptsClass[i] * allWt[j] * kern(allPts[[i]], allPts[[j]])
    }, iEnd = length(allPts), jEnd = length(allPts)
  )

    return(list(bisect = bisect, pts = allPts, wts = c(pNegWt, pPosWt), ptsClass = ptsClass,
              pos = pos, neg = neg, pPosWt = pPosWt, pNegWt = pNegWt))
}

normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
}

#Old b: pNeg+w/2