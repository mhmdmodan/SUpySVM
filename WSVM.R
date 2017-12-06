library(tidyverse)
library(ggplot2)
library(ggthemes)
library(fpCompare)

findVertex <- function(P, s, n) {
  a <- 0; sum <- 0
  
}

checkSameDim <- function(n, P) {
  if(!is.list(P)) {stop('P is not a list')}
  dimension <- length(n)
  for(pt in P) {
    if(length(pt) != dimension) {return(FALSE)}
  }
  
  return(TRUE)
}