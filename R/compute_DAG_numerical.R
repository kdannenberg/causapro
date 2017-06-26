library(pcalg)
library(Rgraphviz)
library(grid)
library(bnlearn)
library(stringr)

source("ci-tests.R")
# source("compute_DAG_functions.R")

setwd("~/Viren/R/Code")

read_data <- function(filename, transpose=FALSE) {
  data = read.csv2(filename, row.names = 1, check.names=FALSE) # if check.names, an X is prepended to numerical column-names
  if (transpose)
    data <- t(data)
  return(data)
}

# eqivalent to estimateDAG for numerical, Gauss-distributed data instead of an MSA
estimate_DAG_from_numerical_data <- function(data, alpha, outpath) {
  if (!length(only_cols) == 0) {
    data = data[as.character(only_cols)]
  }
  
  n <- nrow(data)
  V <- colnames(data)
  
  suffStat <- list(C = cor(data), n=n, adaptDF = FALSE) #dm = dat$x        ### WHY COR?!
  
  sink(paste(outpath, "-pc.txt", sep = ""))
    pc <- pc(suffStat, indepTest = ci_test_cor, alpha = alpha, labels = V, verbose = TRUE) #p=dim(MSA)[2]
  sink()
 
  return(pc)
}
