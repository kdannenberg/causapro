# all the functions that would be useful in an R package; top level computing functions

library(pcalg)
library(Rgraphviz)
library(grid)
library(bnlearn)
library(stringr)

# setwd("~/Viren/R/Code")

# eqivalent to estimateDAG for numerical, Gauss-distributed data instead of an MSA
estimate_DAG_from_numerical_data <- function(data, alpha, cor_FUN = cor, outpath, solve_conflicts = TRUE, u2pd = "retry", 
                                             conservative = FALSE, maj_rule = FALSE) {
  # if (!length(only_cols) == 0) {
  #   data = data[as.character(only_cols)]
  # }
  
  if (cor_FUN == "none") {
    cor_FUN = function(x) {
      if (dim(x)[1] != dim(x)[2]) {
        stop("Data matrix is not quadratic and can thus not be interpreted as a correlation matix.")
      }
      rownames(x) <- colnames(x)
      return(x)
    }
  } else if (cor_FUN == "") {
    cor_FUN <- cor
  } else {
    cor_FUN <- get(cor_FUN)
  }
  
  n <- nrow(data)
  V <- colnames(data)
  
  suffStat <- list(C = cor_FUN(data), n=n, adaptDF = FALSE) #dm = dat$x        ### WHY COR?!
  
  sink(paste(outpath, "-pc.txt", sep = ""))
    # pc <- pc(suffStat, indepTest = ci_test_cor, alpha = alpha, labels = V, verbose = TRUE) #p=dim(MSA)[2]
    # without solve.confl = true, cycles can emerge in the final CPD"A"G.
    pc <- pc(suffStat, indepTest = gaussCItest, alpha = alpha, labels = V, verbose = TRUE, 
             solve.confl = solve_conflicts, u2pd = u2pd, conservative = conservative, maj.rule = maj_rule) #p=dim(MSA)[2]
  sink()
 
  return(pc)
}
