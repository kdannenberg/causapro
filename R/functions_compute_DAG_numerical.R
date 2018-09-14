# all the functions that would be useful in an R package; top level computing functions

library(pcalg)
library(Rgraphviz)
library(grid)
library(bnlearn)
library(stringr)

# setwd("~/Viren/R/Code")

# eqivalent to estimateDAG for numerical, Gauss-distributed data instead of an MSA
estimate_DAG_from_numerical_data <- function(data, alpha, cor_FUN = cor, outpath,
                                             type_of_variables = "continuous",
                                             indepTest, suffStat, # werden aktuell noch nciht genutzt
                                             solve_conflicts = TRUE, u2pd = "retry",
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
  p <- ncol(data)
  V <- colnames(data)

  if (missing(indepTest) || is.null(indepTest)) {
    if (startsWith(type_of_variables, "cont")) {
      if (missing(suffStat) || is.null(suffStat)) {
        suffStat <- list(C = cor_FUN(data), n=n, adaptDF = FALSE) #dm = dat$x        ### WHY COR?!
      }
      indepTest <- gaussCItest   # partial correlation
    } else if (startsWith(type_of_variables, "ord")) {
      indepTest <- "jt" # Jonckhere-Terpstra
      # if (missing(suffStat) || is.null(indepTest)) {
      #   suffStat <- list(dm = data,#scale_data_for_pc_discrete(data),
      #   nlev = apply(data, 2, function(row){length(unique(row))}),
      #   adaptDF = FALSE)
      # }
      # indepTest <- ci_test_pc("jt")   # Jonckhere-Terpstra
    } else if (startsWith(type_of_variables, "nom")) {
      if (missing(suffStat) || is.null(suffStat)) {
        suffStat <- list(dm = data, #scale_data_for_pc_discrete(data),
                         nlev = apply(data, 2, function(row){length(unique(row))}),
                         adaptDF = FALSE)
      }
      indepTest <- disCItest   # G^2
    }
  } #else {
    if (typeof(indepTest) == "character") {
      if (indepTest %in% c("jt", "mc-jt", "smc-jt")) {
        if (missing(suffStat) || is.null(suffStat)) {
          suffStat <- list(dm = data,#scale_data_for_pc_discrete(data),
                           nlev = apply(data, 2, function(row){length(unique(row))}),
                           adaptDF = FALSE)
        }
        indepTest <- ci_test_pc(indepTest)
      } else if (tolower(indepTest) %in% c("copula", "cpc", "copc", "copulapc", "copula-pc")) {
        cat("Computing Copula...")
        ## fill in some (none!) missing values
        # the expected percentage of missing values
        beta = 0
        # MCAR
        Y = apply(data, 2, function(x){x[sample(n, round(n*runif(1,0,2*beta)))] = NA; x})
        # MAR
        # Y = generateMAR(data, beta)

        ## estimate underlying correlation matrix and (effective) sample size
        # copula object
        # nsamp = 1000
        nsamp = 1000
        cop.obj <- copula.estimate(Y, nsamp = nsamp, S0 = diag(p)/100, plugin.threshold = 50, verb = F)
        C_samples <- cop.obj$C.psamp[,, 501:1000]
        # average correlation matrix
        corr.cop = apply(C_samples, c(1,2), mean)
        # # local effective sample size
        # less.cop = ((1-corr.cop^2)^2)/apply(C_samples,c(1,2), var)
        # ## global effective sample size
        # gess.cop = mean(less.cop[upper.tri(less.cop)])

        suffStat = list(C = corr.cop, n = n)
        indepTest <- gaussCItest
        # ## call the PC algorithm for causal discovery
        # # CoPC + SS
        # graph.cpc.ss = pc(suffStat = list(C = corr.cop, n = n),
        #                   indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
        # # CoPC + GESS
        # graph.cpc.gess = pc(suffStat = list(C = corr.cop, n = gess.cop),
        #                     indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
        # # CoPC + LESS
        # graph.cpc.less = pc(suffStat = list(C = corr.cop, n = n, ESS.mat = less.cop),
        #                     indepTest = gaussCItest.local, labels = var.names, alpha = alpha, conservative = T)
        cat("done.\n")
      }
    # }
  }


  sink(paste(outpath, ".txt", sep = ""))
    # pc <- pc(suffStat, indepTest = ci_test_cor, alpha = alpha, labels = V, verbose = TRUE) #p=dim(MSA)[2]
    # without solve.confl = true, cycles can emerge in the final CPD"A"G.
    pc <- pc(suffStat = suffStat, indepTest = indepTest, alpha = alpha, labels = V, verbose = TRUE,
             solve.confl = solve_conflicts, u2pd = u2pd, conservative = conservative, maj.rule = maj_rule) #p=dim(MSA)[2]
  sink()

  return(pc)
}

scale_data_for_pc_discrete <- function(data) {
  scaled_data <- 2 * data
  scaled_data <- scaled_data - min(scaled_data)
  return(scaled_data)
}
