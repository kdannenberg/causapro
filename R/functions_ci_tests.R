disCItest_given_nmin <- function(n_min) {
  return(function(x, y, S, stat) {
    stat$dm <- replace_with_numbers(stat$dm)
    gSquareDis(x, y, S, dm = stat$dm, nlev = stat$nlev,
                                            adaptDF = stat$adaptDF, n.min = n_min, verbose = TRUE)})
}
# staggered <- function(measure, ...){
#   return( function(x) applyStaggered(x, measure, ...) )
# }
#### disCItest_nmin <- disCItest_given_nmin(1000)
# disCItest_nmin_1000 <- function(x, y, S, stat) {
#   gSquareDis(x, y, S, dm = stat$dm, nlev = stat$nlev, adaptDF = stat$adaptDF, n.min = 1000, verbose = TRUE)
# }

# ## TODO: Wrapper für ci.test schreiben und hier verwenden
# ci_test_wrapper <- function() {
# }

# ci_test_pc <- function(test, ...) {
#   if (missing(test)) {
#     return(function(x, y, S, stat, ...) {
#       colnames(stat$dm) <- paste("P", seq(1:dim(stat$dm)[2]), sep = "")
#       if (length(S) == 0) {
#         result <- ci.test(x = paste("P", x, sep = ""), y = paste("P", y, sep = ""), data = data.frame(stat$dm), debug = TRUE, ...)
#       } else {
#         z = paste("P", S, sep = "")
#         result <- ci.test(x = paste("P", x, sep = ""), y = paste("P", y, sep = ""), z = z, data = data.frame(stat$dm), debug = TRUE, ...)
#       }
#       return(result$p.value)
#     })
#   } else {
#     return(function(x, y, S, stat, ...) {
#       # colnames(stat$dm) <- paste("P", seq(1:dim(stat$dm)[2]), sep = "")
#       if (length(S) == 0) {
#         result <- ci.test(x = paste("P", x, sep = ""), y = paste("P", y, sep = ""), data = data.frame(stat$dm), test = test, debug = TRUE, B = 1, ...)
#       } else {
#         z = paste("P", S, sep = "")
#         result <- ci.test(x = paste("P", x, sep = ""), y = paste("P", y, sep = ""), z = z, data = data.frame(stat$dm), test = test, debug = TRUE, ...)
#       }
#       return(result$p.value)
#     })
#   }
# }

ci_test_pc <- function(ci_test, ...) {
  if (missing(ci_test)) {
    return(function(x, y, S, stat, ...) {
      # colnames(stat$dm) <- paste("P", seq(1:dim(stat$dm)[2]), sep = "")
      if (length(S) == 0) {
        result <- ci.test(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]), data = data.frame(stat$dm), debug = TRUE, ...)
      } else {
        # z = paste("P", S, sep = "")
        result <- ci.test(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]), z = as.ordered(stat$dm[,S]), data = data.frame(stat$dm), debug = TRUE, ...)
      }
      return(result$p.value)
    })
  } else {
    independence_test <- function(x, y, S, stat, test, ...) {
      if (length(S) == 0) {
        result <- ci.test(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]), data = data.frame(stat$dm, check.names = FALSE), test = test, debug = TRUE, B = 1, ...)
      } else if (length(S) == 1) {
        with_z_1 <- function_set_parameters(ci.test, parameters = list(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]), z = as.ordered(stat$dm[,S]), data = data.frame(stat$dm, check.names = FALSE), test = test, debug = TRUE))
        result <- with_z_1(...)
        # result <- ci.test(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]), z = as.ordered(stat$dm[,S]), data = data.frame(stat$dm, check.names = FALSE), test = test, debug = TRUE, ...)
      } else {
        # print(S)
        z_list <- list()
        S_data <- stat$dm[,S]
        # print(S_data)
        for (i in seq(1:dim(S_data)[2])) {
          z_list[[i]] <- as.ordered(S_data[,i])
        }
        # print(z_list)
        z <- do.call(data.frame, args = z_list)
        with_z_higher <- function_set_parameters(ci.test, parameters = list(x = as.ordered(stat$dm[,x]), y = as.ordered(stat$dm[,y]),
                                                                            # z = apply(stat$dm[,S], 2, function(v) {as.ordered(as.numeric(v))}),
                                                                            z = z,
                                                                            data = data.frame(stat$dm, check.names = FALSE), test = test, debug = TRUE))

        # debug(with_z_higher)
        result <- with_z_higher(...)
      }
      return(result$p.value)
    }
    # debug(independence_test)
    par <- list(test = ci_test)
    independ_test <- function_set_parameters(FUN = independence_test, parameters = par)
    return(independ_test)
  }
}



## ci_test_pc_chi_square <- ci_test_pc("x2")

## ci_test_cor <- ci_test_pc("cor")

## Code for Copula-PC from Cui et al. (2018).
## For details about the arguements and outputs, refer to function 'sbgcop.mcmc' in R package 'sbgcop',
## https://cran.r-project.org/web/packages/sbgcop/index.html.
library(sbgcop)
copula.estimate <- function (Y, n0 = dim(Y)[2] + 1, S0 = diag(dim(Y)[2])/n0, nsamp = 100,
                             odens = max(1, round(nsamp/1000)), impute = any(is.na(Y)),
                             plugin.threshold = 100, plugin.marginal = (apply(Y, 2, function(x) {
                               length(unique(x))
                             }) > plugin.threshold), seed = NULL, verb = TRUE)
{
  require(sbgcop)
  ok_S0 <- all(eigen(S0)$val > 0) & dim(S0)[1] == dim(Y)[2] &
    dim(S0)[2] == dim(Y)[2]
  ok_n0 <- (n0 >= 0)
  if (!ok_S0) {
    stop("Error: S0 must be a positive definite p x p matrix \n")
  }
  if (!ok_n0) {
    stop("Error: n0 must be positive \n")
  }

  vnames <- colnames(Y)
  Y <- as.matrix(Y)
  colnames(Y) <- vnames
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  set.seed(seed)
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  S <- cov(Z)
  Y.pmean <- Y
  if (impute) {
    Y.pmean <- matrix(0, nrow = n, ncol = p)
  }
  LPC <- NULL
  C.psamp <- array(dim = c(p, p, floor(nsamp/odens)))
  Y.imp <- NULL
  if (impute) {
    Y.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  }
  dimnames(C.psamp) <- list(colnames(Y), colnames(Y),
                            1:floor(nsamp/odens))
  for (ns in 1:nsamp) {
    for (j in sample(1:p)) {
      Sjc <- S[j, -j] %*% solve(S[-j, -j])
      sdj <- sqrt(S[j, j] - S[j, -j] %*% solve(S[-j,
                                                 -j]) %*% S[-j, j])
      muj <- Z[, -j] %*% t(Sjc)
      if (!plugin.marginal[j]) {
        for (r in 1:Rlevels[j]) {
          ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
          lb <- suppressWarnings(max(Z[R[, j] == r -
                                         1, j], na.rm = TRUE))
          ub <- suppressWarnings(min(Z[R[, j] == r +
                                         1, j], na.rm = TRUE))
          Z[ir, j] <- qnorm(runif(length(ir), pnorm(lb,
                                                    muj[ir], sdj), pnorm(ub, muj[ir], sdj)),
                            muj[ir], sdj)
        }
      }
      ir <- (1:n)[is.na(R[, j])]
      Z[ir, j] <- rnorm(length(ir), muj[ir], sdj)
    }
    # relocate the mean to zero
    # added by Ruifei Cui
    Z = t( (t(Z)-apply(Z,2,mean)))

    S <- solve(rwish(solve(S0 * n0 + t(Z) %*% Z), n0 +
                       n))
    if (ns%%odens == 0) {
      C <- S/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      lpc <- ldmvnorm(Z %*% diag(1/sqrt(diag(S))),
                      C)
      LPC <- c(LPC, lpc)
      C.psamp[, , ns/odens] <- C
      if (impute) {
        Y.imp.s <- Y
        for (j in 1:p) {
          Y.imp.s[is.na(Y[, j]), j] <- quantile(Y[,
                                                  j], pnorm(Z[is.na(Y[, j]), j], 0, sqrt(S[j,
                                                                                           j])), na.rm = TRUE, type = 1)
        }
        Y.imp[, , ns/odens] <- Y.imp.s
        Y.pmean <- ((ns/odens - 1)/(ns/odens)) * Y.pmean +
          (1/(ns/odens)) * Y.imp.s
      }
    }
    if (verb == TRUE & (ns%%(odens * 10)) == 0) {
      cat(round(100 * ns/nsamp), "percent done ",
          date(), "\n")
    }
  }
  G.ps <- list(C.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp,
               LPC = LPC)
  class(G.ps) <- "psgc"
  return(G.ps)
}

## Code from Cui et al. (2018).
#####################################################################################
# Goal: test conditional independence (using local effective sample size if given,
#       otherwise using sample size)
# Return: p-value of the current test
####################################################################################
gaussCItest.local <- function (x, y, S, suffStat)
{
  matN = suffStat$ESS.mat
  ##
  n = suffStat$n
  if (!(is.null(matN)))
  {
    sub.mat = matN[c(x,y,S),c(x,y,S)]
    n = mean(sub.mat[upper.tri(sub.mat)], na.rm = T)
  }
  ##
  z <- zStat(x, y, S, C = suffStat$C, n = n)
  2 * pnorm(abs(z), lower.tail = FALSE)
}
