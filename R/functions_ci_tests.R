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



ci_test_pc_chi_square <- ci_test_pc("x2")

ci_test_cor <- ci_test_pc("cor")
