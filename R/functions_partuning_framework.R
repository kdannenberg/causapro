#' Apply fucntion for different values of its parameters alpha and minposvar, select the maximum
#' of the (numerical) results
#' @param objective_fun function that computes which of the values returned by FUN for various alphas and minposvars was best
#' @param best_in_array is the value that the function best returns one of the values in its input (e.g. for min, max)?
#' if the function computes an objective function on the values in its input, it must return the $value, and the $index
#' where this value was obtained
partuning_over_alpha_and_minposvar <- function(alphas, minposvars, pc_fun, objective_fun = max, best_in_array = TRUE,
                                               results_to_value_fun,
                                               plot = TRUE, plot_labels_as_rows_and_cols = (length(alphas) >= 5),
                                               plot_no_isolated_nodes = TRUE,
                                               max_rows_in_plot = 5,...) { #"#000000", ...) {

  # results_m <- matrix(nrow = length(alphas), ncol = length(minposvars))
  # rownames(results_m) <- alphas
  # colnames(results_m) <- minposvars

  # The graph must be in results$pc@graph because it is used later, when the FIRST computation is skipped
  results <- pc_fun(alpha = alphas[1], min_pos_var = minposvars[1])
  first_value <- results_to_value_fun(results)
  dim_3 <- length(first_value)

  results_a <- array(dim = c(length(alphas),length(min_pos_vars), dim_3))
  dimnames(results_a)[[1]] <- as.list(alphas)
  dimnames(results_a)[[2]] <- as.list(minposvars)
  dimnames(results_a)[[3]] <- as.list(names(first_value))

  results_a[1,1,] <- first_value



  # par(mfrow = c(length(alphas),length(minposvars)))

  if (plot) {
    plot.new()
    oma <- c( 2, 0, 2, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
    if (length(alphas) <= max_rows_in_plot) {
      par(mfrow = c(length(alphas), length(minposvars)), oma = oma)
    } else {
      par(mfrow = c(max_rows_in_plot, length(minposvars)), oma = oma)
    }

    # plot row labels
    if (plot_labels_as_rows_and_cols) {
      #inc row and cols
      par(mfrow = par()$mfrow + 1)
      plot.new()
      # legend('center', legend = c(names(edge_types[[1]][[1]][[1]][[1]]), "sum"), col = c('red','green','orange', 'black'), lty = c(1,1,1,1), lwd = 1 )
      for (minposvar in minposvars) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("min_pos_var \n=", minposvar),
             cex = 1.6, col = "black")
      }
    }
  }

  for (alpha in alphas) {
    if (plot) {
      # plot row labels
      if (plot_labels_as_rows_and_cols) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        # text(x = 0.5, y = 0.5, paste("alpha \n=", alpha), cex = 1.6, col = "black")
        text(x = 0.5, y = 0.5, paste("alpha =", alpha), cex = 1.6, col = "black")
      }
    }
    for (minposvar in minposvars) {
      if (any(is.na(results_a[as.character(alpha), as.character(minposvar), ]))) {
        # this must be the first pair results_a[1, 1,], because that is the one filled in earlier.
        # This means that the result (including the graph) is already computed (and has not been overwritten).
        results <- pc_fun(alpha = alpha, min_pos_var = minposvar)
        results_a[as.character(alpha), as.character(minposvar), ] <- results_to_value_fun(results)
      }

      if (plot_no_isolated_nodes) {
        graph <- kernelize_graph(results$pc@graph)
      } else {
        graph <- results$pc@graph
      }
      plot_structure(graph = graph, caption = paste0("alpha = ", alpha, ", minposvar = ", minposvar), ...)

      lines_needed_for_subcaption <- 1
      mgp <- par()$mgp  # get current value
      mgp[1] <- 2 + lines_needed_for_subcaption
      mar <- par()$mar # get current value
      mar[1] <- 3 + lines_needed_for_subcaption
      par(mar = mar, mgp = mgp)

      # if (caption_as_subcaption) {
      #   caption = caption
      #   main_caption = NULL  # soll missing sein
      # } else {
      #   main_caption = caption
      #   caption = NULL
      # }

      # plot_text(paste0("alpha = ", alpha, ", minposvar = ", minposvar), caption = paste0("alpha = ", alpha, ", minposvar = ", minposvar), ...)

      sub <- paste("Value:", paste(results_to_value_fun(results), collapse = " | "), "\n")
      title(sub = sub)
      # dev.off()
    }
  }

  if (best_in_array) {
    best_value <- objective_fun(results_a)
    best_index <-  which(results_a == best_value, arr.ind = TRUE)
  } else {
    best_both <- objective_fun(results_a)
    best_value <- best_both$value
    best_index <- best_both$index
  }

  best_alpha <- rownames(results_a)[best_index[1]]
  best_minposvar <- colnames(results_a)[best_index[2]]
  return(list(best_value = best_value, best_alpha = best_alpha, best_minposvar = best_minposvar, all_results = results_a))
}


#### Functions to extract relevant characteristics from the results of the pc_fun
# EDGES
# maximize for small data sets?
results_to_n_edges <- function(results) {
  value <- results$summary$edges$sum
  return(value)
}

results_to_edges_type_distr <- function(results) {
  value <- unlist(results$summary$edges)[2:4]
  return(value)
}

results_to_conflict_edges <- function(results) {
  value <- results$summary$edges$conflict
  return(value)
}

results_to_avg_degree <- function(results) {
  avg_degree <- (2 * results$summary$edges$sum / results$summary$nodes$sum)
  return(avg_degree)
}



# EMPRIRCAL corELATION
# should be squared or absoluted and added
# best_in_array = FALSE
results_to_sep_empir_cor <- function(results) {
  est <- results$orig$localTests$r$estimate
  neg <- sum(est[est < 0])
  pos <- sum(est[est > 0])
  value <- c(neg = neg, pos = pos)
  return(value)
}

# absolute value / deviation from zero should be minimized
results_to_mean_empir_cor <- function(results) {
  value <- mean((results$orig$localTests$r$estimate))
  return(value)
}

# should be minimized
# requires results to contain results of structure evalutation (evaluation = TRUE in call of pc)
results_to_square_empir_cor <- function(results) {
  value <- sum((results$orig$localTests$r$estimate)^2)
  return(value)
}


# EMPRICAL corELATION AND EDGES
# should be squared or absoluted and added
# best_in_array = FALSE
results_to_sep_empir_cor_n_edges <- function(results) {
  est <- sum(abs(results$orig$localTests$r$estimate))
  n_edges <- results$summary$edges$sum
  value <- c(estimate = est, n_edges = n_edges)
  return(value)
}

results_to_sep_empir_cor_edge_types <- function(results) {
  est <- sum(abs(results$orig$localTests$r$estimate))
  edge_types <- unlist(results$summary$edges)[2:4]
  value <- c(estimate = est, edge_types = edge_types)
  return(value)
}


obj_function_est_per_edges_times_mean_estimate_no_conflict <- function(vector) {
  est <- vector[estimate]
  edge_types <- vector[edge_types]


}


# INCLUDING FUNCTION OF THE CHARACTERISTICS
# should be minimized
# TODO forbid conflict edges to get a general partuning fct
# preivionsly: results_to_empir_cor_dev_from_zero_per_edge
results_to_compute_empir_cor_dev_from_zero_per_edge <- function(results) {
  n_edges <- results$summary$edges$sum
  estimate <- results$orig$localTests$r$estimate
  if (!(is.null(estimate))) {
    estimate_dev_from_zero <- sum(abs(estimate))
  } else {
    estimate_dev_from_zero <- NA
  }
  value <- estimate_dev_from_zero / sqrt(n_edges)
  return(value)
}


#### Functions to determine the best among the results of the pc_fun,
#### by means of relevant characteristics extracted before
basic_obj_function <- function(array, fun_over_components_of_value, fun_over_results_of_other_fun) {
  v <- apply(array, c(1,2), fun_over_components_of_value)
  best_v <- fun_over_results_of_other_fun(v)
  return(list(value = best_v, index = (which(v ==  best_v, arr.ind = TRUE))))
}

obj_fun_min_dev_from <- function(dev_from_what) {
  dev_fun <- function(vector) {
    return(abs(vector - dev_from_what))
  }
  return(function_set_parameters(basic_obj_function,
         parameters = list(fun_over_components_of_value = dev_fun,
                          fun_over_results_of_other_fun = min)))
}

obj_fun_edge_distr <- function(weight_of_conflict_edges = NULL, difference = NULL) {
  quality_edge_distr <- function_set_parameters(quality_of_edge_distribution,
                        parameters = list(weight_of_conflict_edges = weight_of_conflict_edges,
                                          difference = difference))
  return(function_set_parameters(basic_obj_function,
         parameters = list(fun_over_components_of_value = quality_edge_distr,
                           fun_over_results_of_other_fun = max)))
}



obj_fun_est_per_edge <- function_set_parameters(basic_obj_function,
                        parameters = list(fun_over_components_of_value =  function(v) {v[1]/v[2]},
                                          fun_over_results_of_other_fun = min))

obj_fun_est_per_squred_edges <- function_set_parameters(basic_obj_function,
                        parameters = list(fun_over_components_of_value =  function(v) {v[1]/(v[2]^2)},
                                          fun_over_results_of_other_fun = min))






# partuning functions with only pc_fun, alphas and min_pos_vars missing
tune_alpha_mpv_dev_of_edges_from_2 <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
            parameters = list(results_to_value_fun = results_to_avg_degree,
                              objective_fun = obj_fun_min_dev_from(2),
                              best_in_array = FALSE))


tune_alpha_mpv_estimate_per_edges <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_value_sep_empir_cor_n_edges,
                                            objective_fun = obj_fun_est_per_edge,
                                            best_in_array = FALSE))

tune_alpha_mpv_mean_empir_cor <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_mean_empir_cor,
                                            objective_fun = obj_fun_min_dev_from(0),
                                            best_in_array = FALSE))

# tune_alpha_mpv_est_per_edges_times_mean_estimate_no_conflict <-
#   function_set_parameters(partuning_over_alpha_and_minposvar,
#                           parameters = list(results_to_value_fun = ,
#                                             objective_fun = min,
#                                             best_in_array = FALSE))



plot_different_measures_of_all_results <- function(all_results, measures, print = FALSE) {
  par(mfrow = c(1, length(measures)))
  for (measure in measures) {
    matrix <- apply(all_results, c(1,2), measure)
    if (print) {
      print(matrix)
    }
    matplot(x = rownames(matrix), y = matrix, type = "l")
  }
}

plot_partuning <- function_set_parameters(plot_different_measures_of_all_results,
                                          parameters = list(measures = list(function(x) {return(x)})))
