#' Apply fucntion for different values of its parameters alpha and minposvar, select the maximum
#' of the (numerical) results
#' @param objective_fun function that computes which of the values returned by FUN for various alphas and minposvars was best
#' @param best_in_array is the value that the function best returns one of the values in its input (e.g. for min, max)?
#' if the function computes an objective function on the values in its input, it must return the $value, and the $index
#' where this value was obtained
partuning_over_alpha_and_minposvar <- function(alphas, minposvars, pc_fun,
                                               # objective_fun = max, best_in_array = TRUE,
                                               results_to_value_fun,
                                               measure_fun,
                                               objective_fun,
                                               objective_fun_returns_indices = FALSE,
                                               protein,
                                               plot = TRUE, plot_graphs = TRUE,
                                               plot_no_isolated_nodes = TRUE,
                                               max_rows_in_plot = 5,
                                               plot_labels_as_rows_and_cols = (length(alphas) <= max_rows_in_plot),
                                               # TODO: implementieren
                                               stop_when_more_than_c_conflict_edges = 15,
                                               ...) { #"#000000", ...) {

  # results_m <- matrix(nrow = length(alphas), ncol = length(minposvars))
  # rownames(results_m) <- alphas
  # colnames(results_m) <- minposvars

  # TODO Marcel
  # wenn !plot, hier versuchen, results_a zu laden

  stop <- rep(FALSE, length(minposvars))
  names(stop) <- minposvars

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
    if (all(stop)) {
      break()
    }
    for (minposvar in minposvars) {
      if (any(is.na(results_a[as.character(alpha), as.character(minposvar), ])) && !stop[as.character(minposvar)]) {
        # this must be the first pair results_a[1, 1,], because that is the one filled in earlier.
        # This means that the result (including the graph) is already computed (and has not been overwritten).
        results <- pc_fun(alpha = alpha, min_pos_var = minposvar)
        results_a[as.character(alpha), as.character(minposvar), ] <- results_to_value_fun(results)
        ## TODO: HIER stop setzen, wenn mehr als 15 konfliktkanten
        if (conflict_edges(results$pc@graph)$conflict > 15) {
          stop[as.character(minposvar)] <- TRUE
        }
      } else if (stop[as.character(minposvar)]) {
        results_a[as.character(alpha), as.character(minposvar), ] <- NA
      }

      if (plot) {
        if (plot_graphs) {
          if (plot_no_isolated_nodes) {
            graph <- kernelize_graph(results$pc@graph)
          } else {
            graph <- results$pc@graph
          }
          if (!missing(protein)) {
            plot_graph(graph = graph, protein = protein, coloring = "auto", caption = paste0("alpha = ", alpha, ", minposvar = ", minposvar), ...)
          } else {
            plot_structure(graph = graph, caption = paste0("alpha = ", alpha, ", minposvar = ", minposvar), ...)
          }
        } else {
          plot()  # TODO Marcel: hier den Funktionswert plotten
        }
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
  }


  # TODO Marcel:
  # results_a abspeichern (nur, wenn alles da ist? dann alles laden, um teile daraus zu bekommen?)
  # dafür compute_if_not_existent nutzen

  # basic_obj_function <- function(array, measure_fun, fun_over_results_of_other_fun) {
  results_m <- apply(results_a, c(1,2), measure_fun)

  if (objective_fun_returns_indices) {
    best_index <- objective_fun(results_m)
    best_value <- results_m[best_index[1], best_index[2]]
  } else {
    best_value <- objective_fun(results_m)
    best_index <- (which(results_m ==  best_value, arr.ind = TRUE))
  }

    # return(list(value = best_v, index = (which(v ==  best_v, arr.ind = TRUE))))
  # }

  # if (best_in_array) {
  #   best_value <- objective_fun(results_a)
  #   best_index <-  which(results_a == best_value, arr.ind = TRUE)
  # } else {
  #   best_both <- objective_fun(results_a)
  #   best_value <- best_both$value
  #   best_index <- best_both$index
  # }

  best_alpha <- rownames(results_a)[best_index[1]]
  best_minposvar <- colnames(results_a)[best_index[2]]
  return(list(best_value = best_value, best_alpha = best_alpha, best_minposvar = best_minposvar,
              all_results = results_a, all_values = results_m))
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

results_to_sep_abs_and_mean_empir_cor_edge_types <- function(results) {
  est_abs <- sum(abs(results$orig$localTests$r$estimate))
  est_mean <- mean((results$orig$localTests$r$estimate))
  edge_types <- unlist(results$summary$edges)[2:4]
  value <- c(abs_estimate = est_abs, mean_estimate = est_mean, edge_types = edge_types)
  return(value)
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
# basic_obj_function <- function(array, measure_fun, fun_over_results_of_other_fun) {
#   v <- apply(array, c(1,2), measure_fun)
#   best_v <- fun_over_results_of_other_fun(v)
#   return(list(value = best_v, index = (which(v ==  best_v, arr.ind = TRUE))))
# }


charfun_dev_from <- function(dev_from_what) {
  dev_fun <- function(ist, soll) {
    return(abs(ist - soll))
  }
  return(function_set_parameters(dev_fun,
         parameters = list(soll = dev_from_what)))
}

measure_est_per_edges_times_mean_estimate_no_conflict <- function(vector) {
  est_abs <- vector["abs_estimate"]
  est_mean <- vector["mean_estimate"]
  # edge_types <- vector["edge_types"]
  edge_types <- vector[grepl("edge_types", names(vector))]
  abs_est_per_n_edges <- est_abs / sum(edge_types)
  if (edge_types[grepl("conflict", names(edge_types))] > 0) {
    return(Inf)
  } else {
    return(abs_est_per_n_edges * est_mean)
  }
}

measure_est_per_edges_times_mean_estimate_times_conflict_plus_1 <- function(vector) {
  est_abs <- vector["abs_estimate"]
  est_mean <- vector["mean_estimate"]
  # edge_types <- vector["edge_types"]
  edge_types <- vector[grepl("edge_types", names(vector))]
  abs_est_per_n_edges <- est_abs / sum(edge_types)
  return(abs_est_per_n_edges * est_mean * (edge_types[grepl("conflict", names(edge_types))] + 1))
}

measure_est_per_edges_times_mean_estimate_times_edge_score <- function(vector) {
  est_abs <- vector["abs_estimate"]
  est_mean <- abs(vector["mean_estimate"])
  # edge_types <- vector["edge_types"]
  edge_types <- vector[grepl("edge_types", names(vector))]
  names(edge_types) <- sub("edge_types.", "", names(edge_types))
  q <- quality_of_edge_distribution(edge_types, conflict_to_dir_ratio = TRUE, weight_of_conflict_edges = 10)
  abs_est_per_n_edges <- est_abs / sum(edge_types)
  return(abs_est_per_n_edges * est_mean * q)
  # return(abs_est_per_n_edges)
  # return(est_mean)
  # return(q)
}

measure_edge_score <- function(vector, weight_of_conflict_edges) {
  edge_types <- vector#[grepl("edge_types", names(vector))]
  # names(edge_types) <- sub("edge_types.", "", names(edge_types))
  q <- quality_of_edge_distribution(edge_types, conflict_to_dir_ratio = TRUE, weight_of_conflict_edges = weight_of_conflict_edges)
  return(q)
}

measure_quality_of_edge_distr <- function(weight_of_conflict_edges = NULL, difference = NULL) {
  quality_edge_distr <- function_set_parameters(quality_of_edge_distribution,
                                                parameters = list(weight_of_conflict_edges = weight_of_conflict_edges,
                                                                  difference = difference))
  return(function_set_parameters(basic_obj_function,
                                 parameters = list(measure_fun = quality_edge_distr,
                                                   fun_over_results_of_other_fun = max)))
}



# partuning functions with only pc_fun, alphas and min_pos_vars missing
na_min <- function_set_parameters(min, parameters = list(na.rm = TRUE))

# remember to set:
# objective_fun_returns_indices = TRUE
min_pos_gradient <- function(values) {
  gradient <- apply(values, 2, function(col) {diff(col)})
  rownames(gradient) <- rownames(values)[1:nrow(values)-1]

  values[gradient<=0] <- NA

  best <- which(values == max(values, na.rm = TRUE), arr.ind = TRUE)

  if (is.matrix(best)) {   # more than one optimum
    best <- best[1,]
  }

  return(best)
}

# remember to set:
# objective_fun_returns_indices = TRUE
max_gradient_value_below_t <- function(values, t) {
  gradient <- apply(values, 2, function(col) {diff(col)})
  rownames(gradient) <- rownames(values)[1:nrow(values)-1]

  # drop last line (for which the gradient is unknown)
  values_gradient_format <- values[1:(nrow(gradient)), ]

  gradient[values_gradient_format > t] <- NA

  best <- which(gradient == max(gradient, na.rm = TRUE), arr.ind = TRUE)

  if (is.matrix(best)) {   # more than one optimum
    best <- best[1,]
  }

  return(best)
}

#TODO
first_change_of_sign <- function(values, t) {
}

tune_alpha_mpv_dev_of_edges_from_2 <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
            parameters = list(results_to_value_fun = results_to_avg_degree,
                              measure_fun = charfun_dev_from(2),
                              objective_fun = na_min))

# min value of conflict edges, with an increase directly after
# # remember to set:
# objective_fun_returns_indices = TRUE
tune_alpha_mpv_conflict_edges_local_min <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
            parameters = list(results_to_value_fun = results_to_conflict_edges,
                              measure_fun = identity,
                              objective_fun = min_pos_gradient))


# max gradient in conflict edges, while less than t
# remember to set:
# objective_fun_returns_indices = TRUE
tune_alpha_mpv_max_gradient_of_conflict_edges_below_t_factory <- function(t_max_number_of_conflict_edges) {
  max_gradient_value_below_6 <- function_set_parameters(max_gradient_value_below_t,
                                                        parameters = list(t = t_max_number_of_conflict_edges))
  return(function_set_parameters(partuning_over_alpha_and_minposvar,
                                 parameters = list(results_to_value_fun = results_to_conflict_edges,
                                                   measure_fun = identity,
                                                   objective_fun = max_gradient_value_below_6)))
}



tune_alpha_mpv_estimate_per_edges <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_sep_empir_cor_n_edges,
                                            measure_fun = function(v) {v[1]/v[2]},
                                            objective_fun = na_min))

tune_alpha_mpv_estimate_per_square_edges <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_sep_empir_cor_n_edges,
                                            measure_fun = function(v) {v[1]/(v[2]^2)},
                                            objective_fun = na_min))

tune_alpha_mpv_mean_empir_cor <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_mean_empir_cor,
                                            measure_fun = charfun_dev_from(0),
                                            objective_fun = na_min))

tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_sep_abs_and_mean_empir_cor_edge_types,
                                            measure_fun = measure_est_per_edges_times_mean_estimate_no_conflict,
                                            objective_fun = na_min))

tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict_plus_1 <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_sep_abs_and_mean_empir_cor_edge_types,
                                            measure_fun = measure_est_per_edges_times_mean_estimate_times_conflict_plus_1,
                                            objective_fun = na_min))

tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score <-
  function_set_parameters(partuning_over_alpha_and_minposvar,
                          parameters = list(results_to_value_fun = results_to_sep_abs_and_mean_empir_cor_edge_types,
                                            measure_fun = measure_est_per_edges_times_mean_estimate_times_edge_score,
                                            objective_fun = na_min))

tune_alpha_mpv_edge_score_factory <- function(weight_of_conflict_edges) {
  measure_edge_score_w <- function_set_parameters(measure_edge_score,
        parameters = list(weight_of_conflict_edges = weight_of_conflict_edges))
  return(function_set_parameters(partuning_over_alpha_and_minposvar,
               parameters = list(results_to_value_fun = results_to_edges_type_distr,
                            measure_fun = measure_edge_score_w,
                            objective_fun = na_min)))
}



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
                                          parameters = list(measures = list(identity)))
