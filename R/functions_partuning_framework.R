#' Apply fucntion for different values of its parameters alpha and minposvar, select the maximum
#' of the (numerical) results
#' @param best function that computes which of the values returned by FUN for various alphas and minposvars was best
#' @param best_in_array is the value that the function best returns one of the values in its input (e.g. for min, max)?
#' if the function computes an objective function on the values in its input, it must return the $value, and the $index
#' where this value was obtained
partuning_over_alpha_and_minposvar <- function(FUN, alphas, minposvars, best = max, best_in_array = TRUE,
                                               plot = TRUE, plot_labels_as_rows_and_cols = alphas >=5,
                                               plot_no_isolated_nodes = TRUE,
                                               max_rows_in_plot = 5,...) { #"#000000", ...) {

  # results_m <- matrix(nrow = length(alphas), ncol = length(minposvars))
  # rownames(results_m) <- alphas
  # colnames(results_m) <- minposvars

  # The graph must be in result$graph because it is used later, when the FIRST computation is skipped
  result <- FUN(alphas[1], minposvars[1])
  first_value <- result$value
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
        result <- FUN(alpha, minposvar)
        results_a[as.character(alpha), as.character(minposvar), ] <- result$value
      }

      if (plot_no_isolated_nodes) {
        graph <- kernelize_graph(result$graph)
      } else {
        graph <- result$graph
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

      sub <- paste("Value:", paste(result$value, collapse = " | "), "\n")
      title(sub = sub)
      # dev.off()
    }
  }

  if (best_in_array) {
    best_value <- best(results_a)
    best <-  which(results_a == best_value, arr.ind = TRUE)
  } else {
    best_both <- best(results_a)
    best_value <- best_both$value
    best <- best_both$index
  }

  best_alpha <- rownames(results_a)[best[1]]
  best_minposvar <- colnames(results_a)[best[2]]
  return(list(best_value = best_value, best_alpha = best_alpha, best_minposvar = best_minposvar, all_results = results_a))

}

example_function <- function(alpha, minposvar) {
  m <- matrix(c(4,3,2,1), ncol = 2, byrow = TRUE)
  rownames(m) <- c(0.01, 0.1)
  colnames(m) <- c(0.001, 0.01)
  return(list(graph = graphNEL(nodes = c("a", "b")), value = m[as.character(alpha), as.character(minposvar)]))
}

# should be minimized
partune_alpha_minposvar_square_localTests_estimate <- function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- sum((results$orig$localTests$r$estimate)^2)
  result$graph <- results$pc@graph
  return(result)
}

# absolute value / deviation from zero should be minimized
partune_alpha_minposvar_mean_localTests_estimate <- function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- mean((results$orig$localTests$r$estimate))
  result$graph <- results$pc@graph
  return(result)
}

# should be minimized
# TODO forbid conflict edges to get a general partuning fct
partune_alpha_minposvar_localTests_estimate_dev_from_zero_per_edge <- function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  n_edges <- results$summary$edges$sum
  estimate <- results$orig$localTests$r$estimate
  if (!(is.null(estimate))) {
    estimate_dev_from_zero <- sum(abs(estimate))
  } else {
    estimate_dev_from_zero <- NA
  }
  result$value <- estimate_dev_from_zero / sqrt(n_edges)
  result$graph <- results$pc@graph
  return(result)
}

# should be squared or absoluted and added
# best_in_array = FALSE
partune_alpha_minposvar_sep_localTests_estimate_n_edges <- function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  result <- list()
  est <- sum(abs(results$orig$localTests$r$estimate))
  n_edges <- results$summary$edges$sum
  result$value <- c(estimate = est, n_edges = n_edges)
  result$graph <- results$pc@graph
  return(result)
}

# should be squared or absoluted and added
# best_in_array = FALSE
partune_alpha_minposvar_sep_localTests_estimate <- function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  est <- results$orig$localTests$r$estimate
  neg <- sum(est[est < 0])
  pos <- sum(est[est > 0])
  result$value <- c(neg = neg, pos = pos)
  result$graph <- results$pc@graph
  return(result)
}


# maximize for small data sets?
partune_alpha_minposvar_n_edges <- function(pc_FUN, alpha, minposvar) {
  results <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- results$summary$edges$sum
  result$graph <- results$pc@graph
  return(result)
}

partune_alpha_minposvar_avg_degree_2 <- function(pc_FUN, alpha, minposvar) {
  results <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  result <- list()
  avg_degree <- (2 * results$summary$edges$sum / results$summary$nodes$sum)
  result$value <- abs(avg_degree - 2)
  result$graph <- results$pc@graph
  return(result)
}

partune_alpha_minposvar_conflict_edges <- function(pc_FUN, alpha, minposvar) {
  results <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  result <- list()
  result$value <- results$summary$edges$conflict
  result$graph <- results$pc@graph
  return(result)
}

partune_alpha_minposvar_edges_type_distr <- function(pc_FUN, alpha, minposvar,
                                                     weight_of_conflict_edges = NULL, difference = NULL) {
  results <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  result <- list()
  result$value <- quality_of_edge_distribution(unlist(results$summary$edges)[2:4],
                                               weight_of_conflict_edges = weight_of_conflict_edges,
                                               difference = difference)
  result$graph <- results$pc@graph
  return(result)
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
                                          parameters = list(measures = list(function(x) {return(x)})))
