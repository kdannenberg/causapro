library(stringr) # for str_sub

analyse_edge_types_by_alpha <- function(edge_types, numbers_of_nodes, plot_labels_as_rows_and_cols = TRUE,
                                        plot_logscale_alpha = FALSE,
                                        ylim_for_plots = c(0,200), # or NULL
                                        conflict_edge_weight = 1,
                                        print = TRUE, plot = TRUE) {
  best_alphas <- list()
  
  if (print) {
    graphics.off()
    if (plot_labels_as_rows_and_cols) {
      # par(mfrow = c(length(measures) + 1, length(min_pos_vars) + 1))
      par(mfcol = c(length(min_pos_vars) + 1, length(measures) + 1))
    } else {
      # par(mfrow = c(length(measures), length(min_pos_vars)))
      par(mfcol = c(length(min_pos_vars), length(measures)))
    }
    # for (measure_type_sub in c(# "DDS", 
    #   "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")) {
    
    # plot row labels
    if (plot_labels_as_rows_and_cols) {
      plot.new()
      legend('center', legend = c(names(edge_types[[1]][[1]][[1]][[1]]), "sum"), col = c('red','green','orange', 'black'), lty = c(1,1,1,1), lwd = 1 )
      for (min_pos_var in min_pos_vars) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("min_pos_var \n=", min_pos_var), 
             cex = 1.6, col = "black")
      }
    }
  }
  
  for (measure_type_sub in measures) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data = "-"
    }
    best_alphas[[measure_type_sub]] <- list()
    if (plot) {
      # plot row labels
      if (plot_labels_as_rows_and_cols) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, measure_type_sub, 
             cex = 1.6, col = "black")
      }
    }
    
    for (min_pos_var in min_pos_vars) {
      print(paste(measure_type_sub, " - min_pos_var =", min_pos_var))
      edge_types_by_alpha <- matrix(unlist(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]]), ncol = 3, byrow = TRUE)
      
      rownames(edge_types_by_alpha) <- names(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]])
      colnames(edge_types_by_alpha) <- names(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]][[1]])
      
      sum_of_edges <- apply(edge_types_by_alpha, 1, sum)
      edge_types_by_alpha <- cbind(edge_types_by_alpha, sum = sum_of_edges)
      quality <- apply(edge_types_by_alpha, 1, quality_of_edge_distribution, weight_of_conflict_edges = conflict_edge_weight)
      best_alpha <- as.numeric(names(quality)[which(quality == max(quality))])
      best_alphas[[measure_type_sub]][[as.character(min_pos_var)]] <- best_alpha
      print(paste("Best alpha, according to quality:", best_alpha))
      
      if (print) {
        print(cbind(edge_types_by_alpha, quality = quality))
        cat("\n")
      }
      
      if (plot) {
        if (plot_logscale_alpha) {
          log <- "x"
        } else { 
          log <- ""
        }
        
        n_nodes <- numbers_of_nodes[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]]
        n_edges <- 1/2 * n_nodes * (n_nodes - 1)
        if (length(ylim_for_plots) == 1) {
          current_ylim_for_plots <- c(0, n_edges / ylim_for_plots)
        } else {
          current_ylim_for_plots <- ylim_for_plots
        }
        
        matplot(edge_types_by_alpha, x = rownames(edge_types_by_alpha), xlab = "alpha", ylab = "# edges",
                type = "l", lty = c(1,1,1,1), col = c('red','green','orange', 'black'), ylim = current_ylim_for_plots, log = log)
        abline(h = 15, col = "red", lty = 2)
        # abline(h = n_edges, col = "black", lty = 1)
        abline(v = 0.01, col = "black", lty = 2)
        abline(v = best_alpha, col = "black", lty = 1)
        axis(1, at = best_alpha, labels = best_alpha, mgp = c(10, 2, 0))
        if (!plot_labels_as_rows_and_cols) {
          legend('top', legend = colnames(edge_types_by_alpha), col = c('red','green','orange', 'black'), lty = c(1,1,1,1), lwd = 1 )
          caption <- paste0(measure_type_sub, ", min_pos_var = ", min_pos_var)
          title(caption)
        }
      }
    }
  }
  return(best_alphas)
}

find_best_alphas <- function(all_scores) {
  result <- list()
  for (measure_type_sub in names(all_scores)) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    # type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    # if (grepl("-", measure_type_sub)) {
    #   subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    # } else {
    #   subtype_of_data <- ""
    #   subtype_of_data_list <- "-"
    # }
    result[[measure_type_sub]] <- list()
    for (min_pos_var in as.numeric(names(all_scores[[measure_type_sub]][[1]]))) {
      values <- list()
      for (alpha in as.numeric(names(all_scores[[measure_type_sub]]))) {
        values[[as.character(alpha)]] <- all_scores[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]]
      }
      result[[measure_type_sub]][[as.character(min_pos_var)]] <- names(which(unlist(values) == max(unlist(values))))
    }
  }
  return(result)
}


# distinct_alphas is a list containing for each measure a list (or vector) with a value alpha for each min_pos_var.
# For all measures in names(distinct_alphas) and all min_pos_vars in names(distinct_alphas$<measure>), the corresponding alpha is chosen, effects computed
# and displayed in grid
# if measures / min_pos_vars are given, only those are used (but need to be existent in distinct_alphas)
effects_for_distinct_alphas <- function(distinct_alphas, with_graphs = FALSE, with_effects = TRUE, for_all_alphas = FALSE,
                                        measures, min_pos_vars, cols_for_measures = TRUE, protein_causality_function) {
  effects <- list()
  
  if (with_graphs && with_effects) {
    if (missing(min_pos_vars)) {
      lines_of_plot <- 2 * length(names(distinct_alphas[[1]]))
    } else {
      lines_of_plot <- 2 * length(min_pos_vars)
    }
  } else {
    if (missing(min_pos_vars)) {
      lines_of_plot <- length(names(distinct_alphas[[1]]))
    } else {
      lines_of_plot <- length(min_pos_vars)
    }
  }
  if (for_all_alphas) {
    cols_of_plot = 5  # just sth
  } else {
    if (missing(measures)) {
      if (cols_for_measures) {
        cols_of_plot <- length(names(distinct_alphas)) # Anzahl measures
      } else {
        cols_of_plot <- sum(unlist(lapply(distinct_alphas, length))) # Anzahl measures und min_pos_vars
      }
    }  else {
      cols_of_plot <- length(measures)
    }
  }
  
  par(mfcol = c(lines_of_plot, cols_of_plot))
  
  
  if (missing(measures)) {
    measures <- names(distinct_alphas)
  }
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    
    if (missing(min_pos_vars)) {
      current_min_pos_vars <- as.numeric(names(distinct_alphas[[measure_type_sub]]))
    }
    for (min_pos_var in current_min_pos_vars) {
      alphas <- distinct_alphas[[measure_type_sub]][[as.character(min_pos_var)]]
      
      if (!for_all_alphas && length(alphas) > 1) {
        # alpha <- alpha[1]
        alphas <- alphas[length(alphas)] 
      }
      for (alpha in alphas) {
        if (missing(protein_causality_function)) {
          protein_causality_function = get(paste0("protein_causality_", measure))
        }
        
        pc_function_ <- function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                                                              plot_clusters = FALSE, alpha = alpha, min_pos_var = min_pos_var,
                                                                                              for_combined_plot = TRUE))
        pc_function_effects <- function_set_parameters(pc_function_, parameters = list(mute_all_plots = TRUE))
        if (with_effects) {
          effects[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]][[as.character(alpha)]] <-
            analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                  direction = "mean", measure = measure, pc_function = pc_function_effects, alpha = alpha, min_pos_var = min_pos_var,
                                  for_combined_plot = TRUE, scale_in_the_end = FALSE, new = FALSE, save = TRUE)
          }
        if (with_graphs) {
          pc_function_(pc_maj_rule = TRUE) # mute_all_plots = FALSE
          # protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data, alpha = alpha, 
          # min_pos_var = min_pos_var, for_combined_plot = TRUE, pc_maj_rule = TRUE)
        } 
      }
    }
  }
}

# which_effects gives the slot of the list of effects (e.g. )
apply_to_all_effects_in_nested_list <- function(all_effects, FUN, which_effects = "overAllGraphs_mean_on_of") {
  result <- list()
  for (measure_type_sub in names(all_effects)) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    # type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    # if (grepl("-", measure_type_sub)) {
    #   subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    # } else {
    #   subtype_of_data <- ""
    #   subtype_of_data_list <- "-"
    # }
    result[[measure_type_sub]] <- list()
    for (alpha in as.numeric(names(all_effects[[measure_type_sub]]))) {
      result[[measure_type_sub]][[as.character(alpha)]] <- list()
      for (min_pos_var in as.numeric(names(all_effects[[measure_type_sub]][[as.character(alpha)]]))) {
        if (!is.null(all_effects[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]][[which_effects]])) {
          result[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]] <- 
            FUN(effects = all_effects[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]][[which_effects]])
        }
      }
      if (length(result[[measure_type_sub]][[as.character(alpha)]]) == 0) {
        result[[measure_type_sub]][[as.character(alpha)]] <- NULL
      }
    }
  }
  return(result)
}

# TODO: sicherstellen, dass oma nur gesetzt wird, wenn auch eine main_caption (title) geprintet wird

analyse_graphs_for_alphas_and_minposvars <- function(measure = "S", type_of_data = "DDS", subtype_of_data = "", 
                                                     alphas = c(0.001, 0.005, 0.01, 0.05, 0.1), min_pos_vars = c(0.0001, 0.001, 0.01),
                                                     protein_causality_function = get(paste0("protein_causality_", measure)),
                                                     # TODO: DAS FUNKTIONIERT NICHT!! DIE PARAMETER, DIE NCIHT ÜBERGEBEN WERDEN, WERDEN NCIHT JETZT SCHON BELEGT!
                                                     # ALLES WIRD ERST BEI AUFRUF AUSGEWERTET
                                                     # pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
                                                     #   protein_causality_function(type_of_data = eval(type_of_data), subtype_of_data = eval(subtype_of_data), min_pos_var = min_pos_var, alpha = alpha,
                                                     #                              pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                                                     #                              evaluation = evaluation, analysis = analysis, mute_all_plots = TRUE)
                                                     # pc_function = f_protein_causality_pc_parameters_eval_analysis(measure = measure, type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                     #### min_pos_var = min_pos_var, alpha = alpha, 
                                                     # mute_all_plots = TRUE)
                                                     pc_function = function_set_parameters(protein_causality_function, 
                                                          parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                          mute_all_plots = TRUE)),
                                                     new = FALSE, save = TRUE
) {
  
  # if (missing(pc_function)) {
  #   pc_function <- function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
  #                                                                         mute_all_plots = TRUE))
  # }
  
  effects <- list()
  
  graphics.off()
  oma <- c( 2, 0, 2, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
  if (length(alphas) <= 4) {
    par(mfrow = c(length(alphas), length(min_pos_vars)), oma = oma)
  } else {
    par(mfrow = c(4, length(min_pos_vars)), oma = oma)
  }
  
  for (alpha in alphas) {
    effects[[as.character(alpha)]] <- list()
    for (min_pos_var in min_pos_vars) {
      # for (measure in c("G", "S")) {
      # if (
      #     (
      #       (type_of_data == "DDDG" && subtype_of_data == "5") &&
      #       ((alpha == 0.05 && min_pos_var == 0.001) #||  # More than 15 conflict edges. (22)
      #       || (alpha == 0.1 && min_pos_var == 0.001)  # More than 15 conflict edges.
      #       || (alpha == 0.05 && min_pos_var == 0.0001) # More than 15 conflict edges.
      #       || (alpha == 0.1 && min_pos_var == 0.0001) # More than 15 conflict edges.
      #       )
      #      # ) || (
      #      #   (type_of_data == "DDDG" && subtype_of_data == "10") &&
      #      #   (
      #      #   (alpha == 0.005 && min_pos_var == 0.001) ##  "Fehler in wgt.unique[x, ] : Indizierung außerhalb der Grenzen"
      #      #   # || (alpha == 0.05 && min_pos_var == 0.0001)
      #      #   # || (alpha == 0.1 && min_pos_var == 0.0001)
      #      #   )
      #      ) || (
      #        (type_of_data == "DDG" && subtype_of_data == "10") &&
      #        (
      #          (alpha == 0.1 && min_pos_var == 0.0001) ##  15 conflict edges
      #          || (alpha == 0.1 && min_pos_var == 0.001) ##  13 conflict edges
      #        )
      #      ) || (
      #        (measure == "S") &&
      #        (alpha == 0.1 && min_pos_var == 0.0001)
      #      ) || (
      #        (type_of_data == "DDG" && subtype_of_data == "5") &&
      #        (
      #          (alpha == 0.05 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
      #       || (alpha == 0.05 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
      #       || (alpha == 0.05 && min_pos_var == 0.01) #  More than 15 conflict edges. (18)
      #       || (alpha == 0.1 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
      #       || (alpha == 0.1 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
      #       || (alpha == 0.1 && min_pos_var == 0.01) #  More than 15 conflict edges. (21)
      #        )
      #      )
      #     ) {
      #   plot.new()
      # } else {
      cat("\n")
      cat(paste0("alpha = ", alpha, ", min_pos_var = ", min_pos_var, "\n"))
      pc_function_ <- function_set_parameters(pc_function, parameters = list(alpha = alpha, min_pos_var = min_pos_var))
      effects[[as.character(alpha)]][[as.character(min_pos_var)]] <- analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                                                           direction = "mean", measure = measure, pc_function = pc_function_, alpha = alpha, min_pos_var = min_pos_var,
                                                                                           for_combined_plot = TRUE, scale_in_the_end = FALSE, new = new, save = save)
      
      # }
      # }
    }
  }
  
  # title(main = paste0(type_of_data, "-", subtype_of_data),
  #       sub = "Mean causal effects over all conflict graphs and over the effects on and of position 372", outer = TRUE)
  title(main = paste0("Mean causal effects over all conflict graphs and over the effects on and of position 372",
                      "\n", type_of_data, "-", subtype_of_data), outer = TRUE)
  return(effects)
}




# QUALITY MEASURES
# BY EDGES
quality_of_edge_distribution <- function(edge_types, weight_of_conflict_edges = 1) {
  if (edge_types["conflict"] > 15) {
    return(0)                 # infeasible
  } else if (sum(edge_types) == 0) {  # no edges
    return(0)
  } else {
    return(edge_types["directed"]/(edge_types["undirected"] + weight_of_conflict_edges * edge_types["conflict"]))
  }
}

# BY EFFECTS
quality_of_effects_classifier <- function(effects, int_pos, percentile = 1 - (length(int_pos) / length(effects)), perturbed_position) {
  most_influenced_pos <- get_most_influenced_positions(effects, percentile = percentile)
  most_influenced_pos_interesting <- intersect(int_pos, most_influenced_pos)
  fp = setdiff(most_influenced_pos, most_influenced_pos_interesting)
  fn = setdiff(int_pos, most_influenced_pos_interesting)
  if (length(fp) != length(fn)) {
    # Das kann passieren, wenn es weniger als lenght(int_pos) nicht-null Effekte gibt!
    warning("Not as many false positives as false negatives in quality_of_effects_classifier.") 
  }
  # return((mean(c(length(fp), length(fn)))) / length(setdiff(int_pos, perturbed_position)))
  return((max(c(length(fp), length(fn)))) / length(setdiff(int_pos, perturbed_position)))
}

score_for_effects <- function(effects, int_pos, perturbed_position = 372, effect_quality_height, false_pos_neg, 
                              percentile = 1 - (length(int_pos) / length(effects))) {
  if (missing(effect_quality_height)) {
    effect_quality_height <- quality_of_effects_distibution(effects = effects, int_pos = int_pos, perturbed_position = perturbed_position)
  }
  if (missing(false_pos_neg)) {
    quality_of_effects_classifier <- quality_of_effects_classifier(effects, int_pos, percentile, perturbed_position)
  } else {
    quality_of_effects_classifier <- false_pos_neg / length(setdiff(int_pos, perturbed_position))
  }
  
  s_height <- score_of_effects_height(effect_quality_height)
  s_classifier <- score_of_effects_classifier(quality_of_effects_classifier)
  return(s_height + s_classifier)
}

# give either effect_quality_height directly, or effects, int_pos and perturbed_position,
# so it can be computed internally
# previously: quality_score 
score_of_effects_height <- function(effect_quality_height, effects, int_pos, perturbed_position, range_to_six = TRUE) {
  if (missing(effect_quality_height)) {
    effect_quality_height <- quality_of_effects_distibution(effects = effects, int_pos = int_pos, perturbed_position = perturbed_position)
  }
  if (!range_to_six) {
    if (effect_quality_height < 0.5) {
      return(0)
    } else if (effect_quality_height < 1) {
      return(1)
    } else if (effect_quality_height < 1.5) {
      return(2)
    } else if (effect_quality_height < 2) {
      return(3)
    } else if (effect_quality_height < 2.5) {
      return(4)
    } else if (effect_quality_height < 3) {
      return(5)
    } else {
      return(6)
    }
  } else {
    if (effect_quality_height == 0) {
      return(0)
    } else if (effect_quality_height < 1) {
      return(1)
    } else if (effect_quality_height < 1.5) {
      return(2)
    } else if (effect_quality_height < 2.5) {
      return(3)
    } else if (effect_quality_height < 3.5) {
      return(4)
    } else if (effect_quality_height < 5) {
      return(5)
    } else if (effect_quality_height < 7) {
      return(6)
    # } else if (effect_quality_height < 6.5) {
    #   return(6)
    } else if (effect_quality_height < 10) {
      return(7)
    } else if (effect_quality_height == Inf) {
      return(9)
    } else {
      return(8)
    }
  }
}

# give either false_pos_neg_fraction directly, or effects, int_pos and perturbed_position,
# so it can be computed internally
# previously: false_pos_neg_score
score_of_effects_classifier <- function(false_pos_neg_fraction, effects, int_pos, perturbed_position) {
  if (missing(false_pos_neg_fraction)) {
    false_pos_neg_fraction <- quality_of_effects_classifier(effects, int_pos = int_pos, perturbed_position = perturbed_position) 
  }
  # if (false_pos_neg_fraction >= 1) {
  #   return(0)
  # } else if (false_pos_neg_fraction >= 5/6) {
  #   return(1)
  # } else if (false_pos_neg_fraction >= 4/6) {
  #   return(2)
  # } else if (false_pos_neg_fraction >= 3/6) {
  #   return(3)
  # } else if (false_pos_neg_fraction >= 2/6) {
  #   return(4)
  # } else if (false_pos_neg_fraction >= 1/6) {
  #   return(5)
  # } else {
  #   return(6)
  # }
  if (false_pos_neg_fraction == 0) {
    return(9)
  } else if (false_pos_neg_fraction <= 1/9) {
    return(8)
  } else if (false_pos_neg_fraction <= 2/9) {
    return(7)
  } else if (false_pos_neg_fraction <= 3/9) {
    return(6)
  } else if (false_pos_neg_fraction <= 4/9) {
    return(5)
  } else if (false_pos_neg_fraction <= 5/9) {
    return(4)
  } else if (false_pos_neg_fraction <= 6/9) {
    return(3)
  } else if (false_pos_neg_fraction <= 7/9) {
    return(2)
  } else if (false_pos_neg_fraction <= 8/9) {
    return(1)
  } else {
    return(0)
  }
}

# neg_effects: abs -> absolut value, discard -> drop
# option "sep" fehlt hier gegenüber scale_effects
# TODO: rename: quality_of_effects_height
quality_of_effects_distibution <- function(effects, int_pos, neg_effects = "discard", function_over_effects = mean, perturbed_position = "372") {
  if (neg_effects == "abs") {
    effects <- abs(effects)
  } else if (neg_effects == "discard" || neg_effects == "drop") {
    effects <- effects[!effects < 0]
  } else if (neg_effects != "") {
    warning("Unknown treatment of negative effects in quality_of_effects_distibution")
  }
  effects <- effects[!names(effects) %in% perturbed_position]
  mean_effect_int <- function_over_effects(effects[names(effects) %in% int_pos])
  mean_effect_other <- function_over_effects(effects[!names(effects) %in% int_pos])
  
  quality <- mean_effect_int / mean_effect_other
  
  if (is.na(quality)) {
    quality = 0
  }  
  return(quality)
}