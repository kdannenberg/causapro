determine_set_of_graphs <- function(type_of_graph_set, pc_function, ida_function, s, new, save, outpath,
                                    pc_maj_rule_conflict, pc_conservative_conflict) {
  if (type_of_graph_set == "retry") {
    if (new || !file.exists(file = paste0(outpath, "-pc-retry_results.RData")) || !file.exists(file = paste0(outpath, "-pc-retry_graphs.RData"))) {
      ## all_graphs saved
      all_results <- list()
      all_graphs <- list()
      for (i in 1:s) {
        print(paste("DURCHLAUF", i))
        # source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
        results <- pc_function(pc_solve_conflicts = FALSE, pc_u2pd = "retry", pc_maj_rule = FALSE, 
                               pc_conservative = FALSE, evaluation = FALSE, analysis = FALSE, 
                               mute_all_plots = mute_all_plots)
          # protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
                                       # pc_solve_conflicts = FALSE, pc_u2pd = "retry",
                                       # evaluation = FALSE, analysis = FALSE)
        edges <- conflict_edges(results$pc@graph)
        all_results[[i]] <- results
        all_graphs[[i]] <- results$pc@graph
        if ((edges$conflict == 0) && (edges$bidirected == 0)) {
          break
        }
      }
      if (save) {
        save(all_graphs, file = paste0(outpath, "-pc-retry_graphs.RData"))
      }
      
      
      #### save(all_results, file = paste0(outpath, "-pc-retry_results.RData"))
      
      
      #### all_results saved
      for (i in 1:s) {
        print(paste("DURCHLAUF", i))
        all_results[[i]] <- ida_function(results)
      }
      if (save) {
        save(all_results, file = paste0(outpath, "-pc-retry_results.RData"))
      }
    } else {
      load(file = paste0(outpath, "-pc-retry_graphs.RData"))
      load(file = paste0(outpath, "-pc-retry_results.RData"))
    }
  } else if (type_of_graph_set == "conflict") {
    if (new || !file.exists(file = paste0(outpath, "-all_confl_comb_results.RData")) || !file.exists(file = paste0(outpath, "-all_confl_comb_graphs.RData"))) {
      results <- pc_function(pc_solve_conflicts = TRUE, pc_u2pd = "relaxed", pc_maj_rule = pc_maj_rule_conflict, pc_conservative = pc_conservative_conflict, evaluation = FALSE, analysis = FALSE)
      # protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
      #                                             pc_solve_conflicts = TRUE, pc_u2pd = "relaxed", 
      #                                             pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
      #                                             evaluation = FALSE, analysis = FALSE)
      edges <- conflict_edges(results$pc@graph)
      print(edges)
      
      Sys.sleep(2)
      
      if (edges$conflict > 15) {
        stop("More than 15 conflict edges.")
      }
      cat("Enumerating DAGs...")
      all_graphs <- enumerate_graphs(results$pc@graph) # in Zeile 1 berechnet
      cat(" Done. \n")
      # all_results <- pblapply(all_graphs, graph_to_results, ida_function = ida_function)   ## schneller (?) # library("pbapply")
      all_results <- list()
      for (i in 1:length(all_graphs)) {
        print(paste("DURCHLAUF", i, "VON", length(all_graphs)))            ## mit Angabe des aktuellen Druchlaufs
        all_results[[i]] <- graph_to_results(all_graphs[[i]], ida_function = ida_function)
      }
      if (save) {
        save(all_graphs, file = paste0(outpath, "-all_confl_comb_graphs.RData"))
        save(all_results, file = paste0(outpath, "-all_confl_comb_results.RData"))
      }
    } else {
      load(file = paste0(outpath, "-all_confl_comb_graphs.RData"))
      load(file = paste0(outpath, "-all_confl_comb_results.RData"))
    }
  }
  
  return(list(graphs = all_graphs, results = all_results))
}

find_graphs_with_highest_int_pos <- function(all_results, obj_fct = list, dir = c("on", "of")) {
  
  quality_measure <- function(distribution, percentile) {
    statistics_of_influenced_positions(effects = distribution, percentile = percentile, interesting_positions = int_pos, print = FALSE)
  }
  
  compute_quality_measure_for_results <- function(results, percentile) {
    of_effects <- results$ida$`372`$of$scaled_effects
    of_max <- quality_measure(apply(of_effects, 1, max), percentile = percentile)
    of_min <- quality_measure(apply(of_effects, 1, min), percentile = percentile)
    
    on_max <- quality_measure(results$ida$`372`$on$scaled_effects[, 1], percentile = percentile)
    on_min <- quality_measure(results$ida$`372`$on$scaled_effects[, 2], percentile = percentile)
    
    return(list(of_max = of_max, of_min = of_min, on_max = on_max, on_min = on_min))
  }
  
  
  # int_pos_in_percentile_75 <- lapply(all_results, compute_quality_measure_for_results)
  int_pos_in_percentile_75 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.75)))
  int_pos_in_percentile_85 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.85)))
  int_pos_in_percentile_95 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.95)))
  
  if (dir == "on") {
    int_pos_in_percentile_75 <- lapply(int_pos_in_percentile_75, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
    int_pos_in_percentile_85 <- lapply(int_pos_in_percentile_85, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
    int_pos_in_percentile_95 <- lapply(int_pos_in_percentile_95, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
  } else if (dir == "of") {
    int_pos_in_percentile_75 <- lapply(int_pos_in_percentile_75, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
    int_pos_in_percentile_85 <- lapply(int_pos_in_percentile_85, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
    int_pos_in_percentile_95 <- lapply(int_pos_in_percentile_95, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
  }
  
  total_number_of_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unique(unlist(x)))) # fast immer 11! sonst 10!!
  
  total_number_of_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unique(unlist(x)))) # oft 11, tw. bis zu 6
  
  # int_pos_in_percentile_95 <- lapply(all_results, compute_quality_measure_for_results)
  total_number_of_int_pos_in_percetile_95 <- sapply(int_pos_in_percentile_95, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_95 <- sapply(int_pos_in_percentile_95, function(x) length(unique(unlist(x)))) # 5-8
  
  max_pos_75 <- which(total_number_of_int_pos_in_percetile_75 == max(total_number_of_int_pos_in_percetile_75))
  max_pos_85 <- which(total_number_of_int_pos_in_percetile_85 == max(total_number_of_int_pos_in_percetile_85))
  max_pos_95 <- which(total_number_of_int_pos_in_percetile_95 == max(total_number_of_int_pos_in_percetile_95))
  
  max_diff_pos_75 <- which(total_number_of_diff_int_pos_in_percetile_75 == max(total_number_of_diff_int_pos_in_percetile_75))
  max_diff_pos_85 <- which(total_number_of_diff_int_pos_in_percetile_85 == max(total_number_of_diff_int_pos_in_percetile_85))
  max_diff_pos_95 <- which(total_number_of_diff_int_pos_in_percetile_95 == max(total_number_of_diff_int_pos_in_percetile_95))
  
  
  print("Graphs with highest effects for interesting positions:")
  print(paste("max_pos_75", paste(max_pos_75, collapse = ", "), sep = ": "))
  print(paste("max_pos_85", paste(max_pos_85, collapse = ", "), sep = ": "))
  print(paste("max_pos_95", paste(max_pos_95, collapse = ", "), sep = ": "))
  
  print(paste("max_diff_pos_75", paste(max_diff_pos_75, collapse = ", "), sep = ": "))
  print(paste("max_diff_pos_85", paste(max_diff_pos_85, collapse = ", "), sep = ": "))
  print(paste("max_diff_pos_95", paste(max_diff_pos_95, collapse = ", "), sep = ": "))
  
  
  # cat(paste0("max_pos_75: ", paste(max_pos_75, collapse = ", "), "\n"))
  # cat(paste0("max_pos_85: ", paste(max_pos_85, collapse = ", "), "\n"))
  # cat(paste0("max_pos_95: ", paste(max_pos_95, collapse = ", "), "\n"))
  # 
  # cat(paste0("max_diff_pos_75: ", paste(max_pos_75, collapse = ", "), "\n"))
  # cat(paste0("max_diff_pos_85: ", paste(max_pos_85, collapse = ", "), "\n"))
  # cat(paste0("max_diff_pos_95: ", paste(max_pos_95, collapse = ", "), "\n"))
  
  return(obj_fct(max_pos_75 = max_pos_75, max_pos_85 = max_pos_85, max_pos_95 = max_pos_95, 
                 max_diff_pos_75 = max_diff_pos_75, max_diff_pos_85 = max_diff_pos_85, max_diff_pos_95 = max_diff_pos_95))
}


element_in_most_of_the_6_sets <- function(max_pos_75, max_pos_85, max_pos_95, 
                                          max_diff_pos_75, max_diff_pos_85, max_diff_pos_95) {
  all <- c(max_pos_75, max_pos_85, max_pos_95, 
           max_diff_pos_75, max_diff_pos_85, max_diff_pos_95)
  return(as.integer(names(which(table(all) == max(table(all))))))
}






# select for a results object the mean results on and of position 372, respectively and return both as a list
# mean of min and max?
mean_effects_min_max <- function(results, weight_effects_on_by, scaled_effects = FALSE) {
  on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  
  if (scaled_effects) {
    of_effects <- results$ida$`372`$of$scaled_effects
  } else {
    of_effects <- results$ida$`372`$of$effects
  }
  of_max <- apply(of_effects, 1, max)
  of_min <- apply(of_effects, 1, min)
  
  if (scaled_effects) {
    on_max <- results$ida$`372`[[on]]$scaled_effects[, 1]
    on_min <- results$ida$`372`[[on]]$scaled_effects[, 2]
  } else {
    on_max <- results$ida$`372`[[on]]$effects[, 1]
    on_min <- results$ida$`372`[[on]]$effects[, 2]
  }
  
  ret_list <- list()
  ret_list$of <- apply(cbind(of_max, of_min), 1, mean)
  ret_list[[on]] <- apply(cbind(on_max, on_min), 1, mean)
  return(ret_list)
}

# sum all effects:
# should rather be devided by 100, thus mean
compute_over_all_graphs <- function(all_results, weight_effects_on_by, use_scaled_effects_for_sum = FALSE, scale_in_the_end = FALSE, 
                                    function_over_all_graphs = "mean", direction = c("on", "of"), 
                                    scale_effects_on = "372" %in% rownames(all_results[[1]]$ida$`372`$of$effects)) {
  
  min_max_mean_effects_on_of <- lapply(all_results, mean_effects_min_max, weight_effects_on_by = weight_effects_on_by, scaled_effects = use_scaled_effects_for_sum)
  min_max_mean_effects_of <- do.call(cbind, (lapply(min_max_mean_effects_on_of, function(list) return(list$of))))
  effect_over_all_graphs_of <- apply(min_max_mean_effects_of, 1, function_over_all_graphs)
  on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  # if (grepl("mean", weight_effects_on_by)) { # (weight_effects_on_by == "mean_abs_effect") {
  #   on <- "on-rel-to-mean"
  # } else if (grepl("median", weight_effects_on_by)) { # if (weight_effects_on_by == "median_abs_effect") {
  #   on <- "on-rel-to-median"
  # } else if (weight_effects_on_by == "var" || weight_effects_on_by == "vars") {
  #   on <- "on-rel-to-var"
  # }
  min_max_mean_effects_on <- do.call(cbind, (lapply(min_max_mean_effects_on_of, function(list) return(list[[on]]))))
  effect_over_all_graphs_on <- apply(min_max_mean_effects_on, 1, function_over_all_graphs) 
  if (scale_effects_on) {
    effect_over_all_graphs_on <- effect_over_all_graphs_on / effect_over_all_graphs_on["372"]
  }
  
  effects_over_all_graphs_on_of <- list(overAllGraphs_of = effect_over_all_graphs_of, overAllGraphs_on = effect_over_all_graphs_on)
  
  return(effects_over_all_graphs_on_of) #, all_mean = all_mean_effects))
}

display_effects <- function(effects_on_of, direction, int_pos, scale_in_the_end, weight_effects_on_by, 
                            function_over_all_graphs, ida_percentile = ida_percentile, caption, main_caption,
                            print = TRUE, plot = TRUE, for_combined_plot = FALSE) {
  
  if (plot) {
    if (!for_combined_plot && !is.null(main_caption) && !missing(main_caption)) {
      oma <- c( 0, 0, length(main_caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
      # par(mfrow=c(lines, columns), oma = oma) 
      par(oma = oma)
    }
    
    if ("on" %in% direction && "of" %in% direction) {
      par(mfrow = c(2,1))
    }
    
    if ("of" %in% direction) {
      if (missing(caption) || is.null(caption)) {
        caption <- paste0(function_over_all_graphs, " over all graphs of effects of position 372")
      }
      plot_effects(effects_on_of$overAllGraphs_of, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption)
    }
    if ("on" %in% direction) {
      on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
      if (missing(caption) || is.null(caption)) {
        caption <- paste0(function_over_all_graphs, " over all graphs of effects ", on, " position 372")
      }
      
      plot_effects(effects_on_of$overAllGraphs_on, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption)
    }
    
    if (! (("on" %in% direction) ||  ("of" %in% direction))) {
      effects_matrix <- rbind(effects_on_of$overAllGraphs_of, effects_on_of$overAllGraphs_on)
      on_of_mean_effects <- apply(effects_matrix, 2, get(direction))

      on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
      if (missing(caption) || is.null(caption)) {
        caption <- paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372")
      }
      
      plot_effects(on_of_mean_effects, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption)
      
      # scaled_effects_for_coloring_on_of_mean <- scale_effects(as.matrix(on_of_mean_effects), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
      # colors_by_effect <- color_by_effect(scaled_effects_for_coloring_on_of_mean, int_pos, mode = "#FFFFFF")
      # if (!scale_in_the_end) {
      #   barplot(on_of_mean_effects, 
      #           main = paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372"), 
      #           col = colors_by_effect, las = 2)
      # } else {
      #   barplot(as.vector(scaled_effects_for_coloring_on_of_mean), 
      #           main = paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372"), 
      #           col = colors_by_effect, las = 2, 
      #           names.arg = rownames(scaled_effects_for_coloring_on_of_mean))
      # }
      
      effects_on_of[[paste0("overAllGraphs_", direction, "_on_of")]] <- on_of_mean_effects
    }
    if (!for_combined_plot && !missing(main_caption) && !is.null(main_caption)) {
      title(main_caption, outer = TRUE)
    }
  }
  
  if (print) {
    if ("of" %in% direction) {
      print("SUM EFFECTS OF:")
      stat_sum_on <- statistics_of_influenced_positions(effects_on_of$sum_of, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
    }
    if ("on" %in% direction) {
      print("SUM EFFECTS ON:")
      stat_sum_of <- statistics_of_influenced_positions(effects_on_of$sum_on, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
    }
    if (! (("on" %in% direction) ||  ("of" %in% direction))) {
      print("SUM EFFECTS MEAN(ON, OF):")
      # effects_matrix <- do.call(rbind, effects_over_all_graphs_on_of)
      # on_of_mean_effects <- apply(effects_matrix, 2, get(direction))
      stat_sum_of <- statistics_of_influenced_positions(effects_on_of[[paste0("overAllGraphs_", direction, "_on_of")]], percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
    }
  }
  
  return(effects_on_of)
}

plot_effects <- function(effects, int_pos, scale_in_the_end, caption) {
  scaled_effects_for_coloring <- scale_effects(as.matrix(effects), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
  colors_by_effect <- color_by_effect(scaled_effects_for_coloring, int_pos, mode = "#FFFFFF")
  if (!scale_in_the_end) {
    barplot(effects, 
            main = caption, 
            col = colors_by_effect, las = 2)
  } else {
    barplot(as.vector(scaled_effects_for_coloring), 
            main = caption, 
            col = colors_by_effect, las = 2, 
            names.arg = rownames(scaled_effects_for_coloring))
  }
}

deviation_from_mean <- function(all_results, dir, weight_effects_on_by, plot_graphs) {
  # for (dir in direction) {
  # how much do the distinct distributions deviate from the mean?
  # OF
  
  if (dir == "on") {
    dir <- pastes("on", weight_effects_on_by, sep = "-rel-to-")  
  }
  
  all_mean_effects_all <- lapply(all_results, mean_effects_min_max, weight_effects_on_by = weight_effects_on_by)
  all_mean_effects <- do.call(cbind, (lapply(all_mean_effects_all, function(list) return(list[[dir]]))))
  mean_effect <- apply(all_mean_effects, 1, sum) / length(all_results)
  
  dev_from_mean <- apply(all_mean_effects, 2, function(distribution) {return(distribution - mean_effect)})
  
  graphics.off()
  par(mfrow = c(5,5))
  
  for (i in plot_graphs) {
    barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, 
            main = paste("Graph", i, "\n effects", dir))
  }
  # 25 und 19 sehen GLEICH (gut aus)
  # Graphen aber nicht:
  # graphics.off()
  # par(mfrow = c(1,2))
  # plot(all_graphs[[25]])
  # plot(all_graphs[[19]])
  
  # for (i in 26:50) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }
  
  # for (i in 51:75) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }
  
  # for (i in 76:100) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }
  # 98 sieht gut aus! (min_var = 0.01)
  
  # apply(dev_from_mean, 2, function(x) barplot(x, col = colors_for_barplot))  # geht nur für alle auf einmal, dafür ist aber die Grafik zu klein
  
  # TODO: same_procedure for on!
  # btw TODO: on ohne Gewichtung der einzelnen Positionen mit dem durchschnittlichen Effekt
  # }
  # } else if (plot == "on") {
  #   # ON
  #   substract_mean_on <- function(distribution) {
  #     return(distribution - (sum_effect_on/length(all_results)))
  #   }
  #   dev_from_mean_on <- apply(all_mean_effects_on, 2, substract_mean_on)
  #   
  #   graphics.off()
  #   par(mfrow = c(5,5))
  #   
  #   # for (i in 1:25) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #   
  #   # for (i in 26:50) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #   
  #   
  #   # for (i in 51:75) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #   # min_var = 0.001: Graph 66
  #   
  #   for (i in 76:100) {
  #     barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   }
  #   
  #   # apply(dev_from_mean_on, 2, function(x) barplot(x, col = colors_for_barplot))  # geht nur für alle auf einmal, dafür ist aber die Grafik zu klein
  # }
}

graph_to_results <- function(graph, ida_function) {
  results <- list()
  # if (is.null(graph)) {
  #   
  # }
  pc <- new("pcAlgo", graph = graph)
  results$pc <- pc
  # results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
  #                               protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
  #                               amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
  #                               pymol_bg_color = "grey",
  #                               barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
  results <- ida_function(results)
}

# neg_effects: abs -> absolut value, discard -> drop
quality_of_effects_distibution <- function(effects, int_pos, neg_effects = "drop", function_over_effects = mean, perturbed_position = "372") {
  if (neg_effects == "abs") {
    effects <- abs(effects)
  } else if (neg_effects == "discard") {
    effects <- effects[!effects < 0]
  }
  effects <- effects[!names(effects) %in% perturbed_position]
  mean_effect_int <- function_over_effects(effects[names(effects) %in% int_pos])
  mean_effect_other <- function_over_effects(effects[!names(effects) %in% int_pos])
  
  return(mean_effect_int / mean_effect_other)
}
