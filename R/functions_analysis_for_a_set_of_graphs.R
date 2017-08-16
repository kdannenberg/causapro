determine_set_of_graphs <- function(type_of_graph_set, s, new, save, outpath,
                                    pc_maj_rule_conflict, pc_conservative_conflict) {
  if (type_of_graph_set == "retry") {
    if (new || !file.exists(file = paste0(outpath, "-pc-retry_results.RData")) || !file.exists(file = paste0(outpath, "-pc-retry_graphs.RData"))) {
      ## all_graphs saved
      all_results <- list()
      all_graphs <- list()
      for (i in 1:s) {
        print(paste("DURCHLAUF", i))
        # source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
        results <- pc_function(pc_solve_conflicts = FALSE, pc_u2pd = "retry", pc_maj_rule = FALSE, pc_conservative = FALSE, evaluation = FALSE, analysis = FALSE)
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
      if (edges$conflict > 15) {
        stop("More than 15 conflict edges.")
      }
      all_graphs <- enumerate_graphs(results$pc@graph) # in Zeile 1 berechnet
      # all_results <- pblapply(all_graphs, graph_to_results, ida_function = ida_function)   ## schneller (?) # library("pbapply")
      all_results <- list()
      for (i in 1:length(all_graphs)) {
        print(paste("DURCHLAUF", i))                          ## mit Angabe des aktuellen Druchlaufs
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

find_graphs_with_highest_int_pos <- function(all_results, obj_fct = list) {
  
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
  total_number_of_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unique(unlist(x)))) # fast immer 11! sonst 10!!
  
  int_pos_in_percentile_85 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.85)))
  total_number_of_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unique(unlist(x)))) # oft 11, tw. bis zu 6
  
  # int_pos_in_percentile_95 <- lapply(all_results, compute_quality_measure_for_results)
  int_pos_in_percentile_95 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.95)))
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
mean_effects <- function(results, weight_effects_on_by, scaled_effects = FALSE) {
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
sum_all_effects <- function(all_results, weight_effects_on_by, use_scaled_effects_for_sum = FALSE, print = TRUE) {
  all_mean_effects <- lapply(all_results, mean_effects, weight_effects_on_by = weight_effects_on_by, scaled_effects = use_scaled_effects_for_sum)
  all_mean_effects_of <- do.call(cbind, (lapply(all_mean_effects, function(list) return(list$of))))
  sum_effect_of <- apply(all_mean_effects_of, 1, sum)
  on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  # if (grepl("mean", weight_effects_on_by)) { # (weight_effects_on_by == "mean_abs_effect") {
  #   on <- "on-rel-to-mean"
  # } else if (grepl("median", weight_effects_on_by)) { # if (weight_effects_on_by == "median_abs_effect") {
  #   on <- "on-rel-to-median"
  # } else if (weight_effects_on_by == "var" || weight_effects_on_by == "vars") {
  #   on <- "on-rel-to-var"
  # }
  all_mean_effects_on <- do.call(cbind, (lapply(all_mean_effects, function(list) return(list[[on]]))))
  sum_effect_on <- apply(all_mean_effects_on, 1, sum)   
  
  if (print) {
    par(mfrow = c(2,1))
    scaled_effects_for_coloring_of <- scale_effects(as.matrix(sum_effect_of), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
    colors_by_effect <- color_by_effect(scaled_effects_for_coloring_of, int_pos, mode = "#FFFFFF")
    if (use_scaled_effects_for_sum) {
      barplot(sum_effect_of, main = "sum of effects of position 372", col = colors_by_effect)
    } else {
      barplot(as.vector(scaled_effects_for_coloring_of), main = "sum of effects of position 372", col = colors_by_effect)
    }
    
    scaled_effects_for_coloring_on <- scale_effects(as.matrix(sum_effect_on), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
    colors_by_effect <- color_by_effect(scaled_effects_for_coloring_on, int_pos, mode = "#FFFFFF")
    if (use_scaled_effects_for_sum) {
      barplot(sum_effect_on, main = paste("sum of effects", on, "position 372"), col = colors_by_effect)
    } else {
      barplot(as.vector(scaled_effects_for_coloring_on), main = paste("sum of effects", on, "position 372"), col = colors_by_effect)
    }
  }
  # if (plot == "of") {
  
  return(list(sum_of = sum_effect_of, sum_on = sum_effect_on))#, all_mean = all_mean_effects))
}

deviation_from_mean <- function(all_results, dir, weight_effects_on_by, plot_graphs) {
  # for (dir in direction) {
  # how much do the distinct distributions deviate from the mean?
  # OF
  
  if (dir == "on") {
    dir <- pastes("on", weight_effects_on_by, sep = "-rel-to-")  
  }
  
  all_mean_effects_all <- lapply(all_results, mean_effects, weight_effects_on_by = weight_effects_on_by)
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
  pc <- new("pcAlgo", graph = graph)
  results$pc <- pc
  # results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
  #                               protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
  #                               amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
  #                               pymol_bg_color = "grey",
  #                               barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
  results <- ida_function(results)
}

