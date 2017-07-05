source("configuration_code.R")

source("functions_causal_effects.R")

## all_graphs saved
# all_results <- list()
# all_graphs <- list()
# for (i in 18:100) {
#   source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
#   edges <- conflict_edges(results$pc@graph)
#   all_results[[i]] <- results
#   all_graphs[[i]] <- results$pc@graph
#   if ((edges$conflict == 0) && (edges$bidirected == 0)) {
#     break
#   }
# }

conflicts_sorted <- sapply(all_graphs, conflict_edges)
conflicts_sorted <- matrix(unlist(conflicts_sorted), nrow = 3)
rownames(conflicts_sorted) = c("conflict", "directed", "bidirected")
colnames(conflicts_sorted) <- as.character(1:100) # Graph number
conflicts_sorted <- conflicts_sorted[, order(conflicts_sorted["bidirected",])]

print("Different types of edges in the different graphs (columns)")
print(conflicts_sorted)

## all_results saved
# for (i in 1:100) {
#   all_results[[i]] <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", relatve_effects_on_pos = TRUE,
#                                 protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
#                                 amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
#                                 pymol_bg_color = "grey",
#                                 barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)
# }

# of_max <- quality_measure(apply(of_effects, 1, max))

find_graphs_with_highest_int_pos <- function() {
  protein = "PDZ"
  int_pos <- interesting_positions(protein = protein, coloring = "")
  
  quality_measure <- function(distribution) {
    statistics_of_influenced_positions(effects = distribution, percentile = 0.85, interesting_positions = int_pos, print = FALSE)
  }
  
  compute_quality_measure_for_results <- function(results) {
    of_effects <- results$ida$`372`$of$scaled_effects
    of_max <- quality_measure(apply(of_effects, 1, max))
    of_min <- quality_measure(apply(of_effects, 1, min))
    
    on_max <- quality_measure(results$ida$`372`$on$scaled_effects[, 1])
    on_min <- quality_measure(results$ida$`372`$on$scaled_effects[, 2])
    
    return(list(of_max = of_max, of_min = of_min, on_max = on_max, on_min = on_min))
  }
  
  # qualities_75 <- lapply(all_results, compute_quality_measure_for_results)
  # total_number_of_int_pos_in_percetile_75 <- sapply(qualities_75, function(x) length(unlist(x)))
  # total_number_of_diff_int_pos_in_percetile_75 <- sapply(qualities_75, function(x) length(unique(unlist(x)))) # fast immer 11! sonst 10!!
  
  qualities_85 <- lapply(all_results, compute_quality_measure_for_results)
  total_number_of_int_pos_in_percetile_85 <- sapply(qualities_85, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_85 <- sapply(qualities_85, function(x) length(unique(unlist(x)))) # oft 11, tw. bis zu 6
  
  # qualities_95 <- lapply(all_results, compute_quality_measure_for_results)
  # total_number_of_int_pos_in_percetile_95 <- sapply(qualities_95, function(x) length(unlist(x)))
  # total_number_of_diff_int_pos_in_percetile_95 <- sapply(qualities_95, function(x) length(unique(unlist(x)))) # 5-8
  
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
}

cat("\n")
find_graphs_with_highest_int_pos()

# i <- 76 # ist in max_pos_85 und _95 und in max_diff_pos_75 und _95
# causal_effects_ida(data = data, perturbated_position = "372", direction = "both", relatve_effects_on_pos = TRUE,
#                                                    protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
#                                                    amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
#                                                    pymol_bg_color = "grey",
#                                                    barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)
  

# sum all effects:
sum_all_effects <- function(all_results) {
  mean_effects <- function(results) {
    of_effects <- results$ida$`372`$of$scaled_effects
    of_max <- apply(of_effects, 1, max)
    of_min <- apply(of_effects, 1, min)
  
    on_max <- results$ida$`372`$on$scaled_effects[, 1]
    on_min <- results$ida$`372`$on$scaled_effects[, 2]
    
    return(list(of = apply(cbind(of_max, of_min), 1, mean), on = apply(cbind(on_max, on_min), 1, mean)))
  }
  
  all_mean_effects <- lapply(all_results, mean_effects)
  all_mean_effects_of <- do.call(cbind, (lapply(all_mean_effects, function(list) return(list$of))))
  sum_effect_of <- apply(all_mean_effects_of, 1, sum)
  all_mean_effects_on <- do.call(cbind, (lapply(all_mean_effects, function(list) return(list$on))))
  sum_effect_on <- apply(all_mean_effects_on, 1, sum)   
  
  par(mfrow = c(2,1))
  
  scaled_effects_for_coloring_of <- scale_effects(as.matrix(sum_effect_of), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
  colors_by_effect <- color_by_effect(scaled_effects_for_coloring_of, int_pos, mode = "#FFFFFF")
  barplot(sum_effect_of, main = "sum of effects of position 372", col = colors_by_effect)
  
  scaled_effects_for_coloring_on <- scale_effects(as.matrix(sum_effect_on), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
  colors_by_effect <- color_by_effect(scaled_effects_for_coloring_on, int_pos, mode = "#FFFFFF")
  barplot(sum_effect_on, main = "sum of effects on position 372", col = colors_by_effect)
  
  return(list(sum_effect_of, sum_effect_on))
}

sum_all_effects(all_results)


test_sum_all_effects <- function() {
  scaled_results_1_of <- matrix(c(1,0,0.5,1,1,1,2,2,-0.5), byrow = TRUE, ncol = 3)
  scaled_results_1_on <- matrix(c(1,0,1,1,2,1), byrow = TRUE, ncol = 2)
  rownames(scaled_results_1_of) <- rownames(scaled_results_1_on) <- as.character(1:3)
  fake_result_1 <- list(ida=list(`372` = list(of = list(scaled_effects = scaled_results_1_of), 
                                              on = list(scaled_effects = scaled_results_1_on))))
  
  scaled_results_2_of <- matrix(c(2,3,1), byrow = TRUE, ncol = 1)
  scaled_results_2_on <- matrix(c(0,0,1,1,2,2), byrow = TRUE, ncol = 2)
  rownames(scaled_results_2_of) <- rownames(scaled_results_2_on) <- as.character(1:3)
  fake_result_2 <- list(ida=list(`372` = list(of = list(scaled_effects = scaled_results_2_of), 
                                              on = list(scaled_effects = scaled_results_2_on))))
  fake_results <- list(fake_result_1, fake_result_2)
  
  sum_all_effects(fake_results)
}

# test_sum_all_effects()
