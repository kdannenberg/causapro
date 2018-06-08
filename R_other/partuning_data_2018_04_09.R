source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)
# debug(partuning_over_alpha_and_minposvar)


# pc_fun <- function_set_parameters(protein_causality_NoV, parameters =
#                                     list(filename = "Nov_NMR-Mut_apo",
#                                          mute_all_plots = TRUE))

pc_fun <- function_set_parameters(protein_causality_NoV, parameters =
                                    list(filename = "Nov_NMR-Tit",
                                         mute_all_plots = TRUE))

#TODO: generalisieren: auch für andere Daten
partuning_square_localTests_estimate <- function_set_parameters(function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results_NoV <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- sum((results_NoV$orig$localTests$r$estimate)^2)
  result$graph <- results_NoV$pc@graph
  return(result)
}, parameters = list(pc_FUN = pc_fun))


partuning_n_edges <- function_set_parameters(function(pc_FUN, alpha, minposvar) {
  results_NoV <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- results_NoV$summary$edges$sum
  result$graph <- results_NoV$pc@graph
  return(result)
}, parameters = list(pc_FUN = pc_fun))


partuning_avg_degree_2 <- partuning_n_edges <- function_set_parameters(function(pc_FUN, alpha, minposvar) {
  results_NoV <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  result <- list()
  avg_degree <- (2 * results_NoV$summary$edges$sum / results_NoV$summary$nodes$sum)
  result$value <- abs(avg_degree - 2)
  result$graph <- results_NoV$pc@graph
  return(result)
}, parameters = list(pc_FUN = pc_fun))

# debug(partuning_Nov_NMR_2018_06_01)

# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 1, 0.1))
alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.0001)
min_pos_vars = c(0, 0.06) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059

# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_square_localTests_estimate, alphas = alphas, minposvars = min_pos_vars, best = min)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_n_edges, alphas = alphas, minposvars = min_pos_vars, best = max)
opt <- partuning_over_alpha_and_minposvar(FUN = partuning_avg_degree_2, alphas = alphas, minposvars = min_pos_vars, best = min)

print(opt)



