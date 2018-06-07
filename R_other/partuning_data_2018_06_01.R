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

# NoV_NMR-Tit_B3S-with-unass
type_of_data = "NMR-Tit"
## subtype_of_data = "B3S-with-unass"
# subtype_of_data = "B3S_with-unass_1H" ## options are ..._Euclidean, ..._1H and ..._15N
##subtype_of_data = "BTS"
# subtype_of_data = c("Fuc", "B4S")

subtype_of_data = "" # Daten vom 1.6.18

pc_fun <- function_set_parameters(protein_causality_NoV, parameters =
                                     list(type_of_data = type_of_data,
                                          subtype_of_data = subtype_of_data,
                                          # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                          # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                          # causal_analysis = FALSE, #min_pos_var = 0.01,
                                          # evaluation = TRUE,
                                          # alpha = 0.1, min_pos_var = 0,
                                          # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                          # ranked = TRUE, pc_indepTest = "smc-jt",
                                          mute_all_plots = TRUE))

#TODO: generalisieren: auch für andere Daten
partuning_Nov_NMR_2018_06_01_square_localTests_estimate <- function_set_parameters(function(pc_FUN, alpha, minposvar) {
  pc_fun_ <- function_set_parameters(pc_FUN, parameters = list(evaluation = TRUE))

  results_NoV <- pc_fun_(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- sum((results_NoV$orig$localTests$r$estimate)^2)
  result$graph <- results_NoV$pc@graph
  return(result)
}, parameters = list(pc_FUN = pc_fun))


partuning_Nov_NMR_2018_06_01_n_edges <- function_set_parameters(function(pc_FUN, alpha, minposvar) {
  results_NoV <- pc_fun(alpha = alpha, min_pos_var = minposvar)

  # print(conflict_edges(results_NoV$pc@graph))
  # print(sum((results_NoV$orig$localTests$r$estimate)^2))
  result <- list()
  result$value <- results_NoV$summary$edges$sum
  result$graph <- results_NoV$pc@graph
  return(result)
}, parameters = list(pc_FUN = pc_fun))

# debug(partuning_Nov_NMR_2018_06_01)

alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 1, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.0001)
min_pos_vars = c(0, 0.06) # kleinste vorhandene Varianz: 0.059

# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_Nov_NMR_2018_06_01_square_localTests_estimate, alphas = alphas, minposvars = min_pos_vars, best = min)
opt <- partuning_over_alpha_and_minposvar(FUN = partuning_Nov_NMR_2018_06_01_n_edges, alphas= alphas, minposvars = min_pos_vars, best = max)
print(opt)



