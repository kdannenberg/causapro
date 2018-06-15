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


partuning_avg_degree_2 <- function_set_parameters(partune_alpha_minposvar_avg_degree_2,
                                                  parameters = list(pc_FUN = pc_fun))
partuning_edge_type_distr <- function_set_parameters(partune_alpha_minposvar_edges_type_distr,
                                                     parameters = list(pc_FUN = pc_fun,
                                                                       weight_of_conflict_edges = 2,
                                                                       difference = TRUE))

partuning_localTests_estimate <- function_set_parameters(partune_alpha_minposvar_sep_localTests_estimate,
                                                         parameters = list(pc_FUN = pc_fun))
localTests_obj_fct <- function(array) {
  v <- apply(array, c(1,2), function(x) {sum(abs(x))})
  return(list(value = min(v), index = (which(v ==  min(v), arr.ind = TRUE))))
}



# debug(partuning_Nov_NMR_2018_06_01)

alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 1, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.0001)
min_pos_vars = c(0, 0.08) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059

# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_square_localTests_estimate, alphas = alphas, minposvars = min_pos_vars, best = min)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_n_edges, alphas = alphas, minposvars = min_pos_vars, best = max)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_avg_degree_2, alphas = alphas, minposvars = min_pos_vars, best = min)
opt <- partuning_over_alpha_and_minposvar(FUN = partuning_localTests_estimate,
                                          alphas = alphas, minposvars = min_pos_vars,
                                          best = localTests_obj_fct, best_in_array = FALSE)
print(opt)



