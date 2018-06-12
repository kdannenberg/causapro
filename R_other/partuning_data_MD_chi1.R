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

filename = "NoV_MD_chi1_AB"

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename,
                                         every_n_th_row = 10,
                                         mute_all_plots = TRUE))


partuning_avg_degree_2 <- function_set_parameters(partune_alpha_minposvar_avg_degree_2,
                                                  parameters = list(pc_FUN = pc_fun))
partuning_edge_type_distr <- function_set_parameters(partune_alpha_minposvar_edges_type_distr,
                                                  parameters = list(pc_FUN = pc_fun, weight_of_conflict_edges = 2))

# debug(partuning_Nov_NMR_2018_06_01)

alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, 0.01)
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), 0.01)
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 1, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.0001)
min_pos_vars = c(0, 0.5, 1.5, 2.5)

# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_Nov_NMR_2018_06_01_square_localTests_estimate, alphas = alphas, minposvars = min_pos_vars, best = min)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_avg_degree_2, alphas = alphas, minposvars = min_pos_vars, best = min)
opt <- partuning_over_alpha_and_minposvar(FUN = partuning_edge_type_distr, alphas = alphas, minposvars = min_pos_vars, best = max)
print(opt)



