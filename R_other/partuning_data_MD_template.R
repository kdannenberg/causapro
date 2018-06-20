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

filename = "..."
ranked = FALSE              # data SHOULD be ranked

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename, ranked = ranked,
                                         mute_all_plots = TRUE))

pc_fun_with_eval <- function_set_parameters(protein_causality, parameters =
                                              list(filename = filename, ranked = ranked,
                                                   evaluation = TRUE,
                                                   mute_all_plots = TRUE))


alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.09)
# alphas <- c(0.7)
min_pos_vars = c(0, 1) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059?? f. 0.73??


# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
tuning <- tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
                                                                                       alphas = alphas,
                                                                                       minposvars = min_pos_vars)
print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
plot.new()
plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)

pc_fun(alpha = tuning$best_alpha, min_pos_var = tuning$best_minposvar, causal_analysis = TRUE, intervention_position = "all")




