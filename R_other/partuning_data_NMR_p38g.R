source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(div_cols_by_sd)
# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)
# debug(partuning_over_alpha_and_minposvar)
# debug(cluster_pairwise_effects)


filename = "p38g_NMR-Mut_inact"
protein = "p38g"

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename,
                                         # pre_fun_on_data = "div_cols_by_sd",
                                         mute_all_plots = TRUE,
                                         type_of_variables = "continuous"))

pc_fun_with_eval <- function_set_parameters(pc_fun, parameters =
                                              list(evaluation = TRUE,
                                                   max_edges_for_evaluation = 100))

# angucken und setzen:
var_cutoff = 0
# pc_fun(graph_computation = FALSE, show_variance_cutoff_plot = TRUE, min_pos_var = var_cutoff)
# stop()


# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas <- c(0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 0.9)
# alphas = c(0.001, 0.01, 0.05, 0.1)
alphas <- c(1e-10, 1e-5, 0.09)
# alphas <- c(0.7)
# min_pos_vars = c(0, var_cutoff) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059?? f. 0.73??
min_pos_vars = c(0) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059?? f. 0.73??


# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
# tuning <- tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
#                                                                                        protein = protein,
#                                                                                        # coloring = "S3-pie",
#                                                                                        alphas = alphas,
#                                                                                        minposvars = min_pos_vars)
# print(tuning)
# plot.new()
# plot_partuning(tuning$all_values)

tunes <- list()

# EDGE SCORE
# abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score <-
tune_alpha_mpv_edge_score_w4 <- tune_alpha_mpv_edge_score_factory(4)
tunes[[4]] <- edge_score_w4 <-
  tune_alpha_mpv_edge_score_w4(pc_fun = pc_fun,
                               alphas = alphas,
                               minposvars = min_pos_vars,
                               protein = protein)

pc_fun(alpha = tunes[[4]]$best_alpha, min_pos_var = tunes[[4]]$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all", effects_pre_cluster_fun = "cor",
       effects_cluster_cut_height_h = 3)

protein_causality(filename = file_name, alpha = tuning$best_alpha,
                  min_pos_var = tuning$best_minposvar, causal_analysis = TRUE, effects_cluster_cut_height_h = 3)


