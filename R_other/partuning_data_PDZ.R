source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(read_data)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)
# debug(partuning_over_alpha_and_minposvar)

file_name = "PDZ_DDDG_0_10-0_90"
protein = "PDZ"

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = file_name,
                                         # pre_fun_on_data = "div_cols_by_sd",
                                         # cor_cov_fun="none",
                                         mute_all_plots = TRUE))

pc_fun_with_eval <- function_set_parameters(pc_fun, parameters =
                                              list(evaluation = TRUE,
                                                   max_edges_for_evaluation = 50))

# # angucken und setzen:
# var_cutoff = 0.005  # DDG
# var_cutoff = 0.002 # DDDG
# pc_fun(graph_computation = FALSE, show_variance_cutoff_plot = TRUE, min_pos_var = 0.005)
# stop()

# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.0005, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
alphas <- c(1e-20, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 0.0001, 0.0005, 0.001, 0.01, 0.05, 0.1, 0.2)
# alphas <- c(1e-40, 1e-30, 1e-20, 1e-10, 1e-5, 0.0001, 0.0005, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1)) # for SVD15
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas <- c(0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.09)
# alphas <- c(0.7)
# min_pos_vars = c(0, var_cutoff)
min_pos_vars = c(0)

tunes <- list()



# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
# tunes[[1]] <- tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
                                                                                       # protein = protein,
                                                                                       # alphas = alphas,
                                                                                       # minposvars = min_pos_vars)

# MEAN ESTIMATE CLOSE TO ZERO
# tunes[[2]] <- tune_alpha_mpv_mean_empir_cor (pc_fun = pc_fun_with_eval,
#                                              protein = protein,
#                                              alphas = alphas,
#                                              minposvars = min_pos_vars)


# EDGE SCORE
tune_alpha_mpv_edge_score_w4 <- tune_alpha_mpv_edge_score_factory(4)
tunes[[4]] <- edge_score_w4 <-
  tune_alpha_mpv_edge_score_w4(pc_fun = pc_fun,
                               alphas = alphas,
                               minposvars = min_pos_vars,
                               protein = protein)
#
# # print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
# print(edge_score_w4)
# plot.new()
# # # plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
# plot_partuning(edge_score_w4$all_values)
#
#
# LOWER LOCAL MIN OF NUMBER OF CONFLICT EDGES
tunes[[5]] <- tune_alpha_mpv_conflict_edges_loc_min <-
  tune_alpha_mpv_conflict_edges_local_min(pc_fun = pc_fun,
                               alphas = alphas,
                               minposvars = min_pos_vars,
                               protein = protein,
                               objective_fun_returns_indices = TRUE)

print(tune_alpha_mpv_conflict_edges_loc_min)
plot.new()
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
plot_partuning(tune_alpha_mpv_conflict_edges_loc_min$all_values)

#
# POINT BEFORE HIGHEST INCREASE IN CONFLICT EDGES WHILE still below t=6
tune_alpha_mpv_conflict_edges_gradient <- tune_alpha_mpv_max_gradient_of_conflict_edges_below_t_factory(6)

tunes[[6]] <- tune_alpha_mpv_conflict_edges_grad <-
  tune_alpha_mpv_conflict_edges_gradient(pc_fun = pc_fun,
                                         alphas = alphas,
                                         minposvars = min_pos_vars,
                                         protein = protein,
                                         objective_fun_returns_indices = TRUE,
                                         max_rows_in_plot = 3)

tuning <- tunes[[6]]
print(tuning)
plot.new()
plot_partuning(tuning$all_values)

pc_fun(alpha = tuning$best_alpha, min_pos_var = tuning$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all",
       effects_cluster_cut_height_h = 7)

pc_fun(alpha = tunes[[5]]$best_alpha, min_pos_var = tunes[[5]]$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all",
       effects_cluster_cut_height_h = 7)

pc_fun(alpha = tunes[[4]]$best_alpha, min_pos_var = tunes[[4]]$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all",
       effects_cluster_cut_height_h = 7)

#

# protein_causality(filename = file_name, alpha = tuning$best_alpha,
  # min_pos_var = tuning$best_minposvar, causal_analysis = TRUE, intervention_position = "all")

pc_fun(alpha = 0.01, min_pos_var = tunes[[4]]$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all",
       effects_cluster_cut_height_h = 7)

pc_fun(alpha = 0.05, min_pos_var = tunes[[4]]$best_minposvar,
       causal_analysis = TRUE, intervention_position = "all",
       effects_cluster_cut_height_h = 7)
