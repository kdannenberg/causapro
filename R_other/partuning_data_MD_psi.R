source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(plot_clusters_in_pymol)
# debug(plot_graph)
# debug(protein_graph_clustering)

# debug(partuning_over_alpha_and_minposvar)
# debug(obj_function_est_per_edges_times_mean_estimate_no_conflict)


filename = "pdi_MD_psi_AB"
protein = "pdi"

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename,
                                         every_n_th_row = 100,
                                         # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                         # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                         # causal_analysis = FALSE, #min_pos_var = 0.01,
                                         # evaluation = TRUE,
                                         # alpha = 0.1, min_pos_var = 0,
                                         # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                         # ranked = TRUE, pc_indepTest = "smc-jt",
                                         mute_all_plots = TRUE))

# pc_fun(show_variance_cutoff_plot = TRUE)

pc_fun_with_eval <- function_set_parameters(protein_causality_NoV, parameters =
                                              list(filename = filename,
                                                   every_n_th_row = 100,
                                                   max_edges_for_evaluation = 60,
                                                   # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                                   # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                                   # causal_analysis = FALSE, #min_pos_var = 0.01,
                                                   evaluation = TRUE,
                                                   # alpha = 0.1, min_pos_var = 0,
                                                   # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                                   # ranked = TRUE, pc_indepTest = "smc-jt",
                                                   mute_all_plots = TRUE))


alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.3)
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.05, 0.01))
# alphas = c(0.001, 0.00175, 0.03, 0.05)
# alphas = c(0.001, 0.8)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.8)
min_pos_vars = c(0, 0.5) # kleinste vorhandene Varianz: 0.059

tunes <- list()

# ## CONFLICTS YIELD (MINIMIZED) SCORE INFINITY
# tunes[[1]] <- abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict <-
#   tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict(pc_fun = pc_fun_with_eval,
#                                                                           alphas = alphas,
#                                                                           minposvars = min_pos_vars,
#                                                                           protein = protein)
# print(abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict)
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict$all_values)

# ## CONFLICTS MULTIPLIED WITH (MINIMIZED) SCORE
# tunes[[2]] <- abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict <-
#   tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict_plus_1(pc_fun = pc_fun_with_eval,
#                                                                           alphas = alphas,
#                                                                           minposvars = min_pos_vars,
#                                                                           protein = protein)
# print(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict)
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict$all_values)
#
# # debug(obj_function_est_per_edges_times_mean_estimate_times_edge_score)
#
# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
# tunes[[3]] <- abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score <-
#   tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
#                                                                                alphas = alphas,
#                                                                                minposvars = min_pos_vars,
#                                                                                protein = protein)
# print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
# plot.new()
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
#
# # mean(abs_empir_cor_per_edges$all_values[3:30,], na.rm = TRUE)  # 164.2603
# # mean(mean_empir_cor$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# # mean(edge_score$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# # mean(edge_score_w4$all_values[3:30,], na.rm = TRUE)  # 0.1867485
#
#
# # EDGE SCORE
# # abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score <-
# tune_alpha_mpv_edge_score_w4 <- tune_alpha_mpv_edge_score_factory(4)
# tunes[[4]] <- edge_score_w4 <-
#   tune_alpha_mpv_edge_score_w4(pc_fun = pc_fun,
#                                alphas = alphas,
#                                minposvars = min_pos_vars,
#                                protein = protein)
#
# # print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
# print(edge_score_w4)
# plot.new()
# # plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
# plot_partuning(edge_score_w4$all_values)
#
#
# # LOWER LOCAL MIN OF NUMBER OF CONFLICT EDGES
# tunes[[5]] <- tune_alpha_mpv_conflict_edges_loc_min <-
#   tune_alpha_mpv_conflict_edges_local_min(pc_fun = pc_fun,
#                                alphas = alphas,
#                                minposvars = min_pos_vars,
#                                protein = protein,
#                                objective_fun_returns_indices = TRUE)
#
# print(tune_alpha_mpv_conflict_edges_loc_min)
# plot.new()
# # plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
# plot_partuning(tune_alpha_mpv_conflict_edges_loc_min$all_values)
#
#
# POINT BEFORE HIGHEST INCREASE IN CONFLICT EDGES WHILE still below t=6
tune_alpha_mpv_conflict_edges_gradient <- tune_alpha_mpv_max_gradient_of_conflict_edges_below_t_factory(6)

tunes[[6]] <- tune_alpha_mpv_conflict_edges_grad <-
  tune_alpha_mpv_conflict_edges_gradient(pc_fun = pc_fun,
                                          alphas = alphas,
                                          minposvars = min_pos_vars,
                                          protein = protein,
                                          objective_fun_returns_indices = TRUE)
#
# # mean(abs_empir_cor_per_edges$all_values[3:30,], na.rm = TRUE)  # 164.2603
# # mean(mean_empir_cor$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# # mean(edge_score$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# # mean(edge_score_w4$all_values[3:30,], na.rm = TRUE)  # 0.1867485
#
tuned_pars <- tunes[[6]]
#
# print(tuned_pars)
# plot.new()
# plot_partuning(tuned_pars$all_values)


pc_fun(alpha = tuned_pars$best_alpha, min_pos_var = tuned_pars$best_minposvar, causal_analysis = TRUE, intervention_position = "all")


