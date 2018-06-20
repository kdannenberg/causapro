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

filename = "Nov_NMR-Mut_apo"
ranked = FALSE              # data SHOULD be ranked

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename, ranked = ranked,
                                         mute_all_plots = TRUE))

pc_fun_with_eval <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename, ranked = ranked,
                                         evaluation = TRUE,
                                         mute_all_plots = TRUE))


# debug(partuning_Nov_NMR_2018_06_01)

alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.09)
# alphas <- c(0.7)
min_pos_vars = c(0, 1) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059?? f. 0.73??

# debug(obj_function_est_per_edges_times_mean_estimate_times_edge_score)

# dev_of_edges_from_2 <- tune_alpha_mpv_dev_of_edges_from_2(pc_fun = pc_fun,
#                                                           alphas = alphas,
#                                                           minposvars = min_pos_vars)
# print(dev_of_edges_from_2)

# localTest estimate per (squared/...) edge(s)
# estimate_per_edges <- tune_alpha_mpv_estimate_per_edges(pc_fun = pc_fun_with_eval,
#                                                         alphas = alphas,
#                                                         minposvars = min_pos_vars)
# print(estimate_per_edges)

# mean_estimate <-  tune_alpha_mpv_mean_empir_cor(pc_fun = pc_fun_with_eval,
#                                                 alphas = alphas,
#                                                 minposvars = min_pos_vars)
# print(mean_estimate)

# CONFLICTS YIELD (MINIMIZED) SCORE INFINITY
# abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict <-
#   tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict(pc_fun = pc_fun_with_eval,
#                                                                           alphas = alphas,
#                                                                           minposvars = min_pos_vars)
# print(abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict)
# # TODO: make that possible:
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_no_conflict$all_values)


# CONFLICTS MULTIPLIED WITH (MINIMIZED) SCORE
# abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict <-
#   tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict_plus_1(pc_fun = pc_fun_with_eval,
#                                                                           alphas = alphas,
#                                                                           minposvars = min_pos_vars)
# print(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict)
# plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict$all_values)

# debug(obj_function_est_per_edges_times_mean_estimate_times_edge_score)

# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
tuning <- tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
                                                                                    alphas = alphas,
                                                                                    minposvars = min_pos_vars)
print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
plot.new()
plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)

# mean(abs_empir_cor_per_edges$all_values[3:30,], na.rm = TRUE)  # 164.2603
# mean(mean_empir_cor$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# mean(edge_score$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# mean(edge_score_w4$all_values[3:30,], na.rm = TRUE)  # 0.1867485

# EDGE SCORE
# abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score <-
# tune_alpha_mpv_edge_score_w4 <- tune_alpha_mpv_edge_score_factory(4)
# edge_score_w4 <-
#   tune_alpha_mpv_edge_score_w4(pc_fun = pc_fun,
#                                alphas = alphas,
#                                minposvars = min_pos_vars)
#
# # print(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score)
# print(edge_score_w4)
# plot.new()
# # plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score$all_values)
# plot_partuning(edge_score_w4$all_values)
#
# mean(abs_empir_cor_per_edges$all_values[3:30,], na.rm = TRUE)  # 164.2603
# mean(mean_empir_cor$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# mean(edge_score$all_values[3:30,], na.rm = TRUE)  # 0.1867485
# mean(edge_score_w4$all_values[3:30,], na.rm = TRUE)  # 0.1867485

#TODO: vielleciht stattdessen mit dem Score für die edge-types malnehmen???!


# müsste man die Konfliktkanten vllt exponentiell (2^) eingehen lassen, weil sie so böse sind?

# and den Gewichtungen von mean_est-, abs_est_pro_Kanten- und (Konflikt)-Kantenterm rumspielen


# best partuning so far:
# estimate_per_edges_mat * abs(mean_est_mat)^2 und alles mit konfliktkanten ausschließen (z.B. alpha = 0.7, bei minposvar = )

# best_alpha fixen, dann: (aber mute_all_plots_müsste wieder aus!) (und funktioniner für das clustering nicht!!)
pc_fun(alpha = tuning$best_alpha, min_pos_var = tuning$best_minposvar, causal_analysis = TRUE, intervention_position = "all")

# print(opt)

# TODO: Kann man das wieder fixen?
# for localTest estimate per (squared/...) edge(s)
# debug(plot_different_measures_of_all_results)
# measures <- c(function(v) {v[1]}, function(v) {v[1]/sqrt(v[2])}, function(v) {v[1]/v[2]}, function(v) {v[1]/v[2]^2}, function(v) {v[1]/sqrt(v[2])^v[2]})
# par(mfrow = c(1, length(measures)))
# plot_different_measures_of_all_results(all_results = estimate_per_edges$all_results, measures = measures, print = TRUE)

# TODO: Eine funktion schrieben, die solche plots macht



