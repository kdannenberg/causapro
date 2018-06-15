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
# debug(obj_function_est_per_edges_times_mean_estimate_no_conflict)


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

pc_fun_with_eval <- function_set_parameters(protein_causality_NoV, parameters =
                                    list(type_of_data = type_of_data,
                                         subtype_of_data = subtype_of_data,
                                         # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                         # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                         # causal_analysis = FALSE, #min_pos_var = 0.01,
                                         evaluation = TRUE,
                                         # alpha = 0.1, min_pos_var = 0,
                                         # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                         # ranked = TRUE, pc_indepTest = "smc-jt",
                                         mute_all_plots = TRUE))


alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas = c(0.001, 0.8)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.8)
min_pos_vars = c(0, 0.1) # kleinste vorhandene Varianz: 0.059

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


# # CONFLICTS MULTIPLIED WITH (MINIMIZED) SCORE
abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict <-
  tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict_plus_1(pc_fun = pc_fun_with_eval,
                                                                          alphas = alphas,
                                                                          minposvars = min_pos_vars)
print(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict)
plot_partuning(abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict$all_values)

#TODO: vielleciht stattdessen mit dem Score für die edge-types malnehmen???!

# TODO: ist alpha = 0.4 und min_pos_var= 0.1 wirklich so gut? es gibt 4 Konfliktkanten!

# müsste man die Konfliktkanten vllt exponentiell (2^) eingehen lassen, weil sie so böse sind?

# and den Gewichtungen von mean_est-, abs_est_pro_Kanten- und (Konflikt)-Kantenterm rumspielen


# best partuning so far:
# estimate_per_edges_mat * abs(mean_est_mat)^2 und alles mit konfliktkanten ausschließen (z.B. alpha = 0.7, bei minposvar = )

# best_alpha fixen, dann: (aber mute_all_plots_müsste wieder aus!) (und funktioniner für das clustering nicht!!)
tuned_pars <- abs_empir_cor_per_edges_times_mean_empir_cor_times_conflict$
pc_fun(alpha = tuned_pars$best_alpha, min_pos_var = tuned_pars$best_minposvar, causal_analysis = TRUE, perturbed_positions = "all")

# print(opt)

# TODO: Kann man das wieder fixen?
# for localTest estimate per (squared/...) edge(s)
# debug(plot_different_measures_of_all_results)
# measures <- c(function(v) {v[1]}, function(v) {v[1]/sqrt(v[2])}, function(v) {v[1]/v[2]}, function(v) {v[1]/v[2]^2}, function(v) {v[1]/sqrt(v[2])^v[2]})
# par(mfrow = c(1, length(measures)))
# plot_different_measures_of_all_results(all_results = estimate_per_edges$all_results, measures = measures, print = TRUE)

# TODO: Eine funktion schrieben, die solche plots macht
