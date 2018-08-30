source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(plot_structure)
# debug(protein_graph_clustering)
# debug(partuning_over_alpha_and_minposvar)
# debug(cluster_pairwise_effects)
# debug(plot_different_measures_of_all_results)
debug(max_gradient_value_below_t)

filename = "pdi_MD_contacts_A"
ranked = FALSE              # data SHOULD be ranked

pc_fun <- function_set_parameters(protein_causality, parameters =
                                    list(filename = filename, ranked = ranked,
                                         mute_all_plots = TRUE))

pc_fun_with_eval <- function_set_parameters(protein_causality, parameters =
                                              list(filename = filename, ranked = ranked,
                                                   evaluation = TRUE,
                                                   mute_all_plots = TRUE))

# pc_fun(graph_computation = FALSE, show_variance_cutoff_plot = TRUE, min_pos_var = 0.005)
# stop()
# # angucken und setzen:
# var_cutoff = 0.005  # DDG
var_cutoff = 0.002 # DDDG

# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, 0.001, seq(0.01, 0.09, 0.02), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
alphas <- c(1e-100, 1e-20, 1e-10, 1e-5, 0.01)
# alphas <- c(0.7)
min_pos_vars = c(0)
# min_pos_vars = c(0, var_cutoff) # kleinste vorhandene Varianz Nov_NMR-Mut_apo: 0.059?? f. 0.73??

tuning <- list()

# MAX GRADIENT OF CONFLICT EDGES WHILE STILL LESS THAN 10
tune_fct <- tune_alpha_mpv_max_gradient_of_conflict_edges_below_t_factory(t_max_number_of_conflict_edges = 6)
debug(tune_fct)
tuning[["max_grad_conflict_less_than_10"]] <- tune_fct(pc_fun = pc_fun,
                                                     alphas = alphas,
                                                     minposvars = min_pos_vars,
                                                     objective_fun_returns_indices = TRUE,
                                                     plot_labels_as_rows_and_cols = FALSE)


# EDGE SCORE (KONFLICT AND DRIECTED) MULTIPLIED WITH (MINIMIZED) SCORE
# tuning["edge_score_times_score"] <- tune_alpha_mpv_abs_empir_cor_per_edges_times_mean_empir_cor_times_edge_score(pc_fun = pc_fun_with_eval,
#                                                                                        alphas = alphas,
#                                                                                        minposvars = min_pos_vars,
#                                                                                        plot_labels_as_rows_and_cols = FALSE)

print(tuning[[1]])
plot.new()
plot_partuning(tuning[[1]]$all_values)

pc_fun(alpha = tuning[[1]]$best_alpha, min_pos_var = tuning[[1]]$best_minposvar, causal_analysis = TRUE, intervention_position = "all")




