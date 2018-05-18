source("~/.configuration_code.R")

source_all_function_scripts()

# source("functions_causal_effect_without_IDA.R")
## sets working directory depending on local config file
source("configuration_data.R")

## Data parameters

## available data:
## TODO

# debug(protein_causality)
# debug(plot_graph)
# debug(analyse_set_of_graphs)
# debug(determine_set_of_graphs)
# debug(pymol_mean_effects)
# debug(compute_if_not_existent)

# debug("tune_alpha_bnlearn")

# debug(estimate_DAG_from_numerical_data)

# undebug(pseudo_ida_by_causalEffect)
# debug(causal_effects_ida)
# debug(node_colors_pie)
# debug(idaFast)
# debug(solve)
# debug(scale_effects)
# debug(set_effects_of_unconnected_positions_to_zero)
# debug(compare_effects_per_position)
# debug(int_pos_to_color_vector)
# debug(interesting_positions)

## graphics.off()
# plot.new()
# par(mfrow = c(1,4))
results_S <- protein_causality_S(#alpha = 0.25, min_pos_var = 0, #data_set = "bin_approx",
                                 #alpha = 0.08, min_pos_var = 0.01,
                                 alpha = 0.2, min_pos_var = 0.01, cor_cov_FUN = "", # cor_cov_FUN = "none",
                                 ranked = TRUE, plot_no_isolated_nodes = TRUE,
                                 show_variance_cutoff_plot = FALSE, max_conflict_edges = 15,
                                 causal_analysis = TRUE, pc_maj_rule = TRUE, pc_solve_confl = TRUE,
                                 for_combined_plot = TRUE, data_in_results = TRUE, plot_with_graphviz = TRUE,
                                 #plot_clusters = TRUE,
                                 linkcommunities = TRUE, evaluation = TRUE)
# results_S <- protein_causality_S(type_of_data = "DS", alpha = 0.01, min_pos_var = 0.01, analysis = TRUE, pc_maj_rule = TRUE, pc_solve_confl = TRUE,
#                                  for_combined_plot = TRUE, data_in_results = TRUE)

# Sieht super aus, aber min_pos_var auch unanständig hoch
# results_S <- protein_causality_S(pc_maj_rule = TRUE, pc_solve_confl = TRUE, analysis = FALSE, alpha = 0.1, min_pos_var = 0.06, data_in_results = TRUE)
# sink()
# print(conflict_edges(results_S$pc@graph))

source("~/.configuration_code.R")
