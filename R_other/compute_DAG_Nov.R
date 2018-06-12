source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(plot_structure_evaluation)
# debug(adjust_data)
# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)

# NoV_NMR-Tit_B3S-with-unass
type_of_data = "NMR-Mut"
## subtype_of_data = "B3S-with-unass"
# subtype_of_data = "B3S_with-unass_1H" ## options are ..._Euclidean, ..._1H and ..._15N
##subtype_of_data = "BTS"
# subtype_of_data = c("Fuc", "B4S")

subtype_of_data = "apo" # Daten vom 1.6.18

# TODO: "with-unass" etc should rather ne data_set, shouldn't it?
plot.new()
results_NoV <- protein_causality_NoV(type_of_data = type_of_data,
                                     subtype_of_data = subtype_of_data,
                                     # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                     # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                     evaluation = TRUE,
                                     causal_analysis = FALSE, #min_pos_var = 0.01,
                                     alpha = 0.2, min_pos_var = 0,
                                     # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                     # ranked = TRUE, pc_indepTest = "copula",
                                     plot_no_isolated_nodes = TRUE, plot_with_graphviz = TRUE,
                                     print_connected_components = TRUE, plot_clusters = FALSE)
# sink()

# print(conflict_edges(results_NoV$pc@graph))
# print(sum((results_NoV$orig$localTests$r$estimate)^2))


