source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(plot_clusters_in_pymol)


results_NoV <- protein_causality_NoV(subtype_of_data = "", pc_conservative = FALSE, pc_maj_rule = TRUE, 
                                     pc_u2pd = "relaxed", pc_solve_confl = TRUE, analysis = FALSE,
                                     min_pos_var = 0.01, alpha = 0.06)
sink()
print(conflict_edges(results_NoV$pc@graph))



