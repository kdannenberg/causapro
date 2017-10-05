source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)

# NoV_NMR-Tit_B3S-with-unass
type_of_data = "NMR-tit"
subtype_of_data = "B3S-with-unass"
subtype_of_data = c("Fuc", "B4S")

# TODO: "with-unass" etc should rather ne data_set, shouldn't it?

results_NoV <- protein_causality_NoV(type_of_data = type_of_data,
                                     subtype_of_data = subtype_of_data, pc_conservative = FALSE, pc_maj_rule = TRUE, 
                                     pc_u2pd = "relaxed", pc_solve_confl = TRUE, analysis = FALSE, min_pos_var = 0.01, 
                                     alpha = 0.01)
sink()
print(conflict_edges(results_NoV$pc@graph))



 