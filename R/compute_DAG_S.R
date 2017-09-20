source("~/.configuration_code.R")

source_all_function_scripts()

## sets working directory depending on local config file
source("configuration_data.R")

## Data parameters

## available data:
## TODO

# debug(protein_causality)
# debug(analyse_set_of_graphs)
# debug(determine_set_of_graphs)
# debug(pymol_mean_effects)
# debug(compute_if_not_existent)

graphics.off()
par(mfrow = c(1,4))
results_S <- protein_causality_S(alpha = 0.06, min_pos_var = 0.01, analysis = TRUE, pc_maj_rule = TRUE, pc_solve_confl = TRUE,
                                 for_combined_plot = TRUE, data_in_results = TRUE)
# results_S <- protein_causality_S(type_of_data = "DS", alpha = 0.01, min_pos_var = 0.01, analysis = TRUE, pc_maj_rule = TRUE, pc_solve_confl = TRUE, 
#                                  for_combined_plot = TRUE, data_in_results = TRUE)

# Sieht super aus, aber min_pos_var auch unanständig hoch
# results_S <- protein_causality_S(pc_maj_rule = TRUE, pc_solve_confl = TRUE, analysis = FALSE, alpha = 0.1, min_pos_var = 0.06, data_in_results = TRUE)
# sink()
# print(conflict_edges(results_S$pc@graph))

source("~/.configuration_code.R")