source("~/.configuration_code.R")

source_all_function_scripts()

## sets working directory depending on local config file
source("configuration_data.R")

## Data parameters

## available data:
## TODO

# debug(compute_if_not_existent)
# results_S <- protein_causality_S(alpha = 0.08, min_pos_var = 1e-4, analysis = FALSE, pc_maj_rule = TRUE, pc_solve_confl = TRUE, data_in_results = TRUE)

# Sieht super aus, aber min_pos_var auch unanständig hoch
# results_S <- protein_causality_S(pc_maj_rule = TRUE, pc_solve_confl = TRUE, analysis = FALSE, alpha = 0.1, min_pos_var = 0.06, data_in_results = TRUE)
# sink()
# print(conflict_edges(results_S$pc@graph))

source("~/.configuration_code.R")