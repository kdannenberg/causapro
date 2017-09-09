# previously: analysis_find_good_graph.R
# previously: analysis_for_a_set_of_graphs.R
# next: partuning_alpha_by_mean_effects
source("~/.configuration_code.R")

source("functions_causal_effects.R")
source("functions_general.R")
source("functions_conversions.R")
source("functions_tools.R")

source("functions_analysis_for_a_set_of_graphs.R")
# source("functions_protein_causality.R")

# tests
source("tests_analysis_for_a_set_of_graphs.R")

source("compute_DAG_G.R")
source("~/.configuration_code.R")
source("compute_DAG_S.R")
# source("configuration_data.R")


measures <- c("DDS", "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")
# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.005, 0.01, 0.05, 0.1)
alphas <- c(1e-20, 1e-10)
min_pos_vars = c(0.0001, 0.001, 0.01)

# measures <- c("DDS")
# alphas <- c(0.01)
# min_pos_vars <- c(0.0001)

all_effects <- list()
for (measure_type_sub in measures) {
  measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
  subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
  if (is.na(subtype_of_data)) {
    subtype_of_data <- ""
  }
all_effects[[measure_type_sub]] <- analyse_graphs_for_alphas_and_minposvars(measure = measure, type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
                                                  subtype_of_data = subtype_of_data,
                                                  alphas = alphas, min_pos_vars = min_pos_vars)
}



# do_for_measures <- function(FUN, measures) {
#   res <- list()
#   for (measure_type_sub in measures) {
#     measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
#     subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
#     if (is.na(subtype_of_data)) {
#       subtype_of_data <- ""
#     }
#     res[[measure_type_sub]] <- FUN(measure = measure, type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
#                                    subtype_of_data = subtype_of_data)
#     
#   }
#   return(res)
# }
# 
# 
# analyse_graphs_fct <- set_parameters(analyse_graphs_for_alphas_and_minposvars, parameters = list(alphas = alphas, min_pos_vars = min_pos_vars))
# all_effects_do <- do_for_measures(analyse_graphs_fct, measures)




