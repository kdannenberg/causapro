# previously: analysis_find_good_graph.R
# previously: analysis_for_a_set_of_graphs.R
# next: partuning_alpha_by_mean_effects
source("~/.configuration_code.R")

source_all_function_scripts()

source("functions_analysis_for_a_set_of_graphs.R")
# source("functions_protein_causality.R")

# tests
source("tests_analysis_for_a_set_of_graphs.R")

# source("compute_DAG_G.R")
# source("compute_DAG_S.R")

source("configuration_data.R")

# debug(protein_causality)
# debug(estimate_DAG_from_numerical_data)

# debug(analyse_graphs_for_alphas_and_minposvars)
# debug(analyse_set_of_graphs)
# debug(determine_set_of_graphs)
# debug(graph_to_results)
# debug(compute_over_all_graphs)
# debug(causal_effects_ida)
# debug(score_for_effects)
# debug(mean_effects_min_max)

new_whole = TRUE
save_whole = TRUE
# file <- "RData/temp"
file <- "RData/all_effects_ida-reset_DDS_with_and_wo_cor_FUN.RData"
# OBS!!
cor_cov_FUN_default = ""
# data_set = ""  # "bin_approx"
# mal einbauen, dass man das ändern kann
# ida_function = "IDA-reset"

new_single_runs = TRUE
save_single_runs = TRUE


# weitermaachen: DDG-all alpha = 0.004, min_pos_var = 0.0001

# measures <- c("DDS", "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")
# measures <- c("DDDG-5", "DDDG-all")

# richtige Reihenfolge von -5, -10 und -all:
# measures <- c("DDG-5", "DDG-10", "DDG-all", "DDDG-5", "DDDG-10", "DDDG-all")
measures <- c("DDS", "DDS--none", "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")
# measures <- c("DDS--bin_approx")#, "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")
alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.0001)
min_pos_vars = c(0, 0.0001, 0.001, 0.01)
# min_pos_vars = c(0)



# measures <- c("DDG-5")
# alphas <- c(1e-20)
# min_pos_vars <- c(0.01)

if (new_whole || !file.exists(file)) {
  all_effects <- list()
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
    data_set = strsplit(measure_type_sub, "-")[[1]][3]
    if (is.na(subtype_of_data)) {
      subtype_of_data <- ""
    }
    if (is.na(data_set)) {
      data_set <- ""
    }

    if (grepl("none", data_set)) {
      cor_cov_FUN = "none"
      data_set = gsub("none", "", data_set)
    } else {
      cor_cov_FUN = cor_cov_FUN_default
    }

    if (grepl("bin_approx", data_set)) {
      perturbed_position = "98"
    } else {
      perturbed_position = "372"
    }

    protein_causality_FUN <- function_set_parameters(get(paste0("protein_causality_", measure)),
                                                     parameters = list(mute_all_plots = TRUE,
                                                                       cor_cov_FUN = cor_cov_FUN,
                                                                       type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
                                                                       subtype_of_data = subtype_of_data,
                                                                       data_set = data_set))

    ida_FUN <- function_set_parameters(causal_effects_ida,
                                       parameters = list(perturbed_position = perturbed_position, cov_FUN = cor_cov_FUN))

    # all_effects[[measure_type_sub]] <- tryCatch(analyse_graphs_for_alphas_and_minposvars(measure = measure,
    #                                                 # type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
    #                                                 # subtype_of_data = subtype_of_data,
    #                                                 # data_set = data_set,
    #                                                 alphas = alphas, min_pos_vars = min_pos_vars,
    #                                                 protein_causality_function = protein_causality_FUN,
    #                                                 ida_function = ida_FUN,
    #                                                 new = new_single_runs, save = save_single_runs), finally = plot_text(text = "Error."))
    all_effects[[measure_type_sub]] <- analyse_graphs_for_alphas_and_minposvars(measure = measure,
                                                    # type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
                                                    # subtype_of_data = subtype_of_data,
                                                    # data_set = data_set,
                                                    alphas = alphas, min_pos_vars = min_pos_vars,
                                                    protein_causality_function = protein_causality_FUN,
                                                    ida_function = ida_FUN,
                                                    new = new_single_runs, save = save_single_runs)
  }
  if (save_whole) {
    save(all_effects, file = file)
  }
} else {
  load(file)
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




