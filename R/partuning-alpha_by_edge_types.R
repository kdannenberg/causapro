source("~/.configuration_code.R")

# # source("functions_causal_effects.R")
# # source("functions_general.R")
# # source("functions_conversions.R")
# # source("functions_tools.R")
# 
# # source("analysis_for_a_set_of_graphs.R")
# # source("~/.configuration_code.R")
# source("functions_analysis_for_a_set_of_graphs.R")
# source("functions_partuning.R")

source_all_function_scripts()

# source("compute_DAG_G.R")
# source("compute_DAG_S.R")

source("configuration_data.R")

plot_labels_as_rows_and_cols = TRUE
plot_logscale_alpha = FALSE
# if scalar: height of the y-axis = max. possible number of edges (binom(n_nodes, 2)) / ylim_for_plots
# if vector: = ylim_for_plots
# if NULL: automatically set in each plot
ylim_for_plots = c(0,200) #5 # c(0, 200) # or NULL
opt_alpha_double_weight_conflict = FALSE
compute_anew = FALSE
conflict_edge_weight = 1

measures <- c("DDS",
  "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")

min_pos_vars <- c(0.0001, 0.001, 0.01)

# alphas <- c(1e-10, 1e-5, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05, 0.1)
alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas <- c(0.01)

load(file = "edge_types.RData")

if (compute_anew || !exists("edge_types")) {
  edge_types <- list()
  # TODO: DDS hinzufügen
  # for (type_of_data in c("DDG", "DDDG")) {
  #   edge_types[[type_of_data]] <- list()
  #   for (subtype_of_data in c("5", "10", "all")) {
  #     edge_types[[type_of_data]][[subtype_of_data]] <- list()
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    edge_types[[type_of_data]][[subtype_of_data_list]] <- list()
    
    protein_causality_function = get(paste0("protein_causality_", measure))
    
    for (min_pos_var in min_pos_vars) {
      edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]] <- list()
      data <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                 min_pos_var = min_pos_var, graph_computation = FALSE, analysis = FALSE, 
                                                 evaluation = FALSE, data_in_results = TRUE)$data
      edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]]$max_nodes <- dim(data)[2]
      for (alpha in alphas) {
        protein_causality_function <- get(paste0("protein_causality_", measure))
        results <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                            min_pos_var = min_pos_var, alpha = alpha,
                            pc_solve_conflicts = TRUE, pc_u2pd = "relaxed",
                            pc_maj_rule = TRUE, pc_conservative = FALSE,
                            evaluation = FALSE, analysis = FALSE)
        edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]][[as.character(alpha)]] <- conflict_edges(results$pc@graph)
      }
    }
    }
  #   }
  # }
}


if (compute_anew || !exists("numbers_of_nodes")) {
  numbers_of_nodes <- list()
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    numbers_of_nodes[[type_of_data]][[subtype_of_data_list]] <- list()
    
    protein_causality_function = get(paste0("protein_causality_", measure))
    
    for (min_pos_var in min_pos_vars) {
      data <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                         min_pos_var = min_pos_var, graph_computation = FALSE, analysis = FALSE, 
                                         evaluation = FALSE, data_in_results = TRUE)$data
      numbers_of_nodes[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]] <- dim(data)[2]
    }
  }
}

# best_alphas_1 <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
#                                              plot_logscale_alpha = FALSE, ylim_for_plots = ylim_for_plots, # or NULL
#                                              opt_alpha_double_weight_conflict = FALSE,
#                                              print = TRUE, plot = TRUE)
best_alphas <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
                                             plot_logscale_alpha = FALSE, ylim_for_plots = ylim_for_plots, # or NULL
                                             conflict_edge_weight = conflict_edge_weight,
                                             print = TRUE, plot = TRUE)

# debug(determine_set_of_graphs)
# debug(analyse_set_of_graphs)

effects_best_alphas <- effects_for_distinct_alphas(best_alphas)