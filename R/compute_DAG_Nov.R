source("~/.configuration_code.R")

source("functions_causal_effects.R")
source("functions_ci_tests.R")
source("functions_compute_DAG_categorical.R")
source("functions_compute_DAG_numerical.R")
source("functions_conversions.R")
source("functions_evaluate_DAG.R")
source("functions_general.R")
source("functions_i_o.R")
source("functions_linkcommunities.R")
source("functions_pymol.R")
source("functions_tools.R")

source("configuration_data.R")


protein_causality_NoV <- function(
                                  # data parameters
                                  # available data:
                                  # "NoV_NMR-Tit_BTS"
                                  # "NoV_NMR-Tit_Fuc"
                                  # type_of_data = "NMR-Tit"
                                  # type_of_data = c("NMR_Tit-Fuc", "NMR_Tit-BTS")
                                  # type_of_data = "Fuc-Tit-only-assigned"
                                  # subtype_of_data = "Fuc-old"
                                  # subtype_of_data = "Fuc"
                                  # TODO: wieder ermöglichen
                                  # subtype_of_data = c("Fuc", "BTS")
                                  
                                  numerical = TRUE,
                                  protein = "NoV",
                                  type_of_data = "NMR-Tit",
                                  subtype_of_data = c("Fuc", "BTS"),
                                  data_set = "",
                                  position_numbering = "",
                                  # analysis parameters
                                  min_pos_var = 0,
                                  only_cols = NULL,
                                  only_cols_label = "",
                                  alpha = 0.05,
                                  ranked = FALSE,
                                  pc_solve_conflicts = TRUE,
                                  pc_u2pd = "relaxed",
                                  pc_conservative = FALSE,
                                  pc_maj_rule = FALSE,
                                  weight_effects_on_by = "median", # "var", "mean", ""
                                  # graphical parameters
                                  graph_output_formats = "ps",
                                  graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
                                  coloring = "auto", # "auto", "auto-all", "all"
                                  colors = NULL,
                                  plot_as_subgraphs = FALSE,
                                  plot_only_subgraphs = NULL, # 1 is another option
                                  combined_plot = FALSE,
                                  other = "", # "cov"
                                  # technical parameters
                                  graph_computation = TRUE,
                                  evaluation = FALSE,
                                  analysis = FALSE, # !pc_solve_conflicts
                                  stages = c("orig"), # c("orig", "sub"), "sub"
                                  print_analysis = FALSE,
                                  plot_analysis = TRUE,
                                  plot_types = c("localTests", "graph"),
                                  compute_pc_anew = FALSE,
                                  compute_localTests_anew = FALSE,
                                  unnabbrev_r_to_info = FALSE,
                                  print_r_to_console = TRUE,
                                  lines_in_abbr_of_r = 10,
                                  data_in_results = FALSE,
                                  output_parameters_in_results = FALSE,
                                  ida_percentile = "11",
                                  file_separator = "/"
                                  ) {
  return(protein_causality(numerical = numerical,
                           protein = protein,
                           type_of_data = type_of_data,
                           subtype_of_data = subtype_of_data,
                           data_set = data_set,
                           position_numbering = position_numbering,
                           min_pos_var = min_pos_var,
                           only_cols = only_cols,
                           only_cols_label = only_cols_label,
                           alpha = alpha,
                           ranked = ranked,
                           pc_solve_conflicts = pc_solve_conflicts,
                           pc_u2pd = pc_u2pd,
                           pc_conservative = pc_conservative,
                           pc_maj_rule = pc_maj_rule,
                           weight_effects_on_by = weight_effects_on_by,
                           graph_output_formats = graph_output_formats,
                           graph_layout = graph_layout,
                           coloring = coloring,
                           colors = colors,
                           plot_as_subgraphs = plot_as_subgraphs,
                           plot_only_subgraphs = plot_only_subgraphs,
                           combined_plot = combined_plot,
                           other = other,
                           graph_computation = graph_computation,
                           evaluation = evaluation,
                           analysis = analysis,
                           stages = stages,
                           print_analysis = print_analysis,
                           plot_analysis = plot_analysis,
                           plot_types = plot_types,
                           compute_pc_anew = compute_pc_anew,
                           compute_localTests_anew = compute_localTests_anew,
                           unabbrev_r_to_info = unabbrev_r_to_info,
                           print_r_to_console = print_r_to_console,
                           lines_in_abbr_of_r = lines_in_abbr_of_r,
                           data_in_results = data_in_results,
                           output_parameters_in_results = output_parameters_in_results,
                           ida_percentile = ida_percentile,
                           file_separator = file_separator
                           )
         )
  
}

results_NoV <- protein_causality_NoV(pc_conservative = FALSE, pc_maj_rule = FALSE, pc_u2pd = "relaxed", pc_solve_confl = TRUE, analysis = FALSE, alpha = 0.01)
sink()
print(conflict_edges(results_NoV$pc@graph))


