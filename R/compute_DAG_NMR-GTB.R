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
## set working directory from configuration.R
source("configuration_data.R")

protein_causality_NMRGTB <- function(
  # Data parameters
  # 
  # available data:
  # PDZ_DG
  # PDZ_DDG-all
  # PDZ_DDDG-5     # best results so far
  # PDZ_DDDG-10
  # PDZ_DDDG-all_372
  # PDZ_DDDG-all
  # PDZ_DDDG-all_SVD
  numerical = TRUE,
  protein = "GTB",
  # 
  # type_of_data = "DG",
  # type_of_data = "DDG",
  type_of_data = "NMR",
  # 
  # subtype_of_data = "all",
  subtype_of_data = "2",
  # subtype_of_data = "10",
  # 
  data_set = "don+acc",
  # data_set = "SVD",
  # data_set = "372",
  # 
  position_numbering = "crystal",
  # 
  # Analysis parameters
  # remove_positions_with_low_variance = TRUE,
  min_pos_var = 0,
  only_cols = NULL,
  only_cols_label = "",
  # 
  alpha = 0.01,
  ranked = FALSE,
  # 
  # pc_solve_conflicts = FALSE,
  # pc_u2pd = "retry",
  pc_solve_conflicts = TRUE,
  pc_u2pd = "relaxed",
  pc_conservative = FALSE,
  pc_maj_rule = FALSE,
  # 
  # weight_effects_on_by = "",
  # weight_effects_on_by = "var",
  # weight_effects_on_by = "mean",
  weight_effects_on_by = "median",
  # 
  # 
  # 
  # Graphical parameters
  graph_output_formats = "ps",
  graph_layout = "dot", # "dot", "circo", "fdp", "neato", "osage", "twopi"
  coloring = "auto", #"es",#"auto",  # "auto-all" or "all"
  colors = NULL,
  # 
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 oder NULL
  # TODO Marcel: dafür sorgen dass, wenn diese Option aktiv ist, kein graphics.off() o.ä. ausgeführt wird (und nur der graph geplottet wird)
  for_combined_plot = FALSE,
  mute_plot = FALSE,
  # 
  # description of other settings that should be appended to output-filename
  other = "", # cov", 
  # 
  # 
  # Technical parameters (print, plot, save, analysis)
  # steps = c("evaluation", "analysis"),
  graph_computation = TRUE,
  evaluation = FALSE,
  analysis = FALSE,#!pc_solve_conflicts,
  stages = c("orig"), #c("orig", "sub"), # "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graphs"),
  # 
  compute_pc_anew = FALSE,
  compute_localTests_anew = FALSE,
  # if (compute_everything_anew) {
  #   compute_pc_anew <- TRUE
  # }
  unabbrev_r_to_info = FALSE,
  print_r_to_console = TRUE,
  lines_in_abbr_of_r = 10,
  data_in_results = FALSE,
  output_parameters_in_results = FALSE,
  # 
  ida_percentile = "11", # top 11
  # ida_percentile = 0.75, # top 75%
  #
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
                           for_combined_plot = FALSE,
                           mute_plot = FALSE,
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

results_NMRGTB <- protein_causality_NMRGTB(pc_conservative = FALSE, pc_maj_rule = FALSE, pc_u2pd = "relaxed", pc_solve_confl = TRUE, analysis = FALSE, alpha = 0.01)
