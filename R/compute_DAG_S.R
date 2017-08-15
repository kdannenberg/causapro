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
## sets working directory depending on local config file
source("configuration_data.R")

## Data parameters

## available data:
## TODO

protein_causality_S <- function(
                                # data parameters
                                numerical = TRUE,
                                protein = "PDZ",
                                type_of_data = "DDS",
                                subtype_of_data = "",
                                data_set = "",
                                position_numbering = "crystal",
                                # analysis parameters
                                min_pos_var = 0,
                                only_cols = NULL,
                                only_cols_label = "",
                                alpha = 0.01,
                                ranked = FALSE,
                                pc_solve_conflicts = TRUE,
                                pc_u2pd = "relaxed",
                                pc_conservative = FALSE,
                                pc_maj_rule = FALSE,
                                weight_effects_on_by = "median", # "var", "mean", ""
                                # graphical parameters
                                graph_output_formats = "ps",
                                graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
                                coloring = "es", # "auto", "auto-all", "all"
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
  graphics.off()
  data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data <- subtype_of_data, data_set = data_set)
  data <- read_data(data_description, only_cols = only_cols)
  data <- adjust_data(data = data, rank = ranked, min_var = min_pos_var)
  type_of_data <- type_of_data_after_adjustment(type_of_data = type_of_data, rank = ranked, min_var = min_pos_var)
  if(!is.numeric(ida_percentile)) {
    ida_percentile <- 1 - (as.numeric(ida_percentile) / dim(data)[2])
  }
  outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other,
                         alpha = alpha, min_pos_var = min_pos_var, only_cols_label = only_cols_label, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                         pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule, file_separator = file_separator)
  directories <- strsplit(outpath, file_separator)
  filename <- directories[[1]][length(directories[[1]])]
  output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
  caption <- get_caption(protein = protein, data = data_description, alpha = alpha, min_pos_var = min_pos_var, chars_per_line = 45)
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, 
                                                       only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = file_separator))  
 
}


parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  

results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                                output_dir = output_dir, filename = filename, parameters_for_info_file = parameters_for_info_file,
                                alpha = alpha, caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                                graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                                unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                                print_analysis = print_analysis, plot_analysis = plot_analysis)
graph <- results$orig$graph$NEL
cat(paste("Number of edges", number_of_edges(graph), "\n"))
dg = dagitty(conv_to_r(graph))
cat(dconnected(dg, "372", "330"))

plot_connected_components_in_pymol(protein = protein, position_numbering = position_numbering, graph = results$orig$graph$NEL, 
                                   outpath = paste(output_dir, filename, sep = "/"), only_int_pos = TRUE, color_int_pos = FALSE, 
                                   coloring_for_int_pos = "auto")

