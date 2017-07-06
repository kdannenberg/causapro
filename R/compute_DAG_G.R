source("configuration_code.R")

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

# source("compute_DAG_numerical.R")
# source("general_functions.R")
# source("evaluate_DAG.R")
# source("causal_effects.R")
# source("linkcommunities.R")
# source("ci-tests.R")

# setwd("/Volumes/Causality/Viren/R/")
source("configuration_data.R")

# Data parameters

# available data:
# PDZ_DG
# PDZ_DDG-all
# PDZ_DDDG-5     # best results so far
# PDZ_DDDG-10
# PDZ_DDDG-all_372
# PDZ_DDDG-all
# PDZ_DDDG-all_SVD
numerical = TRUE
protein = "PDZ"

# type_of_data = "DG"
# type_of_data = "DDG"
type_of_data = "DDDG"

# subtype_of_data = "all"
subtype_of_data = "5"
# subtype_of_data = "10"

data_set <- ""
# data_set <- "SVD"
# data_set <- "372"

position_numbering = "crystal"

# Analysis parameters
# remove_positions_with_low_variance = TRUE
min_pos_var = 0.001
only_cols = NULL
only_cols_label = ""

alpha = 0.01
ranked = FALSE

pc_solve_conflicts <- FALSE
pc_u2pd = "retry"

effects_on_weighting <- FALSE

stages <- c("orig") # "sub"
plot_types <- c("localTests", "graphs")


# Graphical parameters
graph_output_formats = "ps"
graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
coloring = "auto"  # "auto-all" or "all"
colors <- NULL

plot_as_subgraphs = FALSE
plot_only_subgraphs = NULL # 1 oder NULL
combined_plot = FALSE

# description of other settings that should be appended to output-filename
other = "" # cov" 


# Technical parameters (print, plot, save, analysis)
analysis = FALSE
print_analysis = FALSE
plot_analysis = FALSE
compute_pc_anew <- TRUE
compute_localTests_anew <- FALSE
# if (compute_everything_anew) {
#   compute_pc_anew <- TRUE
# }
unabbrev_r_to_info <- FALSE
print_r_to_console <- TRUE
lines_in_abbr_of_r <- 10

file_separator = "/"



# INIT
graphics.off()

data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set)

# filename_data <- paste("Data/", source_of_data, ".csv", sep = "")
data <- read_data(data_description, only_cols = only_cols)

data <- adjust_data(data = data, rank = ranked, min_var = min_pos_var)
type_of_data <- type_of_data_after_adjustment(type_of_data = type_of_data, rank = ranked, min_var = min_pos_var)


# colnames(data) <- paste("X", colnames(data), sep = "")

# if (other != "") {
#   data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other)
# }

outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other,
                       alpha = alpha, only_cols_label = only_cols_label, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                       file_separator = file_separator)

directories <- strsplit(outpath, file_separator)
filename <- directories[[1]][length(directories[[1]])]
output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = "/", sep = "/")

# output_dir = paste("../Outputs/", type_of_data, sep = "")
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# }

# TODO: data_description nutzen; nicht source_of_data

# filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
# # output_dir <- paste("~/Documents/Uni/Viren/R/Outputs/", protein, "/", type_of_data, "/", filename, sep = "")
# output_dir <- paste("Outputs", protein, type_of_data, filename, sep = "/")
# outpath <- paste(output_dir, filename, sep = "/")
# 
# filename <- paste(filename, other, sep = "-")

# outpath = paste("/Outputs/", type_of_data, "/", only_cols_label, source_of_data, "-alpha=", alpha, sep="") 

caption <- caption(protein = protein, data = data_description, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, 
                                                     only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  

results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                     output_dir = output_dir, filename = filename, outpath = outpath, parameters_for_info_file = parameters_for_info_file,
                     alpha = alpha, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
                     caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                     graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                     unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                     compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                     print_analysis = print_analysis, plot_analysis = plot_analysis, graph_output_formats = graph_output_formats)


if (!is.null(plot_only_subgraphs)) {
  # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
  node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
  subgraphs <- subgraphs_from_node_clusters(node_clustering, graph = results$orig$graph$NEL, protein = protein)
  graph <- subgraphs[[plot_only_subgraphs]]$graph
} else {
  graph <- results$orig$graph$NEL
}

plot_connected_components_in_pymol(protein = protein, position_numbering = position_numbering, graph = graph, 
                                   outpath = outpath, show_int_pos = TRUE, show_positions = FALSE, file_separator = file_separator)



# paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(314), to = c(383), all_paths = FALSE)
# plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE, 
#                     label = TRUE, show_positions = FALSE, file_separator = file_separator)

results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", relatve_effects_on_pos = TRUE,
               protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
               amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
               pymol_bg_color = "grey",
               barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)

# TODO: wenn nicht ranked die trotzdem alles positiv machen (offset?) bevor die farben vergeben werden

# effects_of_76 <- idaFast(which(colnames(data) == "372"), 1:92, cov(data), results$pc@graph)
# rownames(effects_of_76) <- colnames(data)
# int_pos <- interesting_positions("PDZ", "crystal", coloring = "all")
# colors_by_effect <- hue_by_effect(effects_of_76, int_pos)
# plot_total_effects_in_pymol(colors_by_effect, perturbated_position = "76", protein = protein, outpath = outpath, amplification_exponent = 10)
# td_analysis(graph = results$orig$graph$NEL, outpath = outpath)
