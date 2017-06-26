setwd("~/Viren/R/Code")

source("compute_DAG_numerical.R")
source("general_functions.R")
source("evaluate_DAG.R")

# delta-G - no good results
# type_of_data = "DG"
# source_of_data = "delta-G-pdz-crystal"
# position_numbering = "crystal"

## Technical parameters (print, plot, save, analysis)
## those options are either set to FALSE or to TRUE (unused/used option)
## analyse DAG using dagitty
analysis = TRUE
## do not print the dagitty analysis, but save it somewhere?
print_analysis = FALSE
## plot the dagitty analysis
plot_analysis = FALSE
## compute new dag/analysis (TRUE) or use precomputed one (FALSE)
compute_pc_anew <- TRUE
compute_localTests_anew <- FALSE
## if (compute_everything_anew) {
##   compute_pc_anew <- TRUE
## }
## what is this doing, more information to info file?
unabbrev_r_to_info <- FALSE
## and this?
print_r_to_console <- TRUE
lines_in_abbr_of_r <- 10

## what are stages again?
stages <- c("orig") # "sub" "orig"
plot_types <- c("localTests", "graphs")

## DDDG
## choose data type - DDDG is measure proposed by Lockless and Ranganathan (Science, 1999)
numerical = TRUE
type_of_data = "DDDG"
# source_of_data = "DG_pdz"
source_of_data = "DDDG-5_pdz"       # best results so far
# source_of_data = "DDDG-10_pdz"
# source_of_data = "DDDG-all-372"
# source_of_data = "DDDG-all"
# source_of_data = "DDDG-all_SVD"
data_descr <- source_of_data
## numbering differs in paper and alignment - other options?
position_numbering = "crystal"
protein = "PDZ"
## choose one of the layouts offered by Rgraphviz
graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
## auto-all and auto color more nodes - what exactly do they do?
coloring = "auto"  # "auto", "auto-all" or "all"
## other options?
colors <- NULL

## plot with clustering if set TRUE
plot_as_subgraphs = FALSE
## only plots cluster - shouldn't this be stored in a different file?
plot_only_subgraphs = NULL # 1 oder NULL
## what is this option doing?
combined_plot = TRUE

graphics.off()
## get Data
filename_data <- paste("../Data/", source_of_data, ".csv", sep = "")


data <- read_data(filename_data, transpose = FALSE)
# colnames(data) <- paste("X", colnames(data), sep = "")
## choose significance level
alpha = 0.05

## what are those options for?
only_cols = NULL
only_cols_label = ""

# output_dir = paste("../Outputs/", type_of_data, sep = "")
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# }

filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
## maybe one should pass/set the data directory and everything is put relative to this dir?
## or any other way to get rid of hard coded paths
output_dir <- paste("~/Viren/R/Outputs/", protein, "/", type_of_data, "/", filename, sep = "")
outpath <- paste(output_dir, filename, sep = "/")

# outpath = paste("../Outputs/", type_of_data, "/", only_cols_label, source_of_data, "-alpha=", alpha, sep="") 

caption <- caption(protein = protein, data = data_descr, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  

results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                     output_dir = output_dir, filename = filename, parameters_for_info_file = parameters_for_info_file,
                     alpha = alpha, caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                     graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                     unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                     compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                     print_analysis = print_analysis, plot_analysis = plot_analysis)


if (!is.null(plot_only_subgraphs)) {
  # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
  node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
  subgraphs <- subgraphs_from_node_clusters(node_clustering, graph = results$orig$graph$NEL, protein = protein)
  graph <- subgraphs[[plot_only_subgraphs]]$graph
} else {
  graph <- results$orig$graph$NEL
}

cat(paste("Number of edges", number_of_edges(graph), "\n"))

plot_connected_components_in_pymol(protein = protein, position_numbering = position_numbering, graph = graph, 
                                   outpath = paste(output_dir, filename, sep = "/"), color_int_pos = FALSE, show_positions = FALSE)



paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(314), to = c(383), all_paths = FALSE)
plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE, label = TRUE, show_positions = FALSE)

# td_analysis(graph = results$orig$graph$NEL, outpath = outpath)

results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", relatve_effects_on_pos = TRUE, 
               protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath, 
               amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF", 
               pymol_bg_color = "grey", 
               barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)
