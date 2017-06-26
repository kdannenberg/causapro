setwd("~/Viren/R/Code")
source("compute_DAG_numerical.R")
source("general_functions.R")
source("evaluate_DAG.R")

## Technical parameters (print, plot, save, analysis)
## those options are either set to FALSE or to TRUE (unused/used option)
## analyse DAG using dagitty
analysis = FALSE
## do not print the dagitty analysis, but save it somewhere?
print_analysis = TRUE
## plot the dagitty analysis
plot_analysis = FALSE
## compute new dag/analysis (TRUE) or use precomputed one (FALSE)
compute_pc_anew <- FALSE
compute_localTests_anew <- FALSE
## if (compute_everything_anew) {
##   compute_pc_anew <- TRUE
## }
## what is this doing, more information to info file?
unabbrev_r_to_info <- FALSE
## and this?
print_r_to_console <- FALSE
lines_in_abbr_of_r <- 20

## what are stages again?
stages <- c("orig", "anc") # "sub"
plot_types <- c("localTests", "graphs")
## choose one of the layouts offered by Rgraphviz
graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
## auto-all and auto color more nodes - what exactly do they do?
coloring = "auto" # "auto", "auto-all" or "all"
colors <- NULL

## plot with clustering if set TRUE
plot_as_subgraphs = FALSE
## only plots cluster - shouldn't this be stored in a different file?
plot_only_subgraphs = 1 #1 NULL

## DDS
## choose data type - DDS is measure proposed by Halabi et al. (Cell, 2009)
numerical = TRUE
protein = "PDZ"
type_of_data = "DDS"
source_of_data = "DDS_pdz"
# source_of_data = "DDS_pdz_SVD"
position_numbering = "crystal"

# graphics.off()
# par(mfrow = c(3,2))

filename <- paste("../Data/", source_of_data, ".csv", sep = "")

data <- read_data(filename, transpose = FALSE)
# colnames(data) <- paste("X", colnames(data), sep = "")

## choose level of significance
alpha = 0.125

## what is this for?
only_cols = NULL
only_cols_label = ""

# outpath = paste("../Outputs/", type_of_data, "/", only_cols_label, source_of_data, "-alpha=", alpha, sep = "") 

filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
output_dir <- paste("~/Viren/R/Outputs/", protein, "/", type_of_data, "/", filename, sep = "")

# outpath = paste("../Outputs/", type_of_data, "/", only_cols_label, source_of_data, "-alpha=", alpha, sep="") 
data_descr <- source_of_data

caption <- caption(protein = protein, data = data_descr, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
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

results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", relatve_effects_on_pos = TRUE,
                  protein = protein, results = results, coloring = "all", no_colors = FALSE,
                  outpath = outpath, amplification_exponent = 0.5, amplification_factor = TRUE, rank_effects = FALSE,
                  effect_to_color_mode = "#FFFFFF", pymol_bg_color = "grey", barplot = TRUE, caption = caption,
                  show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE)
