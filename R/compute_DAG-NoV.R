source("compute_DAG_numerical.R")
source("general_functions.R")
source("evaluate_DAG.R")
source("causal_effects.R")
source("linkcommunities.R")

# setwd("/Volumes/Causality/Viren/R/Code")
# setwd("/Volumes/Causality/Viren/R/")
source("configuration.R")


# Data parameters

# available data:
# "NoV_NMR-Tit_BTS"
# "NoV_NMR-Tit_Fuc"

numerical = TRUE
protein = "NoV"
type_of_data = "NMR-Tit"
# type_of_data = c("NMR_Tit-Fuc", "NMR_Tit-BTS")
# type_of_data = "Fuc-Tit-only-assigned"
# subtype_of_data = "Fuc-old"
# subtype_of_data = "Fuc"
# TODO: wieder ermöglichen
subtype_of_data = c("Fuc", "BTS")
data_set = ""

position_numbering = ""

# transpose = FALSE

# Analysis parameters
alpha = 0.05
ranked = FALSE

# TODO: add in other scripts
pc_solve_conflicts <- FALSE

# TODO:
# remove_cols <- c(210) 
only_cols = NULL
only_cols_label = ""
# only_cols = c(389,244,283,293,509,231,317,330)
# only_cols_label = "only_assigned"


# Graphical parameters
graph_output_formats = "ps"
graph_layout = "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
coloring = NULL
colors = NULL

plot_as_subgraphs = FALSE
plot_only_subgraphs = NULL
combined_plot = FALSE

# description of other settings that should be appended to output-filename
other = "" # cov" 

# Technical parameters (print, plot, save, analysis)
analysis = FALSE
print_analysis = TRUE
plot_analysis = FALSE
compute_pc_anew <- TRUE
compute_localTests_anew <- FALSE
unabbrev_r_to_info <- FALSE
print_r_to_console <- FALSE
lines_in_abbr_of_r <- 20
stages <- c("orig") # "sub"
plot_types <- c("localTests", "graphs")

file_separator = "/"

# INIT
sink.reset()
graphics.off()

data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set)

data <- read_data(data_description, only_cols = only_cols)

outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other,
                       alpha = alpha, only_cols_label, pc_solve_conflicts = pc_solve_conflicts, file_separator = file_separator)

directories <- strsplit(outpath, file_separator)
filename <- directories[[1]][length(directories[[1]])]
output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = "/", sep = "/")

# source_of_data = paste(protein, type_of_data, sep = "-")
# data <- read_data(files = source_of_data, transpose = transpose, only_cols = only_cols)
# 
# source_of_data <- paste(source_of_data, collapse = "+")
# type_of_data <- paste(type_of_data, collapse = "+")

# # }
# # colnames(data) <- sapply(strsplit(colnames(data), " "), function(x) x[1])
# 
# # if (!is.null(nuclei) && !(nuclei == "all") && !(nuclei == "")) {
# #   data <- data[grepl(nuclei, rownames(data)), ]
# #   if (dim(data)[1] == 0) {
# #     stop("No data for these nuclei!")
# #   }
# #   state <- paste(state, nuclei, sep = "-")
# #   source_of_data <- paste(source_of_data, nuclei, sep = "-")
# #   # source_of_data <- paste(protein, type_of_data, state, sep = "-")  # sollte das gleiche liefern
# # } else if (nuclei == "all") {
# #   # data[grepl("1H", rownames(data)), ] <- 0.5 * data[grepl("1H", rownames(data)), ]
# # }



# TODO include in read_data
if (ranked) {
  # if (rank) {
  # data <- cbind(apply(data, 2, rank))
  data <- t(apply(data, 1, rank))
  # }
  # if (state == "all") {
  #   stop("Data not available (ranked).")
  # } else {
  type_of_data <- paste(type_of_data, "ranked", sep = "-")
  # }
} 

# filename <- paste(only_cols_label, "_", source_of_data, "-alpha=", alpha, sep = "")
# output_dir <- paste("Outputs/", protein, "/", type_of_data, "/", filename, sep = "")
# # print(paste("Output will be written to ", getwd(), "/", substring(outpath, 0, nchar(outpath)), "...", sep = ""))
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
#   print("Directory created.")
# }
# outpath <- paste(output_dir, filename, sep = "/")

caption <- caption(protein = protein, data = paste(type_of_data, sep = ""), alpha = alpha, chars_per_line = 45)
parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, 
                                                     only_cols = only_cols, coloring = coloring, colors = colors, outpath = outpath) 

results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                                output_dir = output_dir, filename = filename, parameters_for_info_file = parameters_for_info_file,
                                alpha = alpha, pc_solve_conflicts = pc_solve_conflicts,
                                caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                                graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                                unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                                print_analysis = print_analysis, plot_analysis = plot_analysis, graph_output_formats = graph_output_formats)

