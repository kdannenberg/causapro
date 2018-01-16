# setwd("~/Documents/Uni/Viren/R/Code")
source("~/.configuration_code.R")

source("functions_compute_DAG_numerical.R")
source("functions_general.R")
source("evaluate_DAG.R")
source("optimize_parameters.R")

source("configuration_data.R")

# delta-G - no good results
# type_of_data = "DG"
# source_of_data = "delta-G-pdz-crystal"
# position_numbering = "crystal"

# Technical parameters (print, plot, save, analysis)
plot = TRUE
print = TRUE
analysis = FALSE
print_analysis = TRUE
plot_analysis = FALSE
compute_pc_anew <- FALSE
compute_localTests_anew <- FALSE
# if (compute_everything_anew) {
#   compute_pc_anew <- TRUE
# }
unabbrev_r_to_info <- FALSE
print_r_to_console <- FALSE
lines_in_abbr_of_r <- 20

stages <- c("orig") # c("orig", "sub", "anc")
plot_types <- c("graphs")

# DDDG
numerical = TRUE
type_of_data = "DDDG"
# source_of_data = "DG_pdz"
source_of_data = "DDDG-5_pdz"       # best results so far
# source_of_data = "DDDG-10_pdz"
# source_of_data = "DDDG-all-372"
# source_of_data = "DDDG-all"
# source_of_data = "DDDG-all_SVD"
data_descr <- source_of_data
position_numbering = "crystal"
protein = "PDZ"

graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
coloring = "auto-all"  # "auto-all" or "all"
colors <- NULL

plot_as_subgraphs = FALSE
plot_only_subgraphs = FALSE
combined_plot = FALSE

filename_data <- paste("../Data/", source_of_data, ".csv", sep = "")

data <- read_data(filename_data, transpose = FALSE)
# colnames(data) <- paste("X", colnames(data), sep = "")

alpha = 0.1

only_cols = NULL
only_cols_label = ""

# output_dir = paste("../Outputs/", type_of_data, sep = "")
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# }

filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", filename, sep = "")

# outpath = paste("../Outputs/", type_of_data, "/", only_cols_label, source_of_data, "-alpha=", alpha, sep="") 

caption <- caption(protein = protein, data = data_descr, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen

pc_fun <- function(alpha, outpath) {
  return(estimate_DAG_from_numerical_data(data = data, alpha = alpha, outpath = outpath))
}

analysis_after_pc_fun <- function(pc, outpath, caption) {
  return(analysis_after_pc(pc = pc, data, outpath = outpath, protein = protein, position_numbering = position_numbering, 
                           graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, 
                           plot_only_subgraphs = plot_only_subgraphs, coloring = coloring, colors = colors, 
                           stages = stages, plot_types = plot_types, unabbrev_r_to_info = unabbrev_r_to_info, 
                           print_r_to_console, lines_in_abbr_of_r, compute_localTests_anew = compute_localTests_anew, 
                           print = print, plot = plot, caption = caption))
} 


parameters_for_info_file_fun <- function(alpha, outpath) {
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, 
                                                       position_numbering = position_numbering, only_cols = only_cols, 
                                                       coloring = coloring, colors = colors, outpath = outpath)  
}

alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
# alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1) #seq(0.1, 0.2, 0.05))#, 0.99, 0.999, 0.9999, 0.99999)
# alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, seq(0.1, 0.9, 0.05), 0.99, 0.999, 0.9999, 0.99999)


filename_w_o_alpha <- paste(only_cols_label, source_of_data, "-alpha=", sep = "")
output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", sep = "")

# optimize_pc_results(pc_fun, analysis_after_pc_fun, output_dir, filename_w_o_alpha, alpha_start, compute_pc_anew = compute_pc_anew)
alpha_results_G <- linear_optimization_of_pc_results(pc_fun = pc_fun, analysis_after_pc_fun = analysis_after_pc_fun, 
                                                     parameters_for_info_file_fun = parameters_for_info_file_fun, 
                                                     output_dir = output_dir, filename_w_o_alpha = filename_w_o_alpha, 
                                                     alphas_vector = alphas_vector, compute_pc_anew = compute_pc_anew, 
                                                     plot_all_graphs = combined_plot, plot_rows = 4)

dev.copy(postscript, paste(output_dir, "/", filename_w_o_alpha, "several.ps", sep = ""), paper = "special", width = 18, height = 10)
dev.off
