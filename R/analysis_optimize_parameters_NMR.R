setwd("~/Documents/Uni/Viren/R/Code")

source("compute_DAG_numerical.R")
source("general_functions.R")
source("evaluate_DAG.R")
source("optimize_parameters.R")

# General parameters
protein = "p38g"
type_of_data = "NMR"

# state = "allstates"
state = "wild_inact_apo"
# state = "wild_ATP"
# state = "wild_act"
# state = "wild_BIRB796"
# previously given as:
# source_of_data = "p38g-NMR-wild_inact_apo" # Table S1
# source_of_data = "p38g-NMR-wild_ATP" # Table S2
# source_of_data = "p38g-NMR-wild_act" # Table S3
# source_of_data = "p38g-NMR-wild_BIRB796" # Table S4
position_numbering = ""

if ((protein == "p38g") && (type_of_data == "NMR")) {
  transpose=TRUE
} else {
  transpose=FALSE
}

# tuning parameters
alpha_start = 0.15
alpha <- alpha_start


# Output parameters (print, plot, save, analysis)
plot = TRUE
print = TRUE
compute_pc_anew <- FALSE
compute_localTests_anew <- FALSE
unabbrev_r_to_info <- FALSE
print_r_to_console <- FALSE
lines_in_abbr_of_r <- 20
stages <- c("orig") # "sub"
plot_types <- c("graphs")
graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
coloring = "FS3-pie"
colors = "auto"
plot_as_subgraphs = TRUE
plot_only_subgraphs = NULL

# coloring = "FS3-simple"
only_cols = NULL
only_cols_label = ""

combined_plot = FALSE

# Code
if (state == "allstates") {
  data_list <- list()
  for (state in c("wild_inact_apo", "wild_ATP", "wild_act", "wild_BIRB796")) {
    source_of_data = paste(protein, type_of_data, state, sep = "-")
    filename <- paste("../Data/", source_of_data, ".csv", sep = "")
    var <- read_data(filename, transpose = transpose)
    rownames(var) <- paste(rownames(var), state, sep = "-")
    data_list[state][[1]] <- var
  }
  data <- do.call(rbind, data_list)
} else {
  source_of_data = paste(protein, type_of_data, state, sep = "-")
  filename <- paste("../Data/", source_of_data, ".csv", sep = "")
  data <- read_data(filename, transpose = transpose)
}


pc_fun <- function(alpha, outpath) {
  return(estimate_DAG_from_numerical_data(data = data, alpha = alpha, outpath = outpath))
}

analysis_after_pc_fun <- function(pc, outpath, caption) {
  return(analysis_after_pc(pc = pc, data, outpath = outpath, protein = protein, position_numbering = position_numbering, 
                           graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs, 
                           coloring = coloring, colors = colors, stages = stages, plot_types = plot_types, 
                           unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console, lines_in_abbr_of_r, 
                           compute_localTests_anew = compute_localTests_anew, print = print, plot = plot, caption = caption))
} 


parameters_for_info_file_fun <- function(alpha, outpath) {
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, only_cols = only_cols, coloring = coloring, colors = colors, outpath = outpath)  
}

# alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, seq(0.1, 0.9, 0.05))#, 0.99, 0.999, 0.9999, 0.99999)
# alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, seq(0.1, 0.9, 0.05), 0.99, 0.999, 0.9999, 0.99999)
alphas_vector <- c(1e-20, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, seq(0.1, 0.9, 0.1))#, 0.99, 0.999, 0.9999, 0.99999)

filename_w_o_alpha <- paste(only_cols_label, source_of_data, "-alpha=", sep = "")
output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", sep = "")

# optimize_pc_results(pc_fun, analysis_after_pc_fun, output_dir, filename_w_o_alpha, alpha_start, compute_pc_anew = compute_pc_anew)
alpha_results_NMR <- linear_optimization_of_pc_results(pc_fun = pc_fun, analysis_after_pc_fun = analysis_after_pc_fun, 
                                                       parameters_for_info_file_fun = parameters_for_info_file_fun, 
                                                       output_dir = output_dir, filename_w_o_alpha = filename_w_o_alpha, 
                                                       alphas_vector = alphas_vector, compute_pc_anew = compute_pc_anew, 
                                                       mute = FALSE, logscale = FALSE, plot_all_graphs = combined_plot)

dev.copy(postscript, paste(output_dir, "/", filename_w_o_alpha, "several.ps", sep = ""))
dev.off
# untested:
# plot_connected_components_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath)