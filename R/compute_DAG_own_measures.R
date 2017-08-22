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
source("functions_measures_MSA.R")
source("functions_pymol.R")
source("functions_tools.R")
## sets working directory depending on local config file
source("configuration_data.R")


# Technical parameters (print, plot, save, analysis)
compute_pc_anew <- FALSE
compute_localTests_anew <- FALSE
unabbrev_r_to_info = FALSE
stages <- c("orig", "anc") # "sub"
plot_types <- c("localTests", "graphs")
graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"
coloring = "auto"
# colors <- NULL
plot_as_subgraphs = TRUE
plot_only_subgraphs = 1

# General parameters
# gen_params <- get_data_params(source = "al_pdz")
numerical <- TRUE
source_of_data <- "al_pdz"
protein <- "PDZ"
position_numbering <- "alignment"
data_descr <- paste(type_of_data, " (", source_of_data, ")", sep = "")

# Computation parameters
type_of_data <- "C_ij_xy_add"
remove_cols_gaps_threshold = 0.2
k <- 2
alpha = 0.05


output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", sep = "")

if (remove_cols_gaps_threshold == 1) {
  filename_out = paste(source_of_data, "-k=", k, "-alpha=", alpha, sep="") 
  filename_data <- paste("../Data/", source_of_data, "-", type_of_data, "-k=", k, sep = "")
} else {
  filename_out = paste(source_of_data, "-ga_th=", remove_cols_gaps_threshold, "-k=", k, "-alpha=", alpha, sep="") 
  filename_data <- paste("../Data/", source_of_data, "-ga_th=", remove_cols_gaps_threshold, "-", type_of_data, "-k=", k, sep = "")
}

# filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", filename_out, sep = "")

# outpath <- paste(output_dir, filename, sep = "/")
# print(paste("Output will be written to ", getwd(), "/", substring(outpath, 0, nchar(outpath)), "...", sep = ""))
print(paste("Output will be written to ", getwd(), "/", output_dir, "/...", sep = ""))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  print("Directory created.")
}
outpath <- paste(output_dir, filename_out, sep = "/")

filename_MSA <- paste("../Data/", source_of_data, sep = "")

MSA <- readAlignment(filename_MSA)

if (!(remove_cols_gaps_threshold == 1)) {
  MSA <- remove_gaps(MSA, threshold=remove_cols_gaps_threshold, n_lev, allAS, outpath)
}


fun <- function() {
  return(mat_from_measure_i_j_mostfreqi_mostfreqj(MSA, k, get(type_of_data)))
}
data <- compute_if_not_existent(filename_data, FUN = fun)

colnames(data) <- colnames(MSA)
# colnames(data) <- paste("X", colnames(MSA), sep = "")   # necessary for evaluation with localTests


# Cijxy <- C_ij_xy_mat(MSA, k)

# Cijxy <- mat_from_measure_i_j_mostfreqi_mostfreqj(MSA, k, C_ij_xy)
# colnames(Cijxy) <- colnames(MSA)
# type_of_data <- "C_ij_xy"

# Cijxy_add <- mat_from_measure_i_j_mostfreqi_mostfreqj(MSA, k, C_ij_xy_adapted_add)
# colnames(Cijxy_add) <- colnames(MSA)
# type_of_data <- "C_ij_xy_add"

# Cijxy_mul <- mat_from_measure_i_j_mostfreqi_mostfreqj(MSA, k, C_ij_xy_adapted_mul)
# colnames(Cijxy_mul) <- colnames(MSA)
# type_of_data <- "C_ij_xy_mul"

# print("Data computed.")



# output_dir = paste("../Outputs/", type_of_data, sep = "")
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# }

caption <- caption(protein = protein, data = data_descr, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  

# Computation of pc 
# pc_fun <- function(outpath) {
#   return(estimate_DAG_from_numerical_data(data, type_of_data = type_of_data, alpha = alpha, outpath = outpath, protein = protein, position_numbering = position_numbering, colors = "auto"))
# }
# pc <- get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file)

# plot_graph(numerical = numerical, pc = pc, outpath = outpath, caption = caption, protein = protein, position_numbering = position_numbering, layout = graph_layout, coloring = coloring, colors = colors)

results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                                output_dir = output_dir, filename = filename, parameters_for_info_file = parameters_for_info_file,
                                alpha = alpha, caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                                graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                                unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                                print_analysis = print_analysis, plot_analysis = plot_analysis)
