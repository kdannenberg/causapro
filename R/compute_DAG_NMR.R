source("~/.configuration_code.R")

source_all_function_scripts()
# source("functions_causal_effects.R")
# source("functions_ci_tests.R")
# source("functions_compute_DAG_categorical.R")
# source("functions_compute_DAG_numerical.R")
# source("functions_conversions.R")
# source("functions_evaluate_DAG.R")
# source("functions_general.R")
# # source("functions_i_o.R")
# source("functions_linkcommunities.R")
# source("functions_pymol.R")
# source("functions_tools.R")

## set working directory to
source("configuration_data.R")

# data_set = "inact",


results_p38g <- protein_causality_p38g(data_set = "inact",
                                      alpha = 0.1, min_pos_var = 0, # show_variance_cutoff_plot = TRUE,
                                      ranked = FALSE, plot_no_isolated_nodes = FALSE, plot_with_graphviz = TRUE, 
                                      pymol_sort_connected_components_by_length = FALSE, pymol_mix_connected_components = TRUE,
                                      print_connected_components = TRUE, linkcommunities_k = 4,
                                      # data_set = "inact", only_cols = c("26", "78", "89", "109", "112", "170", "116", "119", "161"), only_cols_label = "Fig.7a.1"
                                      # data_set = "inact", only_cols = c("77", "78", "81", "89", "149", "343", "109", "112", "116", "170"), only_cols_label = "Fig.7a.2"
                                      # data_set = "inact", only_cols = c("77", "78", "81", "149", "174", "182", "186", "187", "201", "333"), only_cols_label = "Fig.7a.3"
                                      # data_set = "inact", only_cols = c("77", "78", "81", "170", "174", "182"), only_cols_label = "Fig.7b"
                                      # data_set = "act", only_cols = c("150", "187", "198", "250", "253", "262"), only_cols_label = "Fig.7c"
)
# sink()
# print(conflict_edges(results_p38g$pc@graph))




# only_cols = c("V26", "L76", "L89", "M109", "M112", "L170", "L116", "L119", "V161"), only_cols_label = "Fig.7a.1"
# only_cols = c("L76", "L77", "M81", "L89", "L149", "L343"), only_cols_label = "Fig.7a.2"
# only_cols = c("L77", "L78", "M81", "L149", "L174", "M182", "V186", "V187", "M201", "V333"), only_cols_label = "Fig.7a.3"
# only_cols = c("L77", "L78", "M81", "L170", "L174", "M182"), only_cols_label = "Fig.7b"
# only_cols = c("I150", "V187", "L198", "V250", "L253", "M262"), only_cols_label = "Fig.7c"

# sink.reset()

# ## Data parameters

# ## available data:
# ## TODO

# numerical = TRUE
# protein = "p38g"
# type_of_data = "NMR"

# ## state of the protein
# # state = "allstates"
# state = "wild_inact_apo"
# # state = "wild_ATP"
# # state = "wild_act"
# # state = "wild_BIRB796"
# # previously given as:
# # source_of_data = "p38g-NMR-wild_inact_apo" # Table S1
# # source_of_data = "p38g-NMR-wild_ATP" # Table S2
# # source_of_data = "p38g-NMR-wild_act" # Table S3
#                                         # source_of_data = "p38g-NMR-wild_BIRB796" # Table S4
# ## options for position numbering?
# position_numbering = ""

# if ((protein == "p38g") && (type_of_data == "NMR")) {
#   transpose=TRUE
# } else {
#   transpose=FALSE
# }

# ## Analysis parameters
# ## TODO: explain purpose of those parameters
# only_cols = NULL
# only_cols_label = ""
# ## level of significance
# alpha = 0.5 # 0.999999 infeasable
# ## if rank is TRUE instead of the numerical NMR data simply the ranking of the positions will be used
# rank = TRUE
# ## TODO: explain stages
# stages <- c("orig") # "sub"
# plot_types <- c("localTests", "graphs")

# ## Graphical parameters

# # plot all
# # par(mfrow = c(3,2))

# ## choose one of the layouts offered by Rgraphviz
# graph_layout <- "dot" # "dot", "circo", "fdp", "neato", "osage", "twopi"

# ## plot_as_subgraphs = TRUE
# ## NULL seems to be the only option here
# plot_only_subgraphs = NULL

# ## gives a variety of slightly different plots of the DAG
# ## what do they represent
# several_plots = FALSE

# if (several_plots) {
#   coloring = c("FS4", "FS4", "FS3-pie", "FS3-pie", "FS3-mix", "FS3-mix")
#   colors =  c(rep("", 4), 2, 2)
#   # coloring = c(
#   # "FS3-mix", "S3-mixed-manually-simple",
#   # "FS3-mix", "S3-mixed-manually",
#   # "FS3-pie", "S4")
#   # colors = c(1, "", 4, "", 
#   #   "", "")
#   plot_as_subgraphs = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
# } else {
#   # coloring = "none"
#   # colors = ""
  
#   # coloring = "auto"
#   # coloring = NULL
#   coloring = "FS3-pie" #FS3-mix" # "FS3-pie" #"FS3-mix"  # different Method
#   # coloring = "FS3-pie"
#   # coloring = "FS3-mixed-maunally-simple"
#   # coloring = "FS3-mixed-manually"
#   # coloring = "FS3-mix"
#   # coloring = "modules"
#   colors = 2
#   plot_as_subgraphs = FALSE
# }


# # colors = c("#FFD700", "#1874CD", "#CC0000",  "#69A019") # untested
# # colors = list("1" = "#69A019", # green
# #               "2" = "#1874CD") # blue

# # colors = ""

# ## Technical parameters (print, plot, save, analysis)
# ## those options are either set to FALSE or to TRUE (unused/used option)
# ## analyse DAG using dagitty
# analysis = FALSE
# ## do not print the dagitty analysis, but save it somewhere?
# print_analysis = TRUE
# ## plot the dagitty analysis
# plot_analysis = FALSE
# ## compute new dag/analysis (TRUE) or use precomputed one (FALSE)
# compute_pc_anew <- FALSE
# compute_localTests_anew <- FALSE
# ## what is this doing, more information to info file?
# unabbrev_r_to_info <- FALSE
# ## and this?
# print_r_to_console <- FALSE
# lines_in_abbr_of_r <- 10


# # Computation of Output-Location and Output Infos
# if (state == "allstates") {
#   data_list <- list()
#   for (state in c("wild_inact_apo", "wild_ATP", "wild_act", "wild_BIRB796")) {
#     source_of_data = paste(protein, type_of_data, state, sep = "-")
    
#     filename <- paste("~/Viren/R/Data/", source_of_data, ".csv", sep = "")
#     # filename <- paste("../Data/", source_of_data, ".csv", sep = "")
#     var <- read_data(filename, transpose = transpose)
#     rownames(var) <- paste(rownames(var), state, sep = "-")
#     data_list[state][[1]] <- var
#   }
#   data <- do.call(rbind, data_list)
# } else {
#   source_of_data = paste(protein, type_of_data, state, sep = "-")
#   filename <- paste("~/Viren/R/Data/", source_of_data, ".csv", sep = "")
#   # filename <- paste("../Data/", source_of_data, ".csv", sep = "")
#   data <- read_data(filename, transpose = transpose)
# }


# # TODO Marcel: use adjust_data as in _NMR-GTB
# if (rank) {
#   # data <- cbind(apply(data, 2, rank)) # positionsweise (über alle Obeservationen)
#   # type_of_data <- paste(type_of_data, "ranked-otherway", sep = "-")
#   data <- t(apply(data, 1, rank))  # observationsweise (über alle Positionen)
#   type_of_data <- paste(type_of_data, "ranked", sep = "-")
# }


# filename <- paste(only_cols_label, source_of_data, "-alpha=", alpha, sep = "")
# output_dir <- paste("~/Viren/R/Outputs/", protein, "/", type_of_data, "/", filename, sep = "")
# # output_dir <- paste("../Outputs/", protein, "/", type_of_data, "/", filename, sep = "")
# outpath <- paste(output_dir, filename, sep = "/")
 
# source_of_data <- paste(state, " (", type_of_data, ")", sep = "")

# caption <- caption(protein = protein, data = source_of_data, alpha = alpha, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
# parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  

# garbage <- graphics.off()
# if (several_plots) {
#   par(mfrow = c(3,2))
# }
# results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
#                                 output_dir = output_dir, filename = filename, parameters_for_info_file = parameters_for_info_file,
#                                 alpha = alpha, caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
#                                 graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
#                                 unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
#                                 compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
#                                 print_analysis = print_analysis, plot_analysis = plot_analysis)


# # td_analysis(graph = results$orig$graph$NEL, outpath = outpath)

# paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(16, 81), to = c(300, 337), all_paths = FALSE)
# plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE, label = TRUE)

# plot_connected_components_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath)


# # graphics.off()                                                                                                                                                          #yellow       #blue     #red       #green
# # base_colors <- c("#FFD700", "#1874CD", "#CC0000",  "#69A019")
# # cols <- compute_link_communities(results$orig$graph$NEL, k = 4, plot_bar_plot = FALSE, 
# #                                  classify_nodes = TRUE, pie_nodes = FALSE, color_edges = TRUE,
# #                                  round_categories = 1, base_colors = base_colors , protein = protein, 
# #                                  outpath = outpath)
# # ### cols <- compute_link_communities(results$orig$graph$NEL, k = 4, plot_bar_plot = FALSE, classify_nodes = TRUE, pie_nodes = FALSE, color_edges = TRUE,
# #                                 #round_categories = 2, colors = c("#FFD700", "#1874CD", "#CC0000",  "#69A019"))
# # 
# # 
# # coloring = "FS3-mix"
# # colors <- 1
# # node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
# # # node_clustering <- classify_nodes(cols, round_categories = 1, mix = TRUE, base_colors = base_colors) # das ist gecheated und überprüft nur die "linkComm"-Funktion
# # graph <- results$orig$graph$NEL

# # cluster_edges <- edges_between_clusters(graph = graph, clusters = node_clustering)
# # print(cluster_edges)
# # print(paste(sum(diag(cluster_edges)), ":", sum(cluster_edges) - sum(diag(cluster_edges)), "=", sum(diag(cluster_edges)) / (sum(cluster_edges) - sum(diag(cluster_edges)))))

# # "#FFD700" (yellow), 
# # "#FFFFFF" (white),
# # "#69A019" (green), 
# # "#1874CD" (blue),  
# # "#CC0000" (red), 
# # "#FF9933" (orange)
