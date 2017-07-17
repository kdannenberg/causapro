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


protein_causality_G <- function(
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
  protein = "PDZ",
  # 
  # type_of_data = "DG",
  # type_of_data = "DDG",
  type_of_data = "DDDG",
  # 
  # subtype_of_data = "all",
  subtype_of_data = "5",
  # subtype_of_data = "10",
  # 
  data_set = "",
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
  coloring = "es",#"auto",  # "auto-all" or "all"
  colors = NULL,
  # 
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 oder NULL
  # TODO Marcel: dafür sorgen dass, wenn diese Option aktiv ist, kein graphics.off() o.ä. ausgeführt wird (und nur der graph geplottet wird)
  combined_plot = FALSE,
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
  file_separator = "/"
) {
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
                         alpha = alpha, min_pos_var = min_pos_var, only_cols_label = only_cols_label, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                         pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule, file_separator = file_separator)
  
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
  
  caption <- get_caption(protein = protein, data = data_description, alpha = alpha, min_pos_var = min_pos_var, chars_per_line = 45) #TODO rem_gaps_threshold hinzufügen
  # TODO Marcel: add all the new ones
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, 
                                                       only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = "/"))  
  
  
  graph_computation <- graph_computation || evaluation || analysis
  # Computation of the Graph
  results <- list()
  
  if (graph_computation) {
    results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                         output_dir = output_dir, filename = filename, outpath = outpath, parameters_for_info_file = parameters_for_info_file,
                         alpha = alpha, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule,
                         caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                         graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                         unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                         compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                         graph_output_formats = graph_output_formats, numerical = numerical)
  }
  
  # Evaluation
  # if ("evaluation" %in% steps) {
  if (evaluation) {
    results <- analysis_after_pc(results$pc, data, outpath = outpath, protein = protein, position_numbering = position_numbering, 
                                 graph_layout = graph_layout, coloring = coloring, colors = colors,  stages = stages, 
                                 plot_types = plot_types, unabbrev_r_to_info = unabbrev_r_to_info, 
                                 print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r, 
                                 compute_localTests_anew = compute_localTests_anew, print = print_analysis, plot = plot_analysis, 
                                 caption = caption, graph_output_formats = graph_output_formats, combined_plot = combined_plot)
    # print_analysis <- FALSE
    # plot_analysis <- FALSE
  } 
  
  # Pymol
  if (graph_computation) {
    if (!is.null(plot_only_subgraphs)) {
      # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
      node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
      subgraphs <- subgraphs_from_node_clusters(node_clustering, graph = results$orig$graph$NEL, protein = protein)
      graph <- subgraphs[[plot_only_subgraphs]]$graph
    } else {
      graph <- results$pc@graph
    }
    
    plot_connected_components_in_pymol(protein = protein, position_numbering = position_numbering, graph = graph, 
                                     outpath = outpath, show_int_pos = TRUE, show_positions = FALSE, file_separator = file_separator)
  }
  
  
  # paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(314), to = c(383), all_paths = FALSE)
  # plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE, 
  #                     label = TRUE, show_positions = FALSE, file_separator = file_separator)
  
  # Analysis
  # if ("analysis" %in% steps) {
  if (analysis) {
    results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                   protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
                   amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                   pymol_bg_color = "grey",
                   barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)
  }
  
  if (data_in_results) {
    results$data <- data  
  }
  
  # TODO: alle Eingabeparameter in results$pars$<parametername> schreiben
  # if (parameters_in_results) {
  #   
  # }
  
  if (output_parameters_in_results) {
    results$caption = caption
    results$outpath = outpath
  }
  
  return(results)
}

# results_G <- protein_causality_G(pc_conservative = FALSE, pc_u2pd = "retry", pc_solve_confl = FALSE, analysis = TRUE)

# TODO: warum geht das nur für effects of pos. 372
# results_G <- protein_causality_G(pc_conservative = FALSE, pc_u2pd = "retry", pc_solve_confl = TRUE, analysis = TRUE)
results_G <- protein_causality_G(pc_conservative = FALSE, pc_u2pd = "relaxed", pc_solve_confl = TRUE, analysis = FALSE, weight_effects_on_by = "", min_pos_var = 0)



# results_G <- protein_causality_G(min_pos_var = 0.05)
# g <- results_G$pc@graph
# g_w <- set_edge_weights_for_graph(g, cov(data))
# # ida()
# print(wgtMatrix(g))
# print(wgtMatrix(g_w))
# 
# 
# print(paste("1,2: ", paste(ida(x.pos = 1, y.pos = 2, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("1,2: ", causalEffect(x = 1, y = 2, g = g_w)))
# 
# print(paste("1,3: ", paste(ida(x.pos = 1, y.pos = 3, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("1,3: ", causalEffect(x = 1, y = 3, g = g_w)))
# 
# print(paste("7,2: ", paste(ida(x.pos = 7, y.pos = 2, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("7,2: ", causalEffect(x = 7, y = 2, g = g_w)))
