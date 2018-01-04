protein_causality <- function(
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
  numerical,
  protein,
  # 
  # type_of_data = "DG",
  # type_of_data = "DDG",
  type_of_data,
  # 
  # subtype_of_data = "all",
  subtype_of_data,
  # subtype_of_data = "10",
  # 
  data_set,
  # data_set = "SVD",
  # data_set = "372",
  # 
  position_numbering = "crystal",
  # 
  # Analysis parameters
  # remove_positions_with_low_variance = TRUE,
  min_pos_var = 0.3,
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
  pc_maj_rule = TRUE,
  # 
  # weight_effects_on_by = "",
  # weight_effects_on_by = "var",
  # weight_effects_on_by = "mean",
  weight_effects_on_by = "median",
  # 
  # 
  # 
  # Graphical parameters
  graph_output_formats = c("ps", "svg"),
  ## graph_layout = "dot", # "dot", "circo", "fdp", "neato", "osage", "twopi"
  ## "layout_nicely" uses recommended layouts
  ## "layout_with_sugiyama" plots layered dags
  ## "layout_on_sphere" places the vertices uniformly on the surface of a spere (3d layout)
  ## "layout_with_lgl" is a layout for large graphs
  ## "layout_with_dh" uses simulated annealing for the graph layouting
  ## "layout_with_fr" uses a force-directed algorithm
  ## "layout_with_kk" is a physical model based on springs
  graph_layout_igraph = "layout_nicely",
  ## TODO: use graph_layout_graphviz instead
  graph_layout = "dot",
  coloring = "auto", #"es",#"auto",  # "auto-all" or "all"
  colors = NULL,
  # 
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 oder NULL
  # TODO Marcel: dafür sorgen dass, wenn diese Option aktiv ist, kein graphics.off() o.ä. ausgeführt wird (und nur der graph geplottet wird)
  # nur Graphen plotten
  for_combined_plot = FALSE,
  # gar ncihts plotten
  mute_all_plots = FALSE,
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
  plot_ida = FALSE,                                  # NEW!
  plot_clusters = TRUE,                              # NEW!
  plot_no_isolated_nodes = TRUE,  # TODO: make true possible even for edgeless -> empty graphs
  plot_with_graphviz = FALSE, # NEW!
  # 
  pymol_show_int_pos = TRUE,                        # NEW!
  pymol_sort_connected_components_by_length = TRUE, # NEW!
  pymol_mix_connected_components = FALSE,           # NEW!
  #
  print_connected_components = FALSE,                # NEW!
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
  ida_percentile = "11", # top 11
  # ida_percentile = 0.75, # top 75%
  #
  file_separator = "/",
  cluster_methods = c("edge_betweenness", "infomap"),
  add_cluster_of_conserved_positions = TRUE
) {
  # INIT
  if (!(mute_all_plots || for_combined_plot)) {
    graphics.off()
  
    if (!((plot_analysis && analysis) || (plot_ida && evaluation))) {
      if (plot_clusters) {
        cols <- length(cluster_methods) + 1
      } else {
        cols <- 1
      }
      par(mfrow = c(1, cols))
    }
  }
  
  data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set)
  
  # filename_data <- paste("Data/", source_of_data, ".csv", sep = "")
  data_orig <- read_data(data_description, only_cols = only_cols)
  data <- adjust_data(data = data_orig, rank = ranked, min_var = min_pos_var) #mute = combined_plot)
  data_description <- adjust_data_description(data_description = data_description, ranked = ranked)
  
  removed_cols <- setdiff(colnames(data_orig), colnames(data))
  removed_cols <- apply(data_orig[, removed_cols, drop = F], 2, var)
  
  subtype_of_data <- subtype_of_data_after_adjustment(subtype_of_data = subtype_of_data, rank = ranked, min_var = min_pos_var)
  
  if (!is.numeric(ida_percentile)) {
    ida_percentile <- 1 - (as.numeric(ida_percentile) / dim(data)[2])
  }
  
  # colnames(data) <- paste("X", colnames(data), sep = "")
  
  # if (other != "") {
  #   data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other)
  # }
  
  outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other,
                         alpha = alpha, min_pos_var = min_pos_var, only_cols_label = only_cols_label, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                         pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule, file_separator = file_separator)
  
  directories <- strsplit(outpath, file_separator)
  filename <- directories[[1]][length(directories[[1]])]
  output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
  
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
  
  chars_per_line <- 45
  if (for_combined_plot) {
    chars_per_line <- 35
  }
  caption <- get_caption(protein = protein, data = data_description, alpha = alpha, min_pos_var = min_pos_var, chars_per_line = chars_per_line) #TODO rem_gaps_threshold hinzufügen
  # TODO Marcel: add all the new ones
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering, 
                                                       only_cols = only_cols, coloring = coloring, colors = colors, outpath = paste(output_dir, filename, sep = file_separator))  
  
  
  graph_computation <- graph_computation || evaluation || analysis
  # Computation of the Graph
  results <- list()

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
  
  if (graph_computation) {
    results <- protein_causal_graph(data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering, 
                                    output_dir = output_dir, filename = filename, outpath = outpath, parameters_for_info_file = parameters_for_info_file,
                                    alpha = alpha, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule,
                                    caption = caption, analysis = analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors, 
                                    graph_layout = graph_layout, graph_layout_igraph = graph_layout_igraph, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                                    unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                    compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew, 
                                    graph_output_formats = graph_output_formats, numerical = numerical, mute_all_plots = mute_all_plots,
                                    plot_no_isolated_nodes = plot_no_isolated_nodes, plot_with_graphviz = plot_with_graphviz)
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
    # if (!is.null(clustering)) {
    if (plot_clusters) {
      protein_graph_clustering(results = results, protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors, outpath = outpath, output_formats = graph_output_formats, file_separator = file_separator,
                               caption = caption, mute_all_plots = mute_all_plots, cluster_methods = cluster_methods,
                               add_cluster_of_conserved_positions = add_cluster_of_conserved_positions, removed_cols = removed_cols)
    }
    
    if (!is.null(plot_only_subgraphs)) {
      # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
      node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
      subgraphs <- subgraphs_from_node_clusters(node_clustering, graph = results$orig$graph$NEL, protein = protein)
      graph <- subgraphs[[plot_only_subgraphs]]$graph
    } else {
      graph <- results$pc@graph
    }
    
    plot_connected_components_in_pymol(protein = protein, position_numbering = position_numbering, graph = graph, 
                                       outpath = outpath, show_int_pos = pymol_show_int_pos, show_positions = TRUE, 
                                       file_separator = file_separator, sort_connected_components_by_length = 
                                       pymol_sort_connected_components_by_length, 
                                       mix_connected_components = pymol_mix_connected_components)
    if (print_connected_components) {
      conn_comp <- nonsingular_connected_components(graph)
      print(conn_comp)
    }
  }

  
  
  # paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(314), to = c(383), all_paths = FALSE)
  # plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE, 
  #                     label = TRUE, show_positions = FALSE, file_separator = file_separator)
  
  # Analysis
  # if ("analysis" %in% steps) {
  if (analysis) {
    if (conflict_edges(results$pc@graph)$conflict == 0) {
      results <- causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                                  protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
                                  amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                                  pymol_bg_color = "grey",
                                  mute_all_plots = FALSE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
    } else {
      results_copy <- results
      results_copy$data <- data  
      results_copy$caption = caption
      results_copy$outpath = outpath
      analyse_set_of_graphs(results = results_copy, protein = "PDZ")
    }
  }
  
  return(results)
}

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
  type_of_data = "DDG",
  # type_of_data = "DDDG",
  # 
  # subtype_of_data = "all",
  # subtype_of_data = "5",
  subtype_of_data = "10",
  # 
  data_set = "",
  # data_set = "SVD",
  # data_set = "372",
  # 
  position_numbering = "crystal",
  # 
  # Analysis parameters
  # remove_positions_with_low_variance = TRUE,
  min_pos_var = 0.01,
  only_cols = NULL,
  only_cols_label = "",
  # 
  alpha = 0.03,
  ranked = FALSE,
  # 
  # pc_solve_conflicts = FALSE,
  # pc_u2pd = "retry",
  pc_solve_conflicts = TRUE,
  pc_u2pd = "relaxed",
  pc_conservative = FALSE,
  pc_maj_rule = TRUE,
  # 
  # weight_effects_on_by = "",
  # weight_effects_on_by = "var",
  # weight_effects_on_by = "mean",
  weight_effects_on_by = "median",
  # 
  # 
  # 
  # Graphical parameters
  graph_output_formats = "pdf",
  ## graph_layout = "dot", # "dot", "circo", "fdp", "neato", "osage", "twopi"
  ## "layout_nicely" uses recommended layouts
  ## "layout_with_sugiyama" plots layered dags
  ## "layout_on_sphere" places the vertices uniformly on the surface of a spere (3d layout)
  ## "layout_with_lgl" is a layout for large graphs
  ## "layout_with_dh" uses simulated annealing for the graph layouting
  ## "layout_with_fr" uses a force-directed algorithm
  ## "layout_with_kk" is a physical model based on springs
  graph_layout_igraph = "layout_with_lgl",
  graph_layout = "dot",
  coloring = "auto", #"es",#"auto",  # "auto-all" or "all"
  colors = NULL,
  # 
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 oder NULL
  # TODO Marcel: dafür sorgen dass, wenn diese Option aktiv ist, kein graphics.off() o.ä. ausgeführt wird (und nur der graph geplottet wird)
  for_combined_plot = FALSE,
  mute_all_plots = FALSE,
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
  plot_ida = FALSE,                                  # NEW!
  plot_clusters = TRUE,                              # NEW!
  plot_no_isolated_nodes = TRUE,
  plot_with_graphviz = FALSE,
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
  ida_percentile = "11", # top 11
  # ida_percentile = 0.75, # top 75%
  #
  file_separator = "/"
) {
  return(protein_causality(numerical = numerical,
                           protein = protein,
                           type_of_data = type_of_data,
                           subtype_of_data = subtype_of_data,
                           data_set = data_set,
                           position_numbering = position_numbering,
                           min_pos_var = min_pos_var,
                           only_cols = only_cols,
                           only_cols_label = only_cols_label,
                           alpha = alpha,
                           ranked = ranked,
                           pc_solve_conflicts = pc_solve_conflicts,
                           pc_u2pd = pc_u2pd,
                           pc_conservative = pc_conservative,
                           pc_maj_rule = pc_maj_rule,
                           weight_effects_on_by = weight_effects_on_by,
                           graph_output_formats = graph_output_formats,
                           graph_layout = graph_layout,
                           coloring = coloring,
                           colors = colors,
                           plot_as_subgraphs = plot_as_subgraphs,
                           plot_only_subgraphs = plot_only_subgraphs,
                           for_combined_plot = for_combined_plot,
                           mute_all_plots = mute_all_plots,
                           other = other,
                           graph_computation = graph_computation,
                           evaluation = evaluation,
                           analysis = analysis,
                           stages = stages,
                           print_analysis = print_analysis,
                           plot_analysis = plot_analysis,
                           plot_types = plot_types,
                           plot_ida = plot_ida,
                           plot_clusters = plot_clusters,
                           plot_no_isolated_nodes = plot_no_isolated_nodes,
                           plot_with_graphviz = plot_with_graphviz,
                           graph_layout_igraph = graph_layout_igraph,
                           compute_pc_anew = compute_pc_anew,
                           compute_localTests_anew = compute_localTests_anew,
                           unabbrev_r_to_info = unabbrev_r_to_info,
                           print_r_to_console = print_r_to_console,
                           lines_in_abbr_of_r = lines_in_abbr_of_r,
                           data_in_results = data_in_results,
                           output_parameters_in_results = output_parameters_in_results,
                           ida_percentile = ida_percentile,
                           file_separator = file_separator
  )
  )
}

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
  pc_maj_rule = TRUE,
  weight_effects_on_by = "median", # "var", "mean", ""
  # graphical parameters
  graph_output_formats = "ps",
  ## graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
  ## "layout_nicely" uses recommended layouts
  ## "layout_with_sugiyama" plots layered dags
  ## "layout_on_sphere" places the vertices uniformly on the surface of a spere (3d layout)
  ## "layout_with_lgl" is a layout for large graphs
  ## "layout_with_dh" uses simulated annealing for the graph layouting
  ## "layout_with_fr" uses a force-directed algorithm
  ## "layout_with_kk" is a physical model based on springs
  graph_layout_igraph = "layout_nicely",
  graph_layout = "dot",
  coloring = "auto", # "auto", "auto-all", "all"
  colors = NULL,
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 is another option
  for_combined_plot = FALSE,
  mute_all_plots = FALSE,
  other = "", # "cov"
  # technical parameters
  graph_computation = TRUE,
  evaluation = FALSE,
  analysis = FALSE, # !pc_solve_conflicts
  stages = c("orig"), # c("orig", "sub"), "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graph"),
  plot_ida = FALSE,                                  # NEW!
  plot_clusters = TRUE,                              # NEW!
  plot_with_graphviz = FALSE,
  #
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
  return(protein_causality(numerical = numerical,
                           protein = protein,
                           type_of_data = type_of_data,
                           subtype_of_data = subtype_of_data,
                           data_set = data_set,
                           position_numbering = position_numbering,
                           min_pos_var = min_pos_var,
                           only_cols = only_cols,
                           only_cols_label = only_cols_label,
                           alpha = alpha,
                           ranked = ranked,
                           pc_solve_conflicts = pc_solve_conflicts,
                           pc_u2pd = pc_u2pd,
                           pc_conservative = pc_conservative,
                           pc_maj_rule = pc_maj_rule,
                           weight_effects_on_by = weight_effects_on_by,
                           graph_output_formats = graph_output_formats,
                           graph_layout = graph_layout,
                           graph_layout_igraph = graph_layout_igraph,
                           coloring = coloring,
                           colors = colors,
                           plot_as_subgraphs = plot_as_subgraphs,
                           plot_only_subgraphs = plot_only_subgraphs,
                           for_combined_plot = for_combined_plot,
                           mute_all_plots = mute_all_plots,
                           other = other,
                           graph_computation = graph_computation,
                           evaluation = evaluation,
                           analysis = analysis,
                           stages = stages,
                           print_analysis = print_analysis,
                           plot_analysis = plot_analysis,
                           plot_types = plot_types,
                           plot_ida = plot_ida,
                           plot_clusters = plot_clusters,
                           plot_with_graphviz = plot_with_graphviz,
                           compute_pc_anew = compute_pc_anew,
                           compute_localTests_anew = compute_localTests_anew,
                           unabbrev_r_to_info = unabbrev_r_to_info,
                           print_r_to_console = print_r_to_console,
                           lines_in_abbr_of_r = lines_in_abbr_of_r,
                           data_in_results = data_in_results,
                           output_parameters_in_results = output_parameters_in_results,
                           ida_percentile = ida_percentile,
                           file_separator = file_separator
  )
  )
  
}

protein_causality_p38g <- function(
  # data parameters
  numerical = TRUE,
  protein = "p38g",
  type_of_data = "NMR",
  subtype_of_data = "",
  data_set = "act",
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
  pc_maj_rule = TRUE,
  weight_effects_on_by = "median", # "var", "mean", ""
  # graphical parameters
  graph_output_formats = "ps",
  graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
  coloring = "auto", # "auto", "auto-all", "all"
  colors = NULL,
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 is another option
  for_combined_plot = FALSE,
  mute_plot = FALSE,
  other = "", # "cov"
  # technical parameters
  graph_computation = TRUE,
  evaluation = FALSE,
  analysis = FALSE, # !pc_solve_conflicts
  stages = c("orig"), # c("orig", "sub"), "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graph"),
  plot_no_isolated_nodes = TRUE,
  plot_with_graphviz = FALSE,
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
  return(protein_causality(numerical = numerical,
                           protein = protein,
                           type_of_data = type_of_data,
                           subtype_of_data = subtype_of_data,
                           data_set = data_set,
                           position_numbering = position_numbering,
                           min_pos_var = min_pos_var,
                           only_cols = only_cols,
                           only_cols_label = only_cols_label,
                           alpha = alpha,
                           ranked = ranked,
                           pc_solve_conflicts = pc_solve_conflicts,
                           pc_u2pd = pc_u2pd,
                           pc_conservative = pc_conservative,
                           pc_maj_rule = pc_maj_rule,
                           weight_effects_on_by = weight_effects_on_by,
                           graph_output_formats = graph_output_formats,
                           graph_layout = graph_layout,
                           graph_layout_igraph = graph_layout_igraph,
                           coloring = coloring,
                           colors = colors,
                           plot_as_subgraphs = plot_as_subgraphs,
                           plot_only_subgraphs = plot_only_subgraphs,
                           for_combined_plot = for_combined_plot,
                           mute_all_plots = mute_all_plots,
                           other = other,
                           graph_computation = graph_computation,
                           evaluation = evaluation,
                           analysis = analysis,
                           stages = stages,
                           print_analysis = print_analysis,
                           plot_analysis = plot_analysis,
                           plot_types = plot_types,
                           plot_ida = plot_ida,                                  # NEW!
                           plot_clusters = plot_clusters,                              # NEW!
                           plot_no_isolated_nodes = plot_no_isolated_nodes,  # NEW!
                           plot_with_graphviz = plot_with_graphviz,
                           pymol_show_int_pos = pymol_show_int_pos,compute_pc_anew,    # NEW!
                           pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length, # NEW!
                           pymol_mix_connected_components = pymol_mix_connected_components,  # NEW!
                           print_connected_components = print_connected_components,    # NEW!
                           compute_pc_anew = compute_pc_anew,
                           compute_localTests_anew = compute_localTests_anew,
                           unabbrev_r_to_info = unabbrev_r_to_info,
                           print_r_to_console = print_r_to_console,
                           lines_in_abbr_of_r = lines_in_abbr_of_r,
                           data_in_results = data_in_results,
                           output_parameters_in_results = output_parameters_in_results,
                           ida_percentile = ida_percentile,
                           file_separator = file_separator
  )
  )
  
}

protein_causality_NoV <- function(
  # data parameters
  # available data:
  # "NoV_NMR-Tit_B4S"
  # "NoV_NMR-Tit_Fuc"
  # and
  # "NoV_NMR-Tit_B3S-with-unass"
  # type_of_data = "NMR-Tit"
  # type_of_data = c("NMR_Tit-Fuc", "NMR_Tit-BTS")
  # type_of_data = "Fuc-Tit-only-assigned"
  # subtype_of_data = "Fuc-old"
  # subtype_of_data = "Fuc"
  # TODO: wieder ermöglichen
  # subtype_of_data = c("Fuc", "BTS")
  
  numerical = TRUE,
  protein = "NoV",
  # type_of_data = "NMR-Tit",
  type_of_data = "DDS",
  subtype_of_data = "",
  # subtype_of_data = c("Fuc", "BTS"),
  data_set = "",
  position_numbering = "",
  # analysis parameters
  min_pos_var = 0,
  only_cols = NULL,
  only_cols_label = "",
  alpha = 0.05,
  ranked = FALSE,
  pc_solve_conflicts,
  pc_u2pd,
  pc_conservative,
  pc_maj_rule,
  weight_effects_on_by, # "var", "mean", ""
  # graphical parameters
  graph_output_formats = "pdf",
  graph_layout,
  graph_layout_igraph,
  coloring,
  colors,
  plot_as_subgraphs,
  plot_only_subgraphs, # 1 is another option
  plot_ida,                                  # NEW!
  plot_clusters,                              # NEW!
  plot_no_isolated_nodes,  # TODO: make true possible even for edgeless -> empty graphs #NEW
  for_combined_plot,
  mute_all_plots,
  other, # "cov"
  # technical parameters
  graph_computation,
  evaluation,
  analysis, # !pc_solve_conflicts
  stages, # c("orig", "sub"), "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graph"),
  plot_with_graphviz = FALSE,
  pymol_show_int_pos = FALSE,
  pymol_sort_connected_components_by_length, # NEW!
  pymol_mix_connected_components,  # NEW!
  print_connected_components,
  compute_pc_anew,
  compute_localTests_anew,
  unnabbrev_r_to_info,
  print_r_to_console,
  lines_in_abbr_of_r,
  data_in_results,
  output_parameters_in_results,
  ida_percentile,
  file_separator
) {
  return(protein_causality(numerical = numerical,
                           protein = protein,
                           type_of_data = type_of_data,
                           subtype_of_data = subtype_of_data,
                           data_set = data_set,
                           position_numbering = position_numbering,
                           min_pos_var = min_pos_var,
                           only_cols = only_cols,
                           only_cols_label = only_cols_label,
                           alpha = alpha,
                           ranked = ranked,
                           pc_solve_conflicts = pc_solve_conflicts,
                           pc_u2pd = pc_u2pd,
                           pc_conservative = pc_conservative,
                           pc_maj_rule = pc_maj_rule,
                           weight_effects_on_by = weight_effects_on_by,
                           graph_output_formats = graph_output_formats,
                           graph_layout = graph_layout,
                           graph_layout_igraph = graph_layout_igraph,
                           coloring = coloring,
                           colors = colors,
                           plot_as_subgraphs = plot_as_subgraphs,
                           plot_only_subgraphs = plot_only_subgraphs,
                           for_combined_plot = for_combined_plot,
                           mute_all_plots = mute_all_plots,
                           other = other,
                           graph_computation = graph_computation,
                           evaluation = evaluation,
                           analysis = analysis,
                           stages = stages,
                           print_analysis = print_analysis,
                           plot_analysis = plot_analysis,
                           plot_types = plot_types,
                           plot_ida = plot_ida,                                  # NEW!
                           plot_clusters = plot_clusters,                              # NEW!
                           plot_no_isolated_nodes = plot_no_isolated_nodes,  # NEW!
                           plot_with_graphviz = plot_with_graphviz,
                           pymol_show_int_pos = pymol_show_int_pos,compute_pc_anew,    # NEW!
                           pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length, # NEW!
                           pymol_mix_connected_components = pymol_mix_connected_components,  # NEW!
                           print_connected_components = print_connected_components,    # NEW!
                           compute_pc_anew = compute_pc_anew,
                           compute_localTests_anew = compute_localTests_anew,
                           unabbrev_r_to_info = unabbrev_r_to_info,
                           print_r_to_console = print_r_to_console,
                           lines_in_abbr_of_r = lines_in_abbr_of_r,
                           data_in_results = data_in_results,
                           output_parameters_in_results = output_parameters_in_results,
                           ida_percentile = ida_percentile,
                           file_separator = file_separator
  )
  )
  
}
