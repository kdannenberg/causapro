protein_causal_graph <- function(data, protein, type_of_data, source_of_data, position_numbering, output_dir, filename, outpath,
                                 parameters_for_info_file, alpha, pc_solve_conflicts, pc_u2pd, pc_conservative, pc_maj_rule,
                                 caption, analysis, stages, plot_types, coloring, colors, 
                                 graph_layout = "dot", plot_as_subgraphs = plot_as_subgraphs, 
                                 plot_only_subgraphs = plot_only_subgraphs, unabbrev_r_to_info, print_r_to_console, 
                                 lines_in_abbr_of_r, compute_pc_anew, compute_localTests_anew, graph_output_formats,
                                 numerical, mute_all_plots = FALSE, plot_no_isolated_nodes) {
  print(paste("Output will be written to ", getwd(), "/", output_dir, "/...", sep = ""))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
    print("Directory created.")
  } 
  
  if (missing(outpath)) {
    outpath <- paste(output_dir, filename, sep = "/")
  }
  
  # Computation of pc 
  pc_fun <- function(outpath) {
    return(estimate_DAG_from_numerical_data(data, alpha = alpha, outpath = outpath, 
                                            solve_conflicts = pc_solve_conflicts, u2pd = pc_u2pd, 
                                            conservative = pc_conservative, maj_rule = pc_maj_rule))
  }
  
  results <- list()
  
  parameters_to_info_file(parameters_for_info_file, outpath)
  
  pc_func <- function_set_parameters(pc_fun, parameters = list(outpath = outpath))
  loaded_object_ok_fun <- function(pc) {return(length(pc@graph@nodes) != dim(data)[2])}
  results$pc <- compute_if_not_existent(filename = paste0(outpath, "-pc"), FUN = pc_func, obj_name = "pc", 
                                        compute_anew = compute_pc_anew, loaded_object_ok_fun = loaded_object_ok_fun)
  # results$pc <- get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file, data = data)
  
  # garbage <- graphics.off()
  # if (!mute_all_plots) {
  if (plot_no_isolated_nodes) {
    graph <- kernelize_graph(results$pc@graph)
  } else {
    graph <- results$pc@graph
  }
  plot_graph(graph = graph, caption = caption, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout, 
             coloring = coloring, colors = colors, outpath = outpath, numerical = numerical, plot_as_subgraphs = plot_as_subgraphs, 
             plot_only_subgraphs = plot_only_subgraphs, output_formats = graph_output_formats, mute_all_plots = mute_all_plots)
  # }
  
  return(results)
} 


protein_graph_clustering <- function(results, protein, outpath, file_separator, mute_all_plots, caption,
                                     cluster_methods, add_cluster_of_conserved_positions, 
                                     removed_cols, more_levels_of_conservedness = FALSE) {
  
  if (add_cluster_of_conserved_positions) {
    # node_clustering <- c(node_clustering, "#FFFFFF" = list(removed_cols))
    # names(node_clustering)[length(node_clustering)] <- "#FFFFFF"
    
    
    if (more_levels_of_conservedness) {
      removed_cols <- (removed_cols / max(removed_cols)) / 0.9999
    }
    colors_for_rem_pos <- color_by_effect(effects = removed_cols, int_pos = "", color_for_other_positions = "#000000", mode = "#FFFFFF")
    
    add_clusters <- sapply(names(table(colors_for_rem_pos)), 
                           FUN = function(color) {return(names(colors_for_rem_pos[which(colors_for_rem_pos == color)]))}, 
                           simplify = FALSE, USE.NAMES = TRUE)
  }
  
  for (clustering in cluster_methods) {
    igraph <- graph_from_graphnel(results$pc@graph)
    cluster_fct <- get(paste0("cluster_", clustering))
    cl <- cluster_fct(igraph)
    # TODO: save plot, instaed of plotting
    if (!mute_all_plots) {
      plot(cl, igraph, main = paste0(caption, "\n", clustering))
    }
    node_clustering <- groups(cl)
    node_clustering <- unname(node_clustering)  # otherwise interpreted as colors
    
    length_sort_clusters <- TRUE # TODO: rausziehen
    
    if (length_sort_clusters) {
      # if ((length(which(names(node_clustering) != "")) == 1)
      #     && names(node_clustering)[length(node_clustering)] != "") { # only one color defined, and that for the last cluster
      #   length_except_last <- vapply(node_clustering[1:(length(node_clustering) - 1)], length, 1L)
      #   node_clustering <- node_clustering[c(order(length_except_last, decreasing = TRUE), length(node_clustering))]
      # } else {
      node_clustering <- node_clustering[order(vapply(node_clustering, length, 1L), decreasing = TRUE)]
      # }
    }
    
    
    if (add_cluster_of_conserved_positions) {
      node_clustering <- c(node_clustering, add_clusters)
    }
    
    if (clustering == "edge_betweenness") {
      type <- "eb"
    } else if (clustering == "infomap") {
      type <- "im"
    } else {
      type <- "igraph"
    }
    
    plot_clusters_in_pymol(node_clustering = node_clustering, protein = protein, outpath = outpath, 
                           file_separator = file_separator, type_of_clustering = type)  
  }
}