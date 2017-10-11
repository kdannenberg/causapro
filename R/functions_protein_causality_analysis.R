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
  ## plot_graph(graph = graph, caption = caption, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout, 
  ##          coloring = coloring, colors = colors, outpath = outpath, numerical = numerical, plot_as_subgraphs = plot_as_subgraphs, 
  ##        plot_only_subgraphs = plot_only_subgraphs, output_formats = graph_output_formats, mute_all_plots = mute_all_plots)

  call_plot_igraph(g = graph, protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors, clusters = FALSE, caption = caption, outpath = outpath, output_formats = graph_output_formats, mute_all_plots = FALSE, layout_str = graph_layout)
  
  # }
  
  return(results)
} 


protein_graph_clustering <- function(results, protein, position_numbering, coloring, colors, outpath,output_formats, file_separator, mute_all_plots, caption,
                                     cluster_methods, add_cluster_of_conserved_positions, 
                                     removed_cols, more_levels_of_conservedness = FALSE, sort_clusters = length) {

  ## sort_clusters = "DDS-SVD") {
  if (add_cluster_of_conserved_positions) {
    # node_clustering <- c(node_clustering, "#FFFFFF" = list(removed_cols))
    # names(node_clustering)[length(node_clustering)] <- "#FFFFFF"
    
    
    if (more_levels_of_conservedness) {
      removed_cols <- (removed_cols / max(removed_cols)) / 0.9999
    }
    
    if (!length(removed_cols) == 0) {
      colors_for_rem_pos <- color_by_effect(effects = removed_cols, int_pos = "", color_for_other_positions = "#000000", mode = "#FFFFFF")
      
      add_clusters <- sapply(names(table(colors_for_rem_pos)), 
                             FUN = function(color) {return(names(colors_for_rem_pos[which(colors_for_rem_pos == color)]))}, 
                             simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  
  for (clustering in cluster_methods) {
    ## old conversion method
    ## igraph <- graph_from_graphnel(results$pc@graph)
    igraph <- igraph.from.graphNEL(results$pc@graph)
    cluster_fct <- get(paste0("cluster_", clustering))
    cl <- cluster_fct(igraph)
    ## TODO: save plot, instaed of plotting
    if (!mute_all_plots) {
      ## edge.arrow.size determines size of arrows (1 is default), vertex.size determines size of the vertices (15 is default), edge.width determines width of edges (1 is default)
      ## plot(cl, igraph, main = paste0(caption, "\n", clustering), edge.arrow.size=0.2, vertex.size=8, edge.width=0.7)
      call_plot_igraph(g = results$pc@graph, protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors, clusters = TRUE, clustering = cl, caption = caption, outpath = outpath, output_formats = output_formats, mute_all_plots = mute_all_plots)
      ##this is the old version, just in case my adjustments don't work for you
      ##plot(cl, igraph, main = paste0(caption, "\n", clustering))
    }
    node_clustering <- groups(cl)
    node_clustering <- unname(node_clustering)  # otherwise interpreted as colors
    
    
    
    if (!is.null(sort_clusters)) {
      # if ((length(which(names(node_clustering) != "")) == 1)
      #     && names(node_clustering)[length(node_clustering)] != "") { # only one color defined, and that for the last cluster
      #   length_except_last <- vapply(node_clustering[1:(length(node_clustering) - 1)], length, 1L)
      #   node_clustering <- node_clustering[c(order(length_except_last, decreasing = TRUE), length(node_clustering))]
      # } else {
      if (typeof(sort_clusters) == "closure" || typeof(sort_clusters) == "builtin") {
        node_clustering <- node_clustering[order(vapply(node_clustering, sort_clusters, 1L), decreasing = TRUE)]
      } else if (sort_clusters == "DDS-SVD") {
        DDS_SVD <- read_data(files = paste0(protein, "_DDS-SVD-1"))
        cluster_weights <- lapply(node_clustering, function(nodes) {
            nodes <- nodes[which((nchar(nodes) >= 3) && (nodes %in% colnames(DDS_SVD)))]
            len <- length(nodes)
            if (len > 0) {
            sum <- sum(DDS_SVD[,nodes])
            return(sum / len)
            } else {
              return(0)
            }
          })
        node_clustering <- node_clustering[order(unlist(cluster_weights), decreasing = TRUE)]
      }
    }
    
    
    if (add_cluster_of_conserved_positions && exists("add_clusters")) {
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
