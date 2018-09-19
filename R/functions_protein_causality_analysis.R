protein_causal_graph <- function(results = list(), data, protein, type_of_data, source_of_data, position_numbering,
                                 output_dir, filename, outpath, type_of_variables,
                                 indepTest, suffStat, compute_pc_anew = FALSE,
                                 alpha, cor_cov_FUN = cov, pc_solve_conflicts, pc_u2pd, pc_conservative, pc_maj_rule) {



  # Computation of pc
  # pc_fun <- function(outpath) {
  #   return(estimate_DAG_from_numerical_data(data, alpha = alpha, outpath = outpath,
  #                                           cor_FUN = cor_cov_FUN, type_of_variables = type_of_variables,
  #                                           solve_conflicts = pc_solve_conflicts, u2pd = pc_u2pd,
  #                                           conservative = pc_conservative, maj_rule = pc_maj_rule))
  # }
  # pc_func <- function_set_parameters(pc_fun, parameters = list(outpath = outpath))

  pc_func <- function_set_parameters(estimate_DAG_from_numerical_data,
                                     parameters = list(data = data, alpha = alpha, outpath = outpath,
                                            cor_FUN = cor_cov_FUN, type_of_variables = type_of_variables,
                                            indepTest = indepTest, suffStat = suffStat,
                                            solve_conflicts = pc_solve_conflicts, u2pd = pc_u2pd,
                                            conservative = pc_conservative, maj_rule = pc_maj_rule))

  loaded_object_ok_fun <- function(pc) {return(length(pc@graph@nodes) == dim(data)[2])}
  results$pc <- compute_if_not_existent(filename = outpath,
                                        FUN = pc_func, obj_name = "pc", compute_anew = compute_pc_anew,
                                        fun_loaded_object_ok = loaded_object_ok_fun)

  return(results)
}

plot_pc <- function(graph, caption, outpath, protein, position_numbering, plot_types, coloring, colors,
                    graph_layout = "dot", graph_layout_igraph, plot_as_subgraphs = plot_as_subgraphs,
                    plot_only_subgraphs = plot_only_subgraphs, unabbrev_r_to_info, print_r_to_console,
                    lines_in_abbr_of_r, compute_pc_anew, compute_localTests_anew, graph_output_formats,
                    numerical, mute_all_plots = FALSE, plot_no_isolated_nodes, plot_with_graphviz) {

  if (plot_no_isolated_nodes) {
    if (!sum(unlist(conflict_edges(graph))) == 0) {
      graph <- kernelize_graph(graph)
    } else {
      if (!mute_all_plots) {
        plot_text(text = "No non-isolated nodes. (No Edges.)")
      }
    }
  }

  if(plot_with_graphviz) {
    plot_graph(graph = graph, caption = caption, protein = protein, position_numbering = position_numbering,
               graph_layout = graph_layout, coloring = coloring, colors = colors, outpath = outpath,
               numerical = numerical, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
               output_formats = graph_output_formats, mute_all_plots = mute_all_plots)
  } else {
    call_plot_igraph(g = graph, protein = protein, position_numbering = position_numbering,
                     coloring = coloring, colors = colors, clusters = FALSE, caption = caption,
                     outpath = outpath, output_formats = graph_output_formats,
                     mute_all_plots = mute_all_plots, layout_str = graph_layout_igraph,
                     plot_as_subgraphs = plot_as_subgraphs)
  }
}

