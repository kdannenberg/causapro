insbesondere <- function(
  a,
  b,
  c
) {
  allgemein(a = a, b = ifelse(missing(b), missing, b), c = ifelse(missing(c), NA, c))
  # allgemein(a = a, ifelse(missing(b), NULL, b = b), c = ifelse(missing(c), NA, c))
}
  
allgemein <- function(
  a = 1,
  b = 2,
  c = 3
) {
  print(a)
  print(b)
  print(c)
}
  
# insbesondere(a=2)
  
  
insbesondere2 <- function(a = 'A', b = NULL) {
  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$a <- a
  argList$b <- b
  do.call(allgemein, argList)
}

# insbesndere2(a = 5)

protein_causality_test <- function(
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
  min_pos_var = NULL,
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
  # return(protein_causality(numerical = numerical,
  #                          protein = protein,
  #                          type_of_data = type_of_data,
  #                          subtype_of_data = subtype_of_data,
  #                          data_set = data_set,
  #                          position_numbering = position_numbering,
  #                          min_pos_var = min_pos_var,
  #                          only_cols = only_cols,
  #                          only_cols_label = only_cols_label,
  #                          alpha = alpha,
  #                          ranked = ranked,
  #                          pc_solve_conflicts = pc_solve_conflicts,
  #                          pc_u2pd = pc_u2pd,
  #                          pc_conservative = pc_conservative,
  #                          pc_maj_rule = pc_maj_rule,
  #                          weight_effects_on_by = weight_effects_on_by,
  #                          graph_output_formats = graph_output_formats,
  #                          graph_layout = graph_layout,
  #                          graph_layout_igraph = graph_layout_igraph,
  #                          coloring = coloring,
  #                          colors = colors,
  #                          plot_as_subgraphs = plot_as_subgraphs,
  #                          plot_only_subgraphs = plot_only_subgraphs,
  #                          for_combined_plot = for_combined_plot,
  #                          mute_all_plots = mute_all_plots,
  #                          other = other,
  #                          graph_computation = graph_computation,
  #                          evaluation = evaluation,
  #                          analysis = analysis,
  #                          stages = stages,
  #                          print_analysis = print_analysis,
  #                          plot_analysis = plot_analysis,
  #                          plot_types = plot_types,
  #                          plot_ida = plot_ida,                                  # NEW!
  #                          plot_clusters = plot_clusters,                              # NEW!
  #                          plot_no_isolated_nodes = plot_no_isolated_nodes,  # NEW!
  #                          plot_with_graphviz = plot_with_graphviz,
  #                          pymol_show_int_pos = pymol_show_int_pos,compute_pc_anew,    # NEW!
  #                          pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length, # NEW!
  #                          pymol_mix_connected_components = pymol_mix_connected_components,  # NEW!
  #                          print_connected_components = print_connected_components,    # NEW!
  #                          compute_pc_anew = compute_pc_anew,
  #                          compute_localTests_anew = compute_localTests_anew,
  #                          unabbrev_r_to_info = unabbrev_r_to_info,
  #                          print_r_to_console = print_r_to_console,
  #                          lines_in_abbr_of_r = lines_in_abbr_of_r,
  #                          data_in_results = data_in_results,
  #                          output_parameters_in_results = output_parameters_in_results,
  #                          ida_percentile = ida_percentile,
  #                          file_separator = file_separator
  # )
  # )
  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$numerical <- numerical
  argList$protein = protein
  argList$type_of_data = type_of_data
  argList$subtype_of_data = subtype_of_data
  argList$data_set = data_set
  argList$position_numbering = position_numbering
  argList$min_pos_var = min_pos_var
  argList$only_cols = only_cols
  argList$only_cols_label = only_cols_label
  argList$alpha = alpha
  argList$ranked = ranked
  
  do.call(protein_causality, argList)
}

protein_causality_test(analysis = FALSE, alpha = 0.01, # min_pos_var = 0.01,
                                       # ranked = FALSE, plot_no_isolated_nodes = TRUE, plot_with_graphviz = TRUE,
                                       # pymol_sort_connected_components_by_length = FALSE, pymol_mix_connected_components = TRUE,
                                       print_connected_components = TRUE)