library(pvclust)

protein_causality <- function(
  filename = NULL,
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
  type_of_variables = NULL,
  protein,
  #
  # type_of_data = "DG",
  # type_of_data = "DDG",
  type_of_data,
  #
  # subtype_of_data = "all",
  subtype_of_data = NULL,
  # subtype_of_data = "10",
  #
  data_set,
  # data_set = "SVD",
  # data_set = "372",
  transpose_data = FALSE,
  #
  position_numbering = "crystal",
  #
  # Analysis parameters
  # remove_positions_with_low_variance = TRUE,
  min_pos_var = 0,
  show_variance_cutoff_plot = FALSE,
  only_cols = NULL,
  only_cols_label = "",
  every_n_th_row = 1,
  #
  alpha = 0.01,
  ranked = FALSE,
  rank_obs_per_pos = FALSE,
  pc_indepTest = NULL,
  pc_suffStat = NULL,
  #
  # pc_solve_conflicts = FALSE,
  # pc_u2pd = "retry",
  cor_cov_FUN = "",
  # a STRNG; default:
  # cor_cov_FUN = "", i.e. cor in pc, cov in ida
  # cor_cov_FUN = "cov" -> cov everywhere
  # cor_cov_FUN = "cor" -> cor everywhere
  # cor_cov_FUN = "none"
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
  # graph_output_formats = c("ps", "svg", "pdf"),
  graph_output_formats = "pdf",
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
  # TODO Marcel: dafuer sorgen dass, wenn diese Option aktiv ist, kein graphics.off() o.ä. ausgefuehrt wird (und nur der graph geplottet wird)
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
  intervention_position = "all",
  causal_effects_function = "IDA-reset", # "IDA", "causalEffects"
  ida_direction = "of",
  ida_percentile = "11", # top 11
  # ida_percentile = 0.75, # top 75%
  #
  causal_analysis = FALSE,#!pc_solve_conflicts,
  max_conflict_edges = 11,
  linkcommunities = FALSE,
  linkcommunities_k = NULL,
  # linkcommunities_base_colors = ifelse(k==4, c("#FFD700", "#1874CD", "#CC0000",  "#69A019"), rainbow(linkcommunities_k)),
  linkcommunities_base_colors = NULL,
  effects_cluster_k = 3,
  effects_cluster_method = "pv",
  effects_hclust_method = "ward.D2",  #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
  effects_dist_method = "euclidean",
  effects_pv_nboot = 1000,
  effects_cluster_alpha = 0.95,
  stages = c("orig"), # c("orig", "sub"), # "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graphs"),
  plot_ida = FALSE,                                  # NEW!
  graph_clustering = FALSE,                          # NEW!
  plot_no_isolated_nodes = TRUE,  # TODO: make true possible even for edgeless -> empty graphs
  plot_with_graphviz = TRUE, # NEW!
  #
  pymol_show_int_pos = TRUE,                        # NEW!
  pymol_sort_connected_components_by_length = FALSE, # NEW!
  pymol_mix_connected_components = TRUE,           # NEW!
  #
  print_connected_components = FALSE,               # NEW!
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
  # output_parameters_in_results = FALSE,
  #
  file_separator = "/",
  graph_cluster_methods = c("edge_betweenness", "infomap"),
  add_cluster_of_conserved_positions = TRUE,
  compare_effects = FALSE
  ) {
  # INIT
  if (!(mute_all_plots || for_combined_plot)) {
    graphics.off()

    # TOunDO: graph_clustering_funktion bauen, da rein
    if (!((plot_analysis && causal_analysis) || (plot_ida && evaluation))) {
      if (graph_clustering) {
        cols <- length(graph_cluster_methods) + 1
      } else {
        cols <- 1
      }
      par(mfrow = c(1, cols))
    }
  }

  # Data file description
  if (missing(filename)) {
    filename <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set)
  }

  print(paste("Data:", filename))

  # if (!missing(filename)) {
  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  tryCatch(source(paste0("Data/", filename, ".R"), local = TRUE), error = function(cond) {
    message("Conf file could not be read, maybe it does not exist?")
  }, warning = function(cond) {
    message("Conf file could not be read, maybe it does not exist?")
  })
  list2env(argList, env = environment())
  # }
  # print(graph_layout)
  ## filename_data <- paste("Data/", source_of_data, ".csv", sep = "")
  ## include option to read_data from alignment
  start_with_alignment <- FALSE
  if (data_set == "bin_approx") {
    start_with_alignment <- TRUE
  }

  if (start_with_alignment) {
    ## further specification of alignments needed
    if (protein == "PDZ") {
      position_numbering <- "alignment"
    }
    alignment <- readAlignment(filename = paste0("Data", file_separator, protein, "_ALN"))
    data_orig <- compute_data_from_alignment(alignment = alignment, filename = filename)
  } else {
    data_orig <- read_data(filename, transpose = transpose_data, every_n_th_row = every_n_th_row)
  }



  data <- adjust_data(data = data_orig, rank = ranked, rank_obs_per_pos = rank_obs_per_pos, only_cols = only_cols,
                      min_var = min_pos_var,
                      keep_quadratic = (cor_cov_FUN == "none"),
                      mute_plot = !show_variance_cutoff_plot)

  #TODO:
  # data_set <- adjust_data_set(...)
  if (!(missing(every_n_th_row)) && !is.null(every_n_th_row) && !(every_n_th_row == 1)) {
    data_set <- pastes(data_set, paste("every", every_n_th_row, "th-row", sep = "-"), sep = "_")
  }
  # print(colnames(data))

  # Parameter checks and adjustment
  if (cor_cov_FUN == "none" && (dim(data)[1] != dim(data)[2])) {
    warning("Data matrix is not quadratic and can thus not be interpreted as a correlation/covariance matix.
            Option cor_cov_FUN == \"none\" is ignored.")
    cor_cov_FUN = ""
  }

  if (missing(type_of_variables)) {
    if (grepl("ALN", filename)) {
      type_of_variables <- "nominal"
    } else if (grepl("rank", filename) || ranked) {
      type_of_variables <- "ordinal"
    } else {
      type_of_variables <- "continuous"
    }
  }

  if (ranked) {
    type_of_variables <- "ordinal"
  }


  filename <- adjust_data_description(data_description = filename, ranked = ranked)


  removed_cols <- setdiff(colnames(data_orig), colnames(data))
  removed_cols <- apply(data_orig[, removed_cols, drop = FALSE], 2, var)

  subtype_of_data <- subtype_of_data_after_adjustment(subtype_of_data = subtype_of_data, rank = ranked, min_var = min_pos_var)

  if (!is.numeric(ida_percentile)) {
    ida_percentile <- 1 - (as.numeric(ida_percentile) / dim(data)[2])
  }

  # colnames(data) <- paste("X", colnames(data), sep = "")

  # if (other != "") {
  #   data_description <- get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set, suffix = other)
  # }

  outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                         data_set = data_set, suffix = other,
                         alpha = alpha, min_pos_var = min_pos_var, only_cols_label = only_cols_label,
                         pc_indepTest = pc_indepTest, cor_cov_FUN = cor_cov_FUN,
                         pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
                         pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule,
                         file_separator = file_separator)



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
  caption <- get_caption(protein = protein, data = filename, alpha = alpha, min_pos_var = min_pos_var, chars_per_line = chars_per_line) #TODO rem_gaps_threshold hinzufuegen
  # TODO Marcel: add all the new ones
  parameters_for_info_file <- parameters_for_info_file(protein = protein, type_of_data = type_of_data, alpha = alpha, position_numbering = position_numbering,
                                                       only_cols = only_cols, coloring = coloring, colors = colors, #outpath = paste(output_dir, filename, sep = file_separator)
                                                       outpath = outpath)

  parameters_to_info_file(parameters_for_info_file, outpath)

  graph_computation <- graph_computation || evaluation || causal_analysis
  # Computation of the Graph
  results <- list()
  results$summary <- list()

  if (data_in_results) {
    results$data <- data
  }
  results$summary$data_dim <- dim(data)

  # TODO: alle Eingabeparameter in results$pars$<parametername> schreiben
  # if (parameters_in_results) {
  #
  # }

  # if (output_parameters_in_results) {
    results$summary$caption <- caption
    results$summary$outpath <- outpath

    results$general$int_pos <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
  # }

  # CAUSAL STRUCTURE LEARNING
  if (graph_computation) {
    directories <- strsplit(outpath, file_separator)
    # filename <- directories[[1]][length(directories[[1]])]
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
    print(paste("Output will be written to ", getwd(), "/", output_dir, "/...", sep = ""))

    create_parent_directory_if_necessary(outpath)
    # if (!dir.exists(output_dir)) {
    #   dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
    #   print("Directory created.")
    # }
    #
    # if (missing(outpath)) {
    #   outpath <- paste(output_dir, filename, sep = "/")
    # }

    pc_func_w_o_alpha <- function_set_parameters(estimate_DAG_from_numerical_data,
                                       parameters = list(data = data, outpath = outpath,
                                                         cor_FUN = cor_cov_FUN, type_of_variables = type_of_variables,
                                                         solve_conflicts = pc_solve_conflicts, u2pd = pc_u2pd,
                                                         conservative = pc_conservative, maj_rule = pc_maj_rule))

    # TODO: rausziehen
    # alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
    # alphas <- c(1e-5)
    # alpha <- tune_alpha_bnlearn(data = data, pc_func_w_o_alpha = pc_func_w_o_alpha, alphas = alphas, outpath = outpath, compute_pc_anew = compute_pc_anew)

    pc_func <- function_set_parameters(pc_func_w_o_alpha, parameters = list(alpha = alpha))

    # TODO: gleich pc_func übergeben, alle anderen parameter rausnehmen
    results <- protein_causal_graph(results = results, data = data, protein = protein, type_of_data = type_of_data, source_of_data = source_of_data, position_numbering = position_numbering,
                                    output_dir = output_dir, filename = filename,
                                    outpath = outpath, type_of_variables = type_of_variables,
                                    indepTest = pc_indepTest, suffStat = pc_suffStat,
                                    alpha = alpha, cor_cov_FUN = cor_cov_FUN, pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule,
                                    # caption = caption,
                                    # analysis = causal_analysis, stages = stages, plot_types = plot_types, coloring = coloring, colors = colors,
                                    # graph_layout = graph_layout, graph_layout_igraph = graph_layout_igraph, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
                                    # unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                    compute_pc_anew = compute_pc_anew
                                    # , compute_localTests_anew = compute_localTests_anew,
                                    # graph_output_formats = graph_output_formats, numerical = numerical, mute_all_plots = mute_all_plots,
                                    # plot_no_isolated_nodes = plot_no_isolated_nodes, plot_with_graphviz = plot_with_graphviz
                                    )

    # numerical <- (type_of_variables == "continuous")
    # if (!mute_all_plots) {
      plot_pc(graph = results$pc@graph, caption = caption, outpath = outpath, protein = protein, position_numbering = position_numbering,
              plot_types = plot_types, coloring = coloring, colors = colors,
              graph_layout = graph_layout, graph_layout_igraph = graph_layout_igraph, plot_as_subgraphs = plot_as_subgraphs,
              plot_only_subgraphs = plot_only_subgraphs,unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console,
              lines_in_abbr_of_r = lines_in_abbr_of_r, compute_pc_anew = compute_pc_anew, compute_localTests_anew = compute_localTests_anew,
              graph_output_formats = graph_output_formats, numerical = TRUE, mute_all_plots = mute_all_plots,
              plot_no_isolated_nodes = plot_no_isolated_nodes, plot_with_graphviz = plot_with_graphviz)
    # }

    #### some statistics
    graph <- results$pc@graph

    number_of_edges <- sum(unlist(conflict_edges(graph)))
    number_of_nodes <- length(graph@nodes)

    results$summary$nodes$sum <- number_of_nodes
    results$summary$nodes$avg_deg <- number_of_edges/number_of_nodes
    results$summary$nodes$labels <- graph@nodes

    results$summary$edges$sum <- number_of_edges
    results$summary$edges <- c(results$summary$edges, conflict_edges(graph))
    ### end
  }

  # Evaluation
                                        # if ("evaluation" %in% steps) {
  if (evaluation) {
    results <- analysis_after_pc(results, data, outpath = outpath, protein = protein, position_numbering = position_numbering,
                                 stages = stages, #max_number_of_edges_in_graph = 70,
                                 unabbrev_r_to_info = unabbrev_r_to_info,
                                 print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r,
                                 compute_localTests_anew = compute_localTests_anew, print = print_analysis)
    ## print_analysis <- FALSE
    ## plot_analysis <- FALSE
    if (plot_analysis && !mute_all_plots) {
      if (is.null(results$orig$graph$dagitty)) { # wird in analysis_after_pc gesetzt, wenn der Graph nicht zu groß war
        plot_text("Too many edges in the graph. Graph evaluation (probably) infeasible.")
      } else {
          plot_structure_evaluation(results, stages, plot_types, graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
              coloring = coloring, colors = colors, caption = caption, outpath = outpath, graph_output_formats = graph_output_formats,
              combined_plot = for_combined_plot, position_numbering = position_numbering, protein = protein)
          plot_structure_evaluation(results, stages, plot_types, graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs,
              coloring = coloring, colors = colors, caption = caption, outpath = "", graph_output_formats = graph_output_formats,
              combined_plot = for_combined_plot, position_numbering = position_numbering, protein = protein)
      }
    }
  }

  # Pymol
  if (graph_computation) {
    # if (!is.null(clustering)) {
    if (graph_clustering) {
      protein_graph_clustering(results = results, protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors, outpath = outpath, output_formats = graph_output_formats, file_separator = file_separator,
                               caption = caption, mute_all_plots = mute_all_plots, cluster_methods = graph_cluster_methods,
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
      if (length(conn_comp) > 0) {
        cat("Connected components:")
        print(conn_comp)
      }
    }

    if (linkcommunities && sum(unlist(conflict_edges(results$pc@graph))) > 0) {
      cols <- compute_link_communities(results$pc@graph, k = linkcommunities_k, plot_bar_plot = FALSE,
                                       classify_nodes = TRUE, pie_nodes = FALSE, color_edges = TRUE,
                                       round_categories = 1, base_colors = linkcommunities_base_colors , protein = protein,
                                       outpath = outpath)

    }
  }



  # paths <- paths_between_nodes(graph = results$orig$graph$NEL, from = c(314), to = c(383), all_paths = FALSE)
  # plot_paths_in_pymol(protein = protein, graph = results$orig$graph$NEL, outpath = outpath, paths = paths, no_colors = FALSE,
  #                     label = TRUE, show_positions = FALSE, file_separator = file_separator)

  # Analysis
  # if ("analysis" %in% steps) {
  if (causal_analysis) {
    ida_function_w_o_pos_and_results <-
      function_set_parameters(causal_effects_ida,
                              parameters = list(data = data, direction = ida_direction,
                                                weight_effects_on_by = weight_effects_on_by,
                                                protein = protein, coloring = "all", no_colors = FALSE,
                                                outpath = outpath, amplification_exponent = 1,
                                                amplification_factor = TRUE, rank_effects = FALSE,
                                                effect_to_color_mode = "#FFFFFF", pymol_bg_color = "grey",
                                                mute_all_plots = mute_all_plots, caption = caption,
                                                show_neg_causation = TRUE, neg_effects = "sep",
                                                analysis = TRUE, percentile = ida_percentile,
                                                causal_effects_function = "IDA-reset", cov_FUN = cor_cov_FUN))
    # ida_function_w_o_pos_and_results
    ida_function_w_o_pos <- function_set_parameters(ida_function_w_o_pos_and_results,
                                                    parameters = list(results = results))

    # outpath = paste0(outpath, "_", causal_effects_function, "_interv-at-", intervention_position)
    outpath = paste0(outpath, "_", causal_effects_function)

    # plot.new()  # TODO: Muss das sein?!
    if (intervention_position == "all") {
      # all_pairwise_effects_FUN <- function_set_parameters(compute_if_not_existent, parameters = c(filename = paste(outpath, "pairwise_effects", sep="-"),
      #                                                                             # FUN = function_set_parameters(FUN, parameters = c(results = results)),
      #                                                                             FUN = FUN,
      #                                                                             obj_name = "all_pairwise_effects",
      #                                                                             fun_loaded_object_ok = function(all_pairwise_effects)
      #                                                                             {return(dim(all_pairwise_effects)[1] == dim(data)[2])}))
      #
      if (conflict_edges(results$pc@graph)$conflict == 0) {
        loaded_object_ok = function(all_pairwise_effects, data)
        {return(
          # dim(all_pairwise_effects)[1] == dim(data)[2] &&
          dim(all_pairwise_effects)[2] == dim(data)[2]
        )}

        fun_loaded_object_ok <- function_set_parameters(loaded_object_ok, parameters = list(data = data))


        all_pairwise_effects <- compute_if_not_existent(filename = paste(outpath, "pairwise_effects", sep="-"),
                                                        # FUN = function_set_parameters(FUN, parameters = c(results = results)),
                                                        FUN = function_set_parameters(compute_all_pairwise_effects,
                                                                                      parameters = list(data = data, ida_function_w_o_pos = ida_function_w_o_pos)),
                                                        obj_name = "all_pairwise_effects",
                                                        fun_loaded_object_ok = fun_loaded_object_ok)
                                                        # {return(FALSE)})
        results$all_pairwise_effects <- all_pairwise_effects

        results <- cluster_pairwise_effects(results = results, pairwise_effects = all_pairwise_effects, k = effects_cluster_k,
                                 cluster_method = effects_cluster_method, hclust_method = effects_hclust_method,
                                 dist_measure = effects_dist_method, iterations_pv = effects_pv_nboot, alpha = effects_cluster_alpha,
                                 protein = protein, outpath = outpath, file_separator = file_separator)

      } else {

        set_of_graphs <- determine_set_of_graphs(results = results, data = data, type_of_graph_set = "conflict",
                                                 pc_function = pc_function, ida_function = NULL,
                                                 s = s, new = FALSE, save = TRUE, outpath = outpath,
                                                 pc_maj_rule_conflict = pc_maj_rule_conflict,
                                                 pc_conservative_conflict = pc_conservative_conflict,
                                                 direction = ida_direction, # darf ws nicht "both" sein!
                                                 suffix_effects_type = paste0("pos-", intervention_position),
                                                 max_conflict_edges = max_conflict_edges, no_results = TRUE)

        all_graphs <- set_of_graphs#$graphs
        # all_results <- set_of_graphs$results


        # TODO

        # all_results <- pblapply(all_graphs, graph_to_results, ida_function = ida_function_w_o_results)   ## schneller (?) # library("pbapply")
        if (length(all_graphs) > 0) {
          all_results <- list()
          for (i in 1:length(all_graphs)) {
          # for (i in 1:2) {
            print(paste("GRAPH", i, "VON", length(all_graphs)))            ## mit Angabe des aktuellen Durchlaufs
            # all_results[[i]] <- list()
            # all_results[[i]]$pc <- list()
            # all_results[[i]]$pc@graph <- all_graphs[[i]]
            ida_function_w_o_pos <- function_set_parameters(ida_function_w_o_pos_and_results, parameters = list(results = all_graphs[[i]]))
            # all_pairwise_effects_G <-
            all_results[[i]] <- list()
            # all_results[[i]]$`ida` <- list()
            # all_results[[i]]$`ida`[[intervention_position]] <- list()
            # all_results[[i]]$`ida`[[intervention_position]][[ida_direction]] <- list()
            #
            # all_results[[i]]$`ida`[[intervention_position]][[ida_direction]]$`effects` <-
            all_results[[i]] <-
            compute_if_not_existent(filename = paste(outpath, "pairwise_effects-G", i, sep="-"),
                                                            # FUN = function_set_parameters(FUN, parameters = c(results = results)),
                                                            FUN = function_set_parameters(compute_all_pairwise_effects,
                                                              parameters = list(data = data, ida_function_w_o_pos = ida_function_w_o_pos,
                                                                                # apply_FUN = sapply,
                                                                                results_format = TRUE)),
                                                            obj_name = "all_pairwise_effects",
                                                            fun_loaded_object_ok = function(all_pairwise_effects)
                                                            # {return(dim(all_pairwise_effects)[1] == dim(data)[2])})
                                                            {return(TRUE)}, compute_anew = FALSE)

          }


          # all_pairwise_effects_over_all_graphs <- list()
          all_pairwise_effects_over_all_graphs <- matrix(ncol = 0, nrow = length(all_results[[1]]$ida))
          for (position in names(all_results[[1]]$ida)) {
            all_pairwise_effects_over_all_graphs <- cbind(all_pairwise_effects_over_all_graphs,
              unlist(compute_over_all_graphs(all_results = all_results, position = position,
                                      weight_effects_on_by = weight_effects_on_by,
                                      use_scaled_effects_for_sum = FALSE,
                                      #function_over_all_graphs = all_pairwise_effects_FUN,
                                      direction = "of",
                                      scale_effects_on = scale_effects_on_so_372_equals_1)))
            # mean_effects_min_max_FUN <- function_set_parameters(mean_effects_min_max, parameters = list(position = position, dir = direction))
            # min_max_mean_effects_on_of <- lapply(all_results, mean_effects_min_max_FUN, weight_effects_on_by = weight_effects_on_by, scaled_effects = use_scaled_effects_for_sum)

          }
          colnames(all_pairwise_effects_over_all_graphs) <- names(all_results[[1]]$ida)
          rownames(all_pairwise_effects_over_all_graphs) <- names(all_results[[1]]$ida)

          results$all_pairwise_effects <- all_pairwise_effects_over_all_graphs

          results <- cluster_pairwise_effects(results = results, pairwise_effects = all_pairwise_effects_over_all_graphs, k = effects_cluster_k,
                                   cluster_method = effects_cluster_method, hclust_method = effects_hclust_method,
                                   dist_measure = effects_dist_method, iterations_pv = effects_pv_nboot, alpha = effects_cluster_alpha,
                                   protein = protein, outpath = outpath, file_separator = file_separator)
        } else {
          if (!mute_all_plots) {
            plot_infeasible(caption = caption)
          }
        }
        # und jetzt?
        # results_copy <- results
        # results_copy$data <- data
        # results_copy$caption = caption
        # results_copy$outpath = outpath
        # analyse_set_of_graphs(results = results_copy, protein = "PDZ")
      }
      if (compare_effects) {
        compare_effects_per_position(results, hclust_method = effects_hclust_method,
                                     dist_measure = effects_dist_method, nboot = effects_pv_nboot)
      }
    } else {
      if (conflict_edges(results$pc@graph)$conflict == 0) {
        results <- ida_function_w_o_pos(perturbed_position = intervention_position)
      } else {
        results_copy <- results
        results_copy$data <- data
        # results_copy$summary$caption = caption
        # results_copy$summary$outpath = outpath
        ida_function_w_o_results <- function_set_parameters(ida_function_w_o_pos_and_results,
                                      parameters = list(perturbed_position = intervention_position))
                                      #direction = "both" nötig, damit mean berechnet werden kann
        # TODO: das funktioniert nicht mehr!
        # results$effects <- analyse_set_of_graphs(results = results_copy, type_of_data = type_of_data,
        #                                          subtype_of_data = subtype_of_data,
        #                                          measure = substr(type_of_data, nchar(type_of_data), nchar(type_of_data)),
        #                                          protein = "PDZ", perturbed_position = intervention_position,
        #                                          weight_effects_on_by = weight_effects_on_by,
        #                                          ida_function = ida_function_w_o_results,
        #                                          outpath = outpath, caption = caption)
        results$effects <- analyse_set_of_graphs(results = results_copy,
                                                 # type_of_data = type_of_data,
                                                 # subtype_of_data = subtype_of_data,
                                                 direction = ida_direction,
                                                 measure = substr(type_of_data, nchar(type_of_data), nchar(type_of_data)),
                                                 protein = "PDZ", perturbed_position = intervention_position,
                                                 weight_effects_on_by = weight_effects_on_by,
                                                 ida_function = ida_function_w_o_results,
                                                 outpath = outpath, caption = caption)
      }
    }
  }

  return(results)
}

protein_causality_G <- function(
  filename = NULL,
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
  type_of_variables = "continuous",
  protein = "PDZ",
  type_of_data = "DDG",
  subtype_of_data = "10",
  data_set = "",
  # data_set = "SVD",
  # data_set = "372",
  #
  position_numbering = "crystal",
  ## int_pos specific parameters
  ida_percentile = "11",
  # data-related parameters
  only_cols = NULL,
  only_cols_label = "",
  every_n_th_row = 1,
  # data-dependent graphical parameters
  graph_output_formats = NULL,
  plot_with_graphviz = TRUE,
  graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
  graph_layout_igraph = NULL,
  coloring = "auto", # "auto", "auto-all", "all"
  colors = NULL,
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 is another option
  # analysis parameters
  min_pos_var = 0,
  show_variance_cutoff_plot = NULL,
  ranked = FALSE,
  rank_obs_per_pos = FALSE,
  # analysis parameters: pc
  alpha = NULL,
  pc_indepTest = NULL,
  pc_suffStat = NULL,
  cor_cov_FUN = NULL,
  pc_solve_conflicts = TRUE,
  pc_u2pd = NULL,
  pc_conservative = NULL,
  pc_maj_rule = NULL,
  # analysis parameters: ida
  intervention_position = "372",
  causal_effects_function = NULL,
  ida_direction = NULL,
  plot_ida = FALSE,
  weight_effects_on_by = NULL, # "var", "mean", ""
  # analysis parameters: clustering of effects
  effects_cluster_k = NULL,
  effects_cluster_method = NULL,
  effects_hclust_method = NULL,  #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
  effects_dist_method = NULL,
  effects_pv_nboot = NULL,
  effects_cluster_alpha = NULL,
  # general graphical parameters
  for_combined_plot = NULL,
  mute_all_plots = NULL,
  other = NULL,
  # pymol
  pymol_show_int_pos = FALSE,
  pymol_sort_connected_components_by_length = NULL, # NEW!
  pymol_mix_connected_components = NULL,  # NEW!
  print_connected_components = NULL,
  # technical parameters
  graph_computation = NULL,
  evaluation = NULL,
  causal_analysis = NULL,
  max_conflict_edges = NULL,
  linkcommunities = NULL,
  linkcommunities_k = NULL,
  linkcommunities_base_colors = NULL,
  stages = NULL,
  print_analysis = NULL,
  plot_analysis = NULL,
  plot_types = NULL,
  plot_no_isolated_nodes = TRUE,
  graph_clustering = FALSE,
  compute_pc_anew = NULL,
  compute_localTests_anew = NULL,
  unabbrev_r_to_info = NULL,
  print_r_to_console = NULL,
  lines_in_abbr_of_r = NULL,
  data_in_results = NULL,
  output_parameters_in_results = NULL,
  file_separator = NULL
  ) {
  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$filename <- filename
  argList$type_of_variables <- type_of_variables
  argList$protein = protein
  argList$type_of_data = type_of_data
  argList$subtype_of_data = subtype_of_data
  argList$data_set = data_set
  argList$position_numbering = position_numbering
  argList$min_pos_var = min_pos_var
  argList$show_variance_cutoff_plot = show_variance_cutoff_plot
  argList$only_cols = only_cols
  argList$only_cols_label = only_cols_label
  argList$every_n_th_row = every_n_th_row
  argList$alpha = alpha
  argList$pc_indepTest = pc_indepTest
  argList$pc_suffStat = pc_suffStat
  argList$ranked = ranked
  argList$rank_obs_per_pos = rank_obs_per_pos
  argList$cor_cov_FUN = cor_cov_FUN
  argList$pc_solve_conflicts = pc_solve_conflicts
  argList$pc_u2pd = pc_u2pd
  argList$pc_conservative = pc_conservative
  argList$pc_maj_rule = pc_maj_rule
  argList$weight_effects_on_by = weight_effects_on_by
  argList$effects_cluster_k = effects_cluster_k
  argList$effects_cluster_method = effects_cluster_method
  argList$effects_hclust_method = effects_hclust_method
  argList$effects_dist_method = effects_dist_method
  argList$effects_pv_nboot = effects_pv_nboot
  argList$effects_cluster_alpha = effects_cluster_alpha
  argList$graph_output_formats = graph_output_formats
  argList$graph_layout = graph_layout
  argList$graph_layout_igraph = graph_layout_igraph
  argList$coloring = coloring
  argList$colors = colors
  argList$plot_as_subgraphs = plot_as_subgraphs
  argList$plot_only_subgraphs = plot_only_subgraphs
  argList$for_combined_plot = for_combined_plot
  argList$mute_all_plots = mute_all_plots
  argList$other = other
  argList$graph_computation = graph_computation
  argList$evaluation = evaluation
  argList$causal_analysis = causal_analysis
  argList$max_conflict_edges = max_conflict_edges
  argList$intervention_position = intervention_position
  argList$causal_effects_function = causal_effects_function
  argList$ida_direction = ida_direction
  argList$linkcommunities = linkcommunities
  argList$linkcommunities_k = linkcommunities_k
  argList$linkcommunities_base_colors = linkcommunities_base_colors
  argList$stages = stages
  argList$print_analysis = print_analysis
  argList$plot_analysis = plot_analysis
  argList$plot_types = plot_types
  argList$plot_ida = plot_ida                                 # NEW!
  argList$graph_clustering = graph_clustering                             # NEW!
  argList$plot_no_isolated_nodes = plot_no_isolated_nodes  # NEW!
  argList$plot_with_graphviz = plot_with_graphviz
  argList$pymol_show_int_pos = pymol_show_int_pos    # NEW!
  argList$pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length # NEW!
  argList$pymol_mix_connected_components = pymol_mix_connected_components  # NEW!
  argList$print_connected_components = print_connected_components    # NEW!
  argList$compute_pc_anew = compute_pc_anew
  argList$compute_localTests_anew = compute_localTests_anew
  argList$unabbrev_r_to_info = unabbrev_r_to_info
  argList$print_r_to_console = print_r_to_console
  argList$lines_in_abbr_of_r = lines_in_abbr_of_r
  argList$data_in_results = data_in_results
  argList$output_parameters_in_results = output_parameters_in_results
  argList$ida_percentile = ida_percentile
  argList$file_separator = file_separator

  do.call(protein_causality, argList)
}

protein_causality_S <- function(
  filename = NULL,
  # data parameters
  type_of_variables = "continuous",
  protein = "PDZ",
  type_of_data = "DDS",
  subtype_of_data = "",
  data_set = "",
  position_numbering = "crystal",
  ## int_pos specific parameters
  ida_percentile = "11",
  # data-related parameters
  only_cols = NULL,
  only_cols_label = "",
  every_n_th_row = 1,
  ## data-dependent graphical parameters
  graph_output_formats = NULL,
  plot_with_graphviz = TRUE,
  graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
  graph_layout_igraph = NULL,
  coloring = "auto", # "auto", "auto-all", "all"
  colors = NULL,
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 is another option
  ## analysis parameters
  min_pos_var = 0,
  show_variance_cutoff_plot = NULL,
  ranked = FALSE,
  rank_obs_per_pos = FALSE,
  ## analysis parameters: pc
  alpha = NULL,
  pc_indepTest = NULL,
  pc_suffStat = NULL,
  cor_cov_FUN = NULL,
  pc_solve_conflicts = TRUE,
  pc_u2pd = NULL,
  pc_conservative = NULL,
  pc_maj_rule = NULL,
  ## analysis parameters: ida
  causal_analysis = NULL,
  max_conflict_edges = NULL,
  plot_ida = FALSE,
  weight_effects_on_by = NULL, # "var", "mean", ""
  # analysis parameters: clustering of effects
  effects_cluster_k = NULL,
  effects_cluster_method = NULL,
  effects_hclust_method = NULL,  #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
  effects_dist_method = NULL,
  effects_pv_nboot = NULL,
  effects_cluster_alpha = NULL,
  ## general graphical parameters
  for_combined_plot = NULL,
  mute_all_plots = NULL,
  other = NULL,
  ## pymol
  pymol_show_int_pos = FALSE,
  pymol_sort_connected_components_by_length = NULL, # NEW!
  pymol_mix_connected_components = NULL, # NEW!
  print_connected_components = NULL,
  # technical parameters
  graph_computation = NULL,
  evaluation = NULL,
  intervention_position = "372",
  causal_effects_function = NULL,
  ida_direction = NULL,
  linkcommunities = NULL,
  linkcommunities_k = NULL,
  linkcommunities_base_colors = NULL,
  stages = NULL,
  print_analysis = NULL,
  plot_analysis = NULL,
  plot_types = NULL,
  plot_no_isolated_nodes = TRUE,
  graph_clustering = FALSE,
  compute_pc_anew = NULL,
  compute_localTests_anew = NULL,
  unabbrev_r_to_info = NULL,
  print_r_to_console = NULL,
  lines_in_abbr_of_r = NULL,
  data_in_results = NULL,
  output_parameters_in_results = NULL,
  file_separator = NULL
  ) {

  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$filename <- filename
  argList$type_of_variables <- type_of_variables
  argList$protein = protein
  argList$type_of_data = type_of_data
  argList$subtype_of_data = subtype_of_data
  argList$data_set = data_set
  argList$position_numbering = position_numbering
  argList$min_pos_var = min_pos_var
  argList$show_variance_cutoff_plot = show_variance_cutoff_plot
  argList$only_cols = only_cols
  argList$only_cols_label = only_cols_label
  argList$every_n_th_row = every_n_th_row
  argList$alpha = alpha
  argList$pc_indepTest = pc_indepTest
  argList$pc_suffStat = pc_suffStat
  argList$ranked = ranked
  argList$rank_obs_per_pos = rank_obs_per_pos
  argList$cor_cov_FUN = cor_cov_FUN
  argList$pc_solve_conflicts = pc_solve_conflicts
  argList$pc_u2pd = pc_u2pd
  argList$pc_conservative = pc_conservative
  argList$pc_maj_rule = pc_maj_rule
  argList$weight_effects_on_by = weight_effects_on_by
  argList$effects_cluster_k = effects_cluster_k
  argList$effects_cluster_method = effects_cluster_method
  argList$effects_hclust_method = effects_hclust_method
  argList$effects_dist_method = effects_dist_method
  argList$effects_pv_nboot = effects_pv_nboot
  argList$effects_cluster_alpha = effects_cluster_alpha
  argList$graph_output_formats = graph_output_formats
  argList$graph_layout = graph_layout
  argList$graph_layout_igraph = graph_layout_igraph
  argList$coloring = coloring
  argList$colors = colors
  argList$plot_as_subgraphs = plot_as_subgraphs
  argList$plot_only_subgraphs = plot_only_subgraphs
  argList$for_combined_plot = for_combined_plot
  argList$mute_all_plots = mute_all_plots
  argList$other = other
  argList$graph_computation = graph_computation
  argList$evaluation = evaluation
  argList$causal_analysis = causal_analysis
  argList$max_conflict_edges = max_conflict_edges
  argList$intervention_position = intervention_position
  argList$causal_effects_function = causal_effects_function
  argList$ida_direction = ida_direction
  argList$linkcommunities = linkcommunities
  argList$linkcommunities_k = linkcommunities_k
  argList$linkcommunities_base_colors = linkcommunities_base_colors
  argList$stages = stages
  argList$print_analysis = print_analysis
  argList$plot_analysis = plot_analysis
  argList$plot_types = plot_types
  argList$plot_ida = plot_ida                                 # NEW!
  argList$graph_clustering = graph_clustering                             # NEW!
  argList$plot_no_isolated_nodes = plot_no_isolated_nodes  # NEW!
  argList$plot_with_graphviz = plot_with_graphviz
  argList$pymol_show_int_pos = pymol_show_int_pos    # NEW!
  argList$pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length # NEW!
  argList$pymol_mix_connected_components = pymol_mix_connected_components  # NEW!
  argList$print_connected_components = print_connected_components    # NEW!
  argList$compute_pc_anew = compute_pc_anew
  argList$compute_localTests_anew = compute_localTests_anew
  argList$unabbrev_r_to_info = unabbrev_r_to_info
  argList$print_r_to_console = print_r_to_console
  argList$lines_in_abbr_of_r = lines_in_abbr_of_r
  argList$data_in_results = data_in_results
  argList$output_parameters_in_results = output_parameters_in_results
  argList$ida_percentile = ida_percentile
  argList$file_separator = file_separator

  do.call(protein_causality, argList)
}

protein_causality_p38g <- function(
  filename = NULL,
  # data parameters
  type_of_variables = "continuous",
  protein = "p38g",
  type_of_data = "NMR",
  subtype_of_data = "",
  data_set = "inact",
  transpose_data = TRUE,
  position_numbering = "crystal",
  # int_pos specific parameters
  ida_percentile = "11",
  # data-related parameters
  only_cols = NULL,
  only_cols_label = "",
  every_n_th_row = 1,
  # data-dependent graphical parameters
  graph_output_formats = NULL,
  plot_with_graphviz = TRUE,
  graph_layout = "dot", # "circo", "fdp", "neato", "osage", "twopi"
  graph_layout_igraph = NULL,
  coloring = "auto", # "auto", "auto-all", "all"
  colors = NULL,
  plot_as_subgraphs = FALSE,
  plot_only_subgraphs = NULL, # 1 is another option
  # analysis parameters
  min_pos_var = 0,
  show_variance_cutoff_plot = NULL,
  ranked = TRUE,
  rank_obs_per_pos = TRUE,  # Alvaro's way: TRUE
  # analysis parameters: pc
  alpha = NULL,
  pc_indepTest = NULL,
  pc_suffStat = NULL,
  cor_cov_FUN = NULL,
  pc_solve_conflicts = TRUE,
  pc_u2pd = NULL,
  pc_conservative = NULL,
  pc_maj_rule = NULL,
  # analysis parameters: ida
  causal_analysis = NULL,
  max_conflict_edges = NULL,
  plot_ida = FALSE,
  weight_effects_on_by = NULL, # "var", "mean", ""
  # analysis parameters: clustering of effects
  effects_cluster_k = NULL,
  effects_cluster_method = NULL,
  effects_hclust_method = NULL, #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
  effects_dist_method = NULL,
  effects_pv_nboot = NULL,
  effects_cluster_alpha = NULL,
  # general graphical parameters
  for_combined_plot = NULL,
  mute_all_plots = NULL,
  other = NULL,
  # pymol
  pymol_show_int_pos = FALSE,
  pymol_sort_connected_components_by_length = NULL, # NEW!
  pymol_mix_connected_components = NULL,  # NEW!
  print_connected_components = NULL,
  # technical parameters
  graph_computation = NULL,
  evaluation = NULL,
  intervention_position = "all",
  causal_effects_function = NULL,
  ida_direction = NULL,
  linkcommunities = TRUE,
  linkcommunities_k = 4,
  linkcommunities_base_colors = c("#FFD700", "#1874CD", "#CC0000",  "#69A019"),
  stages = NULL,
  print_analysis = NULL,
  plot_analysis = NULL,
  plot_types = NULL,
  plot_no_isolated_nodes = TRUE,
  graph_clustering = FALSE,
  compute_pc_anew = NULL,
  compute_localTests_anew = NULL,
  unabbrev_r_to_info = NULL,
  print_r_to_console = NULL,
  lines_in_abbr_of_r = NULL,
  data_in_results = NULL,
  output_parameters_in_results = NULL,
  file_separator = NULL
) {

  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$filename <- filename
  argList$type_of_variables <- type_of_variables
  argList$protein = protein
  argList$type_of_data = type_of_data
  argList$subtype_of_data = subtype_of_data
  argList$data_set = data_set
  argList$transpose_data = transpose_data
  argList$position_numbering = position_numbering
  argList$min_pos_var = min_pos_var
  argList$show_variance_cutoff_plot = show_variance_cutoff_plot
  argList$only_cols = only_cols
  argList$only_cols_label = only_cols_label
  argList$every_n_th_row = every_n_th_row
  argList$alpha = alpha
  argList$pc_indepTest = pc_indepTest
  argList$pc_suffStat = pc_suffStat
  argList$ranked = ranked
  argList$rank_obs_per_pos = rank_obs_per_pos
  argList$cor_cov_FUN = cor_cov_FUN
  argList$pc_solve_conflicts = pc_solve_conflicts
  argList$pc_u2pd = pc_u2pd
  argList$pc_conservative = pc_conservative
  argList$pc_maj_rule = pc_maj_rule
  argList$weight_effects_on_by = weight_effects_on_by
  argList$effects_cluster_k = effects_cluster_k
  argList$effects_cluster_method = effects_cluster_method
  argList$effects_hclust_method = effects_hclust_method
  argList$effects_dist_method = effects_dist_method
  argList$effects_pv_nboot = effects_pv_nboot
  argList$effects_cluster_alpha = effects_cluster_alpha
  argList$graph_output_formats = graph_output_formats
  argList$graph_layout = graph_layout
  argList$graph_layout_igraph = graph_layout_igraph
  argList$coloring = coloring
  argList$colors = colors
  argList$plot_as_subgraphs = plot_as_subgraphs
  argList$plot_only_subgraphs = plot_only_subgraphs
  argList$for_combined_plot = for_combined_plot
  argList$mute_all_plots = mute_all_plots
  argList$other = other
  argList$graph_computation = graph_computation
  argList$evaluation = evaluation
  argList$causal_analysis = causal_analysis
  argList$max_conflict_edges = max_conflict_edges
  argList$intervention_position = intervention_position
  argList$causal_effects_function = causal_effects_function
  argList$ida_direction = ida_direction
  argList$linkcommunities = linkcommunities
  argList$linkcommunities_k = linkcommunities_k
  argList$linkcommunities_base_colors = linkcommunities_base_colors
  argList$stages = stages
  argList$print_analysis = print_analysis
  argList$plot_analysis = plot_analysis
  argList$plot_types = plot_types
  argList$plot_ida = plot_ida                                 # NEW!
  argList$graph_clustering = graph_clustering                             # NEW!
  argList$plot_no_isolated_nodes = plot_no_isolated_nodes  # NEW!
  argList$plot_with_graphviz = plot_with_graphviz
  argList$pymol_show_int_pos = pymol_show_int_pos    # NEW!
  argList$pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length # NEW!
  argList$pymol_mix_connected_components = pymol_mix_connected_components  # NEW!
  argList$print_connected_components = print_connected_components    # NEW!
  argList$compute_pc_anew = compute_pc_anew
  argList$compute_localTests_anew = compute_localTests_anew
  argList$unabbrev_r_to_info = unabbrev_r_to_info
  argList$print_r_to_console = print_r_to_console
  argList$lines_in_abbr_of_r = lines_in_abbr_of_r
  argList$data_in_results = data_in_results
  argList$output_parameters_in_results = output_parameters_in_results
  argList$ida_percentile = ida_percentile
  argList$file_separator = file_separator

  do.call(protein_causality, argList)
}

protein_causality_NoV <- function(
  filename = NULL,
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

  type_of_variables = "continuous",
  protein = "NoV",
  # type_of_data = "NMR-Tit",
  type_of_data = "DDS",
  subtype_of_data = "",
  # subtype_of_data = c("Fuc", "BTS"),
  data_set = "",
  position_numbering = "",
  # analysis parameters
  min_pos_var = NULL,
  show_variance_cutoff_plot = NULL,
  only_cols = NULL,
  only_cols_label = "",
  every_n_th_row = 1,
  alpha = 0.05,
  pc_indepTest = NULL,
  pc_suffStat = NULL,
  ranked = FALSE,
  rank_obs_per_pos = FALSE,
  cor_cov_FUN = NULL,
  pc_solve_conflicts = TRUE,
  pc_u2pd = NULL,
  pc_conservative = NULL,
  pc_maj_rule = NULL,
  weight_effects_on_by = NULL, # "var", "mean", ""
  # analysis parameters: clustering of effects
  effects_cluster_k = NULL,
  effects_cluster_method = NULL,
  effects_hclust_method = NULL,  #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
  effects_dist_method = NULL,
  effects_pv_nboot = NULL,
  effects_cluster_alpha = NULL,
  # graphical parameters
  graph_output_formats = "pdf",
  graph_layout = NULL,
  graph_layout_igraph = NULL,
  coloring = NULL,
  colors = NULL,
  plot_as_subgraphs = NULL,
  plot_only_subgraphs = NULL, # 1 is another option
  plot_ida = NULL,                                  # NEW!
  graph_clustering = NULL,                              # NEW!
  plot_no_isolated_nodes = NULL,  # TODO: make true possible even for edgeless -> empty graphs #NEW
  for_combined_plot = NULL,
  mute_all_plots = NULL,
  other = NULL, # "cov"
  # technical parameters
  graph_computation = NULL,
  evaluation = NULL,
  causal_analysis = NULL,
  max_conflict_edges = NULL,
  intervention_position = "all",
  causal_effects_function = NULL,
  ida_direction = NULL,
  linkcommunities = NULL,
  linkcommunities_k = NULL,
  linkcommunities_base_colors = NULL,
  stages = NULL, # c("orig", "sub"), "sub"
  print_analysis = FALSE,
  plot_analysis = TRUE,
  plot_types = c("localTests", "graph"),
  plot_with_graphviz = FALSE,
  pymol_show_int_pos = FALSE,
  pymol_sort_connected_components_by_length = NULL, # NEW!
  pymol_mix_connected_components = NULL,  # NEW!
  print_connected_components = NULL,
  compute_pc_anew = NULL,
  compute_localTests_anew = NULL,
  unabbrev_r_to_info = NULL,
  print_r_to_console = NULL,
  lines_in_abbr_of_r = NULL,
  data_in_results = NULL,
  output_parameters_in_results = NULL,
  ida_percentile = "11",
  file_separator = NULL
) {
  argList <-  as.list(match.call(expand.dots = TRUE)[-1])
  # Enforce inclusion of non-optional arguments
  argList$filename <- filename
  argList$type_of_variables <- type_of_variables
  argList$protein = protein
  argList$type_of_data = type_of_data
  argList$subtype_of_data = subtype_of_data
  argList$data_set = data_set
  argList$position_numbering = position_numbering
  argList$min_pos_var = min_pos_var
  argList$show_variance_cutoff_plot = show_variance_cutoff_plot
  argList$only_cols = only_cols
  argList$only_cols_label = only_cols_label
  argList$every_n_th_row = every_n_th_row
  argList$alpha = alpha
  argList$pc_indepTest = pc_indepTest
  argList$pc_suffStat = pc_suffStat
  argList$ranked = ranked
  argList$rank_obs_per_pos = rank_obs_per_pos
  argList$cor_cov_FUN = cor_cov_FUN
  argList$pc_solve_conflicts = pc_solve_conflicts
  argList$pc_u2pd = pc_u2pd
  argList$pc_conservative = pc_conservative
  argList$pc_maj_rule = pc_maj_rule
  argList$weight_effects_on_by = weight_effects_on_by
  argList$effects_cluster_k = effects_cluster_k
  argList$effects_cluster_method = effects_cluster_method
  argList$effects_hclust_method = effects_hclust_method
  argList$effects_dist_method = effects_dist_method
  argList$effects_pv_nboot = effects_pv_nboot
  argList$effects_cluster_alpha = effects_cluster_alpha
  argList$graph_output_formats = graph_output_formats
  argList$graph_layout = graph_layout
  argList$graph_layout_igraph = graph_layout_igraph
  argList$coloring = coloring
  argList$colors = colors
  argList$plot_as_subgraphs = plot_as_subgraphs
  argList$plot_only_subgraphs = plot_only_subgraphs
  argList$for_combined_plot = for_combined_plot
  argList$mute_all_plots = mute_all_plots
  argList$other = other
  argList$graph_computation = graph_computation
  argList$evaluation = evaluation
  argList$causal_analysis = causal_analysis
  argList$max_conflict_edges = max_conflict_edges
  argList$intervention_position = intervention_position
  argList$linkcommunities = linkcommunities
  argList$linkcommunities_k = linkcommunities_k
  argList$linkcommunities_base_colors = linkcommunities_base_colors
  argList$stages = stages
  argList$print_analysis = print_analysis
  argList$plot_analysis = plot_analysis
  argList$plot_types = plot_types
  argList$plot_ida = plot_ida                                 # NEW!
  argList$graph_clustering = graph_clustering                             # NEW!
  argList$plot_no_isolated_nodes = plot_no_isolated_nodes  # NEW!
  argList$plot_with_graphviz = plot_with_graphviz
  argList$pymol_show_int_pos = pymol_show_int_pos    # NEW!
  argList$pymol_sort_connected_components_by_length = pymol_sort_connected_components_by_length # NEW!
  argList$pymol_mix_connected_components = pymol_mix_connected_components  # NEW!
  argList$print_connected_components = print_connected_components    # NEW!
  argList$compute_pc_anew = compute_pc_anew
  argList$compute_localTests_anew = compute_localTests_anew
  argList$unabbrev_r_to_info = unabbrev_r_to_info
  argList$print_r_to_console = print_r_to_console
  argList$lines_in_abbr_of_r = lines_in_abbr_of_r
  argList$data_in_results = data_in_results
  argList$output_parameters_in_results = output_parameters_in_results
  argList$ida_percentile = ida_percentile
  argList$file_separator = file_separator

  do.call(protein_causality, argList)
}
