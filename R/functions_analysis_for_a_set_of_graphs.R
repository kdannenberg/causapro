
#
analyse_set_of_graphs <- function(
  type_of_graph_set = "conflict", # "retry" or "conflict"
  results,
  protein = "PDZ",
  measure, # = "G",
  # type_of_data = "DDDG",
  # subtype_of_data = "5",
  # data_set = "",
  # int_pos = interesting_positions(protein = protein, coloring = measure),
  coloring = "",
  int_pos = interesting_positions(protein = protein, coloring = coloring),
  # use_DDG = FALSE,  # einfach den default im Skript umstellen
  # protein_causality_function = get(paste0("protein_causality_", measure)),
  alpha = 0.1,
  min_pos_var = 0.01, # TODO: does not work for 0.01 (alpha = 0.01) (idafast for determination of median (effects on))
  # TODO MARCEL: Warum bekomme ich für alpha = 0.01, min_pos_var = 0.01 einen Graphen mit 6 conflict-Kanten, aber VOR remove_dummies 62 (!) Graphen?!
  new = FALSE,
  save = TRUE,
  check_graph_equality = FALSE, # takes forever
  pc_solve_conflicts = TRUE,
  # used when type_of_graph_set = "conflict"
  pc_maj_rule_conflict = TRUE,
  pc_conservative_conflict = FALSE,
  use_scaled_effects_for_each_graph = FALSE,
  scale_in_the_end = FALSE,
  plot_effect_quality = TRUE,
  plot_false_pos_neg = TRUE,
  plot_effect_score = TRUE,
  # effect_hue_by = "effects",
  effect_hue_by = get_conservation(measure = measure, protein = protein) / max(get_conservation(measure = measure, protein = protein)), # DG/DS
  # effect_hue_by = apply(data, 2, var),
  # results <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
  #                           pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
  #                           graph_computation = FALSE, evaluation = FALSE, causal_analysis = FALSE,
  #                           data_in_results = TRUE, output_parameters_in_results = TRUE)
  ida_percentile = 1 - (length(int_pos) / dim(data)[2]), # top 11
  # weight_effects_on_by = "",          # in der Summe ganz schlecht
  # weight_effects_on_by = "var",
  # weight_effects_on_by = "mean",
  weight_effects_on_by = "median",  # sieht (in der Summe) am besten aus
  # perturbed_position = "372", # ist in ida_FUN schon gesetzt
  causal_effects_function = "IDA-reset",
  scale_effects_on_so_372_equals_1 = TRUE,
  # if dir == "on": only effects on
  # if dir == "of": only effects of
  # if dir == c("on", "of") or  dir == "both": both speparately
  # if dir == "mean" : mean of effects on and of (first effects on are scaeld such that 372 has value 1)
  # for best graphs, dir = "both" and dir = "mean" yield the same
  direction = "mean",
  # direction = "both",
  # which part of the analysis should be plotted?
  ## IDA für alle s Graphen berechen und summieren (besser wäre vllt: mitteln, also nochmal durch s teilen. Kannst du gerne machen, Marcel.)
  plot = "over_all_graphs",
  function_over_all_graphs = "mean",
  ## IDA für alle s Graphen brechenen und denjenigen bestimmen, der am besten mit den gewünschten Ergebnissen übereinstimmt
  # plot = "best graph"
  ## für alle Graphen mit Nummern in <plot> die Abweichung von der Summe (dem zukünftigen Mittelwert) über alle Graphen bestimmen,
  ## also gewissermaßen wie repräsentativ der Graph jeweils für die Menge ist
  # plot = 1:25     ## if (is.numeric(plot) && length(plot) > 1) --> deviation from mean for graph(s) nr. <plot>
  # pc_function <- function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
  #   protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
  #                       pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
  #                       evaluation = evaluation, causal_analysis = analysis)
  # }

  # pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
  #   protein_causality_function(min_pos_var = min_pos_var, alpha = alpha,
  #                              pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
  #                              evaluation = evaluation, analysis = analysis, mute_all_plots = for_combined_plot)
  # },
  # pc_function = f_protein_causality_pc_parameters_eval_analysis(measure = measure,
  # pc_function = f_protein_causality_pc_parameters_eval_analysis(measure = measure, type_of_data = type_of_data, subtype_of_data = subtype_of_data,
  #                                                               min_pos_var = min_pos_var, alpha = alpha, mute_all_plots = TRUE),
  pc_function = function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                                                      min_pos_var = min_pos_var, alpha = alpha, mute_all_plots = TRUE,
                                                                                      for_combined_plot = TRUE, plot_clusters = FALSE)),
  # ida_function = function(results) {
  #   causal_effects_ida(data = data, perturbed_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
  #                      protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
  #                      amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
  #                      pymol_bg_color = "grey", caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE,
  #                      percentile = ida_percentile, mute_all_plots = for_combined_plot)
  # },
  # ida_function = f_causal_effects_ida_results(data = data, perturbed_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
  #                                             protein = protein, coloring = "all", no_colors = FALSE, outpath = outpath,
  #                                             amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
  #                                             pymol_bg_color = "grey", caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE,
  #                                             percentile = ida_percentile, mute_all_plots = for_combined_plot),
  effect_to_color_mode = "opacity", #"#FFFFFF",
  pymol_bg_color = "grey",
  show_neg_causation = TRUE,
  neg_effects_in_scaling = "sep",
  amplification_exponent = 1,
  amplification_factor = TRUE,
  rank_effects = FALSE,
  no_colors_in_pymol_ida = FALSE,
  ida_function = causal_effects_ida,
  s = 10,    # sample size
  for_combined_plot = FALSE,
  caption_as_subcaption = for_combined_plot,
  barplot_contour_black = TRUE,
  outpath,
  caption
) {
  ida_func <- function_set_parameters(ida_function,
                                          parameters = list(data = data, direction = "both", weight_effects_on_by = weight_effects_on_by,
                                                            protein = protein, coloring = coloring, no_colors = no_colors_in_pymol_ida, outpath = outpath,
                                                            amplification_exponent = amplification_exponent, amplification_factor = amplification_factor,
                                                            rank_effects = rank_effects, effect_to_color_mode = effect_to_color_mode,
                                                            pymol = FALSE, pymol_bg_color = pymol_bg_color, caption = caption, show_neg_causation = show_neg_causation,
                                                            neg_effects = neg_effects_in_scaling, analysis = TRUE, causal_effects_function = "IDA-reset",
                                                            percentile = ida_percentile, mute_all_plots = for_combined_plot))

  if (!missing(results)) {
    type_of_graph_set == "conflict"
    pc_u2pd = "relaxed"
    data <- results$data
    if (missing(caption)) {
      caption <- results$summary$caption
    }
    if (missing(outpath)) {
      outpath <- results$summary$outpath
    }
  } else {
    if (type_of_graph_set == "conflict") {
      pc_u2pd = "relaxed"
    } else if (type_of_graph_set == "retry") {
      pc_u2pd = "retry"
    }

    pc_function_dummy_version <- function_set_parameters(pc_function, parameters = list(
                                      # type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                      # min_pos_var = min_pos_var, alpha = alpha,
                                      pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
                                      pc_maj_rule = ifelse(type_of_graph_set == "retry", FALSE, pc_maj_rule_conflict),
                                      pc_conservative = ifelse(type_of_graph_set == "retry", FALSE, pc_conservative_conflict),
                                      graph_computation = FALSE,
                                      evaluation = FALSE, causal_analysis = FALSE, data_in_results = TRUE,
                                      # mute_all_plots = TRUE, # sollte schon vorher gesetzt sein
                                      for_combined_plot = TRUE))
    temp_results <- pc_function_dummy_version()
    # temp_results <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data, data_set = data_set,
    #                                     min_pos_var = min_pos_var, alpha = alpha,
    #                                     pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
    #                                     pc_maj_rule = ifelse(type_of_graph_set == "retry", FALSE, pc_maj_rule_conflict),
    #                                     pc_conservative = ifelse(type_of_graph_set == "retry", FALSE, pc_conservative_conflict),
    #                                     graph_computation = FALSE, evaluation = FALSE, causal_analysis = FALSE,
    #                                     data_in_results = TRUE,
    #                                     mute_all_plots = TRUE, for_combined_plot = TRUE)
    data <- temp_results$data
    caption <- temp_results$summary$caption
    outpath <- temp_results$summary$outpath
  }

  # TODO: besser, oder?
  # data <- read_data(get_data_description(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data))



  if (direction == "both") {
    direction <- c("on", "of")
  }

  if (!(("on" %in% direction) ||  ("of" %in% direction)) && !scale_effects_on_so_372_equals_1) {
    warning(paste("Effects on are not scaled. Effects on and of might not be comparably high. Computing the", direction, "might not be meaningful."))
    Sys.sleep(2)
  }




  set_of_graphs <- determine_set_of_graphs(results = results, data = data, type_of_graph_set = type_of_graph_set, pc_function = pc_function, ida_function = ida_func,
                                           s = s, new = new, save = save, outpath = outpath,
                                           pc_maj_rule_conflict = pc_maj_rule_conflict, pc_conservative_conflict = pc_conservative_conflict)
  if (is.null(set_of_graphs) && for_combined_plot) {
    if (caption_as_subcaption) {
      caption = caption
      # main_caption = NULL  # soll missing sein
    } else {
      # main_caption = caption
      caption = NULL
    }
    plot_infeasible(caption = caption)
    return(NULL)
  }
  all_graphs <- set_of_graphs$graphs
  all_results <- set_of_graphs$results


  # determine default colors for barplot
  all_one_effects <- as.matrix(all_results[[1]]$ida[[perturbed_position]][[1]]$effects[,1])
  all_one_effects[,] <- 1
  colors_for_barplot <- color_by_effect(all_one_effects, int_pos, mode = "#FFFFFF")

  # check wether any two graphs in the set are equal
  if (check_graph_equality) {
    equal <- compare_all_graphs(all_graphs)
    main_diagonal <- matrix(as.logical(diag(nrow = dim(equal)[1])), ncol = dim(equal)[2])
    print(which(equal & !main_diagonal, arr.ind = TRUE))
  }

  ## compare all the graphs w.r.t. how many conflict, bidirected and directed edges the have
  # conflicts_sorted <- sapply(all_graphs, conflict_edges)
  # conflicts_sorted <- matrix(unlist(conflicts_sorted), nrow = 3)
  # rownames(conflicts_sorted) = c("conflict", "directed", "bidirected")
  # colnames(conflicts_sorted) <- as.character(1:dim(conflicts_sorted)[2]) # Graph number
  # conflicts_sorted <- conflicts_sorted[, order(conflicts_sorted["bidirected",])]

  # print("Different types of edges in the different graphs (columns)")
  # print(conflicts_sorted)

  # cat("\n")
  # cat(paste0("alpha = ", alpha, ", min_pos_var = ", min_pos_var, "\n"))

  if (plot == "best_graph") {
    # cat("\n")
    best_graphs <- find_graphs_with_highest_int_pos(all_results = all_results, obj_fct = element_in_most_of_the_6_sets, dir = direction)
    cat("\n")

    for (i in best_graphs) {
      # i= 28
      cat(paste0("BEST GRAPH #", which(best_graphs == i), ": ", i))
      causal_effects_ida(data = data, perturbed_position = perturbed_position, direction = "both", weight_effects_on_by = weight_effects_on_by,
                         protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
                         amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                         pymol_bg_color = "grey", barplot = TRUE,
                         caption = c(paste("Graph", i), caption), show_neg_causation = TRUE, neg_effects = "sep", causal_analysis = TRUE, percentile = ida_percentile)

    }
  } else if (plot == "over_all_graphs") {
    # cat("\n")
    effects_over_all_graphs_on_of <- compute_over_all_graphs(all_results = all_results, position = perturbed_position,
                                                             weight_effects_on_by = weight_effects_on_by,
                                                             use_scaled_effects_for_sum = FALSE,
                                                             function_over_all_graphs = function_over_all_graphs, direction_type = direction,
                                                             scale_effects_on = scale_effects_on_so_372_equals_1)

    if (caption_as_subcaption) {
      caption = caption
      main_caption = NULL  # soll missing sein
    } else {
      main_caption = caption
      caption = NULL
    }

    effects_over_all_graphs_on_of <- display_effects(effects = effects_over_all_graphs_on_of, effect_hue_by = effect_hue_by,
                                                     direction = direction, int_pos = int_pos,
                                                     perturbed_position = as.numeric(perturbed_position), scale_in_the_end = scale_in_the_end,
                                                     weight_effects_on_by = weight_effects_on_by,
                                                     function_over_all_graphs = function_over_all_graphs,
                                                     ida_percentile = ida_percentile, caption = caption, main_caption = main_caption,
                                                     print = TRUE, plot = TRUE, plot_effect_quality = plot_effect_quality,
                                                     plot_false_pos_neg = plot_false_pos_neg, plot_effect_score = plot_effect_score,
                                                     for_combined_plot = for_combined_plot, barplot_contour_black = barplot_contour_black)

    pymol_mean_effects(effects_over_all_graphs_on_of, protein = protein, int_pos = int_pos, outpath = outpath, perturbed_position = perturbed_position,
                       amplification_exponent = amplification_exponent, amplification_factor = amplification_factor, rank_effects = rank_effects,
                       effect_hue_by = "effect", #"effect", # alles außer effect sinnlos
                       effect_to_color_mode = "opacity", # effect_to_color_mode, # opacity is best to see
                       pymol_bg_color = ifelse(effect_to_color_mode == "opacity", "grey", "black"), # pymol_bg_color, # if effect_to_color_mode == "opacity",
                       # grey is better; if effect_to_color_mode == "#FFFFFF" black makes more sense, since
                       # completely black nodes should not be seen
                       show_neg_causation = show_neg_causation,
                       neg_effects_in_scaling = neg_effects_in_scaling, no_colors = no_colors_in_pymol_ida)

    return(effects_over_all_graphs_on_of)

  } else if (is.numeric(plot) && length(plot) > 1) {
    deviation_from_mean(all_results = all_results, dir = "of", weight_effects_on_by = weight_effects_on_by, plot_graphs = plot)
  }
}

pymol_mean_effects <- function(effects_over_all_graphs_on_of, protein, int_pos, outpath, perturbed_position = "372", amplification_exponent,
                               amplification_factor, rank_effects = FALSE, effect_hue_by = "effect", effect_to_color_mode = "opacity",
                               pymol_bg_color = "black", show_neg_causation = TRUE, neg_effects_in_scaling = "sep", no_colors) {

  for (slot in names(effects_over_all_graphs_on_of)) {
    descr_split <- str_split(string = slot, pattern = "_")[[1]]
    dir <- str_replace(slot, paste0(descr_split[1], "_"), "")

    effects <- as.matrix(effects_over_all_graphs_on_of[[slot]])

    pymol_outpath <- outpath_for_ida(outpath = outpath, direction = dir, option_nr = "", neg_effects = neg_effects_in_scaling,
                                       perturbed_position = perturbed_position, amplification_exponent = amplification_exponent,
                                       amplification_factor = amplification_factor,
                                       no_colors = no_colors, rank_effects = rank_effects, effect_to_color_mode = effect_to_color_mode)



    # if (effect_hue_by == "effect") {
    #   colors_by_effect <- color_by_effect(scaled_effects, int_pos, mode = effect_to_color_mode)
    # } else if (effect_hue_by == "variance" || effect_hue_by == "var") {
    #   vars <- apply(data, 2, var)
    #   colors_by_effect <- color_by_effect(vars, int_pos, mode = effect_to_color_mode)
    # }


    if (is.character(effect_hue_by)) {
      if (effect_hue_by == "effect") {
        current_effect_hue_by <- effects
      } else {
        warning(paste0("The pymol file ", pymol_outpath, " does not show the effects but ", effect_hue_by, "!"))
        current_effect_hue_by <- get(effect_hue_by)
      }
    } else {
      current_effect_hue_by <- effect_hue_by
      warning(paste0("The pymol file ", pymol_outpath, " does not show the effects but what was explicitly given to the function in effect_hue_by!"))
    }


    # scaled_effects <- scale_effects(effect_hue_by, rank = rank_effects, amplification_factor = amplification_factor, neg_effects = neg_effects)
    colors_by_effect <- color_by_effect(as.matrix(current_effect_hue_by), int_pos, mode = effect_to_color_mode)


    plot_total_effects_in_pymol(positions_with_colors_by_effect = colors_by_effect, perturbed_position = perturbed_position,
                                protein = protein, outpath = pymol_outpath, amplification_exponent = amplification_exponent,
                                amplification_factor = amplification_factor, ranked = rank_effects,
                                index = "", no_colors = no_colors, bg_color = pymol_bg_color, orig_effects = effects)
  }
}

# data is only necessary to check wether a loades graph has the right number of nodes
determine_set_of_graphs <- function(results, data, type_of_graph_set, pc_function, ida_function, s, new, save, outpath,
                                    pc_maj_rule_conflict, pc_conservative_conflict, suffix_effects_type = "",
                                    suffix_graphs = "graphs", suffix_results = "results-ida-reset",
                                    max_conflict_edges = 11, no_results = FALSE) {
  start_new <- new
  if (type_of_graph_set == "retry") {
    suffix_retry_conflict = "-pc-retry_"
    suffix_results = pastes(suffix_effects_type, suffix_results, sep = "_")
    # infix <- pastes("-pc-retry", suffix_effects_type, sep = "_")
    # infix <- paste0(infix, "_")
    filename_graphs <- paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData")
    filename_results <- filename_results
    # if (suffix_effects_type == "" || suffix_effects_type == "372") {
      outpath_where_graphs_exist <- get_old_outpath(outpath, suffix = suffix_retry_conflict, suffix_graphs, ".RData")
      outpath_where_results_exist <- get_old_outpath(outpath, suffix = suffix_retry_conflict, suffix_results, ".RData")
    # } else {
    #   outpath_where_graphs_exist <- NULL
    #   outpath_where_results_exist <- NULL
    # }
    if (is.null(outpath_where_graphs_exist) || is.null(outpath_where_results_exist)) {
      start_new <- TRUE
    }
    if (start_new) {
      ## all_graphs saved
      all_results <- list()
      all_graphs <- list()
      for (i in 1:s) {
        print(paste("DURCHLAUF", i))
        # source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
        results <- pc_function(pc_solve_conflicts = FALSE, pc_u2pd = "retry", pc_maj_rule = FALSE,
                               pc_conservative = FALSE, evaluation = FALSE, causal_analysis = FALSE,
                               mute_all_plots = mute_all_plots)
          # protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
                                       # pc_solve_conflicts = FALSE, pc_u2pd = "retry",
                                       # evaluation = FALSE, analysis = FALSE)
        edges <- conflict_edges(results$pc@graph)
        all_results[[i]] <- results
        all_graphs[[i]] <- results$pc@graph
        if ((edges$conflict == 0) && (edges$bidirected == 0)) {
          break
        }
      }
      if (save) {
        save(all_graphs, file = filename_graphs)
      }

      if (no_results) {
        return(all_graphs)
      }

      #### save(all_results, file = filename_results)


      #### all_results saved
      for (i in 1:s) {
        print(paste("DURCHLAUF", i))
        all_results[[i]] <- ida_function(results)
      }
      if (save) {
        save(all_results, file = filename_results)
      }
    } else {
      load(file = outpath_where_graphs_exist)
      if (!file.exists(filename_graphs)) {
        save(all_graphs, file = filename_graphs)
      }

      if (no_results) {
        return(all_graphs)
      }

      load(file = outpath_where_results_exist)
      if (!file.exists(filename_results)) {
        save(all_results, file = filename_results)
      }

      # if (file.exists(file = filename_results)
      #     && file.exists(file = paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData"))) {
      #   load(file = paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData"))
      #   load(file = filename_results)
      # } else {
      #   load(file = paste0(get_old_outpath(outpath), suffix_retry_conflict, suffix_graphs, ".RData"))
      #   save(all_graphs, file = paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData"))
      #   load(file = paste0(get_old_outpath(outpath), suffix_retry_conflict, suffix_results, ".RData"))
      #   save(all_results, file = filename_results)
      # }
    }
    # TODO: rename: -rel_conflict_graph_set
  } else if (type_of_graph_set == "conflict") {
    suffix_retry_conflict = "-all_confl_comb_"
    suffix_results = pastes(suffix_effects_type, suffix_results, sep = "_")
    # infix <- pastes("-pc-retry", suffix_effects_type, sep = "_")
    # infix <- paste0(infix, "_")
    filename_graphs <- paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData")
    filename_results <- paste0(outpath, suffix_retry_conflict, suffix_results, ".RData")
    # if (suffix_effects_type == "" || suffix_effects_type == "372") {
      outpath_where_graphs_exist <- get_old_outpath(outpath, suffix = paste0(suffix_retry_conflict, suffix_graphs, ".RData"))
      outpath_where_results_exist <- get_old_outpath(outpath, suffix = paste0(suffix_retry_conflict, suffix_results, ".RData"))
    # } else {
    #   outpath_where_graphs_exist <- NULL # TODO: so wird es ws. NIE geladen
    #   outpath_where_results_exist <- NULL
    # }
    if (is.null(outpath_where_graphs_exist) || is.null(outpath_where_results_exist)) {
      start_new <- TRUE
    }

    if (start_new) {
      if (!new && !is.null(outpath_where_graphs_exist)) {
        try(load(file = outpath_where_graphs_exist))
        # if (outpath_where_graphs_exist != paste0(outpath, suffix_retry_conflict, suffix_graphs, ".RData")) oder
        if (!exists("all_graphs")) {
          warning(paste("File not loadable or did not contain an object of name all_graphs!"))
          outpath_where_graphs_exist <- NULL
        } else if (!is.null(all_graphs) && length(all_graphs[[1]]@nodes) != dim(data)[2]) {
          warning(paste0("Loaded graph did not have the right number of nodes: ", outpath_where_graphs_exist, "."))
          outpath_where_graphs_exist <- NULL
        } else {
          if (!file.exists(filename_graphs)) {
            save(all_graphs, file = filename_graphs)
          }
          if (no_results) {
            return(all_graphs)
          }
        }
      }
      if (new || is.null(outpath_where_graphs_exist)) {
        if (missing(results)) {
          results <- pc_function(pc_solve_conflicts = TRUE, pc_u2pd = "relaxed", pc_maj_rule = pc_maj_rule_conflict, pc_conservative = pc_conservative_conflict, evaluation = FALSE, causal_analysis = FALSE)
        }

        edges <- conflict_edges(results$pc@graph)
        print(edges)

        # Sys.sleep(2)

        if (edges$conflict > max_conflict_edges) {  # previously: >=15!!
          # all_graphs = NULL
          # save(all_graphs, file = filename_graphs)
          warning(paste("More than", max_conflict_edges, "conflict edges. Regarded infeasible."))
          return(NULL)
          # stop("More than 15 conflict edges.")
        }

        cat("Enumerating DAGs...")
        all_graphs <- enumerate_graphs(results$pc@graph) # in Zeile 1 berechnet
        cat(" Done. \n")
        if (save) {
          save(all_graphs, file = filename_graphs)
        }
        if (no_results) {
          return(all_graphs)
        }
      }

      if (is.null(all_graphs)) {
        return(NULL)
      }

      # all_results <- pblapply(all_graphs, graph_to_results, ida_function = ida_function)   ## schneller (?) # library("pbapply")
      all_results <- list()
      for (i in 1:length(all_graphs)) {
        print(paste("DURCHLAUF", i, "VON", length(all_graphs)))            ## mit Angabe des aktuellen Druchlaufs
        all_results[[i]] <- graph_to_results(all_graphs[[i]], ida_function = ida_function)
      }
      if (save) {
        save(all_results, file = filename_results)
      }
    } else {
      load(file = outpath_where_graphs_exist)
      # if (outpath_where_graphs_exist != filename_graphs) oder
      if (!file.exists(filename_graphs)) {
        save(all_graphs, file = filefilename_graphsname)
      }

      load(file = outpath_where_results_exist)
      # if (outpath_where_graphs_exist != filename) oder
      if (!file.exists(filename_results)) {
        save(all_results, file = filename_results)
      }
    }
  }
  if (no_results) {
    return(all_graphs)
  } else {
    return(list(graphs = all_graphs, results = all_results))
  }
}

find_graphs_with_highest_int_pos <- function(all_results, obj_fct = list, dir = c("on", "of")) {

  quality_measure <- function(distribution, percentile) {
    statistics_of_influenced_positions(effects = distribution, percentile = percentile, interesting_positions = int_pos, print = FALSE)
  }

  compute_quality_measure_for_results <- function(results, percentile) {
    of_effects <- results$ida$`372`$of$scaled_effects
    of_max <- quality_measure(apply(of_effects, 1, max), percentile = percentile)
    of_min <- quality_measure(apply(of_effects, 1, min), percentile = percentile)

    on_max <- quality_measure(results$ida$`372`$on$scaled_effects[, 1], percentile = percentile)
    on_min <- quality_measure(results$ida$`372`$on$scaled_effects[, 2], percentile = percentile)

    return(list(of_max = of_max, of_min = of_min, on_max = on_max, on_min = on_min))
  }


  # int_pos_in_percentile_75 <- lapply(all_results, compute_quality_measure_for_results)
  int_pos_in_percentile_75 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.75)))
  int_pos_in_percentile_85 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.85)))
  int_pos_in_percentile_95 <- lapply(all_results, function(x) return(compute_quality_measure_for_results(x, percentile = 0.95)))

  if (dir == "on") {
    int_pos_in_percentile_75 <- lapply(int_pos_in_percentile_75, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
    int_pos_in_percentile_85 <- lapply(int_pos_in_percentile_85, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
    int_pos_in_percentile_95 <- lapply(int_pos_in_percentile_95, function(x) return(list(on_max = x$on_max, on_min = x$on_min)))
  } else if (dir == "of") {
    int_pos_in_percentile_75 <- lapply(int_pos_in_percentile_75, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
    int_pos_in_percentile_85 <- lapply(int_pos_in_percentile_85, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
    int_pos_in_percentile_95 <- lapply(int_pos_in_percentile_95, function(x) return(list(of_max = x$of_max, of_min = x$of_min)))
  }

  total_number_of_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_75 <- sapply(int_pos_in_percentile_75, function(x) length(unique(unlist(x)))) # fast immer 11! sonst 10!!

  total_number_of_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_85 <- sapply(int_pos_in_percentile_85, function(x) length(unique(unlist(x)))) # oft 11, tw. bis zu 6

  # int_pos_in_percentile_95 <- lapply(all_results, compute_quality_measure_for_results)
  total_number_of_int_pos_in_percetile_95 <- sapply(int_pos_in_percentile_95, function(x) length(unlist(x)))
  total_number_of_diff_int_pos_in_percetile_95 <- sapply(int_pos_in_percentile_95, function(x) length(unique(unlist(x)))) # 5-8

  max_pos_75 <- which(total_number_of_int_pos_in_percetile_75 == max(total_number_of_int_pos_in_percetile_75))
  max_pos_85 <- which(total_number_of_int_pos_in_percetile_85 == max(total_number_of_int_pos_in_percetile_85))
  max_pos_95 <- which(total_number_of_int_pos_in_percetile_95 == max(total_number_of_int_pos_in_percetile_95))

  max_diff_pos_75 <- which(total_number_of_diff_int_pos_in_percetile_75 == max(total_number_of_diff_int_pos_in_percetile_75))
  max_diff_pos_85 <- which(total_number_of_diff_int_pos_in_percetile_85 == max(total_number_of_diff_int_pos_in_percetile_85))
  max_diff_pos_95 <- which(total_number_of_diff_int_pos_in_percetile_95 == max(total_number_of_diff_int_pos_in_percetile_95))


  print("Graphs with highest effects for interesting positions:")
  print(paste("max_pos_75", paste(max_pos_75, collapse = ", "), sep = ": "))
  print(paste("max_pos_85", paste(max_pos_85, collapse = ", "), sep = ": "))
  print(paste("max_pos_95", paste(max_pos_95, collapse = ", "), sep = ": "))

  print(paste("max_diff_pos_75", paste(max_diff_pos_75, collapse = ", "), sep = ": "))
  print(paste("max_diff_pos_85", paste(max_diff_pos_85, collapse = ", "), sep = ": "))
  print(paste("max_diff_pos_95", paste(max_diff_pos_95, collapse = ", "), sep = ": "))


  # cat(paste0("max_pos_75: ", paste(max_pos_75, collapse = ", "), "\n"))
  # cat(paste0("max_pos_85: ", paste(max_pos_85, collapse = ", "), "\n"))
  # cat(paste0("max_pos_95: ", paste(max_pos_95, collapse = ", "), "\n"))
  #
  # cat(paste0("max_diff_pos_75: ", paste(max_pos_75, collapse = ", "), "\n"))
  # cat(paste0("max_diff_pos_85: ", paste(max_pos_85, collapse = ", "), "\n"))
  # cat(paste0("max_diff_pos_95: ", paste(max_pos_95, collapse = ", "), "\n"))

  return(obj_fct(max_pos_75 = max_pos_75, max_pos_85 = max_pos_85, max_pos_95 = max_pos_95,
                 max_diff_pos_75 = max_diff_pos_75, max_diff_pos_85 = max_diff_pos_85, max_diff_pos_95 = max_diff_pos_95))
}


element_in_most_of_the_6_sets <- function(max_pos_75, max_pos_85, max_pos_95,
                                          max_diff_pos_75, max_diff_pos_85, max_diff_pos_95) {
  all <- c(max_pos_75, max_pos_85, max_pos_95,
           max_diff_pos_75, max_diff_pos_85, max_diff_pos_95)
  return(as.integer(names(which(table(all) == max(table(all))))))
}


# select for a results object the mean results on and of position 372, respectively
# and return both as a list
# mean of min and max?
# TODO: remove weight_effects_on_by, read from dir instead
mean_effects_min_max <- function(results, position, weight_effects_on_by, scaled_effects = FALSE, dir = c("on", "of")) {
  ret_list <- list()
  if ("of" %in% dir) {
    if (scaled_effects) {
      of_effects <- results$ida[[position]]$of$scaled_effects
    } else {
      of_effects <- results$ida[[position]]$of$effects
    }
    if (is.null(of_effects)) {
      warning("Effects OF not found (in mean_effects_min_max.")
    }
    of_max <- apply(of_effects, 1, max)
    of_min <- apply(of_effects, 1, min)
    ret_list$of <- apply(cbind(of_max, of_min), 1, mean)
  }

  if ("on" %in% dir) { # bzw. etwas, das "on " enthält -> grepl
    on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
    if (scaled_effects) {
      on_max <- results$ida[[position]][[on]]$scaled_effects[, 1]
      on_min <- results$ida[[position]][[on]]$scaled_effects[, 2]
    } else {
      on_max <- results$ida[[position]][[on]]$effects[, 1]
      on_min <- results$ida[[position]][[on]]$effects[, 2]
    }
    if (is.null(on_min) || is.null(on_max)) {
      warning("Effects ON not found (in mean_effects_min_max.")
    }
    ret_list[[on]] <- apply(cbind(on_max, on_min), 1, mean)
  }

  return(ret_list)
}

# sum all effects:
# should rather be devided by 100, thus mean
# computes the mean of the effects computed for each graph
compute_over_all_graphs <- function(all_results, position, weight_effects_on_by, use_scaled_effects_for_sum = FALSE, scale_in_the_end = FALSE,
                                    function_over_all_graphs = "mean", direction_type = c("on", "of"),
                                    scale_effects_on = "372" %in% rownames(all_results[[1]]$ida$`372`$of$effects)) {
  # debug(mean_effects_min_max)
  if (length(direction_type) == 1 && !grepl("on|of", direction_type)) {
    direction_type = c("on", "of")
  }
  effects_over_all_graphs_on_of <- list()
  if (!missing(position) &&  !is.na(as.numeric(position))) {
    # mean_effects_min_max_FUN <- function_set_parameters(mean_effects_min_max, parameters = list(position = position))
    min_max_mean_effects_on_of <- lapply(all_results, mean_effects_min_max, weight_effects_on_by = weight_effects_on_by, scaled_effects = use_scaled_effects_for_sum, dir = direction_type, position = position)
  } else {
    warning("No position given in compute_over_all_graphs!")
  }
  if ("of" %in% direction_type) {
    min_max_mean_effects_of <- do.call(cbind, (lapply(min_max_mean_effects_on_of, function(list) return(list$of))))
    effect_over_all_graphs_of <- apply(min_max_mean_effects_of, 1, function_over_all_graphs)
    effects_over_all_graphs_on_of <- c(effects_over_all_graphs_on_of, list(overAllGraphs_of = effect_over_all_graphs_of))
  }
  # if (grepl("mean", weight_effects_on_by)) { # (weight_effects_on_by == "mean_abs_effect") {
  #   on <- "on-rel-to-mean"
  # } else if (grepl("median", weight_effects_on_by)) { # if (weight_effects_on_by == "median_abs_effect") {
  #   on <- "on-rel-to-median"
  # } else if (weight_effects_on_by == "var" || weight_effects_on_by == "vars") {
  #   on <- "on-rel-to-var"
  # }
  if ("on" %in% direction_type) { # bzw. etwas, das "on " enthält -> grepl
    on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
    min_max_mean_effects_on <- do.call(cbind, (lapply(min_max_mean_effects_on_of, function(list) return(list[[on]]))))
    effect_over_all_graphs_on <- apply(min_max_mean_effects_on, 1, function_over_all_graphs)
    if (scale_effects_on) {
      effect_over_all_graphs_on <- effect_over_all_graphs_on / effect_over_all_graphs_on[as.character(position)]
    }
    effects_over_all_graphs_on_of <- c(effects_over_all_graphs_on_of, list(overAllGraphs_on = effect_over_all_graphs_on))
  }
  # effects_over_all_graphs_on_of <- list(overAllGraphs_of = effect_over_all_graphs_of, overAllGraphs_on = effect_over_all_graphs_on)

  return(effects_over_all_graphs_on_of) #, all_mean = all_mean_effects))
}

display_effects <- function(effects, effect_hue_by = "effect", direction = "mean", int_pos, perturbed_position = int_pos[1],
                            scale_in_the_end = FALSE, weight_effects_on_by = "median", function_over_all_graphs = "mean",
                            ida_percentile, caption = "", main_caption = "",
                            print = TRUE, plot = TRUE, for_combined_plot = FALSE, plot_effect_quality = TRUE,
                            plot_false_pos_neg = TRUE, plot_effect_score = TRUE, barplot_contour_black = TRUE,
                            plot_to_canvas = TRUE, outpath, output_formats = c()) {

  if (missing(ida_percentile)) {
    if (!is.null(dim(effects))) {
      ida_percentile <- 1 - (length(int_pos) / dim(effects)[1])
    } else {
      ida_percentile <- 1 - (length(int_pos) / length(effects))
    }
  }

  if (plot) {
    if (!for_combined_plot && !is.null(main_caption) && !missing(main_caption)) {
      oma <- c( 0, 0, length(main_caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
      # par(mfrow=c(lines, columns), oma = oma)
      par(oma = oma)
    }

    # if ("on" %in% direction && "of" %in% direction) {
    #   par(mfrow = c(2,1))
    # }
    if (length(direction) > 1 && !for_combined_plot) {
      par(mfrow = c(length(direction), 1))
    }
  }

  # aus ida_geklaut; ergibt heir aber keinen Sinn (?)
  # if (direction == "both") {
  #   direction <- c("of", "on")
  # }

  for (dir in direction) {
    if (! is.list(effects)) {
      effects_dir <- effects
    } else {
      if (! ((dir == "on") ||  (dir == "of"))) {
        effects_matrix <- rbind(effects$overAllGraphs_of, effects$overAllGraphs_on)
        on_of_mean_effects <- apply(effects_matrix, 2, get(direction))
        effects[[paste0("overAllGraphs_", direction, "_on_of")]] <- on_of_mean_effects
      }

      if (! ((dir == "on") ||  (dir == "of"))) {
        effects_dir <- effects[[paste0("overAllGraphs_", direction, "_on_of")]]
      } else {
        effects_dir <- effects[[paste0("overAllGraphs_", tolower(dir))]]
      }
    }

    if (print || plot_effect_quality) {
      false_pos_neg <- statistics_of_influenced_positions(effects_dir, percentile = ida_percentile,
                                                          interesting_positions = int_pos, print = FALSE, return_list = TRUE)
      effect_quality <- quality_of_effects_distibution(effects = effects_dir, int_pos = int_pos)
      score <- score_for_effects(effects = effects_dir, int_pos = int_pos, perturbed_position = perturbed_position,
                                 effect_quality, length(false_pos_neg$fn))
    }

    if (print) {
      if (! ((dir == "on") ||  (dir == "of"))) {
        print(paste0("SUM EFFECTS ", toupper(dir), "(ON, OF):"))
      } else {
        print(paste0("SUM EFFECTS ", toupper(dir), ":"))
      }
      statistics_of_influenced_positions(effects_dir, percentile = ida_percentile,
                                         interesting_positions = int_pos, print = TRUE)

      writeLines(paste("Quality of effects:", effect_quality))
                       # quality_of_effects_distibution(effects = effects_over_all_graphs_on_of$overAllGraphs_mean_on_of,
                                                      # int_pos = int_pos)))
      writeLines(paste("Score:", score))
    }


    if (plot) {
      if (dir == "of") {
        dir_print <- dir
      } else if (dir == "on") {
        dir_print <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
      } else {
        dir_print <- paste0(" (", direction, " of: ", pastes("on", weight_effects_on_by, sep = "-rel-to-"), " and of) ")
      }
      if (missing(caption) || is.null(caption)) {
        caption <- paste0(function_over_all_graphs, " over all graphs of effects", dir_print, "position 372")

        # caption <- paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372")
      }
      # plot_effects(effects_on_of$overAllGraphs_of, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption, plot_effect_quality = plot_effect_quality, false_pos_neg = stat_sum_dir)
      if (!plot_effect_score) {
        score = NULL
      }
      if (!effect_quality) {
        effect_quality = NULL
      }
      if (!plot_false_pos_neg) {
        false_pos_neg = NULL
      } else {
        false_pos_neg <- statistics_of_influenced_positions(effects_dir, percentile = ida_percentile,
                                                            interesting_positions = int_pos, print = FALSE, return_list = FALSE) #neu belegen, diesmal mit String
      }
      plot_effects(effects_dir, effect_hue_by = effect_hue_by, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption,
                   effect_quality = effect_quality, false_pos_neg = false_pos_neg, score = score, barplot_contour_black = barplot_contour_black,
                   plot_to_canvas = plot_to_canvas, outpath = outpath, output_formats = output_formats)
    }

      #############
  # if (print) {
  #   if ("of" %in% direction) {
  #     print("SUM EFFECTS OF:")
  #     stat_sum_of<- statistics_of_influenced_positions(effects_on_of$sum_of, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
  #   }
  #   if ("on" %in% direction) {
  #     print("SUM EFFECTS ON:")
  #     stat_sum_on <- statistics_of_influenced_positions(effects_on_of$sum_on, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
  #   }
  #   if (! (("on" %in% direction) ||  ("of" %in% direction))) {
  #     print("SUM EFFECTS MEAN(ON, OF):")
  #     # effects_matrix <- do.call(rbind, effects_over_all_graphs_on_of)
  #     # on_of_mean_effects <- apply(effects_matrix, 2, get(direction))
  #     stat_sum_on_of <- statistics_of_influenced_positions(effects_on_of[[paste0("overAllGraphs_", direction, "_on_of")]], percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
  #   }
  # }
  # if (plot) {
  #   if (!for_combined_plot && !is.null(main_caption) && !missing(main_caption)) {
  #     oma <- c( 0, 0, length(main_caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
  #     # par(mfrow=c(lines, columns), oma = oma)
  #     par(oma = oma)
  #   }
  #
  #   # if ("on" %in% direction && "of" %in% direction) {
  #   #   par(mfrow = c(2,1))
  #   # }
  #
  #   if ("of" %in% direction) {
  #     if (missing(caption) || is.null(caption)) {
  #       caption <- paste0(function_over_all_graphs, " over all graphs of effects of position 372")
  #     }
  #     plot_effects(effects_on_of$overAllGraphs_of, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption, plot_effect_quality = plot_effect_quality, false_pos_neg = stat_sum_of)
  #   }
  #   if ("on" %in% direction) {
  #     on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  #     if (missing(caption) || is.null(caption)) {
  #       caption <- paste0(function_over_all_graphs, " over all graphs of effects ", on, " position 372")
  #     }
  #     # effects_to_plot <- effects_on_of$overAllGraphs_on
  #     plot_effects(effects_on_of$overAllGraphs_on, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption, plot_effect_quality = plot_effect_quality, false_pos_neg = stat_sum_on)
  #   }
  #
  #   if (! (("on" %in% direction) ||  ("of" %in% direction))) {
  #
  #     on <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  #     if (missing(caption) || is.null(caption)) {
  #       caption <- paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372")
  #     }
  #     # effects_to_plot <- on_of_mean_effects
  #     plot_effects(on_of_mean_effects, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption, plot_effect_quality = plot_effect_quality, false_pos_neg = stat_sum_on_of)
  #
  #     # scaled_effects_for_coloring_on_of_mean <- scale_effects(as.matrix(on_of_mean_effects), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
  #     # colors_by_effect <- color_by_effect(scaled_effects_for_coloring_on_of_mean, int_pos, mode = "#FFFFFF")
  #     # if (!scale_in_the_end) {
  #     #   barplot(on_of_mean_effects,
  #     #           main = paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372"),
  #     #           col = colors_by_effect, las = 2)
  #     # } else {
  #     #   barplot(as.vector(scaled_effects_for_coloring_on_of_mean),
  #     #           main = paste0(function_over_all_graphs, " over all graphs of effects (", direction, " of: ", on, " and of) position 372"),
  #     #           col = colors_by_effect, las = 2,
  #     #           names.arg = rownames(scaled_effects_for_coloring_on_of_mean))
  #     # }
  #
  #
  #   }

    # plot_effects(effects_to_plot, int_pos = int_pos, scale_in_the_end = scale_in_the_end, caption = caption, plot_effect_quality = plot_effect_quality)
    if (!for_combined_plot && !missing(main_caption) && !is.null(main_caption)) {
        title(main_caption, outer = TRUE)
    }
  }

  return(effects)
}


plot_effects <- function(effects, effect_hue_by = effects, int_pos, scale_in_the_end, caption, effect_quality,
                         false_pos_neg, score, effect_to_color_mode = "#FFFFFF", barplot_contour_black,
                         plot_to_canvas = TRUE, outpath, output_formats = c()) {

  lines_needed_for_subcaption = 0
  if (!is.null(false_pos_neg)) {
    lines_needed_for_subcaption = lines_needed_for_subcaption + 2
  }
  if (!is.null(effect_quality)) {
    lines_needed_for_subcaption = lines_needed_for_subcaption + 1
  }
  if (!is.null(score)) {
    lines_needed_for_subcaption = lines_needed_for_subcaption + 1
  }

  # TODO Marcel: Platz für die subcaption schaffen (mit mar oder so)
  mgp <- par()$mgp  # get current value
  mgp[1] <- 2 + lines_needed_for_subcaption
  mar <- par()$mar # get current value
  mar[1] <- 3 + lines_needed_for_subcaption
  par(mar = mar, mgp = mgp)

  if (is.character(effect_hue_by)) {
    if (effect_hue_by == "effect") {
      effect_hue_by <- effects
    } else {
      effect_hue_by <- get(effect_hue_by)
    }
  }

  scaled_effect_hue_by <- scale_effects(as.matrix(effect_hue_by), rank = FALSE, amplification_factor = 1, neg_effects = "sep")

  # if (effect_hue_by == "effect") {
  #   colors_by_effect <- color_by_effect(current_scaled_effects, int_pos, mode = effect_to_color_mode)
  # } else if (effect_hue_by == "variance" || effect_hue_by == "var") {
  #   vars <- apply(data, 2, var)
  #   colors_by_effect <- color_by_effect(vars, int_pos, mode = effect_to_color_mode)
  # }

  colors_by_effect <- color_by_effect(scaled_effect_hue_by, int_pos, mode = effect_to_color_mode)

  colors <- colors_by_effect[names(effects)]

  if (!barplot_contour_black) {
    contour = colors
  } else {
    contour = par("fg")
  }

  if (plot_to_canvas) {
    if (!scale_in_the_end) {
      barplot(effects,
              main = caption,
              col = colors, border = contour, las = 2)
    } else {
      scaled_effects <- scale_effects(as.matrix(effects), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
      barplot(as.vector(scaled_effects),
              main = caption,
              col = colors, las = 2,
              names.arg = rownames(scaled_effects_for_coloring))
    }

    sub = ""

    if (!missing(false_pos_neg) && !is.null(false_pos_neg)) {
      sub <- paste0(sub, paste0(false_pos_neg), "\n")
    }
    if (!missing(effect_quality) && !is.null(effect_quality)) {
      sub <- paste0(sub, paste("Quality of effects: ", round(effect_quality, digits = 2), "\n"))
    }
    if (!missing(score) && !is.null(score)) {
      sub <- paste0(sub, paste("Score: ", score, "\n"))
    }

    title(sub = sub)
  }

  #TODO
  for (format in output_formats) {
    if (!nchar(outpath) == 0) {
      if (format == "pdf") {
        pdf(paste(outpath, "_effects.pdf", sep = ""))
      } else if ((format == "ps") || (format == "postscript")) {
        postscript(paste(outpath, "_effects.ps",  sep = ""), paper="special", width = 10, height = 9)
      } else if (format == "svg") {
        svg(paste(outpath, "_effects.svg", sep = ""))
      }

      if (!scale_in_the_end) {
        barplot(effects,
                main = caption,
                col = colors, border = contour, las = 2)
      } else {
        scaled_effects <- scale_effects(as.matrix(effects), rank = FALSE, amplification_factor = 1, neg_effects = "sep")
        barplot(as.vector(scaled_effects),
                main = caption,
                col = colors, las = 2,
                names.arg = rownames(scaled_effects_for_coloring))
      }

      sub = ""

      if (!missing(false_pos_neg) && !is.null(false_pos_neg)) {
        sub <- paste0(sub, paste0(false_pos_neg), "\n")
      }
      if (!missing(effect_quality) && !is.null(effect_quality)) {
        sub <- paste0(sub, paste("Quality of effects: ", round(effect_quality, digits = 2), "\n"))
      }
      if (!missing(score) && !is.null(score)) {
        sub <- paste0(sub, paste("Score: ", score, "\n"))
      }

      title(sub = sub)

      dev.off()
    }
  }
}

# TODO Marcel: for_combined_plot einfügen
deviation_from_mean <- function(all_results, dir, weight_effects_on_by, plot_graphs) {
  # for (dir in direction) {
  # how much do the distinct distributions deviate from the mean?
  # OF

  if (dir == "on") {
    dir <- pastes("on", weight_effects_on_by, sep = "-rel-to-")
  }

  all_mean_effects_all <- lapply(all_results, mean_effects_min_max, weight_effects_on_by = weight_effects_on_by, dir = dir)
  all_mean_effects <- do.call(cbind, (lapply(all_mean_effects_all, function(list) return(list[[dir]]))))
  mean_effect <- apply(all_mean_effects, 1, sum) / length(all_results)

  dev_from_mean <- apply(all_mean_effects, 2, function(distribution) {return(distribution - mean_effect)})

  graphics.off()
  par(mfrow = c(5,5))

  for (i in plot_graphs) {
    barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot,
            main = paste("Graph", i, "\n effects", dir))
  }
  # 25 und 19 sehen GLEICH (gut aus)
  # Graphen aber nicht:
  # graphics.off()
  # par(mfrow = c(1,2))
  # plot(all_graphs[[25]])
  # plot(all_graphs[[19]])

  # for (i in 26:50) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }

  # for (i in 51:75) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }

  # for (i in 76:100) {
  #   barplot(dev_from_mean[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  # }
  # 98 sieht gut aus! (min_var = 0.01)

  # apply(dev_from_mean, 2, function(x) barplot(x, col = colors_for_barplot))  # geht nur für alle auf einmal, dafür ist aber die Grafik zu klein

  # TODO: same_procedure for on!
  # btw TODO: on ohne Gewichtung der einzelnen Positionen mit dem durchschnittlichen Effekt
  # }
  # } else if (plot == "on") {
  #   # ON
  #   substract_mean_on <- function(distribution) {
  #     return(distribution - (sum_effect_on/length(all_results)))
  #   }
  #   dev_from_mean_on <- apply(all_mean_effects_on, 2, substract_mean_on)
  #
  #   graphics.off()
  #   par(mfrow = c(5,5))
  #
  #   # for (i in 1:25) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #
  #   # for (i in 26:50) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #
  #
  #   # for (i in 51:75) {
  #   #   barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   # }
  #   # min_var = 0.001: Graph 66
  #
  #   for (i in 76:100) {
  #     barplot(dev_from_mean_on[,i], col = colors_for_barplot, border = colors_for_barplot, main = paste("Graph", i))
  #   }
  #
  #   # apply(dev_from_mean_on, 2, function(x) barplot(x, col = colors_for_barplot))  # geht nur für alle auf einmal, dafür ist aber die Grafik zu klein
  # }
}

graph_to_results <- function(graph, ida_function) {
  results <- list()
  # if (is.null(graph)) {
  #
  # }
  pc <- new("pcAlgo", graph = graph)
  results$pc <- pc
  # results <- causal_effects_ida(data = data, perturbed_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
  #                               protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
  #                               amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
  #                               pymol_bg_color = "grey",
  #                               barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", causal_analysis = TRUE, percentile = ida_percentile)
  results <- ida_function(results = results)
}

