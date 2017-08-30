# previously: analysis_find_good_graph.R
source("~/.configuration_code.R")

source("functions_causal_effects.R")
source("functions_general.R")
source("functions_conversions.R")
source("functions_tools.R")

source("functions_analysis_for_a_set_of_graphs.R")

# tests
source("tests_analysis_for_a_set_of_graphs.R")

source("compute_DAG_G.R")
source("~/.configuration_code.R")
source("compute_DAG_S.R")
# source("configuration_data.R")

analyse_set_of_graphs <- function(
  type_of_graph_set = "conflict", # "retry" or "conflict"
  protein = "PDZ",
  int_pos = interesting_positions(protein = protein, coloring = ""),
  measure = "G",
  type_of_data = "DDDG",
  subtype_of_data = "5",
  # use_DDG = FALSE,  # einfach den default im Skript umstellen
  protein_causality_function = get(paste0("protein_causality_", measure)),
  alpha = 0.1,
  min_pos_var = 0.01, # TODO: does not work for 0.01 (alpha = 0.01) (idafast for determination of median (effects on))
  # TODO MARCEL: Warum bekomme ich für alpha = 0.01, min_pos_var = 0.01 einen Graphen mit 6 conflict-Kanten, aber VOR remove_dummies 62 (!) Graphen?!
  new = FALSE,
  save = TRUE,
  check_graph_equality = FALSE, # takes forever
  pc_solve_conflicts = FALSE,
  # used when type_of_graph_set = "conflict"
  pc_maj_rule_conflict = TRUE,
  pc_conservative_conflict = FALSE,
  use_scaled_effects_for_each_graph = FALSE,
  scale_in_the_end = FALSE,
  # results <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
  #                           pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
  #                           graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
  #                           data_in_results = TRUE, output_parameters_in_results = TRUE)
  ida_percentile = 1 - (11 / dim(data)[2]), # top 11
  # weight_effects_on_by = "",          # in der Summe ganz schlecht
  # weight_effects_on_by = "var",
  # weight_effects_on_by = "mean",
  weight_effects_on_by = "median",  # sieht (in der Summe) am besten aus
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
  #                       evaluation = evaluation, analysis = analysis)
  # }
  
  pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
    protein_causality_function(min_pos_var = min_pos_var, alpha = alpha,
                        pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                        evaluation = evaluation, analysis = analysis, mute_all_plots = for_combined_plot)
  },
  ida_function = function(results) {
    causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                       protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
                       amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                       pymol_bg_color = "grey", caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, 
                       percentile = ida_percentile, mute_all_plots = for_combined_plot)
  },
  s = 10,    # sample size
  for_combined_plot = FALSE,
  caption_as_subcaption = for_combined_plot
) {

  if (type_of_graph_set == "conflict") {
    pc_u2pd = "relaxed"
  } else if (type_of_graph_set == "retry") {
    pc_u2pd = "retry"
  }
  
  results <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                        min_pos_var = min_pos_var, alpha = alpha,
                                        pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
                                        graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
                                        data_in_results = TRUE, output_parameters_in_results = TRUE, mute_all_plots = TRUE)
  data <- results$data
  caption <- results$caption
  outpath <- results$outpath
  
  if (direction == "both") {
    direction <- c("on", "of")
  }
  
  if (!(("on" %in% direction) ||  ("of" %in% direction)) && !scale_effects_on_so_372_equals_1) {
    warning(paste("Effects on are not scaled. Effects on and of might not be comparably high. Computing the", direction, "might not be meaningful."))
    Sys.sleep(2)
  }
  
  
  
  
  set_of_graphs <- determine_set_of_graphs(type_of_graph_set = type_of_graph_set, pc_function = pc_function, ida_function = ida_function, 
                                           s = s, new = new, save = save, outpath = outpath,
                                           pc_maj_rule_conflict = pc_maj_rule_conflict, pc_conservative_conflict = pc_conservative_conflict)
  all_graphs <- set_of_graphs$graphs
  all_results <- set_of_graphs$results
  
  
  # determine default colors for barplot
  all_one_effects <- as.matrix(all_results[[1]]$ida$`372`$of$effects[,1])
  all_one_effects[,] <- 1
  colors_for_barplot <- color_by_effect(all_one_effects, int_pos, mode = "#FFFFFF")
  
  # check wether any graphs in the set are equal
  if (check_graph_equality) {
    equal <- compare_all_graphs(all_graphs)
    main_diagonal <- matrix(as.logical(diag(nrow = dim(equal)[1])), ncol = dim(equal)[2])
    print(which(equal & !main_diagonal, arr.ind = TRUE))
  }
  
  # compare all the graphs w.r.t. how many conflict, bidirected and directed edges the have
  conflicts_sorted <- sapply(all_graphs, conflict_edges)
  conflicts_sorted <- matrix(unlist(conflicts_sorted), nrow = 3)
  rownames(conflicts_sorted) = c("conflict", "directed", "bidirected")
  colnames(conflicts_sorted) <- as.character(1:dim(conflicts_sorted)[2]) # Graph number
  conflicts_sorted <- conflicts_sorted[, order(conflicts_sorted["bidirected",])]
  
  # print("Different types of edges in the different graphs (columns)")
  # print(conflicts_sorted)
  
  cat("\n")
  cat(paste0("alpha = ", alpha, ", min_pos_var = ", min_pos_var, "\n"))
  
  if (plot == "best_graph") {
    # cat("\n")
    best_graphs <- find_graphs_with_highest_int_pos(all_results = all_results, obj_fct = element_in_most_of_the_6_sets, dir = direction)
    cat("\n")
    
    for (i in best_graphs) {
    # i= 28
      cat(paste0("BEST GRAPH #", which(best_graphs == i), ": ", i))
      causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                        protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
                        amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                        pymol_bg_color = "grey", barplot = TRUE,
                        caption = c(paste("Graph", i), caption), show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
    
    }
  } else if (plot == "over_all_graphs") {
    # cat("\n")
    effects_over_all_graphs_on_of <- compute_over_all_graphs(all_results = all_results, weight_effects_on_by = weight_effects_on_by, 
                                  use_scaled_effects_for_sum = use_scaled_effects_for_each_graph,   
                                  function_over_all_graphs = function_over_all_graphs, direction = direction, 
                                  scale_effects_on = scale_effects_on_so_372_equals_1)
    
    if (caption_as_subcaption) {
      caption = caption
      main_caption = NULL  # soll missing sein
    } else {
      main_caption = caption
      caption = NULL
    }
    
    effects_over_all_graphs_on_of <- display_effects(effects_on_of = effects_over_all_graphs_on_of, direction = direction, 
                                                     int_pos = int_pos, scale_in_the_end = scale_in_the_end, 
                                                     weight_effects_on_by = weight_effects_on_by, 
                                                     function_over_all_graphs = function_over_all_graphs,
                                                     ida_percentile = ida_percentile, caption = caption, main_caption = NULL,
                                                     print = TRUE, plot = TRUE, for_combined_plot = for_combined_plot)
    
    print(paste("Quality of the effects:", quality_of_effects_distibution(effects = effects_over_all_graphs_on_of$overAllGraphs_mean_on_of, int_pos = int_pos)))
    
    # return(effects_over_all_graphs_on_of)
    
  } else if (is.numeric(plot) && length(plot) > 1) {
    deviation_from_mean(all_results = all_results, dir = "of", weight_effects_on_by = weight_effects_on_by, plot_graphs = plot)
  }
}

# TODO: sicherstellen, dass oma nur gesetzt wird, wenn auch eine main_caption (title) geprintet wird

# effects <- analyse_set_of_graphs(direction = "mean", measure = "G", alpha = 0.001, min_pos_var = 0.0001, new = FALSE)

measure = "G"
type_of_data = "DDDG"
subtype_of_data = "5"

protein_causality_function = get(paste0("protein_causality_", measure))
# TODO: DAS FUNKTIONIERT NICHT!! DIE PARAMETER, DIE NCIHT ÜBERGEBEN WERDEN, WERDEN NCIHT JETZT SCHON BELEGT!
# ALLES WIRD ERST BEI AUFRUF AUSGEWERTET
pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
  protein_causality_function(type_of_data = eval(type_of_data), subtype_of_data = eval(subtype_of_data), min_pos_var = min_pos_var, alpha = alpha,
                             pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                             evaluation = evaluation, analysis = analysis, mute_all_plots = TRUE)
}

graphics.off()
oma <- c( 2, 0, 2, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
# par(mfrow=c(lines, columns), oma = oma)
par(mfrow = c(5,3), oma = oma)
# par(mfrow = c(2,3))
for (alpha in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
# for (alpha in c(0.001, 0.005)) {
  for (min_pos_var in c(0.0001, 0.001, 0.01)) {
    # for (measure in c("G", "S")) {
      if (
          (
            (type_of_data == "DDDG" && subtype_of_data == "5") &&
            ((alpha == 0.05 && min_pos_var == 0.001) #||  # More than 15 conflict edges. (22)
            || (alpha == 0.1 && min_pos_var == 0.001)  # More than 15 conflict edges.
            || (alpha == 0.05 && min_pos_var == 0.0001) # More than 15 conflict edges.
            || (alpha == 0.1 && min_pos_var == 0.0001) # More than 15 conflict edges.
            )
           # ) || (
           #   (type_of_data == "DDDG" && subtype_of_data == "10") &&
           #   (
           #   (alpha == 0.005 && min_pos_var == 0.001) ##  "Fehler in wgt.unique[x, ] : Indizierung außerhalb der Grenzen"
           #   # || (alpha == 0.05 && min_pos_var == 0.0001)
           #   # || (alpha == 0.1 && min_pos_var == 0.0001)
           #   )
           ) || (
             (type_of_data == "DDG" && subtype_of_data == "10") &&
             (
               (alpha == 0.1 && min_pos_var == 0.0001) ##  15 conflict edges
               || (alpha == 0.1 && min_pos_var == 0.001) ##  13 conflict edges
             )
           ) || (
             (measure == "S") &&
             (alpha == 0.1 && min_pos_var == 0.0001)
           ) || (
             (type_of_data == "DDG" && subtype_of_data == "5") &&
             (
               (alpha == 0.05 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
            || (alpha == 0.05 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
            || (alpha == 0.05 && min_pos_var == 0.01) #  More than 15 conflict edges. (18)
            || (alpha == 0.1 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
            || (alpha == 0.1 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
            || (alpha == 0.1 && min_pos_var == 0.01) #  More than 15 conflict edges. (21)
             )
           )
          ) {
        plot.new()
      } else {
        analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data, direction = "mean", measure = measure, pc_function = pc_function, alpha = alpha, min_pos_var = min_pos_var,
                              for_combined_plot = TRUE, scale_in_the_end = FALSE, new = TRUE, save = TRUE)

      }
      # }
  }
}

# title(main = paste0(type_of_data, "-", subtype_of_data),
#       sub = "Mean causal effects over all conflict graphs and over the effects on and of position 372", outer = TRUE)
title(main = paste0("Mean causal effects over all conflict graphs and over the effects on and of position 372",
                    "\n", type_of_data, "-", subtype_of_data), outer = TRUE)
