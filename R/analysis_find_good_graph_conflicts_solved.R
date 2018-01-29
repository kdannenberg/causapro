# DEPRECATED, use analysis_find_good_graph, with type_of_graph_set = "conflict"

source("~/.configuration_code.R")

source("functions_causal_effects.R")
source("functions_general.R")
source("functions_conversions.R")

source("analysis_find_good_graph.R")

source("~/.configuration_code.R")
source("compute_DAG_G.R")

# source("configuration_data.R")

protein = "PDZ"
int_pos <- interesting_positions(protein = protein, coloring = "")
perturbed_position <- "372"

alpha = 0.01
min_pos_var = 0

new = FALSE
save = TRUE

pc_solve_conflicts = FALSE
pc_u2pd = "retry"

use_scaled_effects_for_sum = FALSE   # otherwise scaling is done in the end, for the sum

results <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
                               pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                               graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
                               data_in_results = TRUE, output_parameters_in_results = TRUE)
data <- results$data
caption <- results$caption
outpath <- results$outpath

top_11_percentile <- 1 - (11 / dim(data)[2])
ida_percentile <- top_11_percentile

# weight_effects_on_by = ""          # in der Summe ganz schlecht
# weight_effects_on_by = "var"
# weight_effects_on_by = "mean"
weight_effects_on_by = "median"  # sieht (in der Summe) am besten aus

# which part of the analysis should be plotted?
# IDA für alle s Graphen berechen und summieren (besser wäre vllt: mitteln, also nochmal durch s teilen. Kannst du gerne machen, Marcel.)
# plot = "sum over all graphs"
# IDA für alle s Graphen brechenen und denjenigen bestimmen, der am besten mit den gewünschten Ergebnissen übereinstimmt
# plot = "best graph"
# für alle Graphen mit Nummern in <plot> die Abweichung von der Summe (dem zukünftigen Mittelwert) über alle Graphen bestimmen, 
# also gewissermaßen wie repräsentativ der Graph jeweils für die Menge ist
plot = 1:25     ## if (is.numeric(plot) && length(plot) > 1) --> deviation from mean for graph(s) nr. <plot>



graphs <- enumerate_graphs(results_G$pc@graph) # in Zeile 1 berechnet
# all_results <- lapply(graphs, graph_to_results)

find_good_graph_among <- function(graphs, plot){
  all_results <- lapply(graphs, graph_to_results)
  
  best_graphs <- find_graphs_with_highest_int_pos(all_results = all_results, obj_fct = element_in_most_of_the_6_sets)
  
  for (i in best_graphs) {
    # i= 28
    if (plot == "best graph") {
      cat(paste0("BEST GRAPH #", which(best_graphs == i), ": ", i))
      causal_effects_ida(data = data, perturbed_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                         protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
                         amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                         pymol_bg_color = "grey", barplot = TRUE,
                         caption = c(paste("Graph", i), caption), show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = top_11_percetile)
    }
  }
  
  return(best_graphs)
}

sum_effects <- function(all_results, caption, weight_effects_on_by, percentile = top_11_percentile, int_pos, plot) {
  if (plot == "sum over all graphs") {
    oma <- c( 0, 0, length(caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
    # par(mfrow=c(lines, columns), oma = oma) 
    par(oma = oma)
    
    cat("\n")
    sum_effect <- sum_all_effects(all_results, weight_effects_on_by = weight_effects_on_by)
    print("SUM EFFECTS OF:")
    stat_sum_on <- statistics_of_influenced_positions(sum_effect$sum_of, percentile = percentile, interesting_positions = int_pos, print = TRUE)
    print("SUM EFFECTS ON:")
    stat_sum_of <- statistics_of_influenced_positions(sum_effect$sum_on, percentile = percentile, interesting_positions = int_pos, print = TRUE)
    
    title(caption, outer = TRUE)
  }
}

deviat_from_mean <- function(all_results, weight_effects_on_by, plot) {
  if (is.numeric(plot) && length(plot) > 1) {
    deviation_from_mean(all_results = all_results, dir = "on", weight_effects_on_by = weight_effects_on_by, plot_graphs = plot)
  }
}

graph_to_results <- function(graph) {
  results <- list()
  pc <- new("pcAlgo", graph = graph)
  results$pc <- pc
  results <- causal_effects_ida(data = data, perturbed_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                     protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
                     amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                     pymol_bg_color = "grey",
                     barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
}

# best_graphs <- find_good_graph_among(graphs)
# summed_effects <- sum_effects(all_results, caption = caption, weight_effects_on_by = weight_effects_on_by, percentile = 0.87, int_pos = int_pos)
dev_from_mean <- deviat_from_mean(all_results, weight_effects_on_by = weight_effects_on_by, plot = plot)

# print("hallo")
# idaFast(which(as.character(colnames(data)) == perturbed_position), 
#         +         1:dim(data)[2], cov(data), graphs[[1]])
