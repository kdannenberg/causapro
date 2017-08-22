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

type_of_graph_set = "conflict" # "retry" or "conflict"

if (type_of_graph_set == "conflict") {
  pc_u2pd = "relaxed"
} else if (type_of_graph_set == "retry") {
  pc_u2pd = "retry"
}

protein = "PDZ"
int_pos <- interesting_positions(protein = protein, coloring = "")

measure = "S"
alpha = 0.01
min_pos_var = 0.001 # TODO: does not work for 0.01 (alpha = 0.01) (idafast for determination of median (effects on))
# TODO MARCEL: Warum bekomme ich für alpha = 0.01, min_pos_var = 0.01 einen Graphen mit 6 conflict-Kanten, aber VOR remove_dummies 62 (!) Graphen?!

new = FALSE  
save = TRUE

check_graph_equality = FALSE # takes forever

pc_solve_conflicts = FALSE

# used when type_of_graph_set = "conflict"
pc_maj_rule_conflict = TRUE
pc_conservative_conflict = FALSE

use_scaled_effects_for_sum = TRUE   # otherwise scaling is done in the end, for the sum


protein_causality_function <- get(paste0("protein_causality_", measure))



# results <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
#                           pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
#                           graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
#                           data_in_results = TRUE, output_parameters_in_results = TRUE)
results <- protein_causality_function(min_pos_var = min_pos_var, alpha = alpha,
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
plot = "sum over all graphs"
# IDA für alle s Graphen brechenen und denjenigen bestimmen, der am besten mit den gewünschten Ergebnissen übereinstimmt
# plot = "best graph"
# für alle Graphen mit Nummern in <plot> die Abweichung von der Summe (dem zukünftigen Mittelwert) über alle Graphen bestimmen, 
# also gewissermaßen wie repräsentativ der Graph jeweils für die Menge ist
# plot = 1:25     ## if (is.numeric(plot) && length(plot) > 1) --> deviation from mean for graph(s) nr. <plot>





# pc_function <- function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
#   protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
#                       pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
#                       evaluation = evaluation, analysis = analysis)
# }

pc_function <- function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
  protein_causality_function(min_pos_var = min_pos_var, alpha = alpha,
                      pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                      evaluation = evaluation, analysis = analysis)
}

ida_function <- function(results) {
  causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                     protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
                     amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                     pymol_bg_color = "grey",
                     barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = ida_percentile)
}

s = 10     # sample size

set_of_graphs <- determine_set_of_graphs(type_of_graph_set = type_of_graph_set, s = s, new = new, save = save, outpath = outpath,
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

print("Different types of edges in the different graphs (columns)")
print(conflicts_sorted)


if (plot == "best graph") {
cat("\n")
best_graphs <- find_graphs_with_highest_int_pos(all_results = all_results, obj_fct = element_in_most_of_the_6_sets)

# i <- 76 # ist in max_pos_85 und _95 und in max_diff_pos_75 und _95
# i <- 98 # für min_var = 0.001

# # plotten:
  # plot_graph(graph = all_results[[28]]$pc@graph, caption = caption, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout,
  #                       coloring = coloring, colors = colors, outpath = outpath, numerical = numerical, plot_as_subgraphs = plot_as_subgraphs,
  #                       plot_only_subgraphs = plot_only_subgraphs, output_formats = graph_output_formats)

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
} else if (plot == "sum over all graphs") {
  oma <- c( 0, 0, length(caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
  # par(mfrow=c(lines, columns), oma = oma) 
  par(oma = oma)
  
  cat("\n")
  sum_effect <- sum_all_effects(all_results, weight_effects_on_by = weight_effects_on_by, use_scaled_effects_for_sum = use_scaled_effects_for_sum)
  print("SUM EFFECTS OF:")
  stat_sum_on <- statistics_of_influenced_positions(sum_effect$sum_of, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
  print("SUM EFFECTS ON:")
  stat_sum_of <- statistics_of_influenced_positions(sum_effect$sum_on, percentile = ida_percentile, interesting_positions = int_pos, print = TRUE)
  
  title(caption, outer = TRUE)
} else if (is.numeric(plot) && length(plot) > 1) {
  deviation_from_mean(all_results = all_results, dir = "of", weight_effects_on_by = weight_effects_on_by, plot_graphs = plot)
}

