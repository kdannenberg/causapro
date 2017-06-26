# source('/Volumes/Causality/Viren/R/Code/compute_DAG_G.R')
source("general_functions.R")

effects_of_76 <- idaFast(which(colnames(data) == "372"), 1:92, cov(data), results$pc@graph)
rownames(effects_of_76) <- colnames(data)
effects_of_76[order(effects_of_76), , drop = FALSE]     # sort and keep rownames

graphics.off()
barplot(t(effects_of_76))  # TODO: darin gelb und grün färben oder über Fig. 2B im SCIENCE-Paper legen

# sort(effects_of_76)
# by_76_most_influenced_positions <- colnames(data)[as.integer(rownames(effects_of_76)[which(effects_of_76 > 0.6)])]
threshold = quantile(effects_of_76, probs = 0.75)
by_76_most_influenced_positions <- colnames(data)[as.integer(rownames(effects_of_76)[which(effects_of_76 > threshold)])]
print(by_76_most_influenced_positions)
print(paste(length(by_76_most_influenced_positions), "positions over the threshold", threshold))
int_pos <- interesting_positions("PDZ", "crystal")
int_pos_strongly_influenced_by_76 <- intersect(int_pos, by_76_most_influenced_positions)
print(paste("Thereof",  length(int_pos_strongly_influenced_by_76), "out of the", length(int_pos), "interesting positions:", paste(sort(int_pos_strongly_influenced_by_76), collapse = ", ")))
print(paste("Missing: ", paste(setdiff(int_pos, by_76_most_influenced_positions), collapse = ", ")))
## für alpha = 1e-20 : [1] "372" "322" "325" "329" "330" "347" "353" "362" "376" "386"
## für alpha = 0.05  : [1] "372" "322" "325" "329" "330" "347" "353" "362" "376" "386"


# wo stehen die int_pos in der Liste der effected_positions?
ranks <- cbind(apply(effects_of_76, 2, rank))
ranks[sapply(int_pos, function(x) which(rownames(ranks) == x))]



ida <- function(data, perturbated_position,  protein, coloring = "all", outpath, amplification_exponent) {
  effects <- idaFast(which(colnames(data) == perturbated_position), 1:92, cov(data), results$pc@graph)
  rownames(effects) <- colnames(data)
  int_pos <- interesting_positions(protein = protein, coloring = coloring)
  colors_by_effect <- hue_by_effect(effects, int_pos)
  plot_total_effects_in_pymol(colors_by_effect, perturbated_position = perturbated_position, protein = protein, outpath = outpath, amplification_exponent = amplification_exponent)
}