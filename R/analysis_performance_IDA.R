s = 10
# < 100 weil all_results (bei mir) nur so viele Objekte hat

t_before <- Sys.time()
for (i in 1:s) {
  print(paste("DURCHLAUF", i))
  causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
                     protein = protein, results = all_results[[i]], coloring = "all", no_colors = FALSE, outpath = outpath,
                     amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
                     pymol_bg_color = "grey",
                     barplot = TRUE, caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, percentile = 0.75)
}
t_after <- Sys.time()
t <- t_after - t_before
print(t)

# s = 10: 13 sec
# s = 100: 2,13 min