pos <- "353"

effect_on_position <- function(pos, direct = TRUE, with_regard_to_372 = TRUE, colors) {
  if (!direct) {
    pos <- which(as.character(colnames(data)) == pos)
  }
  eff <- idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)
  rownames(eff) <- colnames(data)
  
  print(colors)
  
  if (with_regard_to_372) {
    eff_on_372 <- eff[which(as.character(colnames(data)) == "372")]
    print(eff_on_372)
    
    barplot(as.vector(eff), main = paste(colnames(data)[pos], "\n eff_on_372 = ", round(eff_on_372, digits = 2), sep = ""),
            col.main = colors[colnames(data)[pos]], col = colors, border = colors)
  } else {
    barplot(t(eff))
  }
  
}

# interesting:

# effects_of_position <- function(pos) {
#   # barplot(as.vector(pcalg::idaFast(3, 1:92, cov(data), results$pc@graph)))
#   eff <- idaFast(which(as.character(colnames(data)) == pos), 1:dim(data)[2], cov(data), results$pc@graph)
#   
#   rownames(eff) <- colnames(data)
#   
#   barplot(t(eff))
#   
#   print(eff[which(as.character(colnames(data)) == "372")])
# }

par(mfrow=c(5,9))

int_pos <- interesting_positions(protein = "PDZ", coloring = "all")  
pseudo_effects <- as.matrix(rep(1.0, dim(data)[2]))
rownames(pseudo_effects) <- colnames(data)
colors <- color_by_effect(pseudo_effects, int_pos, mode = "#FFFFFF", color_for_other_positions = "#000000")

for (i in 1:45) {
  effect_on_position(i, colors = colors)
}


# for (i in 46:90) {
#   effect_on_position(i, colors = colors)
# }

