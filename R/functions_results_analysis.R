compare_effects_per_position <- function (results, 
                                          # subset_of_positions = c(1:25), 
                                          subset_of_positions = c(26:39),
                                          direction_of_effects) {
  all_pairwise_effects <- results$all_pairwise_effects[,subset_of_positions]
  plot.new()
  n_lines_cols <- ceiling(sqrt(length(subset_of_positions)))
  par(mfrow = c(n_lines_cols,n_lines_cols))
  
  clusters <- results$effects_clustering$pv[[1]]
  names(clusters) <- NULL
  
  # apply(all_pairwise_effects, 2, function(effects, ...) {
  #   print(effects)
  #   caption <- paste("Position ", colnames(effects))
  #   display_effects(effects, caption = caption, ...)}, int_pos = results$int_pos)
  # apply(all_pairwise_effects, 2, barplot)
  
  lapply(colnames(all_pairwise_effects), function(position, ...){
    in_cluster <- which(sapply(clusters, function(list) return(position %in% list)))
    if (length(in_cluster) == 0) {
      in_cluster <- "none"
    }
    caption <- paste("Position ", position, "- Cluster:", in_cluster)
    display_effects(all_pairwise_effects[,position], caption = caption, ...)}, 
    int_pos = results$general$int_pos)
}