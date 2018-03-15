compare_effects_per_position <- function (results, 
                                          # subset_of_positions = c(1:25), 
                                          subset_of_positions = c(1:39),
                                          direction_of_effects,
                                          heatmap = TRUE) {
  all_pairwise_effects <- results$all_pairwise_effects[,subset_of_positions]
  
  outpath = paste0(results$summary$outpath, "_poswise_effects")
  
  output_formats <- c("", "pdf")
  
  for (format in output_formats) {
    if (!nchar(outpath) == 0) {
      turn_off_dev = TRUE
      if (format == "pdf") {
        pdf(paste(outpath, ".pdf", sep = ""), width = 100, height = 90)
      } else if ((format == "ps") || (format == "postscript")) {
        postscript(paste(outpath, ".ps",  sep = ""), paper="special", width = 10, height = 9)
      } else if (format == "svg") {
        svg(paste(outpath, ".svg", sep = ""), width = 1000, height = 900)
      } else {
        turn_off_dev = FALSE
        if (!heatmap && length(subset_of_positions) > 25) {
          print("Plotting only the first ones...")
        # } else {
          plot.new()
        }
      }
      
      if (heatmap) {
        par(mfrow = c(2,2))
        # heatmap(results$all_pairwise_effects, revC = TRUE, symm = TRUE)
        int_pos <- interesting_positions()
        # heatmap(results$all_pairwise_effects, revC = TRUE, symm = FALSE, col.axis = "green")
        heatmap.2(results_G$all_pairwise_effects, revC = TRUE, symm = FALSE, ColSideColors = int_pos[colnames()])
      } else {
        n_lines_cols <- ceiling(sqrt(length(subset_of_positions)))
        par(mfrow = c(n_lines_cols,n_lines_cols))
        
        clusters <- results$effects_clustering$pv[[1]]
        names(clusters) <- NULL
        
        # apply(all_pairwise_effects, 2, function(effects, ...) {
        #   print(effects)
        #   caption <- paste("Position ", colnames(effects))
        #   display_effects(effects, caption = caption, ...)}, int_pos = results$int_pos)
        # apply(all_pairwise_effects, 2, barplot)
        # debug(plot_effects)
        
        lapply(colnames(all_pairwise_effects), function(position, ...) {
          in_cluster <- which(sapply(clusters, function(list) return(position %in% list)))
          if (length(in_cluster) == 0) {
            in_cluster <- "none"
          }
          caption <- paste("Position ", position, "- Cluster:", in_cluster)
          display_effects(all_pairwise_effects[,position], caption = caption, ...)}, 
          int_pos = results$general$int_pos)
      }
      if (turn_off_dev) {
        dev.off()
      }
    }
  }
}
