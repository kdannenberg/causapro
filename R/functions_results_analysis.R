library(gplots)
library(heatmap3)

compare_effects_per_position <- function (results, 
                                          # subset_of_positions = c(1:25), 
                                          # subset_of_positions = c(1:39),
                                          subset_of_positions,
                                          direction_of_effects,
                                          heatmap = TRUE,
                                          hclust_method,
                                          dist_measure,
                                          nboot) {
  plot.new()
  
  all_pairwise_effects <- results$all_pairwise_effects
  
  if (is.null(all_pairwise_effects)) {
    warning("No effects there to compare. Results object has no slot all_pairwise_effects.")
    return(NULL)
  }
  
  if (missing(subset_of_positions)) {
    subset_of_positions = c(1:dim(results$all_pairwise_effects)[1])
  }
  
  all_pairwise_effects <- all_pairwise_effects[, subset_of_positions]
  
  if (is.null(all_pairwise_effects)) {
    warning(paste0("No effects there to compare. This slot has none of the positions in ", paste(subset_of_positions, collapse = ","), "."))
    return(NULL)
  }
  
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
        
        # one way for colors
        # breaks <- seq(from=min(range(results$all_pairwise_effects)), 
        #               to=max(range(results$all_pairwise_effects)), length.out=100)
        # midpoint <- which.min(abs(breaks - 0))
        # rampCol1 <- colorRampPalette(c("forestgreen", "darkgreen", "black"))(midpoint)
        # rampCol2 <- colorRampPalette(c("black", "darkred", "red"))(100-(midpoint+1))
        # cols <- c(rampCol1,rampCol2)
        
        # another way for good colors
        # quantile.range <- quantile(results$all_pairwise_effects, probs = seq(0, 1, 0.01))
        # palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
        
        # maunally set palette breaks
        # palette.breaks <- c(0, 0.05, 0.15,  0.25,  0.35,  0.45, 0.55, 0.65, 0.75, 1)
        # palette.breaks <- c(0, 0.1, 0.2,  0.3,  0.4,  0.5, 0.6, 0.7, 0.8, 0.9, 1)
        palette.breaks <- c(-1, -0.5, 0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6)
        # palette.breaks <- 16^(seq(0, 0.35, by = 0.01)) - 1
        # palette.breaks <- 2^(seq(0, 1.5, by = 0.1)) - 1
        
        # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
        # cols  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
        cols  <- colorRampPalette(c( "#3C8CC8", "#3C7864", "#69A019"))(length(palette.breaks) - 1)
        
        # heatmap(results$all_pairwise_effects, revC = TRUE, symm = TRUE)
        int_pos <- interesting_positions(results)
        int_pos_cols <- int_pos_to_color_vector(pos = colnames(results$all_pairwise_effects), int_pos = int_pos)
        # heatmap(results$all_pairwise_effects, revC = TRUE, symm = FALSE, col.axis = "green")
        # heatmap.2(results$all_pairwise_effects, revC = TRUE, symm = FALSE, trace = "none",
        #           col = cols, breaks = palette.breaks,
        #           # dendrogram = "col", Rowv = FALSE,
        #           ColSideColors = int_pos_cols, RowSideColors = int_pos_cols)
        # 
        # heatmap.2(results$all_pairwise_effects, revC = FALSE, symm = FALSE, trace = "none",
        #           col = cols, breaks = palette.breaks,
        #           Rowv = FALSE, Colv = FALSE,
        #           ColSideColors = int_pos_cols, RowSideColors = int_pos_cols)
        
        
        
        heatmap3(results$all_pairwise_effects, revC = TRUE, symm = FALSE,
                  col = heat.colors(length(palette.breaks) - 1), #breaks = palette.breaks,
                  Rowv = NA, Colv = NA,
                  ColSideColors = int_pos_cols, ColSideLabs = "",
                  RowSideColors = int_pos_cols, RowSideLabs = "",
                  ColAxisColors = 1, RowAxisColors = 1,
                  na.rm = FALSE, balanceColor = FALSE)
      
        heatmap3(results$all_pairwise_effects, revC = TRUE, symm = FALSE, trace = "none",
                 col = heat.colors(length(palette.breaks) - 1), #breaks = palette.breaks,
                 # dendrogram = "col", Rowv = FALSE,
                 method = "ward.D",
                 reorderfun = function(d, w) rev(reorder(d, w)),
                 ColSideColors = int_pos_cols, ColSideLabs = "",
                 RowSideColors = int_pos_cols, RowSideLabs = "",
                 ColAxisColors = 1, RowAxisColors = 1,
                 na.rm = FALSE, balanceColor = FALSE)
        
        #TODO: compute_if_not_existent benutzen
        wrap_FUN_pvclust <- function(method, ...) {
          return(FUN_pv(...)$hclust)
        }
        FUN_pv <- function_set_parameters(pvclust, parameters = list(data = all_pairwise_effects, 
                                                                     method.hclust = hclust_method,
                                                                     method.dist = dist_measure, nboot = nboot))
        #TODO: compute_if_not_existent einbauen
        
        heatmap3(results$all_pairwise_effects, revC = FALSE, symm = FALSE,
                 col = heat.colors(length(palette.breaks) - 1), #breaks = palette.breaks,
                 hclustfun = wrap_FUN_pvclust,
                 # dendrogram = "col", Rowv = FALSE,
                 reorderfun = function(d, w) rev(reorder(d, w)),
                 ColSideColors = int_pos_cols, ColSideLabs = "",
                 RowSideColors = int_pos_cols, RowSideLabs = "",
                 ColAxisColors = 1, RowAxisColors = 1,
                 na.rm = FALSE, balanceColor = FALSE)
        
        heatmap3(results$all_pairwise_effects, revC = TRUE, symm = FALSE,
                 col = heat.colors(length(palette.breaks) - 1), #breaks = palette.breaks,
                 method = "ward.D2",
                 # dendrogram = "col", Rowv = FALSE,
                 reorderfun = function(d, w) rev(reorder(d, w)),
                 ColSideColors = int_pos_cols, ColSideLabs = "",
                 RowSideColors = int_pos_cols, RowSideLabs = "",
                 ColAxisColors = 1, RowAxisColors = 1,
                 na.rm = FALSE, balanceColor = FALSE)
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
