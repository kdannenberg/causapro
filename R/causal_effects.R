# source('/Volumes/Causality/Viren/R/Code/compute_DAG_G.R')
# source("general_functions.R")

# adds results$ida
causal_effects_ida <- function(data, perturbated_position, direction = "both", relatve_effects_on_pos = TRUE,
                results = results, protein, coloring = "all", outpath, 
                amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, 
                effect_to_color_mode = "#FFFFFF", pymol_bg_color = "black", barplot = TRUE, caption, no_colors, 
                show_neg_causation = TRUE, neg_effects = "", analysis = TRUE, percentile = 0.75) {
  
  lines <- 1 # lines in plot_space
  
  if (direction == "both") {
    direction <- c("of", "on")
    lines <- 2
  }
  
  graphics.off()
  lines <- 4 # 1 für of, 2 für on (min, max)
  columns <- 2
  par(mfrow=c(lines, columns))

  for (dir in direction) {
    cat("CAUSAL EFFECTS ")
    cat(toupper(dir))
    cat(" POSITION ")
    cat(perturbated_position)
    cat("\n")
    
    if (dir == "of" || dir == "by" || dir == "from") {
      effects <- idaFast(which(as.character(colnames(data)) == perturbated_position), 1:dim(data)[2], cov(data), results$pc@graph)
    } else if (dir == "on") {
      ida_rev <- function(pos) {
        return(pcalg::ida(pos, which(as.character(colnames(data)) == perturbated_position), cov(data), results$pc@graph))
      }
      effects_list <- lapply(1:dim(data)[2], ida_rev)
      # vllt lieber relativ zu dem durchsnittlichen effekt von Position pos 
      effects_min <- sapply(effects_list, function(list) return(min(list)))
      effects_max <- sapply(effects_list, function(list) return(max(list)))
      effects <- cbind(effects_max, effects_min)
      colnames(effects) <- c("max", "min")
      if (relatve_effects_on_pos) {
        means <- sapply(1:dim(data)[2], function(pos) {mean(abs(idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)))})
        # medians <- sapply(1:dim(data)[2], function(pos) {median(idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph))})
        effects <- effects/means
      }
    }
    rownames(effects) <- colnames(data)
    
    # print(effects)
    
    
    int_pos <- interesting_positions(protein = protein, coloring = "")   # nicht "coloring" übergeben, da das "all" enhalten kann, wodurch auch die negtiven int_pos hier bei den positiven dabei wären
    
    scaled_effects <- matrix(nrow = dim(effects)[1], ncol = 0)
    rownames(scaled_effects) <- rownames(effects)
    
    for (i in 1:dim(effects)[2]) {
      current_effects <- as.matrix(effects[,i])
      rownames(current_effects) <- rownames(effects)
      
      if (dim(effects)[2] <= 1) {
        index <- ""
      } else {
        index <- i
      }
      current_outpath <- outpath_for_ida(outpath = outpath, direction = dir, relatve_effects_on_pos = relatve_effects_on_pos, option_nr = index, neg_effects = neg_effects, perturbated_position = perturbated_position,
                                 amplification_exponent = amplification_exponent, amplification_factor = amplification_factor, 
                                 no_colors = no_colors, rank_effects = rank_effects, effect_to_color_mode = effect_to_color_mode)
      
      current_scaled_effects <- scale_effects(current_effects, rank = rank_effects, amplification_factor = amplification_factor, neg_effects = neg_effects)
      scaled_effects <- cbind(scaled_effects, current_scaled_effects)
      
      colors_by_effect <- color_by_effect(current_scaled_effects, int_pos, mode = effect_to_color_mode)
      
      
      if (!show_neg_causation) {
        current_effects <- NULL
      }
      plot_total_effects_in_pymol(positions_with_colors_by_effect = colors_by_effect, perturbated_position = perturbated_position, protein = protein, 
                                  outpath = current_outpath,
                                  amplification_exponent = amplification_exponent, amplification_factor = amplification_factor, ranked = opacity_ranked, 
                                  index = i, no_colors = no_colors, bg_color = pymol_bg_color, orig_effects = current_effects)
    
      if (barplot) {
        # graphics.off()
        # par(mfrow=c(m,n))
        
        vector_effects <- as.vector(current_effects)
        names(vector_effects) <- rownames(current_effects)
        barplot(vector_effects, main = paste(caption, "\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
        vector_scaled_effects <- as.vector(current_scaled_effects)
        names(vector_scaled_effects) <- rownames(current_scaled_effects)
        barplot(vector_scaled_effects, main = paste(caption, "\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
      }
        
      if (analysis) {
        threshold = quantile(current_effects, probs = percentile)
        most_influenced_positions <- colnames(data[(rownames(current_effects)[which(current_effects > threshold)])])
        print(paste(length(most_influenced_positions), "positions over the threshold", threshold, ": ", paste(most_influenced_positions, collapse = ", ")))
        # int_pos <- interesting_positions("PDZ", "crystal")
        int_pos_strongly_influenced <- intersect(int_pos, most_influenced_positions)
        print(paste("Thereof",  length(int_pos_strongly_influenced), "out of the", length(int_pos), "interesting positions:", paste(sort(int_pos_strongly_influenced), collapse = ", ")))
        print(paste("Missing: ", paste(setdiff(int_pos, most_influenced_positions), collapse = ", ")))
        
        print("FOR SCALED EFFECTS:")
        threshold = quantile(current_scaled_effects, probs = percentile)
        most_influenced_positions <- colnames(data[(rownames(current_scaled_effects)[which(current_scaled_effects > threshold)])])
        print(paste(length(most_influenced_positions), "positions over the threshold", threshold, ": ", paste(most_influenced_positions, collapse = ", ")))
        # int_pos <- interesting_positions("PDZ", "crystal")
        int_pos_strongly_influenced <- intersect(int_pos, most_influenced_positions)
        print(paste("Thereof",  length(int_pos_strongly_influenced), "out of the", length(int_pos), "interesting positions:", paste(sort(int_pos_strongly_influenced), collapse = ", ")))
        print(paste("Missing: ", paste(setdiff(int_pos, most_influenced_positions), collapse = ", ")))
      
      }
    }
  
    # results$ida <- list(list(effects = effects, scaled_effects = scaled_effects, pos_colors = colors_by_effect))
    # names(results$ida) <- perturbated_position
    
    results$ida[[perturbated_position]][[dir]] <- list(effects = effects, scaled_effects = scaled_effects, pos_colors = colors_by_effect)
    
    # print(cbind(effects, current_scaled_effects, colors_by_effect))
  }
  return(results)
}

# only to store code
first_tests <- function(data, results) {
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
}

### relevant measures: conflict_edges/(bidirected_edges + unidirected_edges)
### or conflict_edges/(2 * bidirected_edges + unidirected_edges) (number of directed edges)
conflict_edges <- function(graph) {
  adj_m <- wgtMatrix(graph)
  # m <- adj_m[apply(adj_m!=0, 1, any), , drop=FALSE]
  # adj_m_no_zero_rows_cols <- m[apply(m!=0, 2, any), , drop=FALSE]
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m, 1)
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m_no_zero_rows_cols, 2)
  print(adj_m_no_zero_rows_cols)
  n_conflic_edges <- length(which(adj_m == 2)) / 2
  # = length(which(unlist(edgeData(graph)) == 2))
  n_unidirected_edges <- length(which(adj_m != t(adj_m))) / 2
  
  n_bidirected_edges <- (length(which(adj_m == 1)) - n_unidirected_edges) / 2
  
  print(paste("<-!-> : ", n_conflic_edges))
  print(paste("<---> : ", n_bidirected_edges))
  print(paste("----> : ", n_unidirected_edges))
  
  return(list(conflict = n_conflic_edges, directed = n_unidirected_edges, bidirected = n_bidirected_edges))
}

remove_zero_rows_or_columns <- function(matrix, rows_or_columns) {
  if (rows_or_columns == 1) {
    return(matrix[rowSums(abs(matrix)) != 0, ])
  } else {
    return(matrix[,colSums(abs(matrix)) != 0])
  }
  # return(matrix[apply(matrix != 0, rows_or_columns, any), , drop = TRUE])
}
