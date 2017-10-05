# source('/Volumes/Causality/Viren/R/Code/compute_DAG_G.R')
# source("general_functions.R")

# adds results$ida

library(colorspace)  # for mixcolor, hex

# parameter barplot = TRUE remove, use mute_allplots = FALSE instead
causal_effects_ida <- function(data, perturbated_position, direction = "both", weight_effects_on_by = "mean_abs_effect",
                results = results, protein, coloring = "all", effect_hue_by = "effect", #effect_hue_by = "variance",
                outpath,amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, 
                effect_to_color_mode = "#FFFFFF", pymol = TRUE, pymol_bg_color = "black", caption, no_colors, 
                show_neg_causation = TRUE, neg_effects = "", analysis = TRUE, percentile = 0.75, mute_all_plots = FALSE,
                causal_effects_function = "IDA-reset") {
  
  lines <- 1 # lines in plot_space
  
  if (direction == "both") {
    direction <- c("of", "on")
    lines <- 2
  }
  
  if (!mute_all_plots) {
    graphics.off()
    # lines <- 6 # 1/2 für of, 2 für on (min, max)
    columns <- 2
    oma <- c( 0, 0, length(caption) + 1, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
    # par(mfrow=c(lines, columns), oma = oma) 
    par(oma = oma)
  }

  cat("\n")
  for (dir in direction) {
    if (dir == "of" || dir == "by" || dir == "from") {
      if (grepl(pattern = "ida", tolower(causal_effects_function))) {
        # IDA
        # effects <- idaFast(which(as.character(colnames(data)) == perturbated_position), 1:dim(data)[2], data, results$pc@graph)
        effects <- idaFast(which(as.character(colnames(data)) == perturbated_position), 1:dim(data)[2], cov(data), results$pc@graph)
        if (grepl(pattern = "reset", tolower(causal_effects_function))) {
        # GGF. danach Effekte von nicht verbundenen Knoten auf null setzten
        effects <- set_effects_of_unconnected_positions_to_zero(effects, graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
        }
      } else if (grepl(pattern = "causaleffect", tolower(causal_effects_function)) || grepl(pattern = "causal_effect", tolower(causal_effects_function))) {
        # CAUSALEFFECTS
        effects <- pseudo_idaFast_by_causalEffect(which(as.character(colnames(data)) == perturbated_position), 1:dim(data)[2], 
                                                cov(data), results$pc@graph, outpath = outpath)
      }
    } else if (dir == "on") {
      if (grepl(pattern = "ida", tolower(causal_effects_function))) {
        # IDA
        ida_rev <- function(pos) {
          # eff <- pcalg::ida(pos, which(as.character(colnames(data)) == perturbated_position), data, results$pc@graph))
          eff <- pcalg::ida(pos, which(as.character(colnames(data)) == perturbated_position), cov(data), results$pc@graph)
          return(eff)
        }
        # if (grepl(pattern = "reset", tolower(causal_effects_function))) {
        #   ida_rev <- function(pos) {
        #     # eff <- pcalg::ida(pos, which(as.character(colnames(data)) == perturbated_position), data, results$pc@graph))
        #     eff <- pcalg::ida(pos, which(as.character(colnames(data)) == perturbated_position), cov(data), results$pc@graph)
        #     # GGF. danach Effekte von nicht verbundenen Knoten auf null setzten
        #     eff <- set_effects_of_unconnected_positions_to_zero(eff, graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
        #     return(eff)
        #   }
        # }
      } else if (grepl(pattern = "causaleffect", tolower(causal_effects_function)) || grepl(pattern = "causal_effect", tolower(causal_effects_function))) {
        # CAUSALEFFECTS
        ida_rev <- function(pos) {
          eff <- pseudo_ida_by_causalEffect(pos, which(as.character(colnames(data)) == perturbated_position), 
                                            cov(data), results$pc@graph, outpath = outpath)
          return(eff)
        }
      }
      effects_list <- lapply(1:dim(data)[2], ida_rev)
      # vllt lieber relativ zu dem durchsnittlichen effekt von Position pos 
      effects_max <- sapply(effects_list, function(list) return(max(list)))
      effects_min <- sapply(effects_list, function(list) return(min(list)))
      effects <- cbind(effects_max, effects_min)
      colnames(effects) <- c("max", "min")
      rownames(effects) <- colnames(data)
      
      if (grepl(pattern = "reset", tolower(causal_effects_function))) {
        effects <- apply(X = effects, MARGIN = 2, FUN = set_effects_of_unconnected_positions_to_zero, 
                         graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
      }
      
      if (grepl("mean", weight_effects_on_by)) { # (weight_effects_on_by == "mean_abs_effect") {
        dir <- "on-rel-to-mean"
        
        if (grepl(pattern = "ida", tolower(causal_effects_function))) {
          # IDA
          # means <- sapply(1:dim(data)[2], function(pos) {mean(abs(idaFast(pos, 1:dim(data)[2], data, results$pc@graph)))})
          means <- sapply(1:dim(data)[2], function(pos) {
           eff <- mean(abs(idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)))
           return(eff)
           })
          # if (grepl(pattern = "reset", tolower(causal_effects_function))) {
          #   # IDA
          #   means <- sapply(1:dim(data)[2], function(pos) {
          #    eff <- idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)
          #    # GGF. danach Effekte von nicht verbundenen Knoten auf null setzten
          #    eff <- set_effects_of_unconnected_positions_to_zero(eff, graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
          #    return(mean(abs(eff)))
          #    })
          # }
        } else if (grepl(pattern = "causaleffect", tolower(causal_effects_function)) || grepl(pattern = "causal_effect", tolower(causal_effects_function))) {
          # CAUSALEFFECTS
          means <- sapply(1:dim(data)[2], function(pos) {mean(abs(pseudo_idaFast_by_causalEffect(pos, 1:dim(data)[2], 
                                                                                                 cov(data), results$pc@graph, outpath = outpath)))})
        }
        if (grepl(pattern = "reset", tolower(causal_effects_function))) {
          means <- apply(X = means, MARGIN = 2, FUN = set_effects_of_unconnected_positions_to_zero, 
                           graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
        }
        effects <- effects / means
      } else if (grepl("median", weight_effects_on_by)) { # if (weight_effects_on_by == "median_abs_effect") {
        dir <- "on-rel-to-median"
        if (grepl(pattern = "ida", tolower(causal_effects_function))) {
          # IDA
          # medians <- sapply(1:dim(data)[2], function(pos) {median(abs(idaFast(pos, 1:dim(data)[2], data, results$pc@graph)))})
          medians <- sapply(1:dim(data)[2], function(pos) {
            eff <- median(abs(idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)))
            return(eff)
          })
          # if (grepl(pattern = "reset", tolower(causal_effects_function))) {
          #   # IDA
          #   medians <- sapply(1:dim(data)[2], function(pos) {
          #     eff <- idaFast(pos, 1:dim(data)[2], cov(data), results$pc@graph)
          #     # GGF. danach Effekte von nicht verbundenen Knoten auf null setzten
          #     eff <- set_effects_of_unconnected_positions_to_zero(eff, graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
          #     return(median(abs(eff)))
          #   })
          # }
        } else if (grepl(pattern = "causaleffect", tolower(causal_effects_function)) || grepl(pattern = "causal_effect", tolower(causal_effects_function))) {
          # CAUSALEFFECTS
          medians <- sapply(1:dim(data)[2], function(pos) {median(abs(pseudo_idaFast_by_causalEffect(pos, 1:dim(data)[2], 
                                                                                                 cov(data), results$pc@graph, outpath = outpath)))})
        }
        
        names(medians) <- colnames(data)
        medians = medians / medians[as.character(perturbated_position)]
        
        if (grepl(pattern = "reset", tolower(causal_effects_function))) {
          medians_m <- as.matrix(medians)
          rownames(medians_m) <- colnames(data)
          medians <- apply(X = medians_m, MARGIN = 2, FUN = set_effects_of_unconnected_positions_to_zero, 
                         graph = results$pc@graph, perturbed_position = perturbated_position, dir = dir)
        }
        # effects <- effects / medians
        effects <- apply(effects, 2, function(effects) return(effects/medians))
      } else if (weight_effects_on_by == "var" || weight_effects_on_by == "vars") {
        dir <- "on-rel-to-var"
        vars <- apply(data, 2, var)
        effects <- effects * vars
      }
    }
    rownames(effects) <- colnames(data)
    
    # avoid errors
    # TODO: make less restrictive; maybe only if all values non-numbers
    effects[is.na(effects)] <- 0
    effects[is.infinite(effects)] <- 1
    
    
    # print(effects)
    cat("CAUSAL EFFECTS ")
    cat(toupper(dir))
    cat(" POSITION ")
    cat(perturbated_position)
    cat("\n")
    
    int_pos <- interesting_positions(protein = protein, coloring = "")   # nicht "coloring" übergeben, da das "all" enhalten kann, wodurch auch die negtiven int_pos hier bei den positiven dabei wären
    
    scaled_effects <- matrix(nrow = dim(effects)[1], ncol = 0)
    rownames(scaled_effects) <- rownames(effects)
    
    if (dir == "of" && ("of" %in% direction)) {
      lines <- dim(effects)[2] + 2
    } else if (dir == "of" && !("of" %in% direction)) {
      lines <- dim(effects)[2]
    } else if (grepl("on", dir) && direction == "on") {
      lines <- 2
    } 
    
    if (!(grepl("on", dir) && ("of" %in% direction)) && !mute_all_plots) {
      while (lines > 8) {
        lines = lines / 2
      }
      par(mfrow=c(ceiling(lines), columns))
    }
      
    
    for (i in 1:dim(effects)[2]) {
      cat("OPTION ")
      cat(i)
      cat("\n")
      
      current_effects <- as.matrix(effects[,i])
      rownames(current_effects) <- rownames(effects)
      
      if (dim(effects)[2] <= 1) {
        index <- ""
      } else {
        index <- i
      }
      current_outpath <- outpath_for_ida(outpath = outpath, direction = dir, option_nr = index, neg_effects = neg_effects, perturbated_position = perturbated_position,
                                 amplification_exponent = amplification_exponent, amplification_factor = amplification_factor, 
                                 no_colors = no_colors, rank_effects = rank_effects, effect_to_color_mode = effect_to_color_mode)
      
      current_scaled_effects <- scale_effects(current_effects, rank = rank_effects, amplification_factor = amplification_factor, neg_effects = neg_effects)
      scaled_effects <- cbind(scaled_effects, current_scaled_effects)
      
      if (effect_hue_by == "effect") {
        colors_by_effect <- color_by_effect(current_scaled_effects, int_pos, mode = effect_to_color_mode)
      } else if (effect_hue_by == "variance" || effect_hue_by == "var") {
        vars <- apply(data, 2, var)
        colors_by_effect <- color_by_effect(vars, int_pos, mode = effect_to_color_mode)
      }
      
      if (!show_neg_causation) {
        current_effects <- NULL
      }
      
      if (pymol) {
        plot_total_effects_in_pymol(positions_with_colors_by_effect = colors_by_effect, perturbated_position = perturbated_position, 
                                    protein = protein, outpath = current_outpath, amplification_exponent = amplification_exponent, 
                                    amplification_factor = amplification_factor, ranked = opacity_ranked,
                                    index = i, no_colors = no_colors, bg_color = pymol_bg_color, orig_effects = current_effects)
      }
    
      # if (barplot) {
      if (!mute_all_plots) {
        # graphics.off()
        # par(mfrow=c(m,n))
        # plot.new()
        # title("My 'Title' in a strange place", side = 3, line = -21, outer = TRUE)
        # mtext( "Centered Overall Title", outer = TRUE )
        
        vector_effects <- as.vector(current_effects)
        is.na(vector_effects) <- sapply(vector_effects, is.infinite)   # set infinite values to NA
        names(vector_effects) <- rownames(current_effects)
        # barplot(vector_effects, main = paste(caption, "\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
        barplot(vector_effects, main = paste("\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
        vector_scaled_effects <- as.vector(current_scaled_effects)
        is.na(vector_scaled_effects) <- sapply(vector_scaled_effects, is.infinite)   # set infinite values to NA
        names(vector_scaled_effects) <- rownames(current_scaled_effects)
        # barplot(vector_scaled_effects, main = paste(caption, "\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
        barplot(vector_scaled_effects, main = paste("\n total causal effect", dir, "position", perturbated_position), col = colors_by_effect)
      }
        
      if (analysis) {
        statistics_of_influenced_positions(effects = current_effects, percentile = percentile, interesting_positions = int_pos, print = TRUE)
        
        
        # print("FOR SCALED EFFECTS:")        # sollte es immer gleich sein
        # statistics_of_influenced_positions(effects = current_scaled_effects, percentile = percentile, interesting_positions = int_pos, print = TRUE)
        
        # threshold = quantile(current_scaled_effects, probs = percentile)
        # most_influenced_positions <- colnames(data[(rownames(current_scaled_effects)[which(current_scaled_effects > threshold)])])
        # print(paste(length(most_influenced_positions), "positions over the threshold", threshold, ": ", paste(most_influenced_positions, collapse = ", ")))
        # # int_pos <- interesting_positions("PDZ", "crystal")
        # int_pos_strongly_influenced <- intersect(int_pos, most_influenced_positions)
        # print(paste("Thereof",  length(int_pos_strongly_influenced), "out of the", length(int_pos), "interesting positions:", paste(sort(int_pos_strongly_influenced), collapse = ", ")))
        # print(paste("Missing: ", paste(setdiff(int_pos, most_influenced_positions), collapse = ", ")))
      
      }
    }
  
    # results$ida <- list(list(effects = effects, scaled_effects = scaled_effects, pos_colors = colors_by_effect))
    # names(results$ida) <- perturbated_position
    
    results$ida[[perturbated_position]][[dir]] <- list(effects = effects, scaled_effects = scaled_effects, pos_colors = colors_by_effect)
    
    # print(cbind(effects, current_scaled_effects, colors_by_effect))
  }
  if (!mute_all_plots) {
    title(caption, outer = TRUE)
  }
  return(results)
}

get_most_influenced_positions <- function(effects, threshold, percentile) {
  if (missing(threshold)) {
    if (missing(percentile)) {
      stop("Neither threshold nor quantile given in get_most_influenced_positions.")
    }
    threshold = quantile(effects, probs = percentile, na.rm = TRUE)
  }
  if (!is.null(dim(effects))) {
    # most_influenced_positions <- colnames(data[(rownames(effects)[which(effects > threshold)])])
    most_influenced_positions <- (rownames(effects)[which(effects > threshold)])
  } else {
    # most_influenced_positions <- colnames(data[(names(effects)[which(effects > threshold)])])
    most_influenced_positions <- (names(effects)[which(effects > threshold)])
  }
  # int_pos <- interesting_positions("PDZ", "crystal")
  return(most_influenced_positions)
}

statistics_of_influenced_positions <- function(effects, percentile, interesting_positions, print = FALSE, 
                                               verbose = FALSE, return_list = TRUE) { #TODO: list = TRUE setzen
  threshold = quantile(effects, probs = percentile, na.rm = TRUE)
  most_influenced_positions <- get_most_influenced_positions(effects = effects, threshold = threshold)
 if (print && verbose) {
   writeLines(paste0(length(most_influenced_positions), " positions over the threshold ", threshold, ": ", 
              paste(most_influenced_positions, collapse = ", ")))
 }
  # int_pos <- interesting_positions("PDZ", "crystal")
  
  int_pos_strongly_influenced <- intersect(interesting_positions, most_influenced_positions)
  fp = setdiff(most_influenced_positions, int_pos_strongly_influenced)
  if(length(fp) == 0) fp = "-"
  fn = setdiff(interesting_positions, most_influenced_positions)
  if(length(fn) == 0) fn = "-"
  sub = paste0("FP: ", paste(fp, collapse = " "), "\nFN: ", paste(fn, collapse = " "))
  
  if (print && verbose) {
    writeLines(paste("Thereof",  length(int_pos_strongly_influenced), "out of the", length(interesting_positions), 
                "interesting positions:", paste(sort(int_pos_strongly_influenced), collapse = ", ")))
    writeLines(paste("and also (FALSE POSITIVE): ", paste(setdiff(most_influenced_positions, int_pos_strongly_influenced), collapse = ", ")))
    writeLines(paste("missing  (FALSE NEGATIVE): ", paste(setdiff(interesting_positions, most_influenced_positions), collapse = ", ")))
    
  } else {
    writeLines(sub)
  }
  # cat("\n")
  if (return_list) {
    return(list(fp = fp, fn = fn))
  } else {
    return(sub)
  }
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
# calls edge_information with full adj_matrix if interesting_positions is missing or the rows and columns of that vector if it is given
conflict_edges <- function(graph, print = FALSE, interesting_positions) {
  gm <- wgtMatrix(graph)
  if(missing(interesting_positions)) {
    return(edge_information(gm, print = print))
  } else {
    igm <- gm[as.character(interesting_positions), as.character(interesting_positions)]
    return(edge_information(igm, print = print))
  }
}

### relevant measures: conflict_edges/(bidirected_edges + unidirected_edges)
### or conflict_edges/(2 * bidirected_edges + unidirected_edges) (number of directed edges)
edge_information <- function(adj_m, print = FALSE) {

  # m <- adj_m[apply(adj_m!=0, 1, any), , drop=FALSE]
  # adj_m_no_zero_rows_cols <- m[apply(m!=0, 2, any), , drop=FALSE]
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m, 1)
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m_no_zero_rows_cols, 2)
  # print(adj_m_no_zero_rows_cols)
  n_conflic_edges <- length(which(adj_m == 2)) / 2
  # = length(which(unlist(edgeData(graph)) == 2))
  # this works as there are no unidirected conflict edges
  n_unidirected_edges <- length(which(adj_m != t(adj_m))) / 2
  
  n_bidirected_edges <- (length(which(adj_m == 1)) - n_unidirected_edges) / 2
  
  if (print) {
    print(paste("<-!-> : ", n_conflic_edges))
    print(paste("<---> : ", n_bidirected_edges))
    print(paste("----> : ", n_unidirected_edges))
  }

  return(list(conflict = n_conflic_edges, directed = n_unidirected_edges, undirected = n_bidirected_edges))
}

remove_zero_rows_or_columns <- function(matrix, rows_or_columns) {
  if (rows_or_columns == 1) {
    return(matrix[rowSums(abs(matrix)) != 0, ])
  } else {
    return(matrix[,colSums(abs(matrix)) != 0])
  }
  # return(matrix[apply(matrix != 0, rows_or_columns, any), , drop = TRUE])
}

set_effects_of_unconnected_positions_to_zero <- function(effects, graph, perturbed_position, dir) {
  if (dir == "of") {
    mode <- "out"
  } else if (grepl("on", dir)) {
    mode <- "in"
  }
  
  igraph <- igraph.from.graphNEL(graphNEL = graph)
  dist <- distances(graph = igraph, v = perturbed_position, mode = mode)
  
  effects[dist == Inf] <- 0
  
  return(effects)
}
