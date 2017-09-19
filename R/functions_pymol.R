pymol_header <- function(protein, pdb_file, chain = "all", file_separator = "/") {
  if (missing(pdb_file)) {
    pdb_file <- paste("..", "..", "..", "", sep = file_separator)
    if (protein == "GTB") {
      # pdb_file <- "../../../5bxc.pdb" 
      pdb_file <- paste0(pdb_file, "5bxc.pdb")
    } else if (protein == "PDZ" || protein == "pdz") {
      # pdb_file <- "../../../1BE9.pdb"
      pdb_file <- paste0(pdb_file, "1BE9.pdb")
    } else if (protein == "p38g") {
      # pdb_file <- "../../../1cm8.pdb"
      pdb_file <- paste0(pdb_file, "1cm8.pdb")
      chain = "chain A"
    } else {
      stop("No pdb-file given.")
    }
  }
  cat("delete all\n")
  cat("load", pdb_file,"\n")
  cat("hide all\n")
  cat("show cartoon,", chain, "\n")
  cat("color white\n")
  if (protein == "GTB") {
    # acceptor
    cat("show sticks, resi 401\n")
    # cat("show spheres, resi 401\n") 
    cat("set_color acceptor_green, [", paste(col2rgb("#69A019") / 255, collapse = ","), "] \n", sep = "")
    cat("color acceptor_green, resi 401\n")
    # cat("label  i. 401, \"acc\"\n")
    # donor
    cat("show sticks, resi 403\n")
    # cat("show spheres, resi 403\n") 
    cat("set_color donor_yellow, [", paste(col2rgb("#FFD700") / 255, collapse = ","), "] \n", sep = "")
    cat("color donor_yellow, resi 403\n")
    # cat("label  i. 403, \"don\"\n")
    # Mn
    # cat("show spheres, resi 402\n") 
    # cat("color grey, resi 402\n")
  } else if (protein == "PDZ" || protein == "pdz") {
    cat("show sticks, chain B\n")
    # cat("set_color ligand_yellow, [", paste(col2rgb("#FFD700") / 255, collapse = ","), "] \n", sep = "")
    cat("color gray40, chain B\n")
  }
  cat("\n")
  cat("set label_position,(-2,-2,0)\n")
  cat("\n")
}

# most general; TODO: can the other fucntions (except for the paths) be implemented using this one? 
plot_clusters_in_pymol <- function(node_clustering, protein, outpath, pdb_file, 
                                   label = TRUE, no_colors = FALSE, show_positions = TRUE,
                                   file_separator = "/", type_of_clustering = "") {
  
  out_file <- paste0(outpath, "-", length(node_clustering), pastes("_clusters", type_of_clustering, sep = "-"),".pml")
  # out_file <- pastes(out_file, type_of_clustering, sep = "-")
  # out_file <- paste0(out_file, ".pml")
  # print(out_file)
  
  sink(file = out_file)
  pymol_header(protein = protein, file_separator = file_separator, pdb_file = pdb_file)
  if (!is.null(names(node_clustering))) {
    colors <- names(node_clustering)
  } else {
    colors <- rainbow(length(node_clustering))
  }
  
  for (i in 1:length(node_clustering)) {
    cat("create sector_", i, ", (resi ", 
        paste(node_clustering[[i]], collapse = ","), ") \n", sep = "")
    if (show_positions) {
      cat("show spheres, sector_", i, "\n", sep = "")
    }
    if (!no_colors) {
      color <- col2rgb(colors[i])
      color <- color / 255
      cat("set_color col_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
      cat("color col_", i, ", sector_", i, "\n", sep = "")
    }
    if (label) {
      # label n. CA and sector_1, resi
      cat("label n. CA and sector_", i, ", resi\n", sep = "") # label only CA
      # cat("label sector_", i, ", resi\n", sep = "")
    }
    cat("\n")
  }
  
  cat("zoom\n")
  sink()
}

# TODO: use file_separator consistently
plot_connected_components_in_pymol <- function(protein, position_numbering, graph, outpath, label = TRUE, pdb_file, only_int_pos = FALSE, 
                                               show_int_pos = TRUE, color_int_pos = TRUE, only_color_int_pos = FALSE, coloring_for_int_pos, no_colors = FALSE, only_dist = FALSE, 
                                               show_positions = TRUE, file_separator = "/") {
  print(paste("Outpath for pymol-file:", outpath))
  connected_components <- connComp(graph)
  real_ones_ind <- which(sapply(connected_components, function(x) length(x) > 1))
  connected_components <- connected_components[real_ones_ind]
  if (length(connected_components) == 0) {
    return(NULL)
  }
  if (only_int_pos) {
    int_pos_flat <- interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = FALSE)
    connected_components <- lapply(connected_components, function(positions) return(
      positions[sapply(positions, function(pos) return(
        as.numeric(pos) %in% int_pos_flat))]))
  }
  
  # change this to be able to print for example only the found links between interesting positions, whcih are colored accoringly
  if (only_color_int_pos) {
    # TODO: outpath aus dem übergebenen outpath konstruieren (per strsplit)
    # directories <- strsplit(outpath, "/")
    # outpath <- paste(directories[[1]][1:(length(directories[[1]])-3)], collapse = "/", sep = "/")
    outpath <- paste("../Outputs/", protein, "/", sep = "")
    if (grepl("all", coloring_for_int_pos)) {
      out_file <- paste(outpath, "interesting-all.pml", sep = "")
    } else {
      out_file <- paste(outpath, "interesting.pml", sep = "")
    }
  } else {
    out_file <- paste(outpath, ".pml", sep = "")
  }
  sink(file = out_file)
  pymol_header(protein = protein, file_separator = file_separator, pdb_file = pdb_file)
  colors <- rainbow(length(connected_components))
  if (!only_dist) {
    # if (!show_int_pos) {
    for (i in 1:length(connected_components)) {
      cat("create sector_", i, ", (resi ", 
          paste(connected_components[[i]], collapse = ","), ") \n", sep = "")
      # if (show_positions) {
      cat("show spheres, sector_", i, "\n", sep = "")
      # }
      if (!no_colors) {
        color <- col2rgb(colors[i])
        color <- color / 255
        cat("set_color col_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
        cat("color col_", i, ", sector_", i, "\n", sep = "")
      }
      if (label) {
        # label n. CA and sector_1, resi
        cat("label n. CA and sector_", i, ", resi\n", sep = "") # label only CA
        # cat("label sector_", i, ", resi\n", sep = "")
      }
      cat("\n")
    }
    # } else {
    if (show_int_pos) {
      if (missing(coloring_for_int_pos) && !no_colors) {
        int_pos <- interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = TRUE)
      } else {
        int_pos <- interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = TRUE, coloring = coloring_for_int_pos)
      }
      for (i in 1:length(int_pos)) {
        cat("create sector_interesting_", i, ", (resi ",
            paste(int_pos[[i]], collapse = ","), ") \n", sep = "")
        # if (show_positions) {
        cat("show surface, sector_interesting_", i, "\n", sep = "")
        if (color_int_pos) {
          color <- col2rgb(names(int_pos[i]))
          color <- color / 255
          cat("set_color col_int_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
          cat("color col_int_", i, ", sector_interesting_", i, "\n", sep = "")
        }
        # }
        # if (!no_colors) {
        #   color <- col2rgb(names(int_pos)[[i]])
        #   color <- color / 255
        #   position <- int_pos[]
        #   cat("set_color col_interesting_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
        #   cat("color col_interesting_", i, ", sector_interesting_", i, "\n", sep = "")
        # }
      }
    }
  }
  # cat("show surface \n")
  # cat("set transparency, 0.4\n")
  edge_list <- as_edgelist(graph_from_graphnel(graph))
  if (only_int_pos) {
    # bereits oben berechnet
    # int_pos_flat <- interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = FALSE)
    edge_list_is_interesting <-  apply(edge_list, 1:2, function(x) return(as.numeric(x) %in% int_pos_flat))
    # edge_list <- # jeah, what? 
  }
  draw_edge <- function(nodes) {
    cat("distance i.", nodes[1], "and n. CA, i.", nodes[2], "and n. CA\n")
  }
  apply(edge_list, 1, draw_edge)
  numb_edges <- dim(edge_list)[1]
  number_of_digits_to_pad_to <- ceiling(log(numb_edges, base=10))  # does not seem to happen
  number_of_digits_to_pad_to <- 2
  for (i in 1:numb_edges) {
    cat("hide labels, dist", str_pad(i, number_of_digits_to_pad_to, pad = "0"), "\n", sep = "")
  }
  cat("zoom\n")
  sink()
}


plot_paths_in_pymol <- function(protein, pdb_file, graph, outpath, paths, no_colors = FALSE, label = TRUE, show_positions = TRUE, file_separator = "/") {
  
  out_file <- paste(outpath, "-paths.pml", sep = "")  # welche Pfade - hinzufügen
  
  sink(file = out_file)
  pymol_header(protein = protein, file_separator = file_separator, pdb_file = pdb_file)
  colors <- rainbow(length(paths))
  for (i in 1:length(paths)) {
    if (length(paths[[i]]) > 0) {
      cat("create sector_path_", i, ", (resi ",
          paste(paths[[i]], collapse = ","), ") \n", sep = "")
      if (show_positions) {
        cat("show spheres, sector_path_", i, "\n", sep = "")        #TODO: fix
      }
      if (!no_colors) {
        color <- col2rgb(colors[i])
        color <- color / 255
        cat("set_color col_path_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
        cat("color col_path_", i, ", sector_path_", i, "\n", sep = "")
      }
      if (label) {
        # label n. CA and sector_1, resi
        cat("label n. CA and sector_path_", i, ", resi\n", sep = "") # label only CA
      }
    } 
  }
  
  sapply(paths, function(path) {mapply(function(x,y) {return(cat("distance i.", x, "and n. CA, i.", y, "and n. CA\n"))}, path[-length(path)], path[-1])})
  
  
  # draw_edge <- function(nodes) {
  #   cat("distance i.", nodes[1], "and n. CA, i.", nodes[2], "and n. CA\n")
  # }
  # apply(edge_list, 1, draw_edge)
  
  numb_edges <- do.call(sum, lapply(paths, function(path) length(path) - 1))
  # number_of_digits_to_pad_to <- ceiling(log(numb_edges, base = 10))  # does not seem to happen
  if (numb_edges > 0) {
    number_of_digits_to_pad_to <- 2
    for (i in 1:numb_edges) {
      cat("hide labels, dist", str_pad(i, number_of_digits_to_pad_to, pad = "0"), "\n", sep = "")
    }
  }
  cat("zoom\n")
  sink()
}


# names of positions_with_colors_by_effect are positions, values effect-adjusted colors
# if original_effects given, and no_colors (otherwise), the negatively influenced positions are colored in red.
plot_total_effects_in_pymol <- function(positions_with_colors_by_effect, perturbated_position, protein, outpath, label = TRUE, ranked = TRUE,
                                        amplification_exponent = 10, amplification_factor = FALSE, index = "", no_colors = FALSE, bg_color = "black", orig_effects) {
  # out_file <- paste0(outpath, "-total_effects")
  
  sink(file = outpath)
  pymol_header(protein = protein)
  
  if (no_colors && !is.null(orig_effects)) {     # color for negatively influenced positions
    color_neg <- col2rgb("#CC0000")
    color_neg <- color_neg / 255
    cat("set_color col_neg, [", paste(color_neg, collapse = ","), "] \n", sep = "")
    
    color_pos <- col2rgb("#FFD700")
    color_pos <- color_pos / 255
    cat("set_color col_pos, [", paste(color_pos, collapse = ","), "] \n", sep = "")
  } 
  
  for (i in 1:length(positions_with_colors_by_effect)) {
    pos <- names(positions_with_colors_by_effect)[[i]]
    
    cat("show spheres, resi ", pos, "\n", sep = "")
    
    if (label) {
      # label n. CA and sector_1, resi
      cat("label n. CA and resi ", pos, ", resi\n", sep = "") # label only CA
      # cat("label sector_", i, ", resi\n", sep = "")
    }
    
    if (!no_colors || (pos == as.integer(perturbated_position))) {
      if ((no_colors)) {  ## (&& (pos == as.integer(perturbated_position) )
        color <- col2rgb("#69A019")
      } else {
        color <- col2rgb(positions_with_colors_by_effect[i])
      }
      color <- color / 255
      cat("set_color col_", i, ", [", paste(color, collapse = ","), "] \n", sep = "")
      cat("color col_", i, ", resi ", pos, "\n", sep = "")
    } else if (no_colors && !is.null(orig_effects)) {
      if (orig_effects[i] < 0) {
        cat("color col_neg, resi ", pos, "\n", sep = "")
      } else {
        cat("color col_pos, resi ", pos, "\n", sep = "")
      }
    }
    
    alpha <- col2rgb(positions_with_colors_by_effect[i], alpha = TRUE)[4] / 255
    cat("set sphere_transparency=", (1 - alpha) ^ amplification_exponent, ", resi ", pos, "\n", sep = "")
    # cat("set transparency, 0.4, n. ", pos, "\n", sep = "")
  }
  
  # if (no_colors) {
  #   cat("color black, chain A\n")
  # }
  
  cat("\n")
  cat(paste0("bg_color ", bg_color, "\n"))
  cat("set cartoon_color, gray75\n")
  cat("zoom\n")
  sink()
}

