library(graph)
library(dagitty)
library(pcalg)

# source("evaluate_DAG.R")
# source("linkcommunities.R")

protein_causal_graph <- function(data, protein, type_of_data, source_of_data, position_numbering, output_dir, filename, outpath,
                                 parameters_for_info_file, alpha, caption, analysis, stages, plot_types, coloring, colors, 
                                 graph_layout = "dot", plot_as_subgraphs = plot_as_subgraphs, 
                                 plot_only_subgraphs = plot_only_subgraphs, unabbrev_r_to_info, print_r_to_console, 
                                 lines_in_abbr_of_r, compute_pc_anew, compute_localTests_anew, print_analysis, plot_analysis) {
  print(paste("Output will be written to ", getwd(), "/", output_dir, "/...", sep = ""))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
    print("Directory created.")
  }
  
  if (missing(outpath)) {
    outpath <- paste(output_dir, filename, sep = "/")
  }
  
  # Computation of pc 
  pc_fun <- function(outpath) {
    return(estimate_DAG_from_numerical_data(data, alpha = alpha, outpath = outpath))
  }
  pc <- get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file)
  
  # garbage <- graphics.off()
  plot_graph(graph = pc@graph, caption = caption, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout, coloring = coloring, colors = colors, outpath = outpath, numerical = numerical, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs)
  
  # Analysis
  if (analysis) {
    results <- analysis_after_pc(pc, data, outpath = outpath, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout, coloring = coloring, colors = colors,  stages = stages, plot_types = plot_types, unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r, compute_localTests_anew = compute_localTests_anew, print = print_analysis, plot = plot_analysis, caption = caption)
    
    # print_analysis <- FALSE
    # plot_analysis <- FALSE
    } 
  else {
    results <- list()
    results$pc <- pc
    
    results$orig <- list()
    results$orig$graph$NEL <- pc@graph
    # results <- analysis_after_pc(pc, data, outpath = outpath, protein = protein, position_numbering = position_numbering, layout = graph_layout, coloring = coloring, colors = colors,  stages = c(), plot_types = plot_types, unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r, compute_localTests_anew = compute_localTests_anew, print = FALSE, plot = FALSE)
  }
    
  return(results)
} 

sink.reset <- function() {
  if (sink.number() > 0) {
    for (i in 1:sink.number()) {
      sink()
    }
  }
}

subgraph_of_interesting_positions <- function(graphNEL, graph_dagitty, positions = NULL, protein = NULL, position_numbering = NULL) {
  if (is.null(positions)) {
    positions <- interesting_positions(protein, position_numbering)
  }
  subgraph <- subGraph(positions, graphNEL)
  return(subgraph)
}

ancestorgraph_of_interesting_positions <- function(graph_dagitty, positions = NULL, protein = NULL, position_numbering = NULL, nodename_prefix = "") {
  if (is.null(positions)) {
    positions <- paste(nodename_prefix, interesting_positions(protein, position_numbering), sep = "")
  }
  ancestor_graph <- ancestorGraph(graph_dagitty, v = positions)
  return(ancestor_graph)
}


compute_if_not_existent <- function(filename, FUN) {
  if (file.exists(paste(filename, ".RData", sep = ""))) {
    filename <- paste(filename, ".RData", sep = "")
    load(filename)
    if (!exists("data")) {
      print("The file did not contain an object of name 'data'!")
    } else {
      print(paste("Data loaded from ", filename, ".", sep = ""))
    }
  } else {
    print("Computing data.")
    data <- FUN()
    save(data, file = paste(filename, ".RData", sep = ""))
    print("Data computed and saved.")
  }
  return(data)
}

readAlignment <- function(filename) {
  if (file.exists(paste(filename, ".RData", sep = ""))) {
    filename <- paste(filename, ".RData", sep = "")
    load(filename)
    if (!exists("MSA")) {
      print("The file did not contain an object of name 'MSA'!")
    } else {
      print(paste("Alignment loaded from ", filename, ".", sep = ""))
    }
  } else {
    filename_fasta <- paste(filename, ".fasta", sep = "")
    print(paste("Loading alignment from ", filename_fasta, ".", sep = ""))
    MSA <- readAlignment(filename_fasta)
    colnames(MSA) <- seq(1:dim(MSA)[2])
    save(MSA, file = paste(filename, ".RData", sep = ""))
    print("Loading done; Object saved for later.")
  }
  return(MSA)
}

readAlignment_fasta <- function(from_file) {
  MSA <- read.alignment(file = from_file, format = "fasta", forceToLower = FALSE)  # from seqinr
  MSA_list <- sapply(MSA$seq, strsplit, "")
  MSA_mat <- t(sapply(MSA_list, rbind))
  rownames(MSA_mat) <- MSA$nam
  return(MSA_mat)
}

allAS = c("-", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "X")

remove_gaps <- function(MSA_in, threshold, n_lev, allAS, outpath) {
  gap_freq <- t(apply(MSA_in, 2, function(x) (table(factor(x, levels = allAS)) / length(x))["-"]))
  
  if (!is.null(outpath)) {
    sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
    cat("positions without more than 'remove_cols_gaps_threshold'=")
    cat(remove_cols_gaps_threshold)
    cat(" gaps (removed): ")
    if (!length(gap_freq[gap_freq > threshold]) == 0) {
      cat(colnames(MSA_in)[gap_freq > threshold])
    } else {
      cat("none")
    }
    cat("\n")
    sink()
  }

  MSA_out <- MSA_in[, gap_freq < threshold]
}

rem_cols_by_colname <- function(data, remove) {
  return(data[, !(colnames(data) %in% remove)])
}

caption <- function(protein, data, alpha, chars_per_line = 50) {
  par_string = paste("protein: ", protein, ", data: ", data, ", alpha: ", alpha, sep = "")
  caption = strwrap(par_string, width = chars_per_line)
  return(caption)
}

plus_minus_x <- function(vector, offset) {
  vector <- unlist(lapply(vector, function(n) return(seq(n - offset, n + offset))))
  return(sort(unique(vector)))
}

# for_coloring -> output hierarchical (list with different sorts of interesting positions as vectors), otherwise one vector
# std-Reihenfolge: grün-gelb-rot-blau
interesting_positions <- function(protein, position_numbering, allpositions, for_coloring = FALSE, coloring = "auto", colors = "", counts) {
  if (coloring == "none") {
   list <- list()
  } else {
    interesting_pos <- NULL     ## list?!
    if ((protein == "pdz") || (protein == "PDZ")) {
      if (missing(position_numbering)) {
        position_numbering <- "crystal"
      }
      if (grepl("crystal", position_numbering)) {
        # numbering crystal structure
        main = c(372) 
        high = c(322, 325, 329, 330, 347, 353, 362, 376, 380, 386) # interesting ones
        low = c(328, 340, 371, 385)
        rather_high = c(336, 341, 363, 365) # 352 (missing)
        # if (for_coloring && grepl("all", coloring)) {
        if (grepl("all", coloring)) {
          list <- list(main = main, high = high, low = low, rather_high = rather_high)
        } else {
          list <- list(main = main, high = high)
        }
      }
      if (position_numbering == "alignment") {
        # numbering alignment
        main = c(98)
        high = c(24, 28, 32, 33, 65, 72, 81, 102, 106, 119) # interesting ones
        list <- list(main = main, high = high)
      }
      names(list) <- c("#69A019", "#FFD700", "#CC0000", "#FF9933")[1:length(list)]
    } else if (protein == "GTB") {
      bind_donor <- c(c(121, 123, 126, 213, 346, 352), c(301, 302, 188, 211))       # source: first group:pdb, second group: GTAB-Paper von Friedemann
      bind_acceptor <- c(c(233, 245, 303, 326, 348), 266)  # source: first group: pdb, 266: GTA/B-Paper von Friedemann
      bind_donor_indirect_H2O <- c(124, 125, 351)
      bind_donor_indirct_Mn <- c(211, 213)
      close_to_bind_donor <- c() # plus_minus_x(bind_donor, 5)
      close_to_bind_donor <- c() # plus_minus_x(bind_donor, 5)
      list <- list("#FFD700" = bind_donor,             # yellow
                   "#69A019" = bind_acceptor,          # green
                   "#FF9933" = c(close_to_bind_donor, bind_donor_indirect_H2O, bind_donor_indirct_Mn),    # orange
                   "#6B8E23" = close_to_bind_donor)    # olive
    } else if (protein == "p38g") {
      if (coloring == "FS4") {
        blue_ <- c(53, 161, 215, 116, 137, 159, 154, 120, 130, 167, 219, 283, 125, 220, 209, 216, 212, 291, 287)
        green_left <- c(16, 26, 112, 170, 58, 337, 75, 87, 323, 78, 89, 109, 343) # left cluster in Fig. S4
        green_right <- c(23, 55, 30, 119, 33, 41, 141, 169, 45, 134)              # right cluster in Fig. S4
        yellow_left <- c(66, 182, 187, 186, 86, 77, 81, 149, 144, 174)
        yellow_right <- c(90, 357, 367, 150, 293, 361, 285, 333, 300)
        red_left <- c(197, 198, 201, 253, 268, 250, 262, 265)
        red_right <- c(225, 294, 238, 276, 239, 241, 277, 288, 292)
        
        list <- list("#69A019" = c(green_left, green_right),     # green
                     "#FFD700" = c(yellow_left, yellow_right),   # yellow
                     "#CC0000" = c(red_left, red_right),         # red
                     "#1874CD" = blue_)                          # blue
        
      } else #if (grepl("FS3", coloring)) {
        # if (grepl("FS3", coloring) && grepl("mix", coloring) && !(grepl("manual", coloring))) {
        #   # if (missing(counts)) {
        #   counts <- read.csv2("../Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
        #   counts <- as.matrix(counts, ncol = 4)
        #   # }
        #   if (!is.na(as.numeric(colors))) {
        #     round_categories <- as.numeric(colors)
        #   } else {
        #     round_categories <- 1 
        #   }
        #   list <- classify_nodes(counts, round_categories = round_categories, mix = TRUE, colors = colnames(counts))
        #   # fillcolor <- colors_for_nodes(clustering = node_clusters, colors = names(node_clusters))
        # } else
          if (grepl("FS3", coloring) && grepl("mix", coloring) && grepl("manual", coloring)) { # mixed manually (simple)
          # Mischungen nur, wenn mind. ein Achtel
          red_ <- c(250, 262, 265)
          red_little_blue <- c(241, 294, 276, 238) # < 3/4 rot
          red_some_blue <- c(198, 225, 239, 253, 268) # > 3/4 rot
          red_with_blue <- c(red_little_blue, red_some_blue)
          red_half_blue <- c(201, 197) # ca. 1/2 rot
          red_blue <- c(red_with_blue, red_half_blue)
          blue_ <- c(23, 55, 86, 187, 300, 333, 285, 186, 174, 144, 357, 66, 90, 150, 367, 361, 293, 292, 288, 277, 287, 220, 291, 216, 212, 45, 33, 283, 130, 209, 134, 120, 125)
          green_ <- c(26, 75, 30, 343)
          green_little_blue <- c(81, 87, 78, 16)
          green_some_blue <- c(149, 337, 77, 58, 141)
          green_with_blue <- c(green_little_blue, green_some_blue)
          green_half_blue <- c(182, 323)
          green_blue <- c(green_with_blue, green_half_blue)
          green_little_yellow <- c(89, 109, 170)
          yellow_half_green <- c(116, 112)
          yellow_some_green <- c(53, 161)
          yellow_green <- c(green_little_yellow, yellow_half_green, yellow_some_green)
          yellow_little_blue <- c(119, 169)
          yellow_some_blue <- c(159, 219, 137, 41)
          yellow_with_blue <- c(yellow_some_blue, yellow_little_blue)
          yellow_half_blue <- c(167)
          yellow_blue <- c(yellow_with_blue, yellow_half_blue)
          yellow_ <- c(154, 215)
          if (grepl("simple", coloring)) {
            yellow_ <- c(yellow_, yellow_blue, yellow_green)
            green_ <- c(green_, green_little_yellow, green_blue)
            red_ <- c(red_, red_blue)
            list <- list("#69A019" = sort(green_),       # green
                         "#FFD700" = sort(yellow_),       # yellow
                         "#CC0000" = sort(red_),         # red
                         "#1874CD" = sort(blue_))        # blue
          } else {
            list <- list("#69A019" = sort(green_),       # green
                         "#FFD700" = sort(yellow_),       # yellow
                         "#CC0000" = sort(red_),         # red
                         "#1874CD" = sort(blue_),        # blue
                         "#8B008B" = sort(red_blue),     # lilac
                         "#008B8B" = sort(green_blue),   # bluegreen
                         "#6B8E23" = sort(yellow_blue),  # olive
                         "#C0FF3E" = sort(yellow_green)) # brigthgreen
          }
      } else if (coloring == "modules") {
        module_list <- list()
        module_list[[1]] <- c(16, 26, 30, 33, 41, 53, 58, 75, 77, 78, 81, 87, 89, 109, 112, 116, 119, 137, 141, 144, 149, 161, 169, 170, 182, 215, 323, 337, 343)
        module_list[[2]] <- c(125, 130, 150, 197, 198, 201, 212, 220, 225, 238, 239, 241, 250, 253, 262, 265, 268, 276, 277, 287, 288, 291, 292, 293, 294, 300)
        module_list[[3]] <- c(16, 26, 30, 33, 41, 53, 78, 81, 89, 109, 112, 116, 119, 120, 130, 134, 137, 141, 154, 159, 161, 167, 169, 170, 215, 219, 220, 283, 343)
        module_list[[4]] <- c(26, 58, 66, 75, 77, 78, 81, 86, 87, 89, 90, 109, 141, 144, 149, 169, 170, 174, 182, 186, 187, 323, 343, 357)
        module_list[[5]] <- c(150, 174, 197, 198, 201, 225, 238, 239, 241, 250, 253, 262, 265, 268, 276, 288, 292, 293, 294)
        module_list[[6]] <- c(33, 41, 53, 55, 58, 77, 78, 81, 87, 89, 90, 109, 112, 116, 119, 130, 137, 141, 144, 149, 159, 161, 167, 169, 170, 174, 187, 201, 212, 215, 220, 283, 323, 333)
        module_list[[7]] <- c(30, 41, 58, 75, 77, 78, 81, 86, 87, 89, 109, 119, 141, 144, 149, 170, 174, 220, 323, 343, 357)
        module_list[[8]] <- c(16, 23, 26, 30, 33, 41, 45, 53, 55, 58, 66, 75, 77, 78, 81, 87, 89, 109, 112, 116, 119, 134, 137, 141, 149, 150, 159, 161, 169, 170, 283, 323, 337, 343, 357)
        module_list[[9]] <- c(58, 66, 77, 78, 81, 86, 89, 90, 144, 149, 174, 182, 186, 187, 239, 300, 333, 357, 367)
        module_list[[10]] <- c(116, 119, 120, 130, 137, 141, 144, 154, 159, 161, 167, 169, 212, 215, 219, 283, 292, 293)
        module_list[[11]] <- c(90, 150, 197, 209, 220, 238, 276, 283, 288, 292, 293, 357, 361, 367)
        module_list[[12]] <- c(45, 66, 77, 86, 90, 141, 144, 174, 182, 186, 187, 293, 323, 357, 361, 367)
        module_list[[13]] <- c(26, 58, 75, 81, 87, 149, 182, 343)
        module_list[[14]] <- c(33, 120, 125, 150, 169, 220, 277, 287, 288, 294, 323)
        module_list[[15]] <- c(137, 209, 216, 288, 291, 293)
        module_list[[16]] <- c(78, 87)
        module_list[[17]] <- c(58, 116, 141)
        module_list[[18]] <- c(89, 149)
        module_list[[19]] <- c(77, 109, 170)
        module_list[[20]] <- c(78, 323)
        module_list[[21]] <- c(141, 323)
        module_list[[22]] <- c(137, 170)
        module_list[[23]] <- c(58, 77)
        module_list[[24]] <- c(78, 174)
        module_list[[25]] <- c(141, 144, 149)
        module_list[[26]] <- c(141, 174)
        module_list[[27]] <- c(285)
        
        list <- module_list
        
        nAttrs <- list()
        nAttrs$fillcolor <- fillcolor
      } else { # different coloring for p38g # Figure S3 automatically mixed
        # return(list())
        if (missing(counts)) {
          counts <- read.csv2("../Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
          counts <- as.matrix(counts, ncol = 4)
  
          if (!is.na(as.numeric(colors))) {
            round_categories <- as.numeric(colors)
          } else {
            round_categories <- 1
          }
          
          list <- classify_nodes(counts, round_categories = round_categories, mix = TRUE, base_colors = colnames(counts))
        }
      }
    }
    
    if (!is.null(interesting_pos)) {  # !is.null(list) ?!
      warning("No interesting positions known")
      list <- list()
    } 
  }
  
  if (for_coloring) {
    return(list)
  } else {
    # do not rename duplicate names
    return(setNames(unlist(list, use.names = FALSE), rep(names(list), lengths(list))))
    # return(unlist(list))
  }
}

colors_for_edges <- function(clustering, colors, graph) {
  edge_groups <- clustering
  edge_groups <- lapply(edge_groups, function(positions) edgeNames(graph)[positions])
  print(edge_groups)
  
  if (missing(colors)) {
    colors <- rainbow(length(edge_groups))
  }
  
  if (length(colors) > length(edge_groups)) {
    colors <- colors[1:length(edge_groups)]
  }
  
  edges_with_colors <- c()
  for (i in 1:length(edge_groups)) {
    color <- colors[i]
    color_vector <- rep(color, length(edge_groups[[i]]))
    names(color_vector) <- edge_groups[[i]]
    edges_with_colors <- c(edges_with_colors, color_vector)
  }
  # nAttrs <- list()
  # nAttrs$fillcolor <- edges_with_colors
  # return(nAttrs)
  return(edges_with_colors)
  
}

# returns the list for nAttrs$fillcolor.
# colors_for_nodes <- function(protein, position_numbering, coloring, colors, clustering) {
colors_for_nodes <- function(node_clusters, protein, coloring, colors, clustering = FALSE) {
  # if (!missing(clustering)) {
  #   pos_list <- clustering
  # } else {
  pos_list <- node_clusters
  if (!clustering) {
    if (is.null(coloring)) {
      return(NULL)
      # return(list()) # frueher
    }
    # pos_list <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
    if (missing(colors) && !is.null(names(pos_list))) {
      colors <- names(pos_list)
    }
    if (protein == "PDZ" || protein == "pdz") {
      if (colors == "auto" || is.null(colors) || colors == "") {
        colors <- names(pos_list)
      }
    } else if (protein == "p38g") { #coloring wird hier ws. gar nicht gebraucht
      if (coloring == "auto") {
        print("coloring == 'auto' is interpreted as coloring == 'FS3-pie'.")
        coloring <- "FS3-pie"
      }
      if (length(pos_list) > 0 && !coloring == "modules") {
        colors <- names(pos_list)
      }
    }
    if (length(pos_list) == 0) {
      if (colors == "" || is.null(colors) || colors == "auto") { # e.g. if coloring == "pie"
        # return(list())
        return(NULL)
      } else {
        # assume that a distinct color is given for each of the nodes
        # nAttrs <- list()
        # nAttrs$fillcolor <- colors
        # return(nAttrs)
        return(colors)
      } 
    }
  }
  
  if (missing(colors)) {
    colors <- rainbow(length(pos_list))
  }
  
  if (colors == "auto" || is.null(colors) || colors == "") {
    colors <- names(pos_list)
  }
  
  if (length(pos_list) == 0) {
    return(pos_list)
  } else {
    if (length(colors) > length(pos_list)) {
      colors <- colors[1:length(pos_list)]
    }
    
    nodes_with_colors <- c()
    for (i in 1:length(pos_list)) {
      color <- colors[i]
      color_vector <- rep(color, length(pos_list[[i]]))
      names(color_vector) <- pos_list[[i]]
      nodes_with_colors <- c(nodes_with_colors, color_vector)
    }
    # nAttrs <- list()
    # nAttrs$fillcolor <- nodes_with_colors
    # return(nAttrs)
    return(nodes_with_colors)
  }
}


plot_graph <- function(graph, fillcolor, edgecolor = NULL, drawnode, caption = "", graph_layout = "dot", protein, 
                       position_numbering, coloring, colors, outpath = "", plot_as_subgraphs = FALSE, 
                       plot_only_subgraphs = NULL, subgraphs, numerical = TRUE) {
  if (numerical) {
    if (missing(coloring) || missing(colors)) {
      plot_graph_numerical(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout, protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors, outpath = outpath, caption = caption, plot_as_subgraphs = plot_as_subgraphs, subgraphs = subgraphs)
    } else {
      for (i in 1:length(coloring)) {
        coloring_i <- coloring[i]
        if (length(colors) >= i) {
          colors_i <- colors[i]
        } else {
          colors_i <- colors[1]
        }
        if (length(plot_as_subgraphs) >= i) {
          plot_as_subgraphs_i <- plot_as_subgraphs[i]
        } else {
          plot_as_subgraphs_i <- plot_as_subgraphs[1]
        }
        if (length(graph_layout) >= i) {
          graph_layout_i <- graph_layout[i]
        } else {
          graph_layout_i <- graph_layout[1]
        }
        plot_graph_numerical(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout_i, protein = protein, position_numbering = position_numbering, coloring = coloring_i, colors = colors_i, outpath = outpath, caption = caption, plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs, subgraphs = subgraphs)
      }
    }
  }
}

plot_graph_numerical <- function(graph, fillcolor, edgecolor = NULL, drawnode, caption = "", graph_layout = "dot", protein,
                                 position_numbering, coloring, colors, outpath = "", plot_as_subgraphs = FALSE, 
                                 plot_only_subgraphs = NULL, subgraphs) {
  # if (!(missing(subgraphs))) {
  #   plot_as_subgraphs <- TRUE
  # }
  
  if (missing(fillcolor) || (missing(subgraphs) && (plot_as_subgraphs || !is.null(plot_only_subgraphs)))) {
    node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
  }
  if (missing(fillcolor)) {
    fillcolor <- colors_for_nodes(node_clusters = node_clustering, protein, coloring = coloring, colors = colors)
  }
  
  nAttrs <- list()
  nAttrs$fillcolor <- fillcolor
  
  eAttrs <- list()
  eAttrs$color <- edgecolor
  
  if (missing(subgraphs)) {
    if (plot_as_subgraphs || !is.null(plot_only_subgraphs)) {
      subgraphs <- subgraphs_from_node_clusters(node_clustering, graph, protein = protein)
    } else {
      subgraphs <- NULL
    }
  }
  
  # node shapes (pie)
  if (missing(drawnode)) {
    drawnode <- node_function_for_graph(!is.null(coloring) && (grepl("pie", coloring)))
  }
  
  if (!is.null(plot_only_subgraphs)) {
    # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
    graph <- subgraphs[[plot_only_subgraphs]]$graph
    subgraphs <- NULL
  }
    
  pc_graph <- agopen(graph, layoutType = graph_layout, nodeAttrs = nAttrs, edgeAttrs = eAttrs, name = "pc", subGList = subgraphs) # circle produziert cluster
  plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = paste(caption), subGList = subgraphs)
  
  if (!nchar(outpath) == 0) {
    if (!is.null(coloring) && !(coloring == "")) {
      postscript(paste(outpath, "_", graph_layout, "_colored-", coloring, ".ps",  sep = ""), paper="special", width = 10, height = 9)
    } else {
      postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9)
    }
    plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = caption) 
    dev.off()
  }
}


parameters_for_info_file <- function(protein, type_of_data, alpha, position_numbering, only_cols, coloring, colors, outpath) {
  par_string <- paste("protein: ", protein, ", data: ", type_of_data, ", alpha: ", alpha, sep = "")
  if (is.null(colors)) {
    colorstring <- paste("coloring:", coloring)
  } else {
    colorstring <- paste("coloring: ", coloring, ", colors:", paste(names(colors), colors, collapse="; ", sep = " - "), sep = "")
  }
  info_str <- paste(par_string, ", position_numbering: ", position_numbering, ", only_cols: ", only_cols, ", outpath: ", outpath, sep = "")
  parameters_for_info <- paste(str_replace_all(info_str, pattern = ", ", ",\n"), "\n\n", sep = "")
  return(parameters_for_info)
}

parameters_to_info_file <- function(parameters_for_info, outpath) {
  out_file <- paste(outpath, "-info.txt", sep = "")
  sink(file = out_file)
    cat(parameters_for_info)
  sink()
}

# Compute pc if necessary
get_pc <- function(pc_fun, outpath, compute_pc_anew, parameters_for_info) {
  parameters_to_info_file(parameters_for_info, outpath)
  if ((file.exists(paste(outpath, "-pc.RData", sep = ""))) && !(compute_pc_anew)) {
    filename <- paste(outpath, "-pc.RData", sep = "")
    load(filename)
    if (!exists("pc")) {
      print("The file did not contain an object of name 'pc'!")
    } else {
      print(paste("pcAlgo-object loaded from ", filename, ".", sep = ""))
      print_pc_results_to_info_file(paste(outpath, sep = ""), pc)
    }
    # file.copy(paste(outpath, "-info-pc.txt", sep = ""), paste(outpath, "-info.txt", sep = ""), overwrite = TRUE)
  } else {
    directories <- strsplit(outpath, "/")
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = "/", sep = "/")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print("Directory created.")
    }
    
    pc <- pc_fun(outpath)
    
    save(pc, file = paste(outpath, "-pc.RData", sep = ""))
    print_pc_results_to_info_file(paste(outpath, sep = ""), pc)
  }
  return(pc)
}


print_pc_results_to_info_file <- function(outpath, pc) {
  # print(paste(outpath, "-info.txt", sep = ""))
  sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
    # print("_________________")
    print("RESULT:")
    print(pc @ graph)
    print(conv_to_r(pc@graph, type_of_graph = "pdag"))
    cat("\n")
    print("COMPUTATION TIME:")
    print(proc.time())
    cat("\n")
  sink()
}

# library(dagitty)
# library(pcalg)
# library(graph)

# library(dagitty)
# library(pcalg)
# library(graph)

## gets a pcAlgo object and returns a corresponding dagitty r object
## nodename_prefix" is prepended to each node name
conv_to_r <- function(g, type_of_graph = "pdag", nodename_prefix = "") {
  ## work with graph class
  # g <- pc @ graph
  ## list of nodes
  ln <- names(nodeData(g))
  ## list of edges in format "a|b"
  le <- names(edgeData(g))
  ## construct adjacency matrix
  n <- length(ln)
  am <- matrix(0, nrow = n, ncol = n)
  ## this string contains a specification of a graph in the
  ## dot language
  ## it starts with the graphtype and the description
  ## of the actual graph is then between the curly brackets
  ## this is given as a list of nodes (with a prefix found
  ## in nodename_prefix) followed by an edge list
  ## given as [firstnode edgetype secondnode] where
  ## first- and secondnode are the same as previous nodes
  ## (and thus including the prefix)
  ## the edgetype is either -- for undirected or -> for directed edges
  ## note that it does not matter if the edges/nodes are
  ## seperated by spaces or newlines
  s <- paste(type_of_graph, "{", sep = "")
  ## print nodes
  for(node in ln) {
    s <- paste(s, nodename_prefix)
    s <- paste(s, node, sep = "")
    ## comment in to set variables to exposure
    ## s <- paste(s, "[exposure]")
  }
  ## go through all edges and mark them in adjacency matrix
  for(edge in le) {
    splt <- strsplit(edge, "\\|")
    ## this takes linear time and could be replaced by hashmaps
    from <- which(splt[[1]][1] == ln)[[1]]
    to <- which(splt[[1]][2] == ln)[[1]]
    am[[from, to]] <- 1
  }
  ## go through adjacency matrix and print edges
  ## we distinguish between bidirectional edges and unidirectional ones
  ## note that edges of the first kind will be given to dagitty as edges
  ## whose direction is unknown
  ## this might be changed later on
  for(i in 1:n) {
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        ## edge is bidirectional (second check needed to print edge only once)
        if(am[[j, i]] == 1 && i > j) {
          otp <- paste(nodename_prefix, ln[i], " -- ", nodename_prefix,ln[j], sep = "")
          s <- paste(s, otp)
        } else if(am[[j, i]] != 1) { ## edge is unidirectional
          otp <- paste(nodename_prefix, ln[i], " -> ", nodename_prefix,ln[j], sep = "")
          s <- paste(s, otp)
        }
      }
    }
  }
  s <- paste(s, "}")
  ## if you want to create dagitty object in another session, you could store s
  ## and create object later on with this string
  ## print(s)
  ## create dagitty object and return it
  dg <- dagitty(s)
  return(dg)
}

## gets a pcAlgo object and prints a string in the dagitty web format
conv_to_web <- function(g) {
  ## work with graph class
  # g <- pc @ graph
  ## list of nodes
  ln <- names(nodeData(g))
  ## list of edges in format "a|b"
  le <- names(edgeData(g))
  ## print node names followed by newline
  for(node in ln) {
    cat(paste(node, "1\n"))
  }
  cat("\n")
  ## for each node find outgoing edges
  for(node in ln) {
    ## otp is the string of this nodes edges
    otp <- node
    for(edge in le) {
      splt <- strsplit(edge, "\\|")
      ## if starting node of edge is this node, print end node
      if(splt[[1]][1] == node) {
        otp <- paste(otp, splt[[1]][2])
      }
    }
    cat(paste(otp, "\n"))
  }
}

# library(graph)
# library(dagitty)
## gets a dagitty object g and returns a graphnel object
conv_to_graphnel <- function(g) {
  ## list of nodes which will also be input for graphNEL function
  ln <- names(g)
  ## data structure which is a list of factors or something like that
  le <- dagitty::edges(g)
  ## this list of character vectors stores the edges of the graph
  eList <- list()
  ## go over all nodes
  for (i in 1:length(ln)) {
    node <- ln[i]
    ## edges of this node as character vector
    vc <- c()
    for (j in 1:length(le[[1]])) {
      ## le[[1]] stores the start nodes and le[[2]] the end nodes
      start <- as.character(le[[1]][j])
      end <- as.character(le[[2]][j])
      ## if the node we consider in this iteration is the start of an edge
      ## we append the end node to vc
      if(start == node) {
        ## this copies the vector every time and is therefore possibly slow
        vc <- c(vc, end)
      }
    }
    if(length(vc) > 0) {
      eList[[i]] <- vc
    } else {
      ## if an node has no neighbors put character(0)
      eList[[i]] <- character(0)
    }
    ## name character vector with node name
    names(eList)[i] <- node
  }
  ## create directed graphNEL object with nodes and eList
  ng <- graphNEL(nodes = ln, edgeL = eList, edgemode = "directed")
  ## you can plot the graph directly or return it
  ## plot(ng)
  return(ng)
}

plot_clusters_in_pymol <- function(node_clustering, protein, outpath, pdb_file, 
                                   label = TRUE, no_colors = FALSE, show_positions = TRUE) {

  out_file <- paste(outpath, "-", length(node_clustering), "_clusters.pml", sep = "")
  print(out_file)
  
  sink(file = out_file)
  pymol_header(protein = protein)
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
  

plot_connected_components_in_pymol <- function(protein, position_numbering, graph, outpath, label = TRUE, pdb_file, only_int_pos = FALSE, 
                                               show_int_pos = TRUE, color_int_pos = TRUE, only_color_int_pos = FALSE, coloring_for_int_pos, no_colors = FALSE, only_dist = FALSE, 
                                               show_positions = TRUE) {
  print(paste("Outpath for pymol-file:", outpath))
  connected_components <- connComp(graph)
  real_ones_ind <- which(sapply(connected_components, function(x) length(x) > 1))
  connected_components <- connected_components[real_ones_ind]
  if (only_int_pos) {
    int_pos_flat <- interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = FALSE)
    connected_components <- lapply(connected_components, function(positions) return(
      positions[sapply(positions, function(pos) return(
        as.numeric(pos) %in% int_pos_flat))]))
  }
  
  # change this to be able to print for example only the found links between interesting positions, whcih are colored accoringly
  if (only_color_int_pos) {
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
    pymol_header(protein = protein)
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


outpath_for_ida <- function(outpath, direction, relatve_effects_on_pos, option_nr, neg_effects, perturbated_position, amplification_exponent, 
                            amplification_factor, no_colors, rank_effects, effect_to_color_mode) {
    outpath <- paste0(outpath, "-total_effects_(", neg_effects, ")")

  out_file <- outpath
  if (option_nr != "") {
    out_file <- paste0(out_file, "_#", option_nr)
  }
  if (direction == "on" && relatve_effects_on_pos) {
    rel = "(rel)_"
  } else {
    rel = ""
  }
  out_file <- paste0(out_file, "_", rel, direction, "_pos_", perturbated_position)
  
  if (effect_to_color_mode == "opacity") {
    out_file <- paste0(out_file, "-opac")
  }
  if (rank_effects && !(no_colors)) {
    out_file <- paste0(out_file, "-ranked")
  } else {
    if (amplification_exponent != 1) {
      out_file <- paste0(out_file, "-ampl_exp=", amplification_exponent)
    } 
    if (amplification_factor) {
      out_file <- paste0(out_file, "-ampl_fac")
    }
  }
  if (no_colors) {
    out_file <- paste0(out_file, "-bw")
  }
  out_file <- paste0(out_file, ".pml")
  
  print(out_file)
  return(out_file)
}


scale_effects <- function(effects, rank = FALSE, amplification_factor = FALSE, neg_effects = "pos") {
  if (neg_effects == "discard") {
    effects[,1][effects[,1] < 0] <- 0
  } else if (neg_effects == "abs") {
    effects <- abs(effects)
  }
  if (rank) {
    if (neg_effects == "sep") {
      effects_pos <- as.matrix(effects[,1][effects[,1] >= 0])
      effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
      effects_neg <- -effects_neg
      effects_pos <- cbind(apply(effects_pos, 2, rank))
      effects_neg <- cbind(apply(effects_neg, 2, rank))
      
      effects_pos <- effects_pos - min(effects_pos) + 1
      effects_neg <- effects_neg - min(effects_neg) + 1
      effects_pos <- effects_pos / max(effects_pos)
      effects_neg <- effects_neg / max(effects_neg)
      
      effects_neg <- -effects_neg
      effects_ <- rbind(effects_pos, effects_neg)
      effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
      effects <- effects_
    } else {
      effects <- cbind(apply(effects, 2, rank))
      
      effects <- effects - min(effects) + 1
      effects <- effects / max(effects)
      #effects <- effects/dim(effects)[1]
    }
    amplification_exponent <- 1               # will man das? I think so
  } else {
    if (neg_effects != "sep") {
      min_eff <- sort(effects)[1]
      if (min_eff < 0) {
        effects = (effects - min_eff) / 2
      }
      sorted_effects <- sort(effects, decreasing = TRUE)
      if (amplification_factor) {
        factor <- 0.9 / sorted_effects[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
        effects <- effects * factor
        effects[,1][effects[,1] > 1] <- 1
        # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
      }
    } else {
      effects_pos <- as.matrix(effects[,1][effects[,1] >= 0])
      effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
      effects_neg <- -effects_neg
      
      sorted_effects_pos <- sort(effects_pos, decreasing = TRUE)
      if (amplification_factor) {
        factor <- 0.9 / sorted_effects_pos[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
        effects_pos <- effects_pos * factor
        effects_pos[,1][effects_pos[,1] > 1] <- 1
        # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
      }
      
      
      sorted_effects_neg <- sort(effects_neg, decreasing = TRUE)
      if (amplification_factor && dim(effects_neg)[1] > 1) {
        factor <- 0.9 / sorted_effects_neg[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
        effects_neg <- effects_neg * factor
        effects_neg[,1][effects_neg[,1] > 1] <- 1
        # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
      }
      
      
      effects_neg <- -effects_neg
      effects_ <- rbind(effects_pos, effects_neg)
      effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
      effects <- effects_
    }
  }
  return(effects)
}

# either ranked or ampl_factor (or exponent) possible
# prviously: hue_by_effect
color_by_effect <- function(effects, int_pos, color_for_other_positions = "#1874CD", mode = "mix") {
  base_color <- function(pos) {
    if (pos %in% int_pos) {
      return(names(int_pos)[which(int_pos == pos)])
    } else {
      return(color_for_other_positions)
      # return("#AAAAAA")
    }
  }
  pos_with_colors <- sapply(rownames(effects), base_color)
  pos_with_colors <- cbind(pos_with_colors, effects)
  
  # pos_with_colors <- sapply(colnames(pos_with_colors)),  function(pos) {})
  
  if (mode == "opacity") { 
    color_function <- function(vector) {
      return(adjustcolor(vector[1], alpha.f = vector[2]))
    }
  } else { # mode assumed to be a color # if (mode == "mix") {
    # with white
    color_function <- function(vector) {
      if (vector[2] >= 0) {
        return(hex(mixcolor(alpha = vector[2], color1 = hex2RGB(mode), color2 = hex2RGB(vector[1]))))
      } else {
        compl_color <- rgb(hex2RGB("#FFFFFF")@coords - hex2RGB(mode)@coords)
        # return(hex(mixcolor(alpha = -as.numeric(vector[2]), color1 = hex2RGB(vector[1]), color2 = hex2RGB(compl_color))))
        return(hex(mixcolor(alpha = -as.numeric(vector[2]), color1 = hex2RGB(compl_color), color2 = hex2RGB(vector[1]))))
      }
    }
  }
  
  # pos_with_colors <- apply(pos_with_colors, 1, function(vector) return(substr(adjustcolor(vector[1], alpha.f = vector[2]), start = 1, stop = 7)))
    pos_with_colors <- apply(pos_with_colors, 1, color_function)
  # pos_with_colors[] # sollte nur eine Position sein
  
  return(pos_with_colors)
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


pymol_header <- function(protein, pdb_file, chain = "all") {
  if (missing(pdb_file)) {
    if (protein == "GTB") {
      pdb_file <- "../../5bxc.pdb" 
    } else if (protein == "PDZ" || protein == "pdz") {
      pdb_file <- "../../1BE9.pdb"
    } else if (protein == "p38g") {
      pdb_file <- "../../1cm8.pdb"
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

# all_paths: plot all paths between from and to, instead of only the shortest
paths_between_nodes <- function(graph, from, to, all_paths = FALSE) {
  igraph <- graph_from_graphnel(graph)
  
  from_to <- paste(from, to)
  
  paths <- list()
  
  for (endpoints in from_to) {
    if (all_paths) {
      path_fun <- all_simple_paths
    } else {
      path_fun <- function(graph, from, to, ...) {
        return(shortest_paths(graph, from = from, to = to)$vpath)}
    }
    
    from <- strsplit(endpoints, " ")[[1]][1]
    to <- strsplit(endpoints, " ")[[1]][2]
    
    path <- path_fun(igraph, from = from, to = to)
    path <- lapply(path, names)
    
    paths <- c(paths, path)
    
    # path_1 <- all_simple_paths(igraph, from = "30", to = "337")
    # path_2 <- all_simple_paths(igraph, from = "197", to = "268")
    # 
    # s_path_1 <- shortest_paths(igraph, from = "30", to = "337")$vpath
    # s_path_2 <- shortest_paths(igraph, from = "197", to = "268")$vpath
  }
  return(paths)
}

plot_paths_in_pymol <- function(protein, pdb_file, graph, outpath, paths, no_colors = FALSE, label = TRUE, show_positions = TRUE) {

  out_file <- paste(outpath, "-paths.pml", sep = "")  # welche Pfade - hinzufügen
  
  sink(file = out_file)
    pymol_header(protein = protein)
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





# library(graph)

## source("~/dag/Code/graphNEL_to_tikz.R")

## g is a graphNEL object
conv_to_DIMACS <- function(g, type_of_graph = "undirected", outpath = "") {
  ## store number of nodes
  n <- length(nodes(g))
  el <- edgeL(g)
  ## m will be the number of edges which we have to count
  m <- 0
  am <- matrix(0, nrow = n, ncol = n)
  if(type_of_graph == "undirected") {
    ## filter edges with adjacency matrix
    for(i in 1:n) {
      for(j in el[[i]][[1]]) {
        if(am[[j, i]] == 0) {
          am[[i, j]] = 1
          m = m + 1
        }
      }
    }
  } else {
    for(i in 1:n) {
      for(j in el[[i]][[1]]) {
        am[[i, j]] = 1
        m =  m+1
      }
    }
  }
  ## our output string
  s <- "p tw"
  s <- paste(s, n, m)
  s <- paste(s, "\n", sep = "")
  for(i in 1:n) {
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        s <- paste(s, i, " ", j, sep="")
        s <- paste(s, "\n", sep = "")
      }
    }
  }
  if(outpath == "") {
    return(s)
  } else {
    sink(paste0(outpath, "-graph-DIMACS.txt"))
      cat(s)
    sink()
  }
}

## g is a graphNEL object and path gives the path to a textfile
## in tree decomposition format
td_to_graphNEL <- function(g, path) {
  ## n will be the number of nodes
  n = 0
  node_labels = c()
  g_labels = nodes(g)
  am = matrix(0, nrow = n, ncol = n)
  ## go through file line by line
  con = file(paste0(path, "-td.td"), "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    items = strsplit(line, " ")
    if(items[[1]][1] == "s") {
      n = strtoi(items[[1]][3])
      node_labels = vector(mode="character", length=n)
      am = matrix(0, nrow = n, ncol = n)
    } else if(items[[1]][1] == "b") {
      label = "~"
      for(i in 3:length(items[[1]])) {
        if(length(items[[1]]) < 3) break
        label = paste(label, g_labels[strtoi(items[[1]][i])], "~", sep="")
      }
      node_labels[strtoi(items[[1]][2])] = label
    } else {           
      am[[strtoi(items[[1]][1]), strtoi(items[[1]][2])]] = 1
      am[[strtoi(items[[1]][2]), strtoi(items[[1]][1])]] = 1
    }
  }
  eList <- list()
  for(i in 1:n) {
    vc <- c()
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        vc <- c(vc, node_labels[j])
      }
    }
    if(length(vc) > 0) {
      eList[[i]] <- vc
    } else {
      ## if an node has no neighbors put character(0)
      eList[[i]] <- character(0)
    }
    ## name character vector with node name
    names(eList)[i] <- node_labels[i]
  }
  ng <- graphNEL(nodes = node_labels, edgeL = eList, edgemode = "undirected")
  ## plot(ng)
  ## toFile(agopen(ng, "test"), filename="test.ps", fileType="ps")
  ## return(ng)
  close(con)
  ## graphNEL_to_tikz(ng, outpath, name)
  return(ng)
}

## g is a graphNEL object and path is the directory to which shall be written which includes a img subdirectory
graphNEL_to_tikz <- function(g, path, suffix = "") {
  dir = path
  for(i in nchar(dir):1) {
    if(substr(dir,i,i) == "/") {
      dir = substr(dir, 1, i)
      break
    }
  }
  sink(paste0(path, suffix, ".tex"))
  ## tex body of document
  cat("\\RequirePackage{luatex85,shellesc}\n")
  cat("\\documentclass[]{article}\n")
  cat("\\nofiles\n")
  cat("\n")
  cat("% LuaLaTeX stuff\n")
  cat("\\usepackage[utf8]{luainputenc}\n")
  cat("\\usepackage{luacode, luatextra, luatexbase}\n")
  cat("\n")
  cat("% TikZ\n")
  cat("\\usepackage{tikz}\n")
  cat("\\usetikzlibrary{graphs, graphdrawing}\n")
  cat("\\usegdlibrary{force, trees, circular}\n")
  cat("\n")
  cat("\\usetikzlibrary{external}\n")
  cat("\\tikzexternalize[prefix=./]\n")
  cat("\n")
  cat("\\tikzset{\n")
  cat("  graph/.style = {semithick},\n")
  cat("  vertex/.style = {draw, circle, inner sep = 0.33ex}\n")
  cat("}")
  cat("\n")
  cat("\\begin{document}\n")
  cat("\n")
  cat("\\tikz\\graph[binary tree layout] {\n")
  ## the binary tree
  for(i in 1:length(nodes(g))) {
    cat(i)
    cat("/ $\\{")
    cat(nodes(g)[i])
    cat("\\}$;\n")
  }
  for(i in 1:length(edgeL(g))) {
    for(j in 1:length(edgeL(g)[[i]][[1]])) {
      cat(i)
      cat(" -- ")
      cat(edgeL(g)[[i]][[1]][j])
      cat(";\n")
    }
  }
  cat("};")
  cat("\n")
  cat("\\end{document}\n")
  ## compiling the document
  system(paste0("cd ", dir, " && ", "lualatex --shell-escape ", path, suffix, ".tex"))
  system(paste0("rm ", path, suffix, ".auxlock"))
  system(paste0("rm ", path, suffix, "-figure0.dpth"))
  system(paste0("rm ", path, suffix, "-figure0.log"))
  system(paste0("rm ", path, suffix, "-figure0.md5"))
  system(paste0("rm ", path, suffix, ".log"))
  system(paste0("rm ", path, suffix, ".pdf"))
  system(paste0("mv ", path, suffix, "-figure0.pdf ", path, suffix, ".pdf"))
  sink()
}

