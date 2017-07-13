library(graph)
library(dagitty)
library(pcalg)


protein_causal_graph <- function(data, protein, type_of_data, source_of_data, position_numbering, output_dir, filename, outpath,
                                 parameters_for_info_file, alpha, pc_solve_conflicts, pc_u2pd, pc_conservative, pc_maj_rule,
                                 caption, analysis, stages, plot_types, coloring, colors, 
                                 graph_layout = "dot", plot_as_subgraphs = plot_as_subgraphs, 
                                 plot_only_subgraphs = plot_only_subgraphs, unabbrev_r_to_info, print_r_to_console, 
                                 lines_in_abbr_of_r, compute_pc_anew, compute_localTests_anew, graph_output_formats,
                                 numerical) {
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
    return(estimate_DAG_from_numerical_data(data, alpha = alpha, outpath = outpath, 
                                            solve_conflicts = pc_solve_conflicts, u2pd = pc_u2pd, 
                                            conservative = pc_conservative, maj_rule = pc_maj_rule))
  }
  pc <- get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file)
  
  # garbage <- graphics.off()
  plot_graph(graph = pc@graph, caption = caption, protein = protein, position_numbering = position_numbering, graph_layout = graph_layout, 
             coloring = coloring, colors = colors, outpath = outpath, numerical = numerical, plot_as_subgraphs = plot_as_subgraphs, 
             plot_only_subgraphs = plot_only_subgraphs, output_formats = graph_output_formats)
  
  results <- list()
  results$pc <- pc
  
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
    positions <- intersect(positions, graphNEL@nodes)
  }
  subgraph <- subGraph(as.character(positions), graphNEL)
  return(subgraph)
}

ancestorgraph_of_interesting_positions <- function(graph_dagitty, positions = NULL, protein = NULL, position_numbering = NULL, nodename_prefix = "") {
  if (is.null(positions)) {
    positions <- paste(nodename_prefix, interesting_positions(protein, position_numbering), sep = "")
  }
  ancestor_graph <- ancestorGraph(graph_dagitty, v = positions)
  return(ancestor_graph)
}



# for_coloring -> output hierarchical (list with different sorts of interesting positions as vectors), otherwise one vector
# std-Reihenfolge: grün-gelb-rot-blau
interesting_positions <- function(protein, position_numbering, allpositions, for_coloring = FALSE, coloring = "auto", colors = "", counts) {
  if (is.null(coloring) || coloring == "none") {
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

# TODO: call plot_graphs oder so
plot_graph <- function(graph, fillcolor, edgecolor, drawnode, caption = "", graph_layout = "dot", protein, 
                       position_numbering, coloring, colors, outpath = "", plot_as_subgraphs = FALSE, 
                       plot_only_subgraphs = NULL, subgraphs, numerical = TRUE, output_formats) {
  if (numerical) {
    if (missing(coloring) || missing(colors)) {
      plot_graph_numerical(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout, protein = protein, 
                           position_numbering = position_numbering, coloring = coloring, colors = colors, outpath = outpath, caption = caption, 
                           plot_as_subgraphs = plot_as_subgraphs, subgraphs = subgraphs, output_formats = output_formats)
    } else {
      for (i in 1:(max(c(1,length(coloring))))) {
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
        ## TODO: zusammenfügen:
        ## node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
        ## fillcolor <- colors_for_nodes(node_clusters = node_clustering, protein, coloring = coloring, colors = colors)
        ## zu einer in bel. skript möglichst eindach aufrufbaren Fkt. die für protein, pos_numbering etc (colors mit default wert) fillcolors so zurückgibt, 
        ## dass man sie für diese plot.graph-Fkt nutzen kann
        # plot_graph_numerical(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout_i, protein = protein, 
        #                     position_numbering = position_numbering, coloring = coloring_i, colors = colors_i, outpath = outpath, caption = caption, 
        #                     plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs, subgraphs = subgraphs, output_formats = output_formats)
        ## can not use missing here because those are not the parameters of this function
        node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
        if (missing(subgraphs)) {
          if (plot_as_subgraphs_i || !is.null(plot_only_subgraphs)) {
            subgraphs <- subgraphs_from_node_clusters(node_clustering, graph, protein = protein)
          } else {
            subgraphs <- NULL
          }
        }
        if (missing(fillcolor)) {
          fillcolor <- colors_for_nodes(node_clusters = node_clustering, protein, coloring = coloring, colors = colors)
        }
        if (missing(drawnode)) {
          drawnode <- node_function_for_graph(!is.null(coloring) && (grepl("pie", coloring)))
        }
        if (missing(edgecolor)) {
          edgecolor <- get_eAttrs(graph)
        }
        if (!is.null(coloring) && !(coloring == "")) {
            outpath <- paste(outpath, "_", graph_layout, "_colored-", coloring, sep = "")
        } 
        plot_graph_new(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout_i, outpath = outpath, caption = caption, plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs, subgraphs = subgraphs, output_formats = output_formats)
      }
    }
  }
}

## calculate fillcolor, already done by colors_for_nodes
# TODO: default-Wert für subgraphs
plot_graph_numerical <- function(graph, fillcolor, edgecolor = NULL, drawnode, caption = "", graph_layout = "dot", protein,
                                 position_numbering, coloring, colors, outpath = "", plot_as_subgraphs = FALSE, 
                                 plot_only_subgraphs = NULL, subgraphs, output_formats = "pdf") {
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
  ## message when no subgraphs but plot_as_subgraphs true
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
  
  for (format in output_formats) {
    if (!nchar(outpath) == 0) {
      if (!is.null(coloring) && !(coloring == "")) {
        if (format == "pdf") {
          pdf(paste(outpath, "_", graph_layout, "_colored-", coloring, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, "_", graph_layout, "_colored-", coloring, ".ps",  sep = ""), paper="special", width = 10, height = 9)
        }
      } else {
        if (format == "pdf") {
          pdf(paste(outpath, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9)
        }
      }
      plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = caption) 
      dev.off()
    }
  }
}

## instead of fixing/working on the above function I think it is easier to write a new plot function the way we need it and then integrate it into the program if necessary
## what to do about drawnode, if it would be missing it was previously computed through
## node_function_for_graph which, however needs coloring
## for now I assume that this has been already computed and is NOT missing
plot_graph_new <- function(graph, fillcolor, edgecolor=NULL, drawnode, caption="", graph_layout="dot", outpath="", plot_as_subgraphs= FALSE, plot_only_subgraphs = NULL, subgraphs = NULL, output_formats = "pdf") {
  
  nAttrs <- list()
  nAttrs$fillcolor <- fillcolor

  # what happens if edgecolor is NULL
  eAttrs <- list()
  eAttrs$color <- edgecolor

  if (!is.null(plot_only_subgraphs)) {
    # graph@edgeL <- do.call(c, sapply(subgraphs, function(list) {return(list$graph@edgeL)}))
    graph <- subgraphs[[plot_only_subgraphs]]$graph
    subgraphs <- NULL
  }
  
  # this plots the graph with the given options
  pc_graph <- agopen(graph, layoutType = graph_layout, nodeAttrs = nAttrs, edgeAttrs = eAttrs, name = "pc", subGList = subgraphs) 
  plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = paste(caption), subGList = subgraphs)

  for (format in output_formats) {
    if (!nchar(outpath) == 0) {
      # if (!is.null(coloring) && !(coloring == "")) {
      #   if (format == "pdf") {
      #     pdf(paste(outpath, "_", graph_layout, "_colored-", coloring, ".pdf", sep = ""))
      #   } else if ((format == "ps") || (format == "postscript")) {
      #     postscript(paste(outpath, "_", graph_layout, "_colored-", coloring, ".ps",  sep = ""), paper="special", width = 10, height = 9)
      #   }
      # } else {
        if (format == "pdf") {
          pdf(paste(outpath, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9)
        }
      # }
      plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = caption) 
      dev.off()
    }
  }
}

# function that computes the edgecolors of a given graph
# edges with weight 2 (conflict edges) are colored red
get_eAttrs <- function(graph) {
  # list of nodes
  ln <- nodes(graph)
  n <- length(ln)
  wm <- wgtMatrix(graph)
  # init eAtrrs
  eAttrs <- list()
  eAttrs$color <- c()
  for (i in 1:n) {
    for (j in 1:n) {
      if(wm[i,j] == 2) {
        str = paste0(ln[[i]], "~", ln[[j]])
        eAttrs$color[str] <- "red"
      }
    }
  }
  return(eAttrs$color)
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


scale_effects <- function(effects, rank = FALSE, amplification_factor = FALSE, neg_effects = "pos") {
  
  amplify_with_factor <- function(effects, element_that_should_be_scaled_to = 2, 
                                  value_the_element_should_be_scaled_to = 0.9, cut_values_at = 1) {
    sorted_effects <- sort(effects, decreasing = TRUE)
    factor = value_the_element_should_be_scaled_to / sorted_effects[element_that_should_be_scaled_to]
    # factor <- 0.9 / sorted_effects_pos[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
    effects <- effects * factor
    if (is.numeric(cut_values_at)) {
      effects[,1][effects[,1] > cut_values_at] <- cut_values_at
    }
    # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
    return(effects)
  }
  
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
      if (amplification_factor) {
        effects <- amplify_with_factor(effects)
      }
    } else {
      effects_pos <- as.matrix(effects[,1][effects[,1] >= 0])
      effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
      effects_neg <- -effects_neg
      
      if (amplification_factor) {
        effects_pos <- amplify_with_factor(effects_pos)
      }
      
      if (amplification_factor && dim(effects_neg)[1] > 1) {
        effects_neg <- amplify_with_factor(effects_neg)
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







