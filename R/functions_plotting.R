library(graph)
library(igraph)


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

# TODO: call plot_graph_by_coloring oder so
# TODO: add par(mfrow=...), so everything fits on one plane
plot_graph <- function(graph, fillcolor, edgecolor, drawnode, caption = "", graph_layout = "dot", protein,
                       position_numbering, coloring, colors, outpath = "", plot_as_subgraphs = FALSE,
                       plot_only_subgraphs = NULL, subgraphs, numerical = TRUE, output_formats, mute_all_plots = FALSE) {
  ## numerical serves no purpose right now, but it might in the future
  if (numerical) {
    if (missing(coloring) || missing(colors)) {
      plot_structure(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout, outpath = outpath, caption = caption,
                           plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs, subgraphs = subgraphs, output_formats = output_formats, mute_all_plots = mute_all_plots)
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
        ## zu einer in bel. skript möglichst einfach aufrufbaren Fkt. die für protein, pos_numbering etc (colors mit default wert) fillcolors so zurückgibt,
        ## dass man sie für diese plot.graph-Fkt nutzen kann
        # plot_graph_numerical(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, graph_layout = graph_layout_i, protein = protein,
        #                     position_numbering = position_numbering, coloring = coloring_i, colors = colors_i, outpath = outpath, caption = caption,
        #                     plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs, subgraphs = subgraphs, output_formats = output_formats)
        ## can not use missing here because those are not the parameters of this function
        node_clustering <- interesting_positions(protein, position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors)
        # print(node_clustering)
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
        plot_structure(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode,
                       graph_layout = graph_layout_i, outpath = outpath, caption = caption,
                       plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs,
                       subgraphs = subgraphs, output_formats = output_formats, mute_all_plots = mute_all_plots)
      }
    }
  }
}



#' Plots a given structure using the Rgraphviz package.
#'
#' @param graph The DAG to plot given as graph object (as specified in the R graph package).
#' @param fillcolor List of colors indicating the node coloring. The name of a list element specifies the node (e.g. character string)
#' @param edgecolor List of colors indicating the edge coloring. The name of a list element specifies the edge (the edge names are e.g. given by the edgeNames(graph) function of the R graph package).
#' @param drawnode Draws the node as a pie. Parameter should be a list of drawing functions (details can be found in the Rgraphviz package documentation).
#' @param caption A character string, the caption of the plot.
#' @param graph_layout A character string. All graph layouts from Rgraphviz are available.
#' @param outpath A character string, the path to the location where the plot should be stored. If an empty string is passed, the plot will not be stored.
#' @param plot_as_subgraphs A boolean, indicating if the plot should be divided into multiple subgraphs.
#' @param plot_only_subgraphs Possibly a list of indices, indicating which subgraphs should be plotted? <- TODO
#' @param subgraphs TODO
#' @param output_formats A character string, the desired output format. All standard R plotting options are available.
#' @param mute_all_plots A boolean, if TRUE the structure will not be plotted.
#' @return No return value.

plot_structure <- function(graph, fillcolor=NULL, edgecolor=NULL, drawnode=drawAgNode, caption="", graph_layout="dot", outpath="",
                           plot_as_subgraphs= FALSE, plot_only_subgraphs = NULL, subgraphs = NULL,
                           output_formats = "pdf", mute_all_plots = FALSE) {

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

  ## filter drawnode
  if(length(drawnode) > 1) {
    drawnode = drawnode[nodes(graph)]
  }

  pc_graph <- agopen(graph, layoutType = graph_layout, nodeAttrs = nAttrs, edgeAttrs = eAttrs, name = "pc", subGList = subgraphs)

  # this plots the graph with the given options
  if (!mute_all_plots) {
    plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = paste(caption), subGList = subgraphs)
  }

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
      } else if (format == "svg") {
        svg(paste0(outpath, ".svg"))
      } else {
        warning(paste("Unknown format:", format))
      }
      # }
      plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = caption)
      dev.off()
    }
  }
}

# function that computes the edgecolors of a given graph
# edges with weight 2 (conflict edges) are colored red
get_eAttrs <- function(graph, igraph=FALSE) {
  # list of nodes
  ln <- nodes(graph)
  n <- length(ln)
  wm <- wgtMatrix(graph)
  # init eAtrrs
  eAttrs <- list()
  eAttrs$color <- c()
  if (n == 0) {
    return(eAttrs$color)
  }
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




## this library of functions is under development
## the goal is to create a set of functions that plot and color graphs based on
## the known clustering of corresponding proteins
## the igraph library is used, but options for the graph library will be added later

## wrapper around plot function below
call_plot_igraph <- function(g, protein = "PDZ", position_numbering = "crystal", coloring = "", colors = "", clusters = FALSE, cluster_str = "infomap", clustering, caption = "", outpath = "", output_formats = c(), mute_all_plots = FALSE, layout_str = "layout_nicely", plot_as_subgraphs = "FALSE") {
  plot_structure_igraph(g = g, nodecolor = get_nodecolor_igraph(g, interesting_positions(protein = protein, position_numbering = position_numbering, coloring = coloring)), edgecolor = get_edgecolor_igraph(g),clusters = clusters, cluster_str = cluster_str, clustering = clustering, caption = caption, outpath = outpath, output_formats = output_formats, mute_all_plots = mute_all_plots, layout_str = layout_str, plot_as_subgraphs = plot_as_subgraphs, subgraphs = get_subgraphs_igraph(node_clusters = interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors), protein = protein))
}

#' Plots a given structure using the igraph package. Note that the graph is still passed as an graph object (not an igraph object).
#'
#' @param g A graph object as specified in the R graph package.
#' @param nodecolor List of colors indexed by the node.
#' @param edgecolor List of colors indexed by edge id.
#' @param clusters Boolean indicating if graph should be clustered.
#' @param clustering Character string choosing a clustering from the ones available in the igraph package.
#' @param caption A character string, the caption of the plot.
#' @param outpath A character string, the path to the location where the plot should be stored. If an empty string is passed, the plot will not be stored.
#' @param output_formats A character string, the desired output format. All standard R plotting options are available.
#' @param mute_all_plots A boolean, if TRUE the structure will not be plotted.
#' @param layout_str A character string. All graph layouts from igraph are available.
#' @param plot_as_subgraphs A boolean, indicating if the plot should be divided into multiple subgraphs.
#' @param subgraphs TODO
#' @return No return value.
plot_structure_igraph <- function(g, nodecolor, edgecolor, clusters, cluster_str, clustering, caption, outpath, output_formats, mute_all_plots, layout_str, plot_as_subgraphs, subgraphs) {
  ## par(mfrow = c(2,2))
  # kernel
  ig = igraph.from.graphNEL(g)
  V(ig)$color = nodecolor
  E(ig)$color = edgecolor
  ## we plot a given clustering
  if(clusters) {
    ## edges between clusters are blue, conflict edges red. If both categories apply the edge is colored purple
    E(ig)$color[crossing(clustering, ig)] = "blue"
    E(ig)$color[intersect(which(crossing(clustering,ig)), which(E(ig)$color == "red"))] = "purple"
    if(!mute_all_plots) {
      plot(clustering, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.8, main = paste(caption, cluster_str))
    }
    for(format in output_formats) {
      if (!nchar(outpath) == 0) {
        if (format == "pdf") {
          pdf(paste(outpath, "_", cluster_str, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, "_", cluster_str, ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
        } else if(format == "svg") {
          svg(paste(outpath, "_", cluster_str, ".svg", sep = ""))
        }
        plot(clustering, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.8, main = paste(caption, cluster_str))
        dev.off()
      }
    }
  } else {
    ## plotting a given layout
    layout_fct = get(layout_str)
    layout = layout_fct(ig)
    ## I make my own communities, only way I found to emulate subgraph (from Rgraphviz)
    mem = c()
    mem[nodes(g)] = 1
    mem[subgraphs[[1]]] = 2
    cl = make_clusters(ig, mem)
    if(!mute_all_plots) {
      ## note that sugiyama does not work too well with CPDAGs (it needs actual DAGs)
      ## therefore it will add new nodes to "solve" the cycle etc
      if(layout_str == "layout_with_sugiyama") {
        plot(layout$extd_graph, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
      } else {
        if(plot_as_subgraphs) {
          plot(cl, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = paste(caption, layout_str, "as_subgraphs"))
        } else {
          plot(ig, layout = layout, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
        }
      }
    }
    for(format in output_formats) {
      print(output_formats)
        if (!nchar(outpath) == 0) {
          if (format == "pdf") {
            pdf(paste(outpath, "_", layout_str, "_as_sg", ".pdf", sep = ""))
          } else if ((format == "ps") || (format == "postscript")) {
            postscript(paste(outpath, "_", layout_str, "_as_sg", ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
          }  else if(format == "svg") {
            svg(paste(outpath, "_", layout_str, "_as_sg", ".svg", sep = ""))
          }
          if(layout_str == "layout_with_sugiyama") {
            plot(layout$extd_graph, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
          } else {
            if(plot_as_subgraphs) {
              plot(cl, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
            } else {
              plot(ig, layout = layout, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
            }
          }
          dev.off()
        }
    }
  }
}

## calculates the colors of the edges
## all conflict edges are colored red
## this function gets an igraph object
get_edgecolor_igraph <- function(g) {
  ig = igraph.from.graphNEL(g)
  am = wgtMatrix(g)
  edgecolor = rep("black", numEdges(g))
  edgecolor[get.edge.ids(ig, c(t(which(wgtMatrix(g) == 2, arr.ind = numNodes(g)))))] = "red"
  return(edgecolor)
}

## calculates the colors of the nodes
## the colors are given by the interesting positions function
## (the int_pos parameter)
get_nodecolor_igraph <- function(g, int_pos) {
  nodecolor <- c()
  nodecolor[nodes(g)] <- "white"
  nodecolor[paste(int_pos)] <- names(int_pos)
  return(nodecolor)
}

## a new cleaner function for getting subgraphs
## currently only working with PDZ
get_subgraphs_igraph <- function(node_clusters, protein) {
  if(protein == "PDZ") {
    node_clusters[[1]] <- paste(c(node_clusters[[1]], node_clusters[[2]]))
    node_clusters[[2]] <- NULL
  }
  return(node_clusters)
}

plot_infeasible <- function(caption = "") {
    # plot(c(0, 1), c(0, 1), xlab = "", ylab = "", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = caption)
    # text(x = 0.5, y = 0.5, "infeasible",
    #      cex = 1.6, col = "black")
  plot_text(text = "infeasible", caption = caption)
}

plot_text <- function(text, caption = "", ...) {
  plot(c(0, 1), c(0, 1), xlab = "", ylab = "", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = caption, ...)
  text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
}

# computes a vector of the same length as pos that for each position either contains
# color_for_other_positions, or, if the position is in int_pos, its name in int_pos;
# this vector can e.g. be used as the color vector in plots
int_pos_to_color_vector <- function(pos, int_pos, color_for_other_positions = "#000000") {
  base_color <- function(pos) {
    if (pos %in% int_pos) {
      return(names(int_pos)[which(int_pos == pos)])
    } else {
      return(color_for_other_positions)
      # return("#AAAAAA")
    }
  }
  return(sapply(pos, base_color))
}

# either ranked or ampl_factor (or exponent) possible
# prviously: hue_by_effect
# color_by_effect <- function(effects, int_pos, color_for_other_positions = "#1874CD", mode = "mix") {
color_by_effect <- function(effects, int_pos, color_for_other_positions = "#1874CD", mode = "#FFFFFF") {
  # TODO Marcel: Geht das auch eleganter, so dass effects Zeilen-, Spaltenmatrix oder Vektor sein kann?
  pos <- rownames(effects)
  if (is.null(pos)) {
    if (is.vector(effects)) {
      pos <- names(effects)
    } else {
      pos <- names(effects)
    }
  }

  pos_with_colors <- int_pos_to_color_vector(pos = pos, int_pos = int_pos, color_for_other_positions = "#1874CD")
  # pos_with_colors <- sapply(rownames(effects), base_color)
  pos_with_colors <- cbind(pos_with_colors, effects)

  # pos_with_colors <- sapply(colnames(pos_with_colors)),  function(pos) {})

  if (mode == "opacity") {
    color_function <- function(vector) {
      if (is.na(vector[2])) {
        vector <- c("#FF6347", 1)
      }
      return(adjustcolor(vector[1], alpha.f = vector[2]))
    }
  } else { # mode assumed to be a color # if (mode == "mix") {
    # with white
    color_function <- function(vector) {
      if (is.na(vector[2])) {
        vector <- c("#FF6347", 1)
      }
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
  igraph <- igraph.from.graphNEL(graph)

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


nonsingular_connected_components <- function(graph, remove_non_assigned = TRUE){
  connected_components <- connComp(graph)
  if (remove_non_assigned) {
    connected_components <- lapply(connected_components, function (list) {return(list[!grepl('\\[|\\]', list)])})
  }
  real_ones_ind <- which(sapply(connected_components, function(x) length(x) > 1))
  connected_components <- connected_components[real_ones_ind]


  return(connected_components)
}
