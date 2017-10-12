library(graph)
library(igraph)

## this library of functions is under development
## the goal is to create a set of functions that plot and color graphs based on
## the known clustering of corresponding proteins
## the igraph library is used, but options for the graph library will be added later

## wrapper around plot function below

call_plot_igraph <- function(g, protein, position_numbering, coloring, colors, clusters = FALSE, clustering, caption, outpath, output_formats, mute_all_plots, layout_str, plot_as_subgraphs) {
  plot_graph_igraph(g = g, nodecolor = get_nodecolor_igraph(g, interesting_positions(protein = protein, position_numbering = position_numbering, coloring = coloring)), edgecolor = get_edgecolor_igraph(g),clusters = clusters, clustering = clustering, caption = caption, outpath = outpath, output_formats = output_formats, mute_all_plots = mute_all_plots, layout_str = layout_str, plot_as_subgraphs = plot_as_subgraphs, subgraphs = get_subgraphs_igraph(node_clusters = interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors), protein = protein))
}

## wrapper around the igraph plot function
## further parameters will be added later
## note that the graph that is given is still an graph object
plot_graph_igraph <- function(g, nodecolor, edgecolor, clusters, clustering, caption, outpath, output_formats, mute_all_plots, layout_str, plot_as_subgraphs, subgraphs) {
  ## par(mfrow = c(2,2))
  ig = igraph.from.graphNEL(g)
  V(ig)$color = nodecolor
  E(ig)$color = edgecolor
  if(clusters) {
    E(ig)$color[crossing(clustering, ig)] = "blue"
    E(ig)$color[intersect(which(crossing(clustering,ig)), which(E(ig)$color == "red"))] = "purple"
    if(!mute_all_plots) {
      plot(clustering, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.8, main = caption)
    }
    for(format in output_formats) {
      if (!nchar(outpath) == 0) {
        if (format == "pdf") {
          pdf(paste(outpath, clustering, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
        }
        plot(clustering, ig, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.8, main = caption)
        dev.off()
      }
    }
  } else {
    layout_fct = get(layout_str)
    layout = layout_fct(ig)
    mem = c()
    mem[nodes(g)] = 1
    mem[subgraphs[[1]]] = 2
    cl = make_clusters(ig, mem)
    if(!mute_all_plots) {
      if(layout_str == "layout_with_sugiyama") {
        plot(layout$extd_graph, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
      } else {
        if(plot_as_subgraphs) {
          plot(cl, ig, col = V(ig)$color, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)
        } else {
          plot(ig, layout = layout, edge.arrow.size=0.1, vertex.size=8, edge.width=0.8, main = caption)          
        }
      }
    }
    for(format in output_formats) {
        if (!nchar(outpath) == 0) {
          if (format == "pdf") {
            pdf(paste(outpath, ".pdf", sep = ""))
          } else if ((format == "ps") || (format == "postscript")) {
            postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
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

get_subgraphs_igraph <- function(node_clusters, protein) {
  if(protein == "PDZ") {
    node_clusters[[1]] <- paste(c(node_clusters[[1]], node_clusters[[2]]))
    node_clusters[[2]] <- NULL
  }
  return(node_clusters)
}
