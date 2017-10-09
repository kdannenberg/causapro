library(graph)
library(igraph)

## this library of functions is under development
## the goal is to create a set of functions that plot and color graphs based on
## the known clustering of corresponding proteins
## the igraph library is used, but options for the graph library will be added later

## wrapper around plot function below

call_plot_igraph <- function(g, protein, position_numbering, coloring, colors, clusters = FALSE, caption, outpath, output_formats, mute_all_plots) {
  plot_graph_igraph(g = g, nodecolor = get_nodecolor_igraph(g, interesting_positions(protein = protein, position_numbering = position_numbering, coloring = coloring)), edgecolor = get_edgecolor_igraph(g),clusters = clusters, caption = caption, outpath = outpath, output_formats = output_formats, mute_all_plots = mute_all_plots)
}

## wrapper around the igraph plot function
## further parameters will be added later
## note that the graph that is given is still an graph object
plot_graph_igraph <- function(g, nodecolor, edgecolor, clusters, caption, outpath, output_formats, mute_all_plots) {
  ## par(mfrow = c(2,2))
  ig = igraph.from.graphNEL(g)
  V(ig)$color = nodecolor
  E(ig)$color = edgecolor
  if(clusters) {
    cluster_methods = c("edge_betweenness", "infomap")
    for(clustering in cluster_methods) {
      cluster_fct <- get(paste0("cluster_", clustering))
      cl <- cluster_fct(ig)
      E(ig)$color[crossing(cl, ig)] = "blue"
      if(!mute_all_plots) {
        plot(cl, ig, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.7, main = paste0(caption, clustering))
      }
      for(format in output_formats) {
        if (!nchar(outpath) == 0) {
          if (format == "pdf") {
            pdf(paste(outpath, clustering, ".pdf", sep = ""))
          } else if ((format == "ps") || (format == "postscript")) {
            postscript(paste(outpath, clustering, ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
          }
          plot(cl, ig, edge.color = E(ig)$color, edge.arrow.size=0.1, vertex.size=10, edge.width=0.7, main = paste0(caption, clustering))
          dev.off()
        }
      }
    }
  } else {
    if(!mute_all_plots) {
      plot(ig, edge.arrow.size=0.1, vertex.size=8, edge.width=0.7, main = caption)
    }
    for(format in output_formats) {
        if (!nchar(outpath) == 0) {
          if (format == "pdf") {
            pdf(paste(outpath, ".pdf", sep = ""))
          } else if ((format == "ps") || (format == "postscript")) {
            postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
          }
          plot(ig, edge.arrow.size=0.1, vertex.size=8, edge.width=0.7, main = caption)
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
