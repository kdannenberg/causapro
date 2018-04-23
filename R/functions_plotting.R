library(graph)
library(igraph)

## this library of functions is under development
## the goal is to create a set of functions that plot and color graphs based on
## the known clustering of corresponding proteins
## the igraph library is used, but options for the graph library will be added later

## wrapper around plot function below
call_plot_igraph <- function(g, protein = "PDZ", position_numbering = "crystal", coloring = "", colors = "", clusters = FALSE, cluster_str = "infomap", clustering, caption = "", outpath = "", output_formats = c(), mute_all_plots = FALSE, layout_str = "layout_nicely", plot_as_subgraphs = "FALSE") {
  plot_graph_igraph(g = g, nodecolor = get_nodecolor_igraph(g, interesting_positions(protein = protein, position_numbering = position_numbering, coloring = coloring)), edgecolor = get_edgecolor_igraph(g),clusters = clusters, cluster_str = cluster_str, clustering = clustering, caption = caption, outpath = outpath, output_formats = output_formats, mute_all_plots = mute_all_plots, layout_str = layout_str, plot_as_subgraphs = plot_as_subgraphs, subgraphs = get_subgraphs_igraph(node_clusters = interesting_positions(protein = protein, position_numbering = position_numbering, for_coloring = TRUE, coloring = coloring, colors = colors), protein = protein))
}

## wrapper around the igraph plot function
## further parameters will be added later
## note that the graph that is given is still an graph object
plot_graph_igraph <- function(g, nodecolor, edgecolor, clusters, cluster_str, clustering, caption, outpath, output_formats, mute_all_plots, layout_str, plot_as_subgraphs, subgraphs) {
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

plot_text <- function(text, caption = "") {
  plot(c(0, 1), c(0, 1), xlab = "", ylab = "", bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', main = caption)
  text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
}
