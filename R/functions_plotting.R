library(graph)
library(igraph)

## this library of functions is under development
## the goal is to create a set of functions that plot and color graphs based on
## the known clustering of corresponding proteins
## the igraph library is used, but options for the graph library will be added later

## wrapper around the igraph plot function
## further parameters will be added later
## note that the graph that is given is still an graph object
plot_graph_igraph <- function(g, nodecolor, edgecolor) {
  ig = igraph.from.graphNEL(graph)
  V(ig)$color = nodecolor
  E(ig)$color = edgecolor
  plot(ig)
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
## this is merely a wrapper
get_nodecolor_igraph <- function(g) {
  
}
