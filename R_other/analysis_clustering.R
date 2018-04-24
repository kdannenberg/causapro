graphics.off()
par(mfrow = c(1,3))
source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_S.R')

nelgraph <- results_S$pc@graph
## igraph <- graph_from_graphnel(nelgraph)
igraph <- igraph.from.graphNEL(nelgraph)

cl_im <- cluster_infomap(igraph)
cl_eb <- cluster_edge_betweenness(igraph)

# plot(cl_eb, igraph, main = "edge_betweenness")
# plot(cl_im, igraph, main = "infomap")



as_clustering <- function(c) {
#   number_of_clusters <- max(c$membership)
#   
#   if (!is.null(c$membership)) {
    return(groups(c))
#     cat("+ groups:\n")
#     hp <- function(o) {
#       head_print(o, max_lines = igraph_opt("auto.print.lines"), 
#                  omitted_footer = "+ ... omitted several groups/vertices\n", 
#       )
#     }
#     indent_print(grp, .printer = hp, .indent = "  ")
#   }
}
