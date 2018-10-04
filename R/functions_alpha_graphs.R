collect_graphs_over_alpha <- function(alphas, pc_fun) {
  list_of_graphs <- list()
  names(list_of_graphs) <- alphas
  for (i in 1:length(alphas)) {
    result <- pc_fun(alpha = alphas[i], mute_all_plots = TRUE)
    list_of_graphs[[i]] <- result$pc@graph
  }

  return(list_of_graphs)
}

first_alpha <- function(alpha_graphs, n = length(alpha_graphs[[1]]@nodes)) {
  alpha_igraphs <- lapply(alpha_graphs, igraph.from.graphNEL)
  # complete_graph <- make_full_citation_graph(n, directed = FALSE)
  supergraph <- do.call(union, lapply(alpha_igraphs, as.directed))
  super_adj_mat <- as_adjacency_matrix(supergraph)

  adj_mats <- array(dim = c(dim(super_adj_mat), length(alpha_igraphs)))
  names(adj_mats)[3] <- "alphas"
  dimnames(adj_mats)[[3]] <- names(alpha_igraphs)
  for (i_alpha in 1:length(alpha_graphs)) {
    adj_mats[,,i_alpha] <- as.matrix(as_adjacency_matrix(alpha_igraphs[[i_alpha]]))
  }
  smallest_non_zero_alpha <- function(vector) {
    if (!is.null(names(vector))) {
      return(names(vector)[min(which(vector != 0))])
    } else {
      warning("Vector not named in smallest_non_zero_alpha! Returning index instead of alpha.")
      return(min(which(vector != 0)))
    }
  }
  # debug(smallest_non_zero_alpha)
  first_alphas <- apply(adj_mats, c(1,2), smallest_non_zero_alpha)

  # Jetzt an die Kanten schreiben usw.
  # E(supergraph)$labels <- first_alphas
}

# Returns the alpha from the names of alpha_graphs, for which the pertaining graph contains an edge
#' between the nodes labelled u and v (if characters) or the u'th and v'th node (if integers)
# first_alpha_with_edge(u, v, alpha_graphs) {
#   for (alpha %in% names(alpha_graphs)) {
#
#   }
#
# }

# plot(graph, edge.lty=0, edge.arrow.size=0, ...)
