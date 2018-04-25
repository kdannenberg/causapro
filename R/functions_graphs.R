## this function gets a graph and returns graph which includes only non-isolated nodes
## it also prints a list of the isolated nodes
kernelize_graph <- function(graph, ret_list = FALSE) {
  n = numNodes(graph)
  vis = logical(n)
  ln = character()
  isolated_nodes <- list()
  cnt = 1
  for(u in nodes(graph)) {
    # TODO Marcel : the following line fails when the graph has only undirected edges:
    if(is.atomic(graph::degree(graph, u))) {
      if(graph::degree(graph, u) == 0) {
        ## create list
        if(ret_list) {
          isolated_nodes[[cnt]] = u
          cnt = cnt + 1
        }
      } else {
        ln <- c(ln, u)
      }
    } else {
      if(graph::degree(graph, u)$inDegree == 0 && graph::degree(graph, u)$outDegree == 0) {
        ## create list
        if(ret_list) {
          isolated_nodes[[cnt]] = u
          cnt = cnt + 1
        }
      } else {
        ln <- c(ln, u)
      }
    }
  }
  eList <- list()
  # for(i in ln) {
  #     if(i %in% ln) {
  #         vc <- c()
  #         for(j in edgeL(graph)[[i]][[1]]) {
  #             if(length(j) == 0) next
  #             if(nodes(graph)[j] %in% ln) {
  #                 vc <- c(vc, nodes(graph)[j])
  #             }
  #         }
  #     }
  #     if(length(vc) > 0) {
  #         eList[[i]] <- vc
  #     } else {
  #         eList[[i]] <- character(0)
  #     }
  # }
  lam = wgtMatrix(graph)
  am = matrix(0, length(ln), length(ln))
  rownames(am) <- ln
  colnames(am) <- ln
  for(i in ln) {
    for(j in ln) {
      am[i,j] = lam[i,j]
    }
  }
  l <- list()
  # l[["graph"]] <- graphNEL(nodes = ln, edgeL = eList, edgemode = "directed")
  l[["graph"]] <- as(t(am), "graphNEL")
  l[["isolated_nodes"]] <- isolated_nodes
  ## do something with the graph
  ## plot(g)
  if (length(l$graph@nodes) == 0) {
    l$graph <- graph    ## if new graph is empty, return old graph, to avoid trouble
  }
  if(ret_list) {
    return(l)
  } else {
    return(l[["graph"]])
  }
}

## I assume that the adjmatrix is given where conflict edges have weight 2
## this corresponds to the wgtMatrix function of the pcalg package
## I enumerate the graphs as true adjacency matrices with a 1 indicating an edge a 0 no edge and there is no other value
enumerate_graphs_rec <- function(am, n, i = 1, j = 1) {
  cat(i,j)
  cat("\n")
  print(Cstack_info())
  if (i == n+1) {
    if (j == n) {
      ## do stuff with graph
      print(am)
      return()
    } else {
      i = 1
      j = j+1
    }
  }
  ## I assume if am[i,j] = 2 then also am[j,i] = 2, but that should hold, right?
  if(am[i,j] == 2) {
    am[i,j] = 1
    am[j,i] = 0
    enumerate_graphs_rec(am, n, i+1, j)
    am[i,j] = 0
    am[j,i] = 1
    enumerate_graphs_rec(am, n, i+1, j)
    am[i,j] = 2
    am[j,i] = 2
  } else {
    enumerate_graphs_rec(am, n, i+1, j)
  }
}

dfs <- function(am, n, i, color) {
  color[i] = 1
  res = TRUE
  for(j in 1:n) {
    if(am[i,j] == 1) {
      if(color[j] == 0) {
        res = res && dfs(am, n, j, color)
      } else if(color[j] == 1) {
        res = FALSE
      }
    }
  }
  return(res)
}


## I realised now that this is pretty useless
dag_check <- function(am) {
  n = dim(am)[1]
  color <- vector(mode="double", length=n)
  res = TRUE
  for(i in 1:n) {
    if(color[i] == 0) {
      res = res && dfs(am, n, i, color)
    }
  }
  return(res)
}

direct_unambigous_undirected_edges <- function(am) {
  n = dim(am)[1]
  pos <- c()
  for(i in 1:n) {
    for(j in (i+1):n) {
      if(j > n) next
      if(am[i,j] == 2) {
        pos <- c(pos, (i-1)*n+(j-1))
      }
    }
  }
  while(TRUE) {
    oldpos = length(pos)
    am = solve_conflicts(am, pos)
    if(is.null(am)) {
      return(NULL)
    }
    pos <- c()
    for(i in 1:n) {
      for(j in (i+1):n) {
        if(j > n) next
        if(am[i,j] == 2) {
          pos <- c(pos, (i-1)*n+(j-1))
        }
      }
    }
    if(length(pos) == oldpos) break
  }
  # for(i in 1:n) {
  #   for(j in (i+1):n) {
  #     if(j > n) next
  #     if(am[i,j] == 1 && am[j,i] == 1) {
  #       case = 3 # 3 = both ways possible, 2 = i -> j possible, 1 = j -> i possible, 0 no direction possible
  #       for(a in 1:n) {
  #         if(am[a,i] == 1) {
  #           case = bitwAnd(case, 2)
  #         }
  #         if(am[a,j] == 1) {
  #           case = bitwAnd(case, 1)
  #         }
  #       }
  #       if(case == 0) {
  #         return(NULL)
  #       } else if(case == 1) {
  #         am[i,j] = 0
  #       } else if(case == 2) {
  #         am[j,i] = 0
  #       }
  #     }
  #   }
  # }

  return(am)
}

solve_conflicts <- function(pdag) {
  # directed_edges <- c()
  # poss = TRUE
  # n = dim(am)[1]
  # for(i in 1:length(pos)) {
  #   k = pos[i]
  #   c = k %/% n + 1
  #   b = k %% n + 1
  #   if(am[b,c] == 1) {
  #     tmp = c
  #     c = b
  #     b = tmp
  #   }
  #   for(a in 1:n) {
  #     if(am[b,a] == 1 && am[a,b] == 1 && !(am[c,a] == 1 || am[a,c] == 1)) {
  #       if(((a-1)*n + (b-1)) %in% directed_edges) {
  #         poss = FALSE
  #       } else {
  #         directed_edges <- c(directed_edges, (b-1)*n + (a-1))
  #       }
  #     }
  #   }
  # }
  # for(i in directed_edges) {
  #   b = i %/% n + 1
  #   a = i %% n + 1
  #   am[a,b] = 0
  #   am[b,a] = 1
  # }
  # if(poss) {
  #   return(am)
  # } else {
  #   return(NULL)
  #   # return(matrix(0, n, n))
  # }
  n = dim(pdag)[1]
  old_pdag <- matrix(0, n, n)
  while (!all(old_pdag == pdag)) {
    old_pdag <- pdag
    ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((pdag[b, ] == 1 & pdag[, b] ==
                       1) & (pdag[a, ] == 0 & pdag[, a] == 0))
      if (length(indC) > 0) {
        pdag[b, indC] <- 1
        pdag[indC, b] <- 0
      }
    }
  }
  return(pdag)
}

remove_dummies <- function(graphs) {
  remove <- which(sapply(graphs, is.null))
  if (length(remove) > 0) {
    # TODO Marcel: nicht-quadratsich
    graphs <- graphs[-remove]
  }
  return(graphs)
}

enumerate_graphs <- function(graph, direct_adjacent_undirected_edges = TRUE, direct_unambig_undirected_edges = FALSE) {
  am = t(wgtMatrix(graph))
  n = dim(am)[1]
  pos <- c()
  for(i in 1:n) {
    for(j in (i+1):n) {
      if(j > n) next
      if(am[i,j] == 2) {
        pos <- c(pos, (i-1)*n+(j-1))
      }
    }
  }

  if (length(pos) < 27) {
    graphs <- vector("list", 2^length(pos))  # fixed length; non-quadratic running time
  } else {
    stop("Not all graphs enumerable (for all combinations of directions for the conflict edges)") # 2^27 * 8 Byte = 1GB
  }

  for(i in 0:(2^length(pos)-1)) {
    for(j in 1:length(pos)) {
      k = pos[j]
      x = k %/% n + 1
      y = k %% n + 1
      if(bitwAnd(2^(j-1), i) > 0) {
        am[x,y] = 1
        am[y,x] = 0
      } else {
        am[y,x] = 1
        am[x,y] = 0
      }
    }
    ## do stuff with graph
    #   m <- matrix(0, n, n)
    #   for(index in 1:n) {
    #     for(index2 in 1:n) {
    #       m[index,index2] = am[index,index2]
    #     }
    #   }

    #   if(direct_adjacent_undirected_edges) {
    #     m <- solve_conflicts(m, pos)
    #   }
    #   if(direct_unambig_undirected_edges && !is.null(m)) {
    #     m <- direct_unambigous_undirected_edges(m)
    #   }
    #   if (is.null(m)) {
    #     graphs[[i+1]] <- NULL              # macht das die Laufzeit wieder quadratisch?
    #   } else {
    #     graphs[[i+1]] <- as(m, "graphNEL")
    #   }
    # }
    # graphs <- remove_dummies(graphs)
    m <- matrix(0, n, n)
    for(index in 1:n) {
      for(index2 in 1:n) {
        m[index,index2] = am[index,index2]
      }
    }

    colnames(m) <- colnames(am)
    rownames(m) <- rownames(am)
    graphs[[i+1]] <- as(solve_conflicts(m), "graphNEL")
  }
  return(graphs)
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

#' Edges Types in Graph
#' calls edge_information with full adj_matrix if interesting_positions is missing or the rows and columns of that vector if it is given
conflict_edges <- function(graph, print = FALSE, interesting_positions) {
  gm <- wgtMatrix(graph)
  if(missing(interesting_positions)) {
    return(edge_information(gm, print = print))
  } else {
    igm <- gm[as.character(interesting_positions), as.character(interesting_positions)]
    return(edge_information(igm, print = print))
  }
}

### relevant measures: conflict_edges/(bidirected_edges + unidirected_edges)
### or conflict_edges/(2 * bidirected_edges + unidirected_edges) (number of directed edges)
#' Different Edge Types in a Graph
#'
#' Get the number of edges of different types (conflict, bidirected, directed)
#'
#' @param adj_m the adjacency matrix of the graph that is to be considered (as returned by wgtMartix(...))
#' @param print should the result also be printed?
#'
#' @return A list with names "conflict", "directed" and "undirected" and the respective numbers of such edges.
edge_information <- function(adj_m, print = FALSE) {

  # m <- adj_m[apply(adj_m!=0, 1, any), , drop=FALSE]
  # adj_m_no_zero_rows_cols <- m[apply(m!=0, 2, any), , drop=FALSE]
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m, 1)
  adj_m_no_zero_rows_cols <- remove_zero_rows_or_columns(adj_m_no_zero_rows_cols, 2)
  # print(adj_m_no_zero_rows_cols)
  n_conflic_edges <- length(which(adj_m == 2)) / 2
  # = length(which(unlist(edgeData(graph)) == 2))
  # this works as there are no unidirected conflict edges
  n_unidirected_edges <- length(which(adj_m != t(adj_m))) / 2

  n_bidirected_edges <- (length(which(adj_m == 1)) - n_unidirected_edges) / 2

  if (print) {
    print(paste("<-!-> : ", n_conflic_edges))
    print(paste("<---> : ", n_bidirected_edges))
    print(paste("----> : ", n_unidirected_edges))
  }

  return(list(conflict = n_conflic_edges, directed = n_unidirected_edges, undirected = n_bidirected_edges))
}
