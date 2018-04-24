library(graph)
library(dagitty)
library(pcalg)

## gets an alignment as a matrix of characters
## returns a binary matrix of the same size with the ones corresponding to
## the columnwise most frequent proteins
alignment_to_binary_matrix <- function(alignment) {
  ## n is the number of columns
  n <- dim(alignment)[2]
  maj <- function(t) {
    return(names(which.max(table(t))))
  }
  ## storing the columnwise most frequent proteins in a vector
  most_frequent_aminoacids <- apply(alignment, 2, maj)
  ## return binary matrix
  return(t(1*apply(alignment, 1, `==`, most_frequent_aminoacids)))
}

set_parameters <- function(FUN, parameters) {
  # return(function(...) return(do.call(FUN, parameters, ...)))
  # return(function(...) {return(do.call(FUN, list(parameters, list(...))))})  # previously: return(do.call(FUN, c(parameters, list(...)))) (problematic when arguemts are matrices)
  return(function(...) {return(do.call(FUN, c(parameters, list(...))))})
}

function_set_parameters <- set_parameters

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
# the for_coloring result of this function can be converted to the other one by conv_for_coloring
# (which is what happens internally if for_coloring is not set), so this result is more useful in general
# std-Reihenfolge: grün-gelb-rot-blau
interesting_positions <- function(protein, position_numbering = "crystal", for_coloring = FALSE, coloring = "", colors = "", counts) {
  if (!is.null(names(protein))) {
    # if ("int_pos" %in% names(protein$general)) {
    #   if (!for_coloring) {
    #     return(protein$general$int_pos)
    #   } else {
    #     return(conv_for_coloring(protein$general$int_pos))
    #     # TODO: function for converting between for_coloring and other
    #   }
    # }
    if ("caption" %in% names(protein$summary)) {
      if (grepl("bin_approx", protein$summary$caption)) {
        position_numbering <- "alignment"
      }
      protein <- strsplit(protein$summary$caption, "_")[[1]][1]
    } else if ("outpath" %in% names(protein$summary)) {
      if (grepl("bin_approx", protein$summary$outpath)) {
        position_numbering <- "alignment"
      }
      protein <- strsplit(protein$summary$outpath, "/")[[1]][2]
    }
  }
  if (is.null(coloring) || coloring == "none") {
    list <- list()
  } else if (grepl("pymol", tolower(coloring))) {
    # something like:
    # rainbow_colors <- rainbow(length(which(names(node_clustering) == "")))
    # names(node_clustering)[names(node_clustering) == ""] <- rainbow_colors
    # colors <- names(node_clustering)
  } else {
    interesting_pos <- NULL     ## list?!
    if ((protein == "pdz") || (protein == "PDZ")) {
      if (missing(position_numbering)) {
        position_numbering <- "crystal"
      }
      if (grepl("crystal", position_numbering)) {
        if (coloring == "es") {
          ligand = c(318, 322, 323, 324, 326, 327, 328, 330, 331, 332, 334, 337, 348, 352, 366, 373, 380)
          cluster = c(304, 305, 306, 309, 312, 354, 355, 357, 391, 395)
          list <- list(ligand = ligand, cluster = cluster)
          names(list) <- c("#3b73b1","#b0c7df")
        } else if (tolower(coloring) == "dds" || tolower(coloring) == "s") {
          high = c(322, 323, 325, 327, 329, 330, 336, 347, 351, 353, 359, 362, 363, 364, 372, 375, 376, 379, 386, 388)
          list <- list(high = high)
          names(list) <- c("#FFD700")
        } else {
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
          names(list) <- c("#69A019", "#FFD700", "#CC0000", "#FF9933")[1:length(list)]
        }
      } else if (position_numbering == "alignment") {
        # numbering alignment
        main = c(98)
        high = c(24, 28, 32, 33, 65, 72, 81, 102, 106, 119) # interesting ones
        list <- list(main = main, high = high)
        names(list) <- c("#69A019", "#FFD700", "#CC0000", "#FF9933")[1:length(list)]
      }
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
    } else if (tolower(protein) == "nov") {
      bind_ligand <- c(444, 443, 344, 345, 374, 442)
      others <- c(347, 396, 348, 391, 392, 375)
      list <- list("#69A019" = bind_ligand,            # green
                   "#FFD700" = others)                  # yellow
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
                         "#FFD700" = sort(yellow_),      # yellow
                         "#CC0000" = sort(red_),         # red
                         "#1874CD" = sort(blue_))        # blue
          } else {
            list <- list("#69A019" = sort(green_),       # green
                         "#FFD700" = sort(yellow_),      # yellow
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
          counts <- read.csv2("Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
          counts <- as.matrix(counts, ncol = 4)

          if (!is.na(as.numeric(colors)) && !is.null(colors)) {
            round_categories <- as.numeric(colors)
          } else {
            round_categories <- 1
          }

          list <- classify_nodes(counts, round_categories = round_categories, mix = TRUE, base_colors = colnames(counts))
        }
      }
    } else {
      list <- list()
    }
    if (!is.null(interesting_pos)) {  # !is.null(list) ?!
      warning("No interesting positions known")
      list <- list()
    }
  }

  if (for_coloring) {
    return(list)
    # n <- interesting_positions(protein = "PDZ")
    # list <- lapply(unique(names(n)), function(color) {return(unname(n[which(names(n) == color)]))})
    # names(list) <- unique(names(n))
  } else {
    # do not rename duplicate names
    return(conv_for_coloring(list))
    # return(unlist(list))
  }
}

## it seems much easier to convert from the for_coloring version to the other one,
## so for simplicity's sake I do it that way for now
## that means that the standardized way to pass int_pos would be the for_coloring version
conv_for_coloring <- function(int_pos) {
  # int and n are undefined in this function. I therefore copied the code from interesting_positions and
  # replaced ti by this function there
  # list <- lapply(unique(names(int)), function(color) {return(unname(int[which(names(int) == color)]))})
  # names(list) <- unique(names(n))
  # return(list)
  return(setNames(unlist(int_pos, use.names = FALSE), rep(names(int_pos), lengths(int_pos))))
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
                       plot_only_subgraphs = NULL, subgraphs, numerical = TRUE, output_formats, mute_all_plots = FALSE) {
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
        plot_graph_new(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode,
                       graph_layout = graph_layout_i, outpath = outpath, caption = caption,
                       plot_as_subgraphs = plot_as_subgraphs_i, plot_only_subgraphs = plot_only_subgraphs,
        subgraphs = subgraphs, output_formats = output_formats, mute_all_plots = mute_all_plots)
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
        } else if (format == "svg") {
          svg(paste(outpath, "_", graph_layout, "_colored-", coloring, ".svg", sep = ""))
        }
      } else {
        if (format == "pdf") {
          pdf(paste(outpath, ".pdf", sep = ""))
        } else if ((format == "ps") || (format == "postscript")) {
          postscript(paste(outpath, ".ps", sep = ""), paper = "special", width = 10, height = 9)
        } else if (format == "svg") {
          svg(paste(outpath, ".svg", sep = ""))
        }
      }
      plot(pc_graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, drawNode = drawnode, main = caption)
      dev.off()
    }
  }
}


# TODO: call plot_graphs oder so

## instead of fixing/working on the above function I think it is easier to write a new plot function the way we need it and then integrate it into the program if necessary
## what to do about drawnode, if it would be missing it was previously computed through
## node_function_for_graph which, however needs coloring
## for now I assume that this has been already computed and is NOT missing
plot_graph_new <- function(graph, fillcolor, edgecolor=NULL, drawnode, caption="", graph_layout="dot", outpath="",
                           plot_as_subgraphs= FALSE, plot_only_subgraphs = NULL, subgraphs = NULL,
                           output_formats = "pdf", mute_all_plots = mute_all_plots) {

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

# deprecated. Use
#
# parameters_to_info_file(parameters_for_info_file, outpath)
# # pc_func <- function(outpath) {return(pc_fun(outpath))}
# pc_func <- function_set_parameters(pc_fun, parameters = list(outpath = outpath))
# loaded_object_ok_fun <- function(pc) {return(length(pc@graph@nodes) != dim(data)[2])}
# compute_if_not_existent(filename = paste(outpath, "-pc.RData", sep = ""), FUN = pc_func, obj_name = "pc",
#                         compute_anew = compute_pc_anew, loaded_object_ok_fun = loaded_object_ok_fun)
#
# instead of get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file, data = data).
#
# Compute pc if necessary
get_pc <- function(pc_fun, outpath, compute_pc_anew, parameters_for_info, data) {
  parameters_to_info_file(parameters_for_info, outpath)
  if ((file.exists(paste(outpath, "-pc.RData", sep = ""))) && !(compute_pc_anew)) {
    filename <- paste(outpath, "-pc.RData", sep = "")
    load(filename)
    if (!exists("pc")) {
      warning("The file did not contain an object of name 'pc'!")
    } else if (length(pc@graph@nodes) != dim(data)[2]) {
      warning("The pc-object in the file did not correspond to the data!")
    } else {
      print(paste("pcAlgo-object loaded from ", filename, ".", sep = ""))
      print_pc_results_to_info_file(paste(outpath, sep = ""), pc)
      return(pc)
    }
    # file.copy(paste(outpath, "-info-pc.txt", sep = ""), paste(outpath, "-info.txt", sep = ""), overwrite = TRUE)
  }  else {
    warning("pc object not found.")
    warning(outpath)
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


scale_effects <- function(effects, rank = FALSE, amplification_factor = FALSE, neg_effects = "pos") {  neg_effects = "abs"

  amplify_with_factor <- function(effects, element_that_should_be_scaled_to = 2,
                                  value_the_element_should_be_scaled_to = 0.9, cut_values_at = 1) {
    sorted_effects <- sort(effects, decreasing = TRUE)
    if (sorted_effects[element_that_should_be_scaled_to] == 0) {
      return(effects)
    }
    factor = value_the_element_should_be_scaled_to / sorted_effects[element_that_should_be_scaled_to]
    # factor <- 0.9 / sorted_effects_pos[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
    effects <- effects * factor
    if (is.numeric(cut_values_at)) {
      effects[,1][effects[,1] > cut_values_at] <- cut_values_at
    }
    # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
    return(effects)
  }

  effects_na <- as.matrix(effects[,1][is.na(effects[,1])])
  effects <- as.matrix(effects[,1][!is.na(effects[,1])])

  if (!length(effects) == 0) {
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
        effects <- rbind(effects_pos, effects_neg)
        # effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
        # effects <- effects_
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
        # effects_pos <- as.matrix(effects[,1][effects[,1] >= 0 & !is.na(effects[,1])])
        effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
        # effects_neg <- as.matrix(effects[,1][effects[,1] < 0 & !is.na(effects[,1])])
        # effects_neg <- cbind(effects_neg[!is.na(effects_neg),]) # NAs nicht doppelt haben -- in neg rausschmeißen!
                                                                # TODO: in dan anderen Fällen ggf. auch !!!!!
        effects_neg <- -effects_neg


        if (amplification_factor) {
          effects_pos <- amplify_with_factor(effects_pos)
        }

        if (amplification_factor && dim(effects_neg)[1] > 1) {
          effects_neg <- amplify_with_factor(effects_neg)
        }

        effects_neg <- -effects_neg
        effects <- rbind(effects_pos, effects_neg)
        # effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
        # effects <- effects_
      }
    }
  }

  effects <- rbind(effects, effects_na)

  effects <- effects[order(rownames(effects)), , drop = FALSE]
  # effects <- effects

  return(effects)
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

# This function is used e.g. for plotting connected components in pymol.
# Normally, the connected components are sorted alphabetically.
# Since in proteins the false impression that neighbouring residues are in the
# same connected components could arise when they have subsequent and thus
# indistinguishable rainbow colors, the components can be mixed (mix argment given).
# Then, mix_offset determines how many other positions there are (at least) between two
# originally neighbouring components

# if mix = "every_mix_offset_th", the mixed version is obtained by appending to the
# new list every i'th element from the original list, where i = mix_offset + 1,
# continuing with incremeted offset every time the end of the list is reached
reorder_list_of_lists <- function(list, ordering, mix_mode = "mix_offset_between_original_neighbours", mix_offset = floor(sqrt(length(list))),
                                  sort_mode, sort_descending = TRUE) {
  if (ordering == "mix") {
    if (mix_offset > 0) {
      mixed_connected_components <- list()
      if (mix_mode == "every_mix_offset_th") {
        for (rest in seq(0, mix_offset)) {
          for (multipicity in seq(0, floor(length(list)/mix_offset))) {
            index_in_conn_comp <- multipicity * mix_offset + rest + 1
            if (index_in_conn_comp <= length(list)) {
              mixed_connected_components[[length(mixed_connected_components) + 1]] <- list[[index_in_conn_comp]]
              # mixed_connected_components <- append(mixed_connected_components, connected_components[[index_in_conn_comp]])
            }
          }## for(i in 1:length(list)) {
##              for(j in 1:length(connected_components)) {
##  if(list[i] %in% connected_components[[j]])
## names(list)[i] = colors[j]
##              }
          ##              }
##           for (i in 1:length(connected_components)) {
## for (j in 1:length(connected_components[[i]])) {

##           nAttrs$fillcolor[connected_components[[i]][j]] = colors[i]
## }
## }
        }
      } else {
        # there are at least mix_offest + 1 blcoks (+ one block with the rest of the division).
        # To obtain the permutation, first, the first elements of all the blocks are taken,
        # then the second ones and so on.
        block_size = floor(length(list) / mix_offset)
        for (offset in 0:(block_size-1)) {
          for (block in seq(0, ceiling(length(list) / block_size))) {
            index_in_conn_comp <- block * block_size + offset + 1
            if (index_in_conn_comp <= length(list)) {
              mixed_connected_components[[length(mixed_connected_components) + 1]] <- list[[index_in_conn_comp]]
              # mixed_connected_components <- append(mixed_connected_components, connected_components[[index_in_conn_comp]])
            }
          }
        }
      }
      list <- mixed_connected_components
    }
  } else if (grepl("sort", ordering)) {
    if (sort_mode == "length") {
      if (sort_descending) {
        sorted_connected_components <- list[order(vapply(list, length, 1L), decreasing = TRUE)]
      } else {
        sorted_connected_components <- list[order(vapply(list, length, 1L), decreasing = FALSE)]
      }
    }
    list <- sorted_connected_components
  }
  return(list)
}

position_clustering_from_clustering_with_duplicates <- function(clustering_with_duplicates) {
  positions <- unique(sapply(names(clustering_with_duplicates), function(long_name) return(gsub("-.*","",long_name))))
  k = max(clustering_with_duplicates)


  ## averaged_clusters <- sapply(positions, function(position) {
  ##  position_clustering <- clustering_with_duplicates[which(grepl(position, names(clustering_with_duplicates)))]
  ##  mean_cluster <- mean(position_clustering)
  ## names(mean_cluster) <- position
  ##  return(mean_cluster)
  ##})

  ## averaged_clusters <- round(averaged_clusters)
  rb_cols = substr(rainbow(k), 1, 7)
  averaged_clusters <- sapply(positions, function(position) {
    position_clustering <- clustering_with_duplicates[which(grepl(position, names(clustering_with_duplicates)))]
    cols <- rb_cols[position_clustering]

    if(length(cols) > 1) {
      mean_col <- mixcolors(cols)
    } else {
      mean_col = cols[[1]]
    }
    ## weird stuff
    names(mean_col) = NULL
    ##  names(mean_col) <- position
    return(mean_col)
    })

  cl <- clusterlist_from_membershiplist(averaged_clusters)

  return(cl)
  # TODO!
}



