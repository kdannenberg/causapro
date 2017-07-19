library("pcalg")
library("plyr")   # für alply (apply to array and return list)
library("ggm")    # for pcor  (partial correlation)
library("RBGL")   # for tsort (top. sorting)

set.seed(17)
p <- 10

n <- 10000
myDAG <- pcalg::randomDAG(p, prob = 0.07) ## true DAG; default: prob = 0.2
dat <- rmvDAG(n, myDAG)
cov.d <- cov(dat)

suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)

par(mfrow = c(3,3))

pDAG <- pc.fit@graph
plot(pDAG)  # Warum sind 3 und 4 hier andersrum??


topological_nodesnames <- function(DAG) {
  sorting <- tsort(DAG)
  # indices <- sapply(sorting, function(i) {return(which(sorting == i))})
  # top_node_names <- order(names(indices), indices)
  top_node_names <- sapply(1:7, function(i) {return(which(sorting == i))})
  # nodes(DAG) <- as.character(top_node_names)
  return(as.character(top_node_names))
}

# only the parameters x, y and graphEst are used;  rest is only to mathc function definition of ida()
pseudo_ida_by_causalEffect <- function(x, y, mcov, graphEst, method = "", y.notparent = FALSE, verbose = FALSE, all.dags = NA) {
  cat(paste("Effect from", x, "on", y, "\n"))
  # allDAGS <- pdag2allDags(wgtMatrix(graphEst))$dags
  allDAGs <- set_of_DAGs(graphEst) 
  # 
  for (DAG_as in allDAGs) {
  #   print(paste("DAG", i))
  #   m <- matrix(allDAGS[i,], nrow = p, ncol = p, byrow = TRUE)
  #   
  #   #Why is it necessary to transpose m?!
  #   DAG_as <-  as(t(m), "graphNEL")
    # in einer Zeile: 
    # DAG_as <- as(t(matrix(pdag2allDags(wgtMatrix(pDAG))$dags[i,],p,p, byrow = TRUE)), "graphNEL")
    plot(DAG_as)
    
    DAG_as_wgt <- set_edge_weights_for_graph(DAG_as, mcov)
    
    print(wgtMatrix(DAG_as_wgt))
    
    # nodenams topologisch machen
    orignial_node_names <- nodes(DAG_as_wgt)
    
    sorting <- tsort(DAG_as_wgt)
    top_node_names <- sapply(1:p, function(i) {return(which(sorting == i))})
    # nodes(DAG_as_wgt) <- topological_nodesnames(DAG_as_wgt) 
    nodes(DAG_as_wgt) <- as.character(top_node_names)
    # x_caus_eff <- which(top_node_names == x)
    # y_caus_eff <- which(top_node_names == y)
    x_caus_eff <- top_node_names[x]
    y_caus_eff <- top_node_names[y]
    
    # dann neue nodenames sortieren 
    wgt_mat <- wgtMatrix(DAG_as_wgt)
    wgt_mat <- wgt_mat[order(rownames(wgt_mat)), order(colnames(wgt_mat))] 
    
    DAG_as_wgt <-  as(t(wgt_mat), "graphNEL")
    
    plot(DAG_as_wgt, col = "green")
    
    causal_effect <- causalEffect(DAG_as_wgt, x = x_caus_eff, y = y_caus_eff)
    
    cat(paste("hier:", causal_effect, "\n"))
    # #maybe better: ftM2graphNEL???
    # source: https://support.bioconductor.org/p/90421/
    # DAG_ftM2 <- ftM2graphNEL()
    # plot(DAG)
    cat(paste("ida:", ida(x.pos = x, y.pos = y, graphEst = DAG_as, mcov = mcov)))
    cat("\n")
  }
  # cat(paste("ida:", ida(x.pos = x, y.pos = y, graphEst = pDAG, mcov = mcov)))
}

set_edge_weights_for_graph <- function(graph, cov) {
  edges <- names(graph@edgeData@data)
  
  setweight <- function(edge, graph, cov, mode = "adj_set") {
    l <- list()
    nodes <- unlist(strsplit(edge, "\\|"))
    if (mode == "adj_set") {
      # inst_vars <- instrumentalVariables(conv_to_r(graph, type_of_graph = "dag"), nodes[1], nodes[2])
      adj_set_single_door <- adjustmentSets(conv_to_r(graph, type_of_graph = "dag"), nodes[1], nodes[2], type = "minimal", effect = "direct")
      
      adj_set <- c()
      for (adj_var in adj_set_single_door[[1]]) {
        adj_set <- c(adj_set, adj_var)
      }
      
      l$weight <- pcor(c(nodes[1], nodes[2], adj_set), cov)
    } else {
      # graph@edgeData@data[[edge]] <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
      # l[[edge]] <- list()
      # l[[edge]]$weight <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
      # l[[edge]] <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
      # return(cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]])
      l$weight <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
      # l$weight <- cov[nodes[2],nodes[1]] / cov[nodes[2],nodes[2]]    # andersrum
    }
    
    
    
    return(l)
  }
  
  l <- lapply(edges, setweight, graph, cov)
  names(l) <- edges
  
  graph@edgeData@data <- l
  
  return(graph)
}


set_of_DAGs <- function(pdag) {
  allDAGS_m <- pdag2allDags(wgtMatrix(pdag))$dags
  m <- function (line) {
    return(matrix(line, nrow = p, ncol = p, byrow = TRUE))
  }
  allDAGS_adj <- alply(allDAGS_m, 1, m)
  # allDAGS_adj <- apply(allDAGS_m, 1, matrix, nrow = p, ncol = p, byrow = TRUE)
    # m <- matrix(allDAGS[i,], nrow = p, ncol = p, byrow = TRUE)
    
    # Why is it necessary to transpose m?!
    # DAG_as <-  as(t(m), "graphNEL")
  allDAGs <- lapply(allDAGS_adj, function(m) {return(as(t(m), "graphNEL"))})
  return(allDAGs)
}


# pseudo_ida_by_causalEffect(x = 2, y = 6, graphEst = pDAG, mcov = cov.d) # in Graph 3 nur fast gleich!
# pseudo_ida_by_causalEffect(x = 5, y = 6, graphEst = pDAG, mcov = cov.d) # in Graph 3 nur fast gleich!

pseudo_ida_by_causalEffect(x = 8, y = 6, graphEst = pDAG, mcov = cov.d) # in Graph 3 nur fast gleich!   für p = 10


ida_for_set_of_DAGs <- function(x, y, mcov, graphEst, method = "", y.notparent = FALSE, verbose = FALSE, all.dags = NA) {
  allDAGs <- set_of_DAGs(graphEst) 
  ida_DAG <- function(graphEst) {
    return(ida(x, y, mcov, graphEst, method = "", y.notparent = FALSE, verbose = FALSE, all.dags = NA))
  }
  effects <- lapply(allDAGs, ida_DAG)
}
