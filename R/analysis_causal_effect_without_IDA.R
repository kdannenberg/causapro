set.seed(17)
p <- 7

n <- 10000
myDAG <- pcalg::randomDAG(p, prob = 0.07) ## true DAG; default: prob = 0.2
dat <- rmvDAG(n, myDAG)
cov.d <- cov(dat)

suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)

par(mfrow = c(3,3))

pDAG <- pc.fit@graph
plot(pDAG)  # Warum sind 3 und 4 hier andersrum??


# only the parameters x, y and graphEst are used;  rest is only to mathc function definition of ida()
pseudo_ida_by_causalEffect <- function(x, y, mcov, graphEst, method = "", y.notparent = FALSE, verbose = FALSE, all.dags = NA) {
  cat(paste("Effect from", x, "on", y, "\n"))
  allDAGS <- pdag2allDags(wgtMatrix(pDAG))$dags
  
  for (i in 1:dim(allDAGS)[1]) {
    print(paste("DAG", i))
    m <- matrix(allDAGS[i,], nrow = p, ncol = p, byrow = TRUE)
    
    #Why is it necessary to transpose m?!
    DAG_as <-  as(t(m), "graphNEL")
    # in einer Zeile: 
    # DAG_as <- as(t(matrix(pdag2allDags(wgtMatrix(pDAG))$dags[i,],p,p, byrow = TRUE)), "graphNEL")
    plot(DAG_as)
    
    DAG_as_wgt <- set_edge_weights_for_graph(DAG_as, mcov)
    
    print(wgtMatrix(DAG_as_wgt))
    
    cat(paste("hier:", causalEffect(DAG_as_wgt, x = x, y = y), "\n"))
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
  
  setweight <- function(edge, graph, cov) {
    nodes <- unlist(strsplit(edge, "\\|"))
    # graph@edgeData@data[[edge]] <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
    l <- list()
    # l[[edge]] <- list()
    # l[[edge]]$weight <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
    # l[[edge]] <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
    # return(cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]])
    l$weight <- cov[nodes[1],nodes[2]] / cov[nodes[1],nodes[1]]
    # l$weight <- cov[nodes[2],nodes[1]] / cov[nodes[2],nodes[2]]    # andersrum
    return(l)
  }
  
  l <- lapply(edges, setweight, graph, cov)
  names(l) <- edges
  
  graph@edgeData@data <- l
  
  return(graph)
}


# pseudo_ida_by_causalEffect(x = 2, y = 6, graphEst = pDAG, mcov = cov.d) # in Graph 3 nur fast gleich!
pseudo_ida_by_causalEffect(x = 5, y = 6, graphEst = pDAG, mcov = cov.d) # in Graph 3 nur fast gleich!
