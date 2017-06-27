set.seed(17)
p <- 7

n <- 10000
dat <- rmvDAG(n, myDAG)
cov.d <- cov(dat)

suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)

par(mfrow = c(3,3))

pDAG <- pc.fit@graph
plot(pDAG)  # Warum sind 3 und 4 hier andersrum??

pseudo_ida_by_causalEffect(x = 2, y = 7, graphEst = pDAG, mcov = cov.d)


# only the parameters x, y and graphEst are used;  rest is only to mathc function definition of ida()
pseudo_ida_by_causalEffect <- function(x, y, mcov, graphEst, method = "", y.notparent = FALSE, verbose = FALSE, all.dags = NA) {
  allDAGS <- pdag2allDags(wgtMatrix(pDAG))$dags
  
  for (i in 1:dim(allDAGS)[1]) {
    m <- matrix(allDAGS[i,],p,p, byrow = TRUE)
    
    #Why is it necessary to transpose m?!
    DAG_as <-  as(t(m), "graphNEL")
    plot(DAG_as)
    
    print(causalEffect(DAG_as, x = x, y = y))
    # #maybe better: ftM2graphNEL???
    # source: https://support.bioconductor.org/p/90421/
    # DAG_ftM2 <- ftM2graphNEL()
    # plot(DAG)
  }
  cat(paste("ida:", ida(x.pos = x, y.pos = y, graphEst = pDAG, mcov = mcov)))
}