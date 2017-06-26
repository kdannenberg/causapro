## Simulate the true DAG
set.seed(123)
p <- 7
# p <- 2
myDAG <- pcalg::randomDAG(p, prob = 0.2) ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

## simulate data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)
cov.d <- cov(dat)

## estimate CPDAG (see help on the function "pc")
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)

if(require(Rgraphviz)) {
  op <- par(mfrow=c(1,3))
  plot(myDAG,        main="true DAG")
  plot(myCPDAG,      main="true CPDAG")
  plot(pc.fit@graph, main="pc()-estimated CPDAG")
  par(op)
}

(eff.est1 <- pcalg::ida(2,5, cov.d, pc.fit@graph))## method = "local" is default
(eff.est2 <- pcalg::ida(2,6, cov.d, pc.fit@graph))
(eff.est3 <- pcalg::ida(2,7, cov.d, pc.fit@graph))
## These three computations can be combinded in an efficient way
## by using idaFast :
(eff.estF <- idaFast(2, c(5,6,7), cov.d, pc.fit@graph))

# ida(1,2,cov.d, pc.fit@graph, verbose = TRUE)
# ida(7,3,cov.d, pc.fit@graph, verbose = TRUE)


check_formula <- function(x, y, par_x) {
  cat("IDA:")
  cat(ida(x,y,cov.d,pc.fit@graph))
  
  cat("\n")
  
  if (!missing(par_x)) {
    cat("Formel:")
    formel_value <- solve(C[x, x], C[x, y, drop = FALSE])[1, ]
    cat(formel_value)
  } else {
    cat("Formel:")
    formel_value <- cov.d[x,y] / cov.d[x,x]
    cat(formel_value)
  }
  
  cat("\n")
  
  cat("equal?: ")
  cat(as.numeric(ida(x,y,cov.d,pc.fit@graph)) == formel_value)
  cat("\n")
}

for (x in c(1,2,3,5,6)) {  # immer TRUE
  for (y in seq(1:7)) {
    print(paste0("x = ", x, ", y = ", y))
    check_formula(x,y)
  }
}


