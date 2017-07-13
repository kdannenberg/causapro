library(pcalg)
## Simulate the true DAG
# set.seed(123) # default
set.seed(12) # used below
# set.seed(14) # interesting
# set.seed(16) # ungerichtete Kanten!
# set.seed(17)
p <- 7 # default
# p <- 100
myDAG <- pcalg::randomDAG(p, prob = 0.2) ## true DAG; default: prob = 0.2
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

## simulate data from the true DAG
# n <- 10000 # default
n <- 1000
# dat <- rmvDAG(n, myDAG)
rm(.Random.seed, envir=globalenv())
dat <- rmvDAG(n, myDAG, errDist = "normal")

# print(dat)

cov.d <- cov(dat)

## estimate CPDAG (see help on the function "pc")
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p, solve.confl = FALSE, u2pd = "relaxed") # retry

# print(conflict_edges(pc.fit))

# if(require(Rgraphviz)) {
#   graphics.off()
#   op <- par(mfrow=c(1,2))
#   # plot(myDAG,        main="true DAG")
#   plot(myCPDAG,      main="true CPDAG")
#   plot(pc.fit@graph, main="pc()-estimated CPDAG")
#   par(op)
# }

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
  ida_value <- ida(x,y,cov.d,pc.fit@graph)
  cat(ida_value)
  
  cat("\n")
  
  if (!missing(par_x)) {
    cat("Formel:")
    formel_value <- solve(cov.d[x, x], cov.d[x, y, drop = FALSE])[1, ]
    cat(formel_value)
  } else {
    cat("Formel:")
    formel_value <- cov.d[x,y] / cov.d[x,x]
    cat(formel_value)
  }
  
  cat("\n")
  
  cat("equal?: ")
  cat(abs(as.numeric(ida_value) - formel_value) < 1e-10)
  cat("\n")
}

# Beispiel set.seed(12)
# for (x in c(1,2,3,5,6)) {  # ohne Eltern
#   for (y in seq(1:7)) {
#     print(paste0("x = ", x, ", y = ", y))
#     check_formula(x,y)
#   }
# }
# immer TRUE

# for (y in 1:7) {     # immer TRUE
#   check_formula(4,y,c(4,2,5))
# }

# for (y in 1:7) {      # immer TRUE
#   check_formula(7,y,c(7,2,6))
# }
# Ende des Beispiels

# Beispiel set.seed(2)
# for (x in c(1)) {  # ohne Eltern 
#   for (y in seq(1:7)) {
#     print(paste0("x = ", x, ", y = ", y))
#     check_formula(x,y)
#   }
# }
# immer TRUE

# for (y in 1:7) {     # immer TRUE
#   print(paste0("x = ", 6, ", y = ", y))
#   check_formula(6,y,c(6,4,5))
# }

# for (y in 1:7) {     # immer TRUE
#   print(paste0("x = ", 7, ", y = ", y))
#   check_formula(7,y,c(7,6,3))
# }

# for (y in 1:7) {     # NICHT immer TRUE (wegen eventueller Eltern (pa2))
#   print(paste0("x = ", 2, ", y = ", y))
#   check_formula(2,y,c(5))
# }
# Ende des Beispiels

# causalEffect(pc.fit@graph, 1,2)


# # causalEffect addiert für alle gerichteten Pfade zwischen x und y das Produkt der Kantengewichte entlang des Pfades:
set.seed(123)
p <- 7
myDAG <- pcalg::randomDAG(p, prob = 0.2)
cov_t <- trueCov(myDAG)

plot(myDAG)

causalEffect(x=2, y=7, g=myDAG)
# [1] 0.7425168
myDAG@edgeData@data$`2|3`$weight*myDAG@edgeData@data$`3|5`$weight*myDAG@edgeData@data$`5|7`$weight+myDAG@edgeData@data$`2|7`$weight
# [1] 0.7425168

edgeWeight <- function(x,y) {
  cov.d[x,y] / cov.d[x,x]
}


# Test own function to set weights:
myDAG_wgt <- set_edge_weights_for_graph(myDAG, cov_t)

print(wgtMatrix(myDAG))
print(causalEffect(x=2, y=7, g=myDAG))

print(wgtMatrix(myDAG_wgt))
print(causalEffect(x=2, y=7, g=myDAG_wgt))


# statistic
# p = 100
# n = 1000000 -> 8 Konflikt-Kanten
# n = 100000  -> 14 Konflikt-Kanten
# n = 10000   -> 9 / 20 / 25 / 6 / 11 Konflikt-Kanten
# n = 1000    -> 50 Konflikt-Kanten