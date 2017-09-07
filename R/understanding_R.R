# value = 2
# 
# funfun <- function(x) {
#   return(x + value)
# }
# 
# funfun2 <- function(value) {
#   function(x) {x + value}
# }
# 
# 
# 
# 
# funfun(3)
# ff2 <- funfun2(3)
# ff2(11)
# 
# value = 6
# 
# funfun(3)
# ff2 <- funfun2(3)
# ff2(11)


allgemeineFunktion <- function(a,b,c) {
  print(paste(a,b,c))
}

b = 2
c = 3


spezielleFunktion <- function(a) {
  allgemeineFunktion(a,b,c)
}

as.list(environment(spezielleFunktion))
spezielleFunktion("a")

b = 11

as.list(environment(spezielleFunktion))
spezielleFunktion("a")

######################################################################

besserespezielleFunktion <- function(b,c) {
  return(function(a) {allgemeineFunktion(a,b,c)})
}

b = 2
c = 3

spezielleFunktion <- besserespezielleFunktion(b, c)

as.list(environment(spezielleFunktion))
spezielleFunktion("a")

b = 11

as.list(environment(spezielleFunktion))
spezielleFunktion("a")

