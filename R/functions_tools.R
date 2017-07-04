# small utilities

plus_minus_x <- function(vector, offset) {
  vector <- unlist(lapply(vector, function(n) return(seq(n - offset, n + offset))))
  return(sort(unique(vector)))
}

# behaves like paste, but does not add seps when an element is "",
# in other words, removes the "" before pasting
pastes <- function(..., sep = " ", collapse = NULL) {
  return(do.call("paste", c(list(...)[list(...) != ""], sep = sep, collapse = collapse)))
}

