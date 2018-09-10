div_rows_cols_by_sd <- function(data) {
  sd_rows <- apply(data, 1, sd)
  sd_cols <- apply(data, 2, sd)
  data <- data/sd_rows
  data <- t(t(data)/sd_cols)
  return(data)
}

div_rows_by_sd <- function(data) {
  sd_rows <- apply(data, 1, sd)
  data <- data/sd_rows
  return(data)
}

div_cols_by_sd <- function(data) {
  sd_cols <- apply(data, 2, sd)
  data <- t(t(data)/sd_cols)
  return(data)
}
