library(seqinr)

## setwd("~/Documents/Uni/Viren/R/Code")

source_of_data <- "al_pdz"

compute_q_ <- function(MSA) {
  return(table(MSA) / (dim(MSA)[1] * dim(MSA)[1]))
}

f_i_x <- function(data, i, x) {
  numerator <- table(data[,i])[x]
  names(numerator) <- c()
  return(numerator / nrow(data)) # z.B.  table(toy[,1])["A"]
}

f_ij_xy <- function(data, i, j, x, y) {
  position <- paste(x,y)
  table <- table(paste(data[,i], data[,j]))
  # return(table(paste(data[,i], data[,j]))) #[paste(x,y)]  / nrow(data))
  if (is.na(table[position])) {
    numerator = 0
  } else {
    numerator <- table[position]
  }
  names(numerator) <- c() 
  return(numerator / nrow(data))
}

C_ij_xy <- function(data, i, j, x, y) {
  fij <- f_ij_xy(data, i, j, x, y)
  fi <- f_i_x(data, i, x)
  fj <- f_i_x(data, j, y)
  diff <- fi * fj
  result <- fij - diff
  return(result)
}

C_ij_xy_mul <- function(data, i, j, x, y) {
  fij <- f_ij_xy(data, i, j, x, y)
  fi <- f_i_x(data, i, x)
  fj <- f_i_x(data, j, y)
  diff <- fi * fj
  result <- fij - ((fi-fij) * (fj-fij))
  return(result)
}

C_ij_xy_add <- function(data, i, j, x, y) {
  fij <- f_ij_xy(data, i, j, x, y)
  fi <- f_i_x(data, i, x)
  fj <- f_i_x(data, j, y)
  diff <- fi * fj
  result <- fij - (fi-fij) - (fj-fij)
  return(result)
}



# elegant?
# Cijxy <- array(dim = c(nrow(MSA), nrow(MSA), 21, 21))
# Cijxy[i,j,x,y] <- C_ij_xy(MSA, i, j, x, y)

# find the k most frequent aa at each position in data
# if there i < k unique elements in a column, the last k-i elements in this column will be "NA" in the output matrix 
# result is a n x k matrix
most_freq_k <- function(data, k) {
  most_freq_list <- apply(data, 2, function(m) sort(table(m), decreasing=TRUE))
  most_freq <- lapply(most_freq_list, function(m) m[1:k])
  return(matrix(sapply(most_freq, names), ncol = ncol(data)))
}

# for each column, aa ar sorted by frequency; result is a list
most_freq <- function(data) {
  most_freq_list <- apply(data, 2, function(m) sort(table(m), decreasing=TRUE))
  return(sapply(most_freq_list, names))
}

# regard k most frequent aa at each position
mat_from_measure_i_j_mostfreqi_mostfreqj <- function(data, k, FUN, output_matrix = TRUE) {
# C_ij_xy_mat <- function(data, k, output_matrix = TRUE) {
  # CC_ij_xy <- matrix(dim=c(ncol(data), k^2*ncol(data)))
  
  # C_ij_xy_array <- array(dim=c(ncol(data), ncol(data), k, k))
  result_array <- array(dim=c(ncol(data), ncol(data), k, k))
  
  most_freq_aa <- most_freq_k(data, k)
  
  # return(most_freq_aa)
  
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      for (k_i in 1:k) {
        for (k_j in 1:k) {
          x <- most_freq_aa[k_i,i]
          y <- most_freq_aa[k_j,j]
          if (is.na(x)) {
            print("Achtung, x=NA.")
          } else if (is.na(y)) {
            print("Achtung, y=NA.")
          } else {
            result_array[i,j,k_i,k_j] <- FUN(data, i, j, x, y)
            # C_ij_xy_array[i,j,k_i,k_j] <- C_ij_xy(data, i, j, x, y)
          }
        }
      }
    }
  }
  
  # if (output_matrix) {
  #   return(matrix(C_ij_xy_array, nrow = k^2 * ncol(data), byrow=TRUE))
  # } else {
  #   return(C_ij_xy_array)
  # }
  if (output_matrix) {
    return(matrix(result_array, nrow = k^2 * ncol(data), byrow=TRUE))
  } else {
    return(result_array)
  }
}

 
