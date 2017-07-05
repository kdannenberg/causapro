# correlation between positions with highest variation (over threshold),
# wenn firgure margins too large --> threshold erhöhen
# rot: 1H, grün: 13C 

setwd("~/Viren/R/Code")
source("compute_DAG_numerical.R")

protein <- "GTB"
type_of_data <- "NMR-2"

scatterplot_NMR <- function(data, position_1, position_2, nuclei) {
  # plot(data[,position_1], data[,position_2], col = c("red", "green"))
  data <- data[grepl(nuclei, rownames(data)), ]
  if (nuclei == "1H") {
    colors <- "red"
  } else if (nuclei == "13C") {
    colors <- "green"
  } else {
    colors <- c("red", "green") 
  }
  plot(data[,as.character(position_1)], data[,as.character(position_2)], xlab = as.character(position_1), ylab = as.character(position_2),  col = colors)
}

# Varianz
var_of_pos <- function(data) {
  garbage <- graphics.off()
  vars <- apply(data, 2 ,var) 
  plot(x = colnames(data), y = vars, xlab = "Position", ylab = "Variation")
  return(vars)
}

rem_cols_by_colname <- function(data, remove) {
  return(data[, !(colnames(data) %in% remove)])
  # MSA_out <- MSA_in[, gap_freq < threshold]
}

# var_of_pos(data[,(as.numeric(colnames(data)) > 180 && as.numeric(colnames(data)) < 220)])

# dimension is the string to be looked for in the rownames (eg."1H")
scale_NMR_dimensions <- function(data, dimension, factor) {
  data[grepl(dimension, rownames(data)), ] <- data[grepl(dimension, rownames(data)), ] * factor
  return(data)
}

#rem: positions to be removed before computation
scatterplots_of_most_variating_positions <- function(data, threshold, nuclei, int_pos) {
  vars <- var_of_pos(data)
  
  if (missing(int_pos)) {
    int_pos <- vars[vars > threshold]
  } else {
    int_pos <- vars[int_pos]
  }
  
  print(int_pos)
  
  n <- length(int_pos)
  combinations <- 0.5 * (n-1) * n
  
  dim1 <- round(sqrt(combinations))
  dim2 <- ceiling(combinations / dim1)
  # dim1 <- 5
  # dim2 <- 5
  
  par(mfrow = c(dim1, dim2))
  
  counter = dim1 * dim2
  for (i in names(int_pos)) {
    for (j in names(int_pos)) {
        if (i < j) {
            ## added writing to file temporarily, looking for a nice solution
            pdf(paste0("~/report/test", i, j, ".pdf"))
            scatterplot_NMR(data, i, j, nuclei)
            dev.off()
            counter <- counter - 1
            if (counter == 0) {
                break;
            }
        }
    }
    if (counter == 0) {
      break;
    }
  }
  dev.off()
}

# scatterplot_NMR(data, 127, 177)

# state = "don"
# state = "acc"
state = "don+acc"
# state = "all"
source_of_data = paste(protein, type_of_data, state, sep = "-")
filename <- paste("../Data/", source_of_data, ".csv", sep = "")
data <- read_data(filename, transpose = TRUE)
colnames(data) <- sapply(strsplit(colnames(data), " "), function(x) x[1])

rem_don_acc <- c(143, 175, 184, 189, 266, 329)
rem_don <- c(210)
# rem <- c(rem_apo, rem_other)
if (state == "all") {
  rem <- c(rem_don_acc, rem_don)
} else if (state == "don+acc") {
  rem <- rem_don_acc
} else if (state == "don") {
  rem <- rem_don
} else if (state == "acc") {
  rem <- c()
}
data <- rem_cols_by_colname(data, rem)

# nuclei <- "13C"
# nuclei <- "1H"
nuclei <- ""


# if (nuclei == "1H") {
#   threshold = 3000
# } else if (nuclei == "13C") {
#   threshold = 3000
# } else {
#   threshold = 4000
# }

int_pos <- c("120", "177", "214", "289", "299")
threshold <- 4000
# scatterplots_of_most_variating_positions(data, threshold = threshold, nuclei, int_pos = int_pos)

## mean(data[grepl("1H", rownames(data)), ])
## # [1] -0.1695813
## mean(data[grepl("13C", rownames(data)), ])
## # [1] -3.520347
## mean(data[grepl("13C", rownames(data)), ]) /  mean(data[grepl("1H", rownames(data)), ])
## # [1] 20.75905
# data_mean_adj <- scale_NMR_dimensions(data, "1H", 20.75905)
# scatterplots_of_most_variating_positions(data_mean_adj, threshold = 500000, nuclei)
## dann sind die C-Werte sehr zusammengedrängt
 
## aber Varianz gleich machen ist besser: 
## mean(var_of_pos(data[grepl("13C", rownames(data)), ])) /  mean(var_of_pos(data[grepl("1H", rownames(data)), ]))
## # [1] 1.546044
data_var_adj <- scale_NMR_dimensions(data, "1H", 1.546044)
scatterplots_of_most_variating_positions(data_var_adj, threshold = 6000, nuclei, int_pos = int_pos)

# data_adjusted <- scale_NMR_dimensions(data, "1H", 2)
# scatterplots_of_most_variating_positions(data_adjusted, threshold = 10000, nuclei)
 
