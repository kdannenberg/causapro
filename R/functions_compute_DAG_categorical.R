# all the functions that would be useful in an R package; top level computing functions

#' if read_from_file, source_of_data should be string giving the path to the file, 
#' otherwise source_of_data should be a matrix containing the alignment
#' map_AS specifies wether to pool aminoacids in groups (vectors) given in a list "cluster", 
#' in this case, all amino acids in one group are replaced by the name of the group in the list
#' ...
#' 
#' 

library(stringr)

# source("general_functions.R")

estimate_DAG <- function(source_of_data, map_AS = TRUE, cluster, chop_MSA_from = 0, chop_MSA_to = NULL, chop_rownames = FALSE, 
                        remove_singular_columns = TRUE, remove_cols_gaps_threshold = 1, test_for_ci = "sp-x2", outpath, alpha = 0.01, output_neato = FALSE, output_twopi = FALSE) {

    par_string = paste("source: ", source_of_data, ", ci_test: ", test_for_ci, ", alpha: ", alpha, ", cluster: ", paste(print.cluster(cluster), " (", length(cluster), " clusters)", sep = ""), ", from_pos: ", chop_MSA_from, ", to_pos: ", chop_MSA_to, ", rem_gaps: ", remove_cols_gaps_threshold, sep = "")
    ## ich erstelle erst einen allgemeinen String mit den
    ## Parameterinfos und baue mir daraus dann info fuer die
    ## Info-Datei und main fuer die Graphiken
    info_str = paste(par_string, ", map_AS: ", map_AS, ", chop_rownames: ", chop_rownames, ", remove_singular_columns: ", remove_singular_columns, ", outpath: ", outpath, sep = "")
    info_str = str_replace_all(info_str, pattern = "from_pos", "chop_MSA_from")
    info_str = str_replace_all(info_str, pattern = "to_pos", "chop_MSA_to")
    info_str = str_replace_all(info_str, pattern = "source", "source_of_data")
    info_str = str_replace_all(info_str, pattern = "rem_gaps", "remove_cols_gaps_threshold")
    info = str_replace_all(info_str, pattern = ", ", ",\n")
    info = paste(info, "\n\n", sep = "")
    capture = strwrap(par_string, width = 50)
    sink(file = paste(outpath, "-info.txt", sep = ""))
      cat(info)
    sink()
    
    if (is.numeric(test_for_ci)) {
        ci_test <- disCItest_given_nmin(test_for_ci) # ordering of the values necessary
    } else {
        ## if (is.null(number_of_permutations) || !(grepl("mc", test_for_ci))) {
        ci_test <- ci_test_pc(test_for_ci)
        ## } else {
        ##   ci_test <- ci_test_pc(test_for_ci, B <- number_of_permutations)
        ## }
    }
  filename <- paste("../Data/", source_of_data, sep = "")
  MSA <- readAlignment(filename)
  
  colnames(MSA) <- seq(1:dim(MSA)[2])
  
  if (chop_MSA_from != 1 || !is.null(chop_MSA_to)) {
    if (is.null(chop_MSA_to)) {
      chop_MSA_to <- dim(MSA)[2]
    }
    MSA <- MSA[,chop_MSA_from:chop_MSA_to]
  }
  
  if (map_AS) {
    MSA <- map_aminoacids(MSA, cluster)
  }
  
  n_lev <- length(unique(as.list(MSA))) #22
  sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
  print(paste("Number of levels in the data:", n_lev))
  sink()
  
  if (!(remove_cols_gaps_threshold==1)) {
    MSA <- remove_gaps(MSA, threshold=remove_cols_gaps_threshold, n_lev, allAS, outpath)
  }
  
  if (remove_singular_columns) {
    MSA <- remove_sing_col(MSA, outpath)
  }
  
  if (chop_rownames) {
    rownames(MSA) <- sapply(rownames(MSA), substr, 0, 10)
  }

  n_variables <- dim(MSA)[2]
  suffStat <- list(dm = MSA, nlev = rep(n_lev, n_variables), adaptDF = FALSE) #dm = dat$x
  
  print(paste("Output will be written to ", getwd(), "/", substring(outpath, 0, nchar(outpath)), "...", sep = ""))
  
  special_parameters_of_pc <- ""                                    # record here if call of pc is changed
  outpath <- paste(outpath, special_parameters_of_pc, sep = "")
  
  sink(paste(outpath, "-pc.txt", sep = ""))
  pc.D <- pc(suffStat, indepTest = ci_test, alpha = alpha, labels = colnames(suffStat$dm) , verbose = TRUE) #p=dim(MSA)[2]
  sink()
  
  if (require(Rgraphviz)) {         ## show estimated CPDAG
    ## plot in R
    garbage <- graphics.off() # clear before plotting
    plot(pc.D, main = capture)
    
    ## plot regularly as .pdf and .ps
    pdf(paste(outpath, ".pdf", sep = "")) 
    plot(pc.D, main = capture) 
    dev.off()
    postscript(paste(outpath, ".ps", sep = ""), paper="special", width=10, height=9)
    plot(pc.D, main = capture) 
    dev.off()
    
    # collect all pdf files in one folder
    outpath_folder <- unlist(strsplit(outpath, "/"))[3]
    outpath_name <- unlist(strsplit(outpath, "/"))[4]
    pdf(paste("../Outputs/_all_graphs_pdf/", outpath_folder, "-", outpath_name, ".pdf", sep = ""), paper="special", width=10, height=9)
    plot(pc.D, main = capture)
    dev.off()
    
    # collect all ps files in one folder
    outpath_folder <- unlist(strsplit(outpath, "/"))[3]
    outpath_name <- unlist(strsplit(outpath, "/"))[4]
    postscript(paste("../Outputs/_all_graphs_ps/", outpath_folder, "-", outpath_name, ".ps", sep = ""), paper="special", width=10, height=9)
    plot(pc.D, main = capture)
    dev.off()
    
    if (output_twopi) {   ## plot as agraph, layout: twopi
      graph <- pc.D @ graph
      graph_ag <- agopen(graph, "pc", layoutType = "twopi")
      ## produce pdf of agraph
      pdf(paste(outpath, "-twopi.pdf", sep = ""), width = 6, height = 4) 
      plot(graph_ag, main = capture)
      dev.off()
    }
    
    if (output_neato) {   ## plot as agraph, layout: neato
      graph <- pc.D @ graph
      graph_ag <- agopen(graph, "pc", layoutType = "neato")
      ## produce pdf of agraph
      pdf(paste(outpath, "-neato.pdf", sep = ""), width = 6, height = 4) 
      plot(graph_ag, main = capture)
      dev.off()
    }
  }
  
  return(pc.D)
}

# replaces the nominal values in non_num_array with integers from one (or 0 if min_zero==TRUE)
replace_with_numbers <- function(non_num_array, replace, allAS, min_zero = TRUE) {
  if (missing(replace)) {
     replace = allAS
  }
  save_rownames <- rownames(non_num_array)
  save_colnames <- colnames(non_num_array)
  pumped_matrix <- apply(non_num_array, 2, function(a,b) c(b,a), replace) # append all appearing values to make sure that in every column, every value appears and thus, every values has the same number in every column
  pumped_matrix_num <- data.matrix(as.data.frame(pumped_matrix))          # data.matrix transforms the values to numericals
  return = pumped_matrix_num[-(1:length(replace)),]
  colnames(return) <- save_colnames
  rownames(return) <- save_rownames
  if (min_zero) {
    return <- return - 1;
  }
  return(return)
}

# Anwendungsbeispiel:
# hydrophob_aliph <- c("I", "L", "M", "V")
# arom <- c("F", "W", "Y")
# pos <- c("H", "K", "R")
# polar_I <- c("S", "T")
# polar_II <- c("D", "E", "N", "Q")
# small <- c("A", "C", "G", "P")
# asp_glu_gln <- c("D", "E", "Q")
# 
# lys_arg <- c("K", "R")
# asn_asp_glu_gln <- c("N", "D", "E", "Q", "B")  # B = Asx = Asn/Asp (N/D)
# # Rest erstmal einzeln
# 
# cluster <- list(ILMV=hydrophob_aliph, FWY=arom, KR=lys_arg, DEQ=asn_asp_glu_gln, H="H", A="A", C="C", G="G", P="P", S="S", "T"="T", UNKNOWN="X", "GAP"="-")
# 
# MSA_mapped_4_13 <- map_AS(MSA, cluster)

# Fasst Aminosäuren zu Gruppen zusammen, die als Vektoren in der Liste cluster angegeben werden
map_aminoacids <- function(MSA, cluster) {
  map <- function(x,  cluster) {
    return(names(cluster)[grep(x, cluster)])
  }
  MSA_mapped <- apply(MSA, 1:2, map, cluster)                                      # langsam
  # MSA_mapped <- structure(sapply(MSA, map, cluster), dim=dim(MSA))               # schneller
  # MSA_mapped <- structure(vapply(MSA, map, character(1), cluster), dim=dim(MSA)) # am schnellsten
  return(MSA_mapped)
}

print.cluster <- function(cluster) {
  if (is.null(cluster)) {
    return("--")
  }
  return_if_cluster <- function(x) {
    if ((nchar(x) > 1) && (x != "UNKNOWN") && (x != "GAP")) {
      return(x)
    }
  }
  
  output <- lapply(names(cluster), return_if_cluster)
  output <- output[!sapply(output, is.null)]
  return(paste(unlist(output), collapse = "-"))
}

remove_sing_col <- function(MSA_in, outpath) {
  verschiedene_AS_an_Pos <- apply(MSA_in, 2, function(vector) {return(unlist(unique(vector)))})
  positions_without_variation <- names(verschiedene_AS_an_Pos)[lapply(verschiedene_AS_an_Pos, length) == 1]
  pos_w_o_var <- sapply(positions_without_variation, as.numeric)
  names(pos_w_o_var) <- NULL
  sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
  cat("positions without variation (removed): ")
  if (!length(pos_w_o_var) == 0) {
    cat(pos_w_o_var)
  } else {
    cat("none")
  }
  cat("\n\n")
  sink()
  
  if (length(pos_w_o_var) > 0) {
    # MSA_out <- MSA_in[,-pos_w_o_var]
    MSA_out <- MSA_in[, !colnames(MSA_in) %in% pos_w_o_var]
  } else {
    MSA_out <- MSA_in
  }
  
  return(MSA_out)
}



