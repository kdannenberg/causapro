library(stringr) # for str_replace_all


get_data_description <- function(protein, type_of_data, subtype_of_data = "", data_set = "",
                                 only_rows_cols_label, pre_fun, cor_cov_FUN, type_of_variables,
                                 suffix = "") {

  # if (!missing(only_cols_label) && !is.null(only_cols_label) && only_cols_label != "") {
  #   only_cols_label <- paste0("cols=", pre_fun)
  # } # TODO: auch für rows! ## ODER GING DAS NICHT ANDERS?
  if (!missing(pre_fun) && !is.null(pre_fun) && pre_fun != "") {
    pre_fun <- paste0("preFUN=", pre_fun)
  }
  if (typeof(cor_cov_FUN) == "closure") {
    cor_cov_FUN <- deparse(substitute(cor_cov_FUN))
  } else if (is.null(cor_cov_FUN) || missing(cor_cov_FUN)) {
    cor_cov_FUN <- ""
  }
  if (cor_cov_FUN != "") {
    cor_cov_FUN <- paste0("corFUN=", cor_cov_FUN)
  }

  data_description <- pastes(protein, type_of_data, paste(subtype_of_data, collapse = "+"),
                             data_set, suffix, only_rows_cols_label, pre_fun, cor_cov_FUN,
                             type_of_variables, sep = "-")



  # data_description <- ""
  # if (is.null(subtype_of_data) || subtype_of_data != "") {
  #   type_of_data <- paste(type_of_data, subtype_of_data, sep = "_")
  # }
  # data_description <- paste(protein, type_of_data, sep = "_")
  # if (data_set != "") {
  #   data_description <- paste(data_description, data_set, sep = "_")
  # }
  # if (suffix != "") {
  #   data_description <- paste(data_description, suffix, sep = "_")
  # }
  return(data_description)
}

# parameter filename was redundant?
read_data <- function(files, read_fct, path_to_data = "Data/", extension = ".csv", #filename,
                      transpose = FALSE, every_n_th_row = 1) {
  read <- function (file, every_n_th_row_ = 1) {
    # if (file.exists(paste0(path_to_data, file, ".RData"))) {
    #   load(paste0(path_to_data, file, ".RData"))
    #   print(paste(paste0(path_to_data, file, ".RData"), "loaded instead of", paste0(path_to_data, file, extension)))
    # } else {
      if (!(missing(every_n_th_row_)) && !is.null(every_n_th_row_) && !(every_n_th_row_ == 1)
          && file.exists(filename = paste0(path_to_data, file, "-every_", every_n_th_row_, "th", extension))) {
          # data_i = read.csv2(filename, row.names = 1, check.names = FALSE) # if check.names, an X is prepended to numerical column-names
          data_i = read_fct(filename = filename)
          every_n_th_row_ = 1
      } else {
          filename <- paste0(path_to_data, file, extension)
          # data_i = read.csv2(filename, row.names = 1, check.names = FALSE) # if check.names, an X is prepended to numerical column-names
          data_i = read_fct(filename = filename)
      }
      i <- which(files == file)
      if ((length(transpose) > i && transpose[i]) || (length(transpose) == 1 && transpose[1])) {
        data_i <- t(data_i)
      }
      rownames(data_i) <- paste(file, rownames(data_i), sep = ".")
      # save(data_i, file = paste0(path_to_data, file, ".RData"))

      if (!(missing(every_n_th_row_)) && !is.null(every_n_th_row_) && !(every_n_th_row_ == 1)) {
        data_i <- data_i[seq(from = 1, to = dim(data_i)[1], by = every_n_th_row_),]
      }

    # }
    return(data_i)
  }

  # debug(read)

  read_or_get_data <- function(file, every_n_th_row) {
    func <- function_set_parameters(read, parameters = list(file = file, every_n_th_row = every_n_th_row))
    filename <- paste0(path_to_data, file)
    if (!(missing(every_n_th_row)) && !is.null(every_n_th_row) && !(every_n_th_row == 1)) {
      filename = paste0(filename, "-every_", every_n_th_row, "th")
    }
    data <- compute_if_not_existent(filename = filename,
                                    FUN = func, obj_name = "data_i")
  }

  # debug(read_or_get_data)

  if (length(files) == 1) {   # for very mysterious reasons, the matrixs gets transposed by do.call(rbind, data_i) otherwise (rbind(data_i) does NOT transpose the matrix!)
    # data <- read(files[1])
    data <- read_or_get_data(files[1], every_n_th_row = every_n_th_row)
  } else {
    data <- do.call(rbind, lapply(files, read_or_get_data))
  }


  # data <- matrix()
  # for (file in files) {
  #   filename <- paste(path_to_data, source_of_data, extension, sep = "")
  #   data_i = read.csv2(filename, row.names = 1, check.names=FALSE) # if check.names, an X is prepended to numerical column-names
  #   i <- which(files == file)
  #   if ((length(transpose) > i && transpose[i]) || (length(transpose) == 1 && transpose[1])) {
  #     data_i <- t(data_i)
  #   }
  #   data <- rbind(data, data_i)
  # }
  return(as.matrix(data))
}

analyse_data <- function(data, int_pos, min_var) {
  if (!missing(min_var)) {
    which(apply(data, 2, var)[as.character(int_pos)] > min_var)
  }
  min(apply(data, 2, var)[as.character(int_pos)])
}

# FKT für p38g (z.B. p38g_NMR-Mut_inact.L77V)
intervention_positions_from_rownames <- function(rownames) {
  return(str_extract(lapply(rownames, function(x) {str_split(x,"\\.")[[1]][2]}), "[0-9]+"))
}


# rank_obs_per_pos: should the ranking be done the other way round?
#   That is, per position, over all observations?
adjust_data <- function(data, data_descr = "", type_of_data, rank = FALSE, rank_obs_per_pos = FALSE,
                        only_cols = NULL, remove_cols = NULL,
                        only_rows = NULL, remove_rows = NULL,
                        rows_cols_subset_grep = TRUE,
                        # every_n_th_row = 1,
                        remove_low_variance = FALSE, zero_var_fct, min_var = 0.01,
                        pre_fun_on_data,
                        keep_quadratic = FALSE, mute_plot = TRUE,
                        adjust_colnames = TRUE) {

  if (adjust_colnames) {
    remove_three_letter_aa_prefixes <- function(string) {
      if (grepl('^[ARG|ASN|ASP|CYS|GLN|GLU|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL]', string)) {
        return(substring(string, 4))
      } else {
        return(string)
      }
    }
    colnames(data) <- sapply(colnames(data), remove_three_letter_aa_prefixes, USE.NAMES = FALSE)
  }

  if (length(only_cols) > 0) {
    # grep the right cols
    data <- subset_data(data = data, selection = only_cols, grep_selection = rows_cols_subset_grep,
                        remove = FALSE, cols = TRUE)
  }
  if (length(remove_cols) > 0) {
    data <- subset_data(data = data, selection = remove_cols, grep_selection = rows_cols_subset_grep,
                        remove = TRUE, cols = TRUE)
  }
  if (length(only_rows) > 0) {
    # grep the right cols
    data <- subset_data(data = data, selection = only_rows, grep_selection = rows_cols_subset_grep,
                        remove = FALSE, cols = FALSE)
  }
  if (length(remove_rows) > 0) {
    data <- subset_data(data = data, selection = remove_rows, grep_selection = rows_cols_subset_grep,
                        remove = TRUE, cols = FALSE)
  }

  # rank data i desired
  if (rank) {
    if (!rank_obs_per_pos) {
      if (!missing(data)) {
        data <- t(apply(data, 1, rank))  # observationsweise (über alle Positionen)
      }
      if (!missing(type_of_data)) {
        type_of_data <- paste(type_of_data, "ranked", sep = "-")
        # type_of_data <- paste(type_of_data, "ranked-pos-per-obs", sep = "-")
      }
    } else {
      if (!missing(data)) {
        data <- cbind(apply(data, 2, rank)) # positionsweise (über alle Obeservationen)
      }
      if (!missing(type_of_data)) {
        type_of_data <- paste(type_of_data, "ranked-obs-per-pos", sep = "-")
      }
    }
    singular <- apply(data, 2, function(v) length(unique(v)) == 1)
    if (length(singular) > 0) {
      # cat("\n")
      print(paste("Removed columns:", paste(colnames(data)[singular], collapse = ", ")))
      if (keep_quadratic) {
        data <- data[!singular, !singular]
      } else {
        data <- data[, !singular]
      }
    }
  }

  # apply pre_fun_on_data on data
  if (!missing(pre_fun_on_data) && !(is.null(pre_fun_on_data)) && nchar(pre_fun_on_data) > 0) {
    data <- get(pre_fun_on_data)(data)
  }

  # Remove cols with less than min_var
  # TODO: statistical test for zero variance
  if (typeof(min_var) == "closure") {
    # remove_low_var_cols <-  nearZeroVar(data, freqCut = 15, saveMetrics = FALSE)
    drop <-  min_var(data)
  } else {
    drop <- which(apply(data, 2, var) < as.numeric(min_var)) # fkt nicht zuverlässig bei der Grenze (vllt wegen Rundungsfehlern)
  }
  colors <- rep("#FFFFFF", dim(data)[[2]])
  colors[drop] <- "#000000"
  if (!mute_plot) {
    caption <- pastes(data_descr, ifelse(rank, "ranked", ""),
                     ifelse(rank, ifelse(rank_obs_per_pos, "per_pos", "per_obs"), ""),
                     pre_fun_on_data,
                     ifelse(min_var > 0, paste0("var_cutoff=", min_var), ""),
                     ifelse(keep_quadratic, "kept_quadratic", ""), sep = ", ")
    barplot(apply(data, 2, var), col = colors, las = 2,
            names.arg = colnames(data),
            main = strwrap(caption, width = 50))
  }
  # var unter min_var wegschmeißen
  # data <- data[,-drop]
  # data <- subset(data, select = -drop)
  if (length(drop) > 0) {
    # cat("\n")
    print(paste("Removed columns:", paste(colnames(data)[drop], collapse = ", ")))
    if (keep_quadratic) {
      data <- data[!colnames(data) %in% names(drop), !colnames(data) %in% names(drop)]
    } else {
      data <- data[, !colnames(data) %in% names(drop)]
    }
  } else {
    print(paste("No columns removed."))
  }

  print(paste("Next lowest variance:", min(apply(data, 2, var))))

  return(data)
}


subset_data <- function(data, selection, grep_selection, remove = FALSE, cols = TRUE) {
  if (cols) {
    name_function = colnames
  } else {
    name_function = rownames
  }

  if (grep_selection) {
    selection_ind <- sapply(selection, function(x) which(grepl(as.character(x), name_function(data))))
    # only_cols <- names(sapply(only_cols, function(x) which(grepl(as.character(x), colnames(data)))))
    not_found <- selection[which(lapply(selection_ind, length) == 0)]
    if (length(not_found) > 0) {
      warning(paste("No columns containing", paste(not_found, collapse = ", "), "found!"))
    }
    selection_ind <- selection_ind[which(!(lapply(selection_ind, length) == 0))]
    selection <- colnames(data)[unlist(selection_ind)]
  } else {
    not_found <- setdiff(as.character(selection), colnames(data))
    if (length(not_found) > 0) {
      warning(paste("Column(s)", paste(not_found, collapse = ", "), "not found!"))
      selection <- setdiff(selection, not_found)
    }
  }

  if (remove) {
    lines_to_remove <- as.character(selection)
  } else {
    lines_to_remove <- setdiff(colnames(data), as.character(selection))
  }
  if (length(selection) > 0) {
    if (cols) {
      print(paste("Column(s)", paste(lines_to_remove, collapse = ", "), "removed by only_cols/remove_cols."))
    } else {
      print(paste("Row(s)", paste(lines_to_remove, collapse = ", "), "removed by only_rows/remove_rows."))
    }
    lines_to_remove_ind <- which(colnames(data) %in% as.character(lines_to_remove))
    if (cols) {
      data = data[, -lines_to_remove_ind]
    } else {
      data = data[-lines_to_remove_ind, ]
    }
  }
  return(data)
  # data = data[,as.character(selection)]

  # if (rows_cols_subset_grep) {
  #   remove_cols_ind <- sapply(remove_cols, function(x) which(grepl(as.character(x), colnames(data))))
  #   # only_cols <- names(sapply(only_cols, function(x) which(grepl(as.character(x), colnames(data)))))
  #   cols_not_found <- names(remove_cols_ind[which(lapply(remove_cols_ind, length) == 0)])
  #   if (length(cols_not_found) > 0) {
  #     warning(paste("No columns containing", paste(cols_not_found, collapse = ", "), "found!"))
  #   }
  #   remove_cols_ind <- remove_cols_ind[which(!(lapply(remove_cols_ind, length) == 0))]
  #   remove_cols <- colnames(data)[unlist(remove_cols_ind)]
  # }
  # data = data[, -(which(colnames(data) %in% as.character(remove_cols)))]

}

adjust_data_description <- function(data_description, ranked = FALSE) {
  if (ranked) {
    data_description <- paste0(data_description, "_ranked")
  }
  return(data_description)
}

subtype_of_data_after_adjustment <- function(data, subtype_of_data, rank = FALSE, rank_obs_per_pos = FALSE, remove_low_variance = FALSE,
                                          zero_var_fct, min_var = 0.01) {
  if (rank) {
    if (!rank_obs_per_pos) {
        subtype_of_data <- pastes(subtype_of_data, "ranked", sep = "_")
        # type_of_data <- pastes(type_of_data, "ranked-pos-per-obs", sep = "-")
    } else {
        subtype_of_data <- pastes(subtype_of_data, "ranked-obs-per-pos", sep = "_")
    }
  }

  # TODO: statistical test for zero variance
  if (!typeof(min_var) == "closure") {
    if ((min_var > 0) || (min_var < 0)) {
      # subtype_of_data <- paste0(subtype_of_data, delimiter, "var>", min_var)
      subtype_of_data <- pastes(subtype_of_data, paste0("var>", min_var), sep = "_")
    }
  } else {
    subtype_of_data <- pastes(subtype_of_data, "var>fct", sep = "_")
    # subtype_of_data <- paste(subtype_of_data, delimiter, "var>fct", sep = "-")
  }



  return(subtype_of_data)
}


get_only_cols_label <- function(only_cols, remove_cols, only_rows, remove_rows,
                                rows_cols_subset_grep) {
  only_cols_label <- ""
  # for (param in c("only_cols", "remove_cols", "only_rows", "remove_rows")) {
  #   if (is.null(get(param))) {
  #     assign(param, "")
  #   }
  # }
  cols_str_only <- paste(only_cols, collapse = "+")
  cols_str_rem <- paste(remove_cols, collapse = "-")
  cols_str = ""
  if (nchar(cols_str_only) + nchar(cols_str_rem) > 0) {
    cols_str = pastes(pastes("c", cols_str_only, sep = "+"), cols_str_rem, sep = "-")
  }
  rows_str_only <- paste(only_rows, collapse = "+")
  rows_str_rem <- paste(remove_rows, collapse = "-")
  rows_str = ""
  if (nchar(rows_str_only) + nchar(rows_str_rem) > 0) {
    rows_str = pastes(pastes("r", rows_str_only, sep = "+"), rows_str_rem, sep = "-")
  }
  return(pastes(cols_str, rows_str, sep = "-"))
}


# loaded_object_ok_fun -> fun_loaded_object_ok, dabei einmal invertieren!
# (Jetzt wirklich checken, ob es ok ist, nicht, ob es nicht ok ist!!)
# TODO: old_outpath an geänderte get_old_outpath-Funktion anpassen
compute_if_not_existent <- function(filename, FUN, obj_name = "data", compute_anew = FALSE,
                                    fun_loaded_object_ok = function(obj) return(TRUE),
                                    try_old_outpath = TRUE) {
  # debug(FUN)
  filename <- paste(filename, ".RData", sep = "")
  if (file.exists(filename) && !compute_anew) {
    # rm(list = ls()[which(ls() == obj_name)])
    try(load(filename))
    if (!exists(obj_name)) {
      warning(paste("File not loadable or did not contain an object of name", obj_name, "!"))
    } else if (!fun_loaded_object_ok(get(obj_name))) {
      warning("The loaded object did not fit!")
    } else {
      print(paste(obj_name, " object loaded from ", filename, ".", sep = ""))
      # TODO speichern, falls aus altem Outpath geladen
      return(get(obj_name))
    }
  # } else {
  #   # if (!missing(old_outpath)) {
  #   old_outpath <- get_old_outpath(filename)
  #   if (try_old_outpath && file.exists(old_outpath) && !compute_anew) {
  #     # rm(list = ls()[which(ls() == obj_name)])
  #     load(old_outpath)
  #     if (!exists(obj_name)) {
  #       warning(paste("The file did not contain an object of name", obj_name, "!"))
  #     } else if (loaded_object_ok_fun(get(obj_name))) {
  #       warning("The loaded object did not fit!")
  #     } else {
  #       print(paste("Old ", obj_name, " object loaded from ", old_outpath, ".", sep = ""))
  #       save(list = obj_name, file = filename) # Objekt an neuem Ort speichern
  #       return(get(obj_name))
  #     }
  #   }
  #   # }
  } # else {
  print(paste("Computing", obj_name, "object."))
  assign(obj_name, FUN())
  create_parent_directory_if_necessary(filename)
  save(list = obj_name, file = filename)
  print(paste(obj_name, "object computed and saved."))
  # }
  return(get(obj_name))
}

compute_data_from_alignment <- function(alignment, filename) {
  if(grepl("bin_approx", filename)) {
    return(alignment_to_binary_matrix(alignment = alignment))
  }
}

readAlignment <- function(filename) {
  if (file.exists(paste(filename, ".RData", sep = ""))) {
    filename <- paste(filename, ".RData", sep = "")
    load(filename)
    if (!exists("MSA")) {
      print("The file did not contain an object of name 'MSA'!")
    } else {
      print(paste("Alignment loaded from ", filename, ".", sep = ""))
    }
  } else {
    filename_fasta <- paste(filename, ".fasta", sep = "")
    print(paste("Loading alignment from ", filename_fasta, ".", sep = ""))
    MSA <- readAlignment_fasta(filename_fasta)
    colnames(MSA) <- seq(1:dim(MSA)[2])
    save(MSA, file = paste(filename, ".RData", sep = ""))
    print("Loading done; Object saved for later.")
  }
  return(MSA)
}

readAlignment_fasta <- function(from_file) {
  MSA <- read.alignment(file = from_file, format = "fasta", forceToLower = FALSE)  # from seqinr
  MSA_list <- sapply(MSA$seq, strsplit, "")
  MSA_mat <- t(sapply(MSA_list, rbind))
  rownames(MSA_mat) <- MSA$nam
  return(MSA_mat)
}

allAS = c("-", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "X")

remove_gaps <- function(MSA_in, threshold, n_lev, allAS, outpath) {
  gap_freq <- t(apply(MSA_in, 2, function(x) (table(factor(x, levels = allAS)) / length(x))["-"]))

  if (!is.null(outpath)) {
    sink(file = paste(outpath, "_info.txt", sep = ""), append = TRUE)
    cat("positions without more than 'remove_cols_gaps_threshold'=")
    cat(remove_cols_gaps_threshold)
    cat(" gaps (removed): ")
    if (!length(gap_freq[gap_freq > threshold]) == 0) {
      cat(colnames(MSA_in)[gap_freq > threshold])
    } else {
      cat("none")
    }
    cat("\n")
    sink()
  }

  MSA_out <- MSA_in[, gap_freq < threshold]
}

rem_cols_by_colname <- function(data, remove) {
  return(data[, !(colnames(data) %in% remove)])
}

## previously: "caption"
get_caption <- function(protein, data, alpha, min_pos_var, chars_per_line = 50) {
  # par_string = paste("protein: ", protein, ", data: ", data, ", alpha: ", alpha, ", var > ", min_pos_var, sep = "")
  par_string = paste(data, " - alpha: ", alpha, ", var > ", min_pos_var, sep = "")
  caption = strwrap(par_string, width = chars_per_line)
  return(caption)
}

library(stringr)
parameters_for_info_file <- function(protein, type_of_data, alpha, position_numbering, only_cols, coloring, colors, outpath) {
  par_string <- paste("protein: ", protein, ", data: ", type_of_data, ", alpha: ", alpha, sep = "")
  if (is.null(colors)) {
    colorstring <- paste("coloring:", coloring)
  } else {
    colorstring <- paste("coloring: ", coloring, ", colors:", paste(names(colors), colors, collapse="; ", sep = " - "), sep = "")
  }
  info_str <- paste(par_string, ", position_numbering: ", position_numbering, ", only_cols: ", only_cols, ", outpath: ", outpath, sep = "")
  parameters_for_info <- paste(str_replace_all(info_str, pattern = ", ", ",\n"), "\n\n", sep = "")
  return(parameters_for_info)
}

parameters_to_info_file <- function(parameters_for_info, outpath) {
  out_file <- paste(outpath, "_info.txt", sep = "")
  sink(file = out_file)
  cat(parameters_for_info)
  sink()
}

print_pc_results_to_info_file <- function(outpath, pc) {
  # print(paste(outpath, "-info.txt", sep = ""))
  sink(file = paste(outpath, "_info.txt", sep = ""), append = TRUE)
  # print("_________________")
  print("RESULT:")
  print(pc @ graph)
  print(conv_to_r(pc@graph, type_of_graph = "pdag"))
  cat("\n")
  print("COMPUTATION TIME:")
  print(proc.time())
  cat("\n")
  sink()
}

get_conservation <- function(measure, protein) {
  type_of_data <- paste0("D", measure)
  data_description <- get_data_description(protein = protein, type_of_data = type_of_data)
  data <- read_data(files = data_description)

  if (dim(data)[1] > 0) {
    data <- apply(data, 2, sum)
  }

  return(data)
}

## set_alpha_in_outpath <- function(outpath, alpha) {
##   pattern <- "alpha=[0-9]*.[0-9]+|alpha=[0-9]*e-[0-9]+"
##   return(str_replace_all(outpath, pattern, paste0("alpha=", alpha)))
## }

set_in_outpath <- function(outpath, parameter_name, parameter_value) {
  outpath <- str_replace_all(outpath, paste0(parameter_name, "=[^_]*_"), paste0(parameter_name, "=", parameter_value, "_"))
  return(outpath)
}

overwrite_in_outpath <- function(outpath, pattern, replacement) {
  outpath <- str_replace_all(outpath, pattern, replacement)
  outpath <- str_replace_all(outpath, "[-]{2}", "-")
  outpath <- str_replace_all(outpath, "[_-]{2}", "_")
  return(outpath)
}

create_parent_directory_if_necessary <- function(outpath, file_separator = "/") {
  directories <- strsplit(outpath, file_separator)
  # filename <- directories[[1]][length(directories[[1]])]
  output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
    print("Directory created.")
  }
}

