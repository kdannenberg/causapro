library(stringr) # for str_replace_all

get_data_description <- function(protein, type_of_data, subtype_of_data = "", data_set = "", suffix = "") {
  data_description <- ""
  if (is.null(subtype_of_data) || subtype_of_data != "") {
    type_of_data <- paste(type_of_data, subtype_of_data, sep = "-")
  }
  data_description <- paste(protein, type_of_data, sep = "_")
  if (data_set != "") {
    data_description <- paste(data_description, data_set, sep = "_")
  }
  if (suffix != "") {
    data_description <- paste(data_description, suffix, sep = "_")
  }
  return(data_description)
}


read_data <- function(files, path_to_data = "Data/", extension = ".csv", filename, transpose = FALSE, only_cols = NULL) {
  read <- function (file) {
    filename <- paste(path_to_data, file, extension, sep = "")
    data_i = read.csv2(filename, row.names = 1, check.names = FALSE) # if check.names, an X is prepended to numerical column-names
    if (!length(only_cols) == 0) {
      data_i = data_i[,as.character(only_cols)]
    }
    i <- which(files == file)
    if ((length(transpose) > i && transpose[i]) || (length(transpose) == 1 && transpose[1])) {
      data_i <- t(data_i)
    }
    rownames(data_i) <- paste(file, rownames(data_i), sep = "-")
    return(data_i)
  }
  
  if (length(files) == 1) {   # for very mysterious reasons, the matrixs gets transposed by do.call(rbind, data_i) otherwise (rbind(data_i) does NOT transpose the matrix!)
    data <- read(files[1])
  } else {
    data <- do.call(rbind, lapply(files, read))
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
  return(data)
}


# rank_obs_per_pos: should the ranking be done the other way round? 
#   That is, per position, over all observations?
adjust_data <- function(data, type_of_data, rank = FALSE, rank_obs_per_pos = FALSE, remove_low_variance = FALSE,
                        zero_var_fct, min_var = 0.01, mute_plot = TRUE) {
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
  }
  # TODO: statistical test for zero variance
  if (typeof(min_var) == "closure") {
    # remove_low_var_cols <-  nearZeroVar(data, freqCut = 15, saveMetrics = FALSE)
    drop <-  min_var(data)
  } else {
    drop <- which(apply(data, 2, var) <= as.numeric(min_var))
  }
  colors <- rep("#FFFFFF", dim(data)[[1]])
  colors[drop] <- "#000000"
  if (!mute_plot) {
    barplot(apply(data, 2, var), col = colors, las = 2, 
            names.arg = colnames(data))
  }
  # var unter min_var wegschmeißen
  # data <- data[,-drop]
  # data <- subset(data, select = -drop)
  if (length(drop) > 0) {
    # cat("\n")
    print(paste("Removed columns:", paste(colnames(data)[drop], collapse = ", ")))
    data2 <- data[, !names(data) %in% names(drop)]
  }
  return(data)
}

adjust_data_description <- function(data_description, ranked) {
  if (ranked) {
    data_description <- paste0(data_description, "_ranked")
  }
  return(data_description)
}

subtype_of_data_after_adjustment <- function(data, subtype_of_data, rank = FALSE, rank_obs_per_pos = FALSE, remove_low_variance = FALSE,
                                          zero_var_fct, min_var = 0.01) {
  if (rank) {
    if (!rank_obs_per_pos) {
        subtype_of_data <- paste(subtype_of_data, "ranked", sep = "-")
        # type_of_data <- paste(type_of_data, "ranked-pos-per-obs", sep = "-")
    } else {
        subtype_of_data <- paste(subtype_of_data, "ranked-obs-per-pos", sep = "-")
    }
  }
  
  # TODO: statistical test for zero variance
  if (subtype_of_data == "") {
    delimiter = ""
  } else {
    delimiter = "-"
  }
  if (!typeof(min_var) == "closure") {
    if ((min_var > 0) || (min_var < 0)) {
      subtype_of_data <- paste0(subtype_of_data, delimiter, "var>", min_var)
    }
  } else {
    subtype_of_data <- paste(subtype_of_data, delimiter, "var>fct", sep = "-")
  }
  
  
  
  return(subtype_of_data)
}



get_outpath <- function(protein, type_of_data, subtype_of_data = "", data_set = "", suffix = "", alpha, min_pos_var, only_cols_label, 
                        pc_solve_conflicts, pc_u2pd, pc_conservative, pc_maj_rule, file_separator = "/", 
                        filename_suffix, main_dir = "Outputs") {   ## last two options: only for get_old_outpath
  dir_1 <- protein
  dir_2 <- type_of_data
  # if (subtype_of_data != "")
  #   dir_3 <- paste(type_of_data, subtype_of_data, sep = "-")
  # else {
  #   dir_3 <- type_of_data
  # }
  dir_3 <- pastes(type_of_data, paste(subtype_of_data, collapse = "+"), data_set, sep = "-")
  
  dir_4 <- paste0(get_data_description(protein = protein, type_of_data = type_of_data, 
                                       subtype_of_data = paste(subtype_of_data, collapse = "+"), 
                                       data_set = data_set, suffix = suffix), "_alpha=", alpha)
  
  # if (min_pos_var > 0) {
  #   dir_4 <- paste0(dir_min_pos_var, "_mv=", min_pos_var)
  # }
 
  
  output_dir <- paste(main_dir, dir_1, dir_2, dir_3, dir_4, sep = file_separator) 
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }
  
  filename <- dir_4
  
  if (!missing(filename_suffix)) {
    filename <- paste0(filename, filename_suffix)
  } else {
    if (pc_solve_conflicts) {
      filename <- paste0(filename, "_sc")
    } 
    if (pc_conservative) {
      filename <- paste0(filename, "_cons")
    }
    if (pc_maj_rule) {
      filename <- paste0(filename, "_maj")
    }
    if (!(pc_solve_conflicts || pc_conservative)) {
      filename <- paste(filename, substr(pc_u2pd, 1, 3), sep = "_")
    }
  }
  
  
  
  return(paste(output_dir, filename, sep = file_separator))
}

# gab früher den old_outpath zurück, 
# jetzt das Verzeichnis in dem es die Datei gibt, oder NULL, wenn es sie nicht gibt.
# TODO: rename: get_file (oder so) / get_outpath_where_file_exists
get_old_outpath <- function(outpath, suffix) {
  if (file.exists(paste0(outpath, suffix))) {
    return(paste0(outpath, suffix))
  } 
  
  if (endsWith(outpath, "_sc_maj")) {
    old_outpath <- str_replace(outpath, "_sc_maj", "_rel")
    if (file.exists(paste0(old_outpath, suffix))) {
      return(paste0(old_outpath, suffix))
    }
  } # else {
  dirs <- str_split(outpath, "/", simplify = TRUE)
    
  protein <- dirs[2]
  type_of_data <- dirs[3]
  extension <- sub(type_of_data, "", dirs[4])
  if (is.null(extension) || extension == "NULL") {
    subtype_of_data <- ""
  } else {
    extension <- substr(extension, 2, nchar(extension))
    if (grepl(pattern = "-", extension)) {
      subtype_of_data <- str_split(extension, "-", simplify = TRUE)[1]
    } else {
      subtype_of_data <- ""
    }
    rest_of_extension <- sub(subtype_of_data, "", extension)
    if (substr(rest_of_extension, 1, 1) == "-") {
      rest_of_extension <- substr(rest_of_extension, 2, nchar(extension))
    }
    type_of_data <- pastes(type_of_data, rest_of_extension, sep = "-")
  }
  
  start_of_alpha <- gregexpr(pattern ='alpha', outpath)[[1]][1]
  end_that_starts_with_first_alpha <- substr(outpath, start_of_alpha + 6, nchar(outpath))
  first_slash <- gregexpr(pattern ='/', end_that_starts_with_first_alpha)[[1]][1]
  alpha <- substr(end_that_starts_with_first_alpha, 1, first_slash - 1)
  # rest_of_extension <- str_replace(extension, ", "")
  
  filename <- dirs[length(dirs)]
  start_of_suffix <- gregexpr(pattern ='alpha', filename)[[1]][1] + 6 + nchar(alpha)
  filename_suffix <- substr(filename, start_of_suffix, nchar(filename)) 
  
  old_outpath <- get_outpath(protein = protein, type_of_data = type_of_data, subtype_of_data = subtype_of_data, alpha = alpha,
                             filename_suffix = filename_suffix, main_dir = "Outputs_2017-09-14")
    # old_outpath <- paste(old_outpath, paste(dirs[c(5,6)], collapse = "/"), sep = "/")
  # }
  
  old_outpath <- str_replace(old_outpath, "_sc_maj", "_rel")
  if (file.exists(paste0(old_outpath, suffix))) {
    return(paste0(old_outpath, suffix))
  }
  
  return(NULL)
}

# TODO: old_outpath an geänderte get_old_outpath-Funktion anpassen
compute_if_not_existent <- function(filename, FUN, obj_name = "data", compute_anew = FALSE, 
                                    loaded_object_ok_fun = function(obj) return(TRUE), 
                                    try_old_outpath = TRUE) {
  filename <- paste(filename, ".RData", sep = "")
  if (file.exists(filename) && !compute_anew) {
    # rm(list = ls()[which(ls() == obj_name)])
    try(load(filename))
    if (!exists(obj_name)) {
      warning(paste("File not loadable or did not contain an object of name", obj_name, "!"))
    } else if (loaded_object_ok_fun(get(obj_name))) {
      warning("The loaded object did not fit!")
    } else {
      print(paste(obj_name, " object loaded from ", filename, ".", sep = ""))
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
  save(list = obj_name, file = filename)
  print(paste(obj_name, "object computed and saved."))
  # }
  return(get(obj_name))
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
    # MSA <- readAlignment(filename_fasta)
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
    sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
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
  out_file <- paste(outpath, "-info.txt", sep = "")
  sink(file = out_file)
  cat(parameters_for_info)
  sink()
}

print_pc_results_to_info_file <- function(outpath, pc) {
  # print(paste(outpath, "-info.txt", sep = ""))
  sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
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

outpath_for_ida <- function(outpath, direction, weight_effects_on_by, option_nr, neg_effects, perturbated_position, amplification_exponent, 
                            amplification_factor, no_colors, rank_effects, effect_to_color_mode) {
  outpath <- paste0(outpath, "-total_effects_(", neg_effects, ")")
  
  out_file <- outpath
  if (option_nr != "") {
    out_file <- paste0(out_file, "_#", option_nr)
  }
  
  out_file <- paste0(out_file, "_", direction, "_pos_", perturbated_position)
  
  if (effect_to_color_mode == "opacity") {
    out_file <- paste0(out_file, "-opac")
  }
  if (rank_effects && !(no_colors)) {
    out_file <- paste0(out_file, "-ranked")
  } else {
    if (amplification_exponent != 1) {
      out_file <- paste0(out_file, "-ampl_exp=", amplification_exponent)
    } 
    if (amplification_factor) {
      out_file <- paste0(out_file, "-ampl_fac")
    }
  }
  if (no_colors) {
    out_file <- paste0(out_file, "-bw")
  }
  out_file <- paste0(out_file, ".pml")
  
  # print(out_file)
  # if (!mute_all_plots) {    ### HÄÄÄÄ???
    return(out_file)
  # }
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