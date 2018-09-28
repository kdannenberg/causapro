# paste if non of the arguments is empty, return empty string otherwise
pastes_parvalue <- function(parameter_name, value, sep = "=") {
  if (is.null(parameter_name) || is.null(value) ||
      missing(parameter_name)  || missing(value) ||
      nchar(parameter_name) == 0 || nchar(value) == 0) {
    return("")
  } else {
    return(paste(parameter_name, value, sep = sep))
  }
}

empty_str_fct <- function() {return("")}


# Projekt erstmal wieder aufgegeben.
# get_outpath <- function(protein, type_of_data, subtype_of_data = "", data_set = "", suffix = "",
#                         alpha, min_pos_var, only_rows_cols_label = "", every_n_th_row,
#                         pre_fun_on_data, cor_cov_FUN = "", type_of_variables = "",
#                         file_separator = "/", filename_suffix, main_dir = "Outputs",
#                         plot_only_subgraphs = NULL, coloring = NULL, graph_layout = NULL) {
#   outpath <- ""
#   # quick and dirty:
#   outpath_data <- get_outpath_data(protein = protein, type_of_data = type_of_data, subtype_of_data = "", data_set = "", suffix = "",
#                                    alpha, min_pos_var, only_rows_cols_label = "", every_n_th_row,
#                                    pre_fun_on_data, cor_cov_FUN = "", type_of_variables = "",
#                                    file_separator = "/", filename_suffix, main_dir = "Outputs")
#
#   outpath <- paste0(outpath, outpath_data)
#   if (!all(sapply(list(...), is.null))) { #bisschen doof, die da alle reinzustecken, oder?
#     # Kommt da in dem Fall nicht sowieso "" raus? Nein, "graph"...
#     outpath_pc <- get_outpath_pc_(...)
#     outpath <- c(outpath, outpath_pc)
#   }
#   outpath_pc_graph <- get_outpath_pc_graph(prefix = outpath_pc, plot_only_subgraphs = plot_only_subgraphs,
#                                            coloring = coloring, graph_layout = graph_layout)
#
#   return(with_graph)
# }



get_outpath_data <- function(protein, type_of_data, subtype_of_data = "", data_set = "", suffix = "",
                             min_pos_var, only_rows_cols_label = "", every_n_th_row,
                             pre_fun_on_data, cor_cov_FUN = "", type_of_variables = "",
                             file_separator = "/", filename_suffix, main_dir = "Outputs") {   ## last two options: only for get_old_outpath
  dir_1 <- protein
  dir_2 <- type_of_data
  # if (subtype_of_data != "")
  #   dir_3 <- paste(type_of_data, subtype_of_data, sep = "-")
  # else {
  #   dir_3 <- type_of_data
  # }
  # dir_3 <- pastes(type_of_data, paste(subtype_of_data, collapse = "+"), data_set, pre_fun_on_data, sep = "_")

  dir_4 <- paste0(get_data_description(protein = protein, type_of_data = type_of_data,
                                       subtype_of_data = paste(subtype_of_data, collapse = "+"),
                                       data_set = data_set, only_rows_cols_label = only_rows_cols_label,
                                       pre_fun = pre_fun_on_data, cor_cov_FUN = cor_cov_FUN, type_of_variables = type_of_variables,
                                       suffix = suffix))

  # if (min_pos_var > 0) {
  #   dir_4 <- paste0(dir_min_pos_var, "_mv=", min_pos_var)
  # }




  if (!missing(filename_suffix)) { # vllt so für get_old_outpath?
    filename <- paste0(filename, filename_suffix)
    dir_4 <- paste0(get_data_description(protein = protein, type_of_data = type_of_data,
                                         subtype_of_data = paste(subtype_of_data, collapse = "+"),
                                         data_set = data_set, only_cols_label = only_cols_label,
                                         suffix = suffix))
  } else {
    dir_4 <- paste0(get_data_description(protein = protein, type_of_data = type_of_data,
                                         subtype_of_data = paste(subtype_of_data, collapse = "+"),
                                         data_set = data_set, only_rows_cols_label = only_rows_cols_label,
                                         pre_fun = pre_fun_on_data, cor_cov_FUN = cor_cov_FUN,
                                         type_of_variables = type_of_variables,
                                         suffix = suffix))
  }
  # else {
  #   filename <- pastes(filename, only_cols_label, sep = "-")
  #   if (typeof(cor_cov_FUN) == "closure") {
  #     cor_cov_FUN <- deparse(substitute(cor_cov_FUN))
  #   }
  #   if (is.null(cor_cov_FUN) || missing(cor_cov_FUN)) {
  #     cor_cov_FUN <- ""
  #   }
  #   if (cor_cov_FUN != "") {
  #     filename <- paste0(filename, "-corFUN=", cor_cov_FUN)
  #   }
  # }

  output_dir <- paste(main_dir, dir_1, dir_2, dir_4, sep = file_separator)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }

  filename <- dir_4

  return(paste(output_dir, filename, sep = file_separator))
}


# TODO: Function add_directory(path, dir_name)
add_directory_to_path <- function(path, dir_name, append_to_cur_file_name = TRUE, file_separator = "/") {
  if (append_to_cur_file_name) {
    directories <- strsplit(path, file_separator)
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
    last_dir_name <- pastes(directories[[1]][length(directories[[1]])-1], dir_name, sep = "_")
    outpath <- paste(output_dir, last_dir_name, last_dir_name, sep = file_separator)
  } else {
    directories <- strsplit(path, file_separator)
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
    last_dir_name <- pastes(directories[[1]][length(directories[[1]])-1], dir_name, sep = "-")

    outpath <- paste(output_dir, dir_name, last_dir_name, sep = file_separator)
  }
  return(outpath)
}

get_outpath_suff_pc <- function(alpha, pc_indepTest = "", cor_cov_FUN, pc_solve_conflicts,
                                pc_u2pd, pc_conservative, pc_maj_rule) {
  # ab hier: pc # soll das überhaupt schon dran?
  filename_pc <- paste0("PC-alpha=", alpha)
  if (typeof(pc_indepTest) == "closure") {
    pc_indepTest <- deparse(substitute(pc_indepTest))
  }
  if (!(missing(pc_indepTest)) && !(is.null(pc_indepTest)) && pc_indepTest != "") {
    filename_pc <- pastes(filename_pc, "test=", pc_indepTest, sep = "-")
  }
  if (pc_solve_conflicts) {
    filename_pc <- pastes(filename_pc, "sc", sep = "-")
  }
  if (pc_conservative) {
    filename_pc <- pastes(filename_pc, "cons", sep = "-")
  }
  if (pc_maj_rule) {
    filename_pc <- pastes(filename_pc, "maj", sep = "-")
  }
  if (!(pc_solve_conflicts || pc_conservative)) {
    filename_pc <- pastes(filename_pc, substr(pc_u2pd, 1, 3), sep = "-")
  }
  # filename <- pastes(filename, filename_pc, sep = "_")

  return(filename_pc)
}

get_outpath_pc <- function(outpath_data, file_separator, alpha, pc_indepTest = "", cor_cov_FUN = "",
                           pc_solve_conflicts = "", pc_u2pd = "",
                           pc_conservative = "", pc_maj_rule = "") {

  # directories <- strsplit(outpath_data, file_separator)
  # output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = file_separator, sep = file_separator)
  #
  outpath_suff_pc <- get_outpath_suff_pc(alpha = alpha, pc_indepTest = pc_indepTest, cor_cov_FUN = cor_cov_FUN,
                                         pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
                                         pc_conservative = pc_conservative, pc_maj_rule = pc_maj_rule)
  # last_dir_name <- pastes(directories[[1]][length(directories[[1]])-1], outpath_suff_pc, sep = "_")
  #
  # outpath_pc <- paste(output_dir, last_dir_name, last_dir_name, sep = file_separator)

  outpath_pc <- add_directory_to_path(path = outpath_data, dir_name = outpath_suff_pc, append_to_cur_file_name = TRUE,
                                      file_separator = file_separator)

  return(outpath_pc)
}

get_outpath_pc_graph <- function(prefix = NULL, plot_only_subgraphs = NULL, coloring = NULL,
                                 graph_layout = NULL, weighted = FALSE, kernelized = TRUE, plot = FALSE) {
  str <- ""
  if (plot) {
    str <- "PLOT"
  }
  if ((nchar(coloring) == 0) || (is.null(coloring)) || (coloring == "auto")) {
    coloring <- ""
  }
  graph_str <- "graph"
  if (weighted) {
    graph_str <- paste0(graph_str, "weighted")
  }
  if (kernelized) {
    kern <- paste("kern", graph_str)
  }
  str <- pastes(str, graph_str, graph_layout,
                pastes_parvalue("col", coloring), sep = "-")
  str <- pastes(prefix, str, sep = "_")

  return(str)
}

get_outpath_graph_evaluation <- function(prefix = "", stage = NULL, object = NULL) {
  outpath <- pastes(prefix, "localTests", sep = "_")
  if (grepl("orig", stage)) {
    stage <- ""
  }
  outpath <- pastes(outpath, stage, object, sep = "-")
  return(outpath)
}

get_outpath_graph_clustering <- function(prefix = "", cluster_type = NULL, object = "") {
  return(pastes(prefix, pastes("CLUST", cluster_type, object, sep = "-"), sep = "_"))
}


get_outpath_causal_analysis <- function(prefix = "", first_graph_i = NULL, last_graph_i = NULL,
                            current_graph_i = NULL, object = "") {
  # outpath <- pastes(causal_effects_function, sep = "-")
  outpath = ""
  if (is.null(current_graph_i) && !is.null(first_graph_i) && !is.null(last_graph_i)) {
    regarded_graphs_indices <- first_graph_i:last_graph_i
    if (length(regarded_graphs_indices) == 1) {
      outpath <- paste0(outpath, "-G", regarded_graphs_indices)
    }
    if (length(regarded_graphs_indices) > 1) {
      outpath <- paste0(outpath, "-G", first_graph_i, "-", last_graph_i)
    }
  }
  outpath <- pastes(outpath, object, sep = "-")

  outpath_with_prefix <- pastes(prefix, outpath, sep = "_")
  if (!is.null(current_graph_i)) {
    suffix <- paste0("G", current_graph_i)
    outpath_with_prefix <- add_directory_to_path(path = outpath_with_prefix, dir_name = suffix)
    create_parent_directory_if_necessary(outpath_with_prefix)
  }

  return(outpath_with_prefix)
}

get_outpath_ida <- function(prefix, causal_effects_fun, direction, weight_effects_on_by = "", option_nr = NULL,
                            neg_effects, perturbed_position) {
  # outpath <- paste0(prefix, "-totaleffects(", neg_effects, ")")

  out_file = ""

  if ((grepl(pattern = "reset", tolower(causal_effects_fun)))) {
    out_file <- pastes(out_file, "reset", sep = "-")
  } else if (grepl(pattern = "causaleffect", tolower(causal_effects_fun)) ||
             grepl(pattern = "causal_effect", tolower(causal_effects_fun))) {
    out_file <- pastes(out_file, "pseudo", sep = "-")
  } else {
    stop("No ida-fun for outpath")
  }

  out_file <- pastes(out_file, pastes_parvalue("neg", neg_effects), sep = "-")

  if (!(!missing(direction) && grepl("on", direction))) {
    weight_effects_on_by = ""
  }

  out_file <- pastes(out_file, direction, weight_effects_on_by, pastes_parvalue("pos", perturbed_position, sep = ""), sep = "-")

  out_file <- pastes(out_file, pastes_parvalue("#", option_nr, sep = ""), sep = "-")

  return(pastes(prefix, paste("IDA", out_file, sep = "-"), sep = "_"))
}

get_outpath_ida_plot <- function(prefix, amplification_exponent, amplification_factor, rank_effects, no_colors, effect_to_color_mode) {
  str <- ""

  if (startsWith(effect_to_color_mode, "opac")) {
    str <- paste0(str, "-opac")
  }
  if (rank_effects && !(no_colors)) {
    str <- paste0(str, "-ranked")
  } else {
    if (amplification_exponent != 1) {
      str <- pastes(str, pastes_parvalue("amplExp=", amplification_exponent), sep = "-")
    }
    if (amplification_factor) {
      str <- pastes(str, "amplFac", sep = "-")
    }
  }
  if (no_colors) {
    str <- paste0(str, "-bw")
  }
  return(pastes(prefix, paste("PLOT", str, sep = "-"), sep = "_"))
}


outpath_pairwise_effects <- function(prefix = "", pre_fct = "", k = "", type = "", object = "") {
  outpath <- pastes("CLUST", pastes_parvalue("prefct", pre_fct), type, pastes_parvalue(pastes_parvalue(">", k, sep = ""), "clusters", sep = "-"), object, sep = "-")
  return(pastes(prefix, outpath, sep = "_"))
}


# gab früher den old_outpath zurück,
# jetzt das Verzeichnis in dem es die Datei gibt, oder NULL, wenn es sie nicht gibt.
# TODO: rename: get_file (oder so) / get_outpath_where_file_exists
# TODO: try - <-> _
get_old_outpath <- function(outpath, suffix) {
  if (file.exists(paste0(outpath, suffix))) {
    return(paste0(outpath, suffix))
  }

  if (endsWith(outpath, "-sc-maj")) {
    old_outpath <- str_replace(outpath, "-sc-maj", "-rel")
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
