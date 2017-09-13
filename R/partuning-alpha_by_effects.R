new = FALSE

# get all_effects
if (!exists("all_effects")) {
  if (!new) {
    load(file = "RData/all_effects.RData")
  } else {
    source('~/Documents/Uni/Viren/ProteinCausalPaths/R/analysis_for_a_set_of_graphs.R')
  }
}

int_pos <- interesting_positions("PDZ")
percentile = 1 - (length(int_pos) / dim(data)[2])

function_score_for_effects <- function_set_parameters(score_for_effects, 
                                                      parameters = list(int_pos = int_pos, perturbed_position = 372))
all_q_score <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_score_for_effects)

find_best_alphas <- function(all_scores) {
  result <- list()
  for (measure_type_sub in names(all_scores)) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    # type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    # if (grepl("-", measure_type_sub)) {
    #   subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    # } else {
    #   subtype_of_data <- ""
    #   subtype_of_data_list <- "-"
    # }
    result[[measure_type_sub]] <- list()
    for (min_pos_var in as.numeric(names(all_scores[[measure_type_sub]][[1]]))) {
      values <- list()
      for (alpha in as.numeric(names(all_scores[[measure_type_sub]]))) {
         values[[as.character(alpha)]] <- all_scores[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]]
      }
      result[[measure_type_sub]][[as.character(min_pos_var)]] <- names(which(unlist(values) == max(unlist(values))))
    }
  }
  return(result)
}

best_alphas <- find_best_alphas(all_q_score)

# swap_minposvar_and_alpha <- function(nested_list) {
#   do.call(list, lapply(nested_list, function(list_of_lists) return))
# }

effects_for_distinct_alphas(best_alphas, with_graphs = TRUE, min_pos_vars = c(0.01), for_all_alphas = TRUE)
# effects_for_distinct_alphas(best_alphas, with_graphs = FALSE)

