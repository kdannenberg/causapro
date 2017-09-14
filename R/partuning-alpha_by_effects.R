source("~/.configuration_code.R")

source("functions_partuning.R")

source("configuration_data.R")

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


best_alphas <- find_best_alphas(all_q_score)

# swap_minposvar_and_alpha <- function(nested_list) {
#   do.call(list, lapply(nested_list, function(list_of_lists) return))
# }

effects_for_distinct_alphas(best_alphas, with_graphs = TRUE, min_pos_vars = c(0.01), for_all_alphas = TRUE)
# effects_for_distinct_alphas(best_alphas, with_graphs = FALSE)

