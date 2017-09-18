source("~/.configuration_code.R")

source_all_function_scripts()

source("compute_DAG_G.R")
source("compute_DAG_S.R")

source("configuration_data.R")


new = FALSE
print_best_scenarios = TRUE
with_graphs = TRUE
for_all_best_alphas = TRUE

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

if (print_best_scenarios) {
  all_scores <- unlist(all_q_score)
  best_scenarios <- all_scores[all_scores == max(all_scores)]
  print(best_scenarios)
}


best_alphas <- find_best_alphas(all_q_score)

# swap_minposvar_and_alpha <- function(nested_list) {
#   do.call(list, lapply(nested_list, function(list_of_lists) return))
# }

# debug(determine_set_of_graphs)

effects_for_distinct_alphas(best_alphas, with_graphs = with_graphs, min_pos_vars = c(0.01), for_all_alphas = for_all_best_alphas)
# effects_for_distinct_alphas(best_alphas, with_graphs = FALSE)

