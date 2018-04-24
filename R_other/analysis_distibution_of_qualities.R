source("~/.configuration_code.R")

source("functions_partuning.R")
source("functions_general.R")
source("functions_causal_effects.R")

source("configuration_data.R")

# debug(score_for_effects)

new = FALSE

# get all_effects
# if (!exists("all_effects")) {
  if (!new) {
    # load(file = "RData/all_effects.RData")
    load(file = "RData/all_effects-ida-reset.RData")
  } else {
    source('~/Documents/Uni/Viren/ProteinCausalPaths/R/analysis_for_a_set_of_graphs.R')
  }
# }

int_pos <- interesting_positions("PDZ")
percentile = 1 - (length(int_pos) / dim(data)[2])

# HEIGHT
function_quality_of_effects_distibution <- function_set_parameters(quality_of_effects_distibution, 
                                                                   parameters = list(int_pos = int_pos))
all_q_height <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_quality_of_effects_distibution)


function_score_of_effects_distibution <- function_set_parameters(score_of_effects_height, 
                                                                   parameters = list(int_pos = int_pos, perturbed_position = 372))
all_s_height <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_score_of_effects_distibution)

# CLASSIFIER
function_quality_of_effects_classifier <- function_set_parameters(quality_of_effects_classifier, 
                                                                  parameters = list(int_pos = int_pos,
                                                                                    perturbed_position = 372))
all_q_classifier <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_quality_of_effects_classifier)


function_score_of_effects_classifier <- function_set_parameters(score_of_effects_classifier, 
                                                                  parameters = list(int_pos = int_pos,
                                                                                    perturbed_position = 372))
all_s_classifier <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_score_of_effects_classifier)


# SUM
function_score_for_effects <- function_set_parameters(score_for_effects, 
                                                      parameters = list(int_pos = int_pos, perturbed_position = 372))
all_q_score <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_score_for_effects)

# TODO: mit lapply
# all_q_classifier <- lapply(all_effects, lapply, lapply, quality_of_effects_distibution, int_pos = int_pos)
# all_q_classifier <- lapply(all_effects, 
#                            function(x) {lapply(x,
#                                                function(x) {lapply(x$overAllGraphs_mean_on_of, 
#                                                                   quality_of_effects_distibution, int_pos = int_pos)})})
par(mfrow = c(1,7))
breaks = 10
hist(unlist(all_q_height), breaks = breaks, main = paste("q_height \n", breaks, "breaks"))
# hist(unlist(all_q_height), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), main = "q_height \n breaks from code")
# hist(unlist(all_q_height), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 40), main = "q_height \n breaks from code")
hist(unlist(all_q_height), breaks = c(0, 1, 2, 3, 4, 5, 40), main = "q_height \n breaks from code")
hist(unlist(all_s_height), breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), main = "s_height")
# hist(unlist(all_q_classifier_old))
hist(unlist(all_q_classifier), breaks = breaks, main = paste("q_classifier \n", breaks, "breaks"))
# hist(unlist(all_q_classifier), breaks = c(0, 1/6, 2/6, 3/6, 4/6, 5/6, 1), main = "q_classifier \n breaks from code")
hist(unlist(all_q_classifier), breaks = c(0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1), main = "q_classifier \n breaks from code")
hist(unlist(all_s_classifier), breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), main = "s_classifier")
# barplot(unlist(all_s_classifier), main = "s_classifier")

hist(unlist(all_q_score), breaks = breaks, main = "q_score")
