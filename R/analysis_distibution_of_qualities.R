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

# which_effects gives the slot of the list of effects (e.g. )
apply_to_all_effects_in_nested_list <- function(all_effects, FUN, which_effects = "overAllGraphs_mean_on_of") {
  result <- list()
  for (measure_type_sub in names(all_effects)) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    # type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    # if (grepl("-", measure_type_sub)) {
    #   subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    # } else {
    #   subtype_of_data <- ""
    #   subtype_of_data_list <- "-"
    # }
    result[[measure_type_sub]] <- list()
    for (alpha in as.numeric(names(all_effects[[measure_type_sub]]))) {
      result[[measure_type_sub]][[as.character(alpha)]] <- list()
      for (min_pos_var in as.numeric(names(all_effects[[measure_type_sub]][[as.character(alpha)]]))) {
        if (!is.null(all_effects[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]][[which_effects]])) {
          result[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]] <- 
            FUN(all_effects[[measure_type_sub]][[as.character(alpha)]][[as.character(min_pos_var)]][[which_effects]])
        }
      }
      if (length(result[[measure_type_sub]][[as.character(alpha)]]) == 0) {
        result[[measure_type_sub]][[as.character(alpha)]] <- NULL
      }
    }
  }
  return(result)
}

function_quality_of_effects_distibution <- function_set_parameters(quality_of_effects_distibution, 
                                                                   parameters = list(int_pos = int_pos))
all_q_height <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_quality_of_effects_distibution)


function_quality_of_effects_classifier <- function_set_parameters(quality_of_effects_classifier, 
                                                                  parameters = list(int_pos = int_pos,
                                                                                    perturbed_position = 372))
all_q_classifier <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_quality_of_effects_classifier)


function_score_for_effects <- function_set_parameters(score_for_effects, 
                                                      parameters = list(int_pos = int_pos, perturbed_position = 372))
all_q_score <- apply_to_all_effects_in_nested_list(all_effects = all_effects, FUN = function_score_for_effects)

# TODO: mit lapply
# all_q_classifier <- lapply(all_effects, lapply, lapply, quality_of_effects_distibution, int_pos = int_pos)
# all_q_classifier <- lapply(all_effects, 
#                            function(x) {lapply(x,
#                                                function(x) {lapply(x$overAllGraphs_mean_on_of, 
#                                                                   quality_of_effects_distibution, int_pos = int_pos)})})
par(mfrow = c(1,3))
breaks = 6
hist(unlist(all_q_height), breaks = breaks)
# hist(unlist(all_q_classifier_old))
hist(unlist(all_q_classifier), breaks = breaks)
hist(unlist(all_q_score), breaks = breaks)
