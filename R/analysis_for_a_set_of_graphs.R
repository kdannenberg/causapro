# previously: analysis_find_good_graph.R
source("~/.configuration_code.R")

source("functions_causal_effects.R")
source("functions_general.R")
source("functions_conversions.R")
source("functions_tools.R")

source("functions_analysis_for_a_set_of_graphs.R")
source("functions_protein_causality.R")
# tests
source("tests_analysis_for_a_set_of_graphs.R")

source("compute_DAG_G.R")
source("~/.configuration_code.R")
source("compute_DAG_S.R")
# source("configuration_data.R")

# TODO: sicherstellen, dass oma nur gesetzt wird, wenn auch eine main_caption (title) geprintet wird

# effects <- analyse_set_of_graphs(direction = "mean", measure = "G", alpha = 0.03, min_pos_var = 0.01, new = TRUE)




analyse_graphs <- function(measure = "S", type_of_data = "DDS", subtype_of_data = "", 
               alphas = c(0.001, 0.005, 0.01, 0.05, 0.1), min_pos_vars = c(0.0001, 0.001, 0.01),
  protein_causality_function = get(paste0("protein_causality_", measure)),
  # TODO: DAS FUNKTIONIERT NICHT!! DIE PARAMETER, DIE NCIHT ÜBERGEBEN WERDEN, WERDEN NCIHT JETZT SCHON BELEGT!
  # ALLES WIRD ERST BEI AUFRUF AUSGEWERTET
  # pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
  #   protein_causality_function(type_of_data = eval(type_of_data), subtype_of_data = eval(subtype_of_data), min_pos_var = min_pos_var, alpha = alpha,
  #                              pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
  #                              evaluation = evaluation, analysis = analysis, mute_all_plots = TRUE)
  # pc_function = f_protein_causality_pc_parameters_eval_analysis(measure = measure, type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                                #### min_pos_var = min_pos_var, alpha = alpha, 
                                                                # mute_all_plots = TRUE)
  pc_function = function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                              mute_all_plots = TRUE))
 ) {

  # if (missing(pc_function)) {
  #   pc_function <- function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
  #                                                                         mute_all_plots = TRUE))
  # }
  
  effects <- list()
  
  graphics.off()
  oma <- c( 2, 0, 2, 0 )  # oberer Rand für Caption: eine Zeile mehr als benötigt
  par(mfrow = c(length(alphas), length(min_pos_vars)), oma = oma)
  for (alpha in alphas) {
    effects[[as.character(alpha)]] <- list()
    for (min_pos_var in min_pos_vars) {
      # for (measure in c("G", "S")) {
        # if (
        #     (
        #       (type_of_data == "DDDG" && subtype_of_data == "5") &&
        #       ((alpha == 0.05 && min_pos_var == 0.001) #||  # More than 15 conflict edges. (22)
        #       || (alpha == 0.1 && min_pos_var == 0.001)  # More than 15 conflict edges.
        #       || (alpha == 0.05 && min_pos_var == 0.0001) # More than 15 conflict edges.
        #       || (alpha == 0.1 && min_pos_var == 0.0001) # More than 15 conflict edges.
        #       )
        #      # ) || (
        #      #   (type_of_data == "DDDG" && subtype_of_data == "10") &&
        #      #   (
        #      #   (alpha == 0.005 && min_pos_var == 0.001) ##  "Fehler in wgt.unique[x, ] : Indizierung außerhalb der Grenzen"
        #      #   # || (alpha == 0.05 && min_pos_var == 0.0001)
        #      #   # || (alpha == 0.1 && min_pos_var == 0.0001)
        #      #   )
        #      ) || (
        #        (type_of_data == "DDG" && subtype_of_data == "10") &&
        #        (
        #          (alpha == 0.1 && min_pos_var == 0.0001) ##  15 conflict edges
        #          || (alpha == 0.1 && min_pos_var == 0.001) ##  13 conflict edges
        #        )
        #      ) || (
        #        (measure == "S") &&
        #        (alpha == 0.1 && min_pos_var == 0.0001)
        #      ) || (
        #        (type_of_data == "DDG" && subtype_of_data == "5") &&
        #        (
        #          (alpha == 0.05 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
        #       || (alpha == 0.05 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
        #       || (alpha == 0.05 && min_pos_var == 0.01) #  More than 15 conflict edges. (18)
        #       || (alpha == 0.1 && min_pos_var == 0.0001) #  More than 15 conflict edges. (20)
        #       || (alpha == 0.1 && min_pos_var == 0.001) #  More than 15 conflict edges. (20)
        #       || (alpha == 0.1 && min_pos_var == 0.01) #  More than 15 conflict edges. (21)
        #        )
        #      )
        #     ) {
        #   plot.new()
        # } else {
        pc_function_ <- function_set_parameters(pc_function, parameters = list(alpha = alpha, min_pos_var = min_pos_var))
        effects[[as.character(alpha)]][[as.character(min_pos_var)]] <- analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                direction = "mean", measure = measure, pc_function = pc_function_, alpha = alpha, min_pos_var = min_pos_var,
                                for_combined_plot = TRUE, scale_in_the_end = FALSE, new = TRUE, save = FALSE)
  
        # }
        # }
    }
  }
  
  # title(main = paste0(type_of_data, "-", subtype_of_data),
  #       sub = "Mean causal effects over all conflict graphs and over the effects on and of position 372", outer = TRUE)
  title(main = paste0("Mean causal effects over all conflict graphs and over the effects on and of position 372",
                      "\n", type_of_data, "-", subtype_of_data), outer = TRUE)
  return(effects)
}

all_effects <- list()
for (measure_type_sub in c("DDS", "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")) {
  measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
  subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
  if (is.na(subtype_of_data)) {
    subtype_of_data <- ""
  }
  all_effects[[measure_type_sub]] <- analyse_graphs(measure = measure, type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
                                                    subtype_of_data = subtype_of_data,
                                                    alphas = c(0.001, 0.005, 0.01, 0.05, 0.1), min_pos_vars = c(0.0001, 0.001, 0.01))
}