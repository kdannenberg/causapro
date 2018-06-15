source("~/.configuration_code.R")

source_all_function_scripts()

source("configuration_data.R")

# debug(protein_causality)
# debug(analysis_after_pc)
# debug(evaluate_DAG)
# debug(adjust_data)
# debug(plot_clusters_in_pymol)
# debug(protein_graph_clustering)

# debug(partuning_over_alpha_and_minposvar)
# debug(partune_alpha_minposvar_empir_cor_dev_from_zero_per_edge)


# NoV_NMR-Tit_B3S-with-unass
type_of_data = "NMR-Tit"
## subtype_of_data = "B3S-with-unass"
# subtype_of_data = "B3S_with-unass_1H" ## options are ..._Euclidean, ..._1H and ..._15N
##subtype_of_data = "BTS"
# subtype_of_data = c("Fuc", "B4S")

subtype_of_data = "" # Daten vom 1.6.18

pc_fun <- function_set_parameters(protein_causality_NoV, parameters =
                                     list(type_of_data = type_of_data,
                                          subtype_of_data = subtype_of_data,
                                          # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                          # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                          # causal_analysis = FALSE, #min_pos_var = 0.01,
                                          # evaluation = TRUE,
                                          # alpha = 0.1, min_pos_var = 0,
                                          # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                          # ranked = TRUE, pc_indepTest = "smc-jt",
                                          mute_all_plots = TRUE))

pc_fun_with_eval <- function_set_parameters(protein_causality_NoV, parameters =
                                    list(type_of_data = type_of_data,
                                         subtype_of_data = subtype_of_data,
                                         # pc_conservative = FALSE, pc_maj_rule = TRUE,
                                         # pc_u2pd = "relaxed", pc_solve_confl = TRUE,
                                         # causal_analysis = FALSE, #min_pos_var = 0.01,
                                         evaluation = TRUE,
                                         # alpha = 0.1, min_pos_var = 0,
                                         # only_cols = c(272, 278, 281, 301, 330, 416, 418, 419, 431, 472, 486, 527),
                                         # ranked = TRUE, pc_indepTest = "smc-jt",
                                         mute_all_plots = TRUE))


# alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, seq(0.2, 0.9, 0.1))
# alphas <- c(0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas = c(0.001, 0.01, 0.05, 0.1)
alphas = c(0.001, 0.8)
# alphas <- c(1e-10, 1e-5, 0.0001)
# alphas <- c(0.8)
min_pos_vars = c(0, 0.08) # kleinste vorhandene Varianz: 0.059

# dev_of_edges_from_2 <- tune_alpha_mpv_dev_of_edges_from_2(pc_fun = pc_fun,
#                                                           alphas = alphas,
#                                                           minposvars = min_pos_vars)
# # localTest estimate per (squared/...) edge(s)
# estimate_per_edges <- tune_alpha_mpv_estimate_per_edges(pc_fun = pc_fun,
#                                                         alphas = alphas,
#                                                         minposvars = min_pos_vars)


mean_estimate <-  tune_alpha_mpv_mean_empir_cor(pc_fun = pc_fun_with_eval,
                                                          alphas = alphas,
                                                          minposvars = min_pos_vars)

stop("fertig")

partuning_mean_empir_cor <- function_set_parameters(partune_alpha_minposvar_mean_empir_cor,
                                                                    parameters = list(pc_FUN = pc_fun))
print(partune_alpha_minposvar_mean_empir_cor)



# localTests_n_edges_obj_fct <- function(array) {
#   v <- apply(array, c(1,2), function(x) {sum(abs(x))})
#   return(list(value = min(v), index = (which(v ==  min(v), arr.ind = TRUE))))
# }

# localTests_sep_estimate_obj_fct <- function_set_parameters(array_obj_function,
#                     parameters = list(fun_over_components_of_value = function(x) {sum(abs(x))},
#                                       fun_over_results_of_other_fun = min))

array_obj_function <- function(array, fun_over_components_of_value, fun_over_results_of_other_fun) {
  v <- apply(array, c(1,2), fun_over_components_of_value)
  best_v <- fun_over_results_of_other_fun(v)
  return(list(value = best_v, index = (which(v ==  best_v, arr.ind = TRUE))))
}

localTests_sep_estimate_n_edges_obj_fct <-
  function_set_parameters(array_obj_function,
  parameters = list(fun_over_components_of_value = function(v) {v[1]/v[2]},
                    fun_over_results_of_other_fun = min))




# debug(partuning_Nov_NMR_2018_06_01)





# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_Nov_NMR_2018_06_01_square_empir_cor, alphas = alphas, minposvars = min_pos_vars, best = min)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_Nov_NMR_2018_06_01_n_edges, alphas = alphas, minposvars = min_pos_vars, best = max)
# opt <- partuning_over_alpha_and_minposvar(FUN = partuning_empir_cor_per_edge, alphas = alphas, minposvars = min_pos_vars, best = min)

# localTest estimate per (squared/...) edge(s)
estimate_per_edges <- partuning_over_alpha_and_minposvar(FUN = partuning_sep_empir_cor_n_edge,
                                          alphas = alphas, minposvars = min_pos_vars,
                                          best = localTests_sep_estimate_n_edges_obj_fct, best_in_array = FALSE)

# zero mean localTests estimate
mean_estimate <- partuning_over_alpha_and_minposvar(FUN = partuning_mean_empir_cor,
                                          alphas = alphas, minposvars = min_pos_vars,
                                          best = function(x) {min(abs(x), na.rm = TRUE)})

# best partuning so far:
# estimate_per_edges_mat * abs(mean_est_mat)^2 und alles mit konfliktkanten ausschließen (z.B. alpha = 0.7, bei minposvar = )

# print(opt)

# for localTest estimate per (squared/...) edge(s)
# debug(plot_different_measures_of_all_results)
measures <- c(function(v) {v[1]}, function(v) {v[1]/sqrt(v[2])}, function(v) {v[1]/v[2]}, function(v) {v[1]/v[2]^2}, function(v) {v[1]/sqrt(v[2])^v[2]})
par(mfrow = c(1, length(measures)))
plot_different_measures_of_all_results(all_results = estimate_per_edges$all_results, measures = measures, print = TRUE)
