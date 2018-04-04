source("~/.configuration_code.R")

source_all_function_scripts()

# source("compute_DAG_numerical.R")
# source("general_functions.R")
# source("evaluate_DAG.R")
# source("causal_effects.R")
# source("linkcommunities.R")
# source("ci-tests.R")

# setwd("/Volumes/Causality/Viren/R/")
source("configuration_data.R")


# protein_causality_G(pc_maj_rule = TRUE)
# results_G <- protein_causality_G(alpha = 0.1, pc_conservative = FALSE, pc_u2pd = "relaxed", pc_maj_rule = TRUE, pc_solve_confl = TRUE, analysis = FALSE)

# TODO: warum geht das nur für effects of pos. 372
                                        # results_G <- protein_causality_G(pc_conservative = FALSE, pc_u2pd = "retry", pc_solve_confl = TRUE, analysis = TRUE)
## results_G <- protein_causality_G(pc_conservative = FALSE, pc_maj_rule = TRUE, pc_u2pd = "relaxed", pc_solve_confl = TRUE, 
##                                  analysis = FALSE, weight_effects_on_by = "mean", min_pos_var = 0, alpha = 0.01)
## sink()
## print(conflict_edges(results_G$pc@graph))

# DDDG-10.0.08.1e-04
results_G <- protein_causality_G(type_of_data = "DDG", subtype_of_data = "", min_pos_var = 0.01, 
                                 alpha = 0.000001234, analysis = TRUE, pc_maj_rule = TRUE, 
                                 intervention_position = "all",#"372","all",
                                 plot_no_isolated_nodes = TRUE, plot_clusters = FALSE, plot_ida = FALSE,
                                 effects_cluster_method = "ward.D",#"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
                                 effects_pv_nboot = 10000, effects_cluster_alpha = 0.7)
# g <- results_G$pc@graph
# g_w <- set_edge_weights_for_graph(g, cov(data))
# # ida()
# print(wgtMatrix(g))
# print(wgtMatrix(g_w))
# 
# 
# print(paste("1,2: ", paste(ida(x.pos = 1, y.pos = 2, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("1,2: ", causalEffect(x = 1, y = 2, g = g_w)))
# 
# print(paste("1,3: ", paste(ida(x.pos = 1, y.pos = 3, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("1,3: ", causalEffect(x = 1, y = 3, g = g_w)))
# 
# print(paste("7,2: ", paste(ida(x.pos = 7, y.pos = 2, mcov = cov(data), graphEst = g), collapse = " / ")))
# print(paste("7,2: ", causalEffect(x = 7, y = 2, g = g_w)))

# source("~/.configuration_code.R")


