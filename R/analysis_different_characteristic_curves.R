source("configuration_code.R")



#necessary??
source("functions_causal_effects.R")
source("functions_ci_tests.R")
source("functions_compute_DAG_categorical.R")
source("functions_compute_DAG_numerical.R")
source("functions_conversions.R")
source("functions_evaluate_DAG.R")
source("functions_general.R")
source("functions_i_o.R")
source("functions_linkcommunities.R")
source("functions_pymol.R")
source("functions_tools.R")

source("compute_DAG_G.R")

source("configuration_code.R")
source("configuration_data.R")

# alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2)   # no more is fesible with localTests 
alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2)

s <- 5 # relaxed leifert sowieso immer den gleichen Graphen

new = FALSE
save = TRUE

min_pos_var = 0


pc_solve_conflicts = FALSE
pc_u2pd = "retry"

results <- protein_causality_G(min_pos_var = min_pos_var,
                               pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, 
                               graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
                               data_in_results = TRUE, output_parameters_in_results = TRUE)
data <- results$data
caption <- results$summary$caption
outpath <- results$summary$outpath
directories <- strsplit(outpath, file_separator)
outpath <- paste(directories[[1]][1:(length(directories[[1]])-2)], collapse = "/", sep = "/")
#TODO: letztes Verzeichnis entfernen


number_of_conflict_edges <- function(results) {
  return(conflict_edges(results$pc@graph)$conflict)
}

fraction_of_conflict_edges <- function(results) {
  return(conflict_edges(results$pc@graph)$conflict / do.call(sum, conflict_edges(results$pc@graph)))
}

edges <- function(results) {
  print(conflict_edges(results$pc@graph))
  print(do.call(sum, conflict_edges(results$pc@graph)))
}

measure <- fraction_of_conflict_edges

if (new) {
  results_by_alpha <- list()
  for (alpha in alphas_vector) {
    print(paste0("alpha = ", alpha))
    results_for_this_alpha <- c()
    for (i in 1:s) {
      results_for_this_alpha[[i]] <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha, pc_solve_conflicts = TRUE, pc_u2pd = "relaxed",
                                                         graph_computation = TRUE, evaluation = FALSE)
    }
    results_by_alpha[[as.character(alpha)]] <- results_for_this_alpha
  }
  if (save) {
    save(results_by_alpha, file = paste0(outpath, "/", s, "-results-per-alpha.RData"))
  }
  
} else {
  load(file = paste0(outpath, "/", s, "-results-per-alpha.RData"))
}


# means_for_alphas <- sapply(values, mean)

graphics.off()

values_per_alpha <- sapply(results_by_alpha, function(results) sapply(results, fraction_of_confict_edges))
mean_per_alpha <- apply(values_per_alpha, 2, mean)

# confl_edges_per_alpha <- sapply(results_by_alpha, function(results) sapply(results, number_of_confict_edges))
# frac_of_confl_edges_per_alpha <- sapply(results_by_alpha, function(results) sapply(results, fraction_of_confict_edges))

# mean_confl_edges_per_alpha <- apply(confl_edges_per_alpha, 2, mean)
# mean_frac_of_confl_edges_per_alpha <- apply(frac_of_confl_edges_per_alpha, 2, mean)

plot_names(mean_per_alpha)
lines(mean_per_alpha)

# plot_names(mean_frac_of_confl_edges_per_alpha)
# lines(mean_frac_of_confl_edges_per_alpha)

plot_names <- function(x) {
  plot(x = names(x), y = x, xlab = "alpha", ylab = deparse(substitute(measure)))
}
