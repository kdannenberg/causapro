source("configuration_code.R")

source("compute_DAG_G.R")

source("configuration_code.R")
source("configuration_data.R")

alphas_vector <- c(1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2)

s <- 2 # relaxed leifert sowieso immer den gleichen Graphen

measure <- number_of_conflict_edges

number_of_conflict_edges <- function(results) {
  return(conflict_edges(results$pc@graph)$conflict)
}

fraction_of_confict_edges <- function(results) {
  return(conflict_edges(results$pc@graph)$conflict / do.call(sum, conflict_edges(results$pc@graph)))
}

edges <- function(results) {
  print(conflict_edges(results$pc@graph))
  print(do.call(sum, conflict_edges(results$pc@graph)))
}

# protein_causality_G(alpha = alpha, pc_solve_conflicts = TRUE, pc_u2pd = "relaxed")

results_by_alpha <- list()
for (alpha in alphas_vector) {
  print(paste0("alpha = ", alpha))
  results_for_this_alpha <- c()
  for (i in 1:s) {
     results_for_this_alpha[[i]] <- protein_causality_G(alpha = alpha, pc_solve_conflicts = TRUE, pc_u2pd = "relaxed", 
                                                        graph_computation = TRUE, evaluation = TRUE)
  }
  results_by_alpha[[as.character(alpha)]] <- results_for_this_alpha
}

# means_for_alphas <- sapply(values, mean)

graphics.off()

confl_edges_per_alpha <- sapply(results_by_alpha, function(results) sapply(results, number_of_confict_edges))
frac_of_confl_edges_per_alpha <- sapply(results_by_alpha, function(results) sapply(results, fraction_of_confict_edges))

mean_confl_edges_per_alpha <- apply(confl_edges_per_alpha, 2, mean)
mean_frac_of_confl_edges_per_alpha <- apply(frac_of_confl_edges_per_alpha, 2, mean)

plot_names(mean_confl_edges_per_alpha)
lines(mean_confl_edges_per_alpha)

# plot_names(mean_frac_of_confl_edges_per_alpha)
# lines(mean_frac_of_confl_edges_per_alpha)

plot_names <- function(x) {
  plot(x = names(x), y = x, xlab = "alpha", ylab = deparse(substitute(x)))
}
