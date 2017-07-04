# all_results <- list()
# all_graphs <- list()
for (i in 18:100) {
  source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
  edges <- conflict_edges(results$pc@graph)
  all_results[[i]] <- results
  all_graphs[[i]] <- results$pc@graph
  if ((edges$conflict == 0) && (edges$bidirected == 0)) {
    break
  }
}