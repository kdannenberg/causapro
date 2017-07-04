for (i in 1:10000) {
  source('~/Documents/Uni/Viren/ProteinCausalPaths/R/compute_DAG_G.R')
  edges <- conflict_edges(results$pc@graph)
  if ((edges$conflict == 0) && (edges$bidirected == 0)) {
    break
  }
}