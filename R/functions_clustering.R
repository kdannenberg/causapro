#' Cluster by Causal Structure (Graph Clustering)
#' @param monochromatic_removed_cols all removed_cols are given (the same shade of) white as color
#' @param more_levels_of_conservedness (bacause of too little variacne) removed positions are
#' clored in shades of white with increasing amounts of blue with increasing variance
protein_graph_clustering <- function(graph, protein, position_numbering, coloring, colors,
                                     outpath, output_formats, file_separator,
                                     mute_all_plots, caption, cluster_methods,
                                     remove_singular_clusters = TRUE,
                                     merge_singular_clusters = FALSE, add_cluster_of_conserved_positions = FALSE,
                                     removed_cols, monochromatic_removed_cols = TRUE,
                                     more_levels_of_conservedness = FALSE,
                                     sort_clusters = length) {

  igraph_plot <- graph$plot
  igraph_cluster <- igraph.from.graphNEL(graph$cluster)
  ## sort_clusters = "DDS-SVD") {
  if (add_cluster_of_conserved_positions) {
    # node_clustering <- c(node_clustering, "#FFFFFF" = list(removed_cols))
    # names(node_clustering)[length(node_clustering)] <- "#FFFFFF"


    if (more_levels_of_conservedness) {
      removed_cols <- (removed_cols / max(removed_cols)) / 0.9999
    } else if (monochromatic_removed_cols) {
      removed_cols <- removed_cols * 0
    }

    if (!length(removed_cols) == 0) {

      colors_for_rem_pos <- color_by_effect(effects = removed_cols, int_pos = "", color_for_other_positions = "#000000", mode = "#FFFFFF")

      add_clusters <- sapply(names(table(colors_for_rem_pos)),
                             FUN = function(color) {return(names(colors_for_rem_pos[which(colors_for_rem_pos == color)]))},
                             simplify = FALSE, USE.NAMES = TRUE)

      # TODO ergibt hier keinen Sinn, oder?
      # if (remove_singular_clusters) {
      #   node_clustering <- remove_singular_clusters(node_clustering)
      # }

    }
  }

  for (clustering in cluster_methods) {
    if (clustering == "edge_betweenness") {
      type <- "eb"
    } else if (clustering == "infomap") {
      type <- "im"
    } else {
      type <- "igraph"
    }
    ## old conversion method
    ## igraph <- graph_from_graphnel(results$pc@graph)
    # graph <- results$pc@graph
    cluster_fct <- get(paste0("cluster_", clustering))
    ## obs: absolute weights
    if (clustering == "infomap") {
      cl <- cluster_fct(igraph_cluster, e.weights = abs(E(igraph_cluster)$weight), modularity = FALSE)
    } else {
      cl <- cluster_fct(igraph_cluster, weights = abs(E(igraph_cluster)$weight), modularity = FALSE)
    }

    node_clustering <- groups(cl)
    node_clustering <- unname(node_clustering)  # otherwise interpreted as colors

    if (remove_singular_clusters) {
      node_clustering <- remove_singular_clusters(node_clustering)
    } else if (merge_singular_clusters) {
      node_clustering <- merge_singular_clusters(node_clustering)
    }

    clustering_colors <- rainbow(length(node_clustering))
    names(node_clustering) <- clustering_colors

    ## TODO: save plot, instaed of plotting
    if (!mute_all_plots) {
      ## edge.arrow.size determines size of arrows (1 is default), vertex.size determines size of the vertices (15 is default), edge.width determines width of edges (1 is default)
      ## plot(cl, igraph, main = paste0(caption, "\n", clustering), edge.arrow.size=0.2, vertex.size=8, edge.width=0.7)
      call_plot_igraph(g = igraph_plot, protein = protein, position_numbering = position_numbering,
                       coloring = coloring, colors = colors, clusters = TRUE, cluster_str = type,
                       clustering = cl, clustering_colors = clustering_colors, caption = caption,
                       outpath = paste0(outpath,"-", length(node_clustering), "_", "clusters"),
                       output_formats = output_formats, mute_all_plots = mute_all_plots)
      ##this is the old version, just in case my adjustments don't work for you
      ##plot(cl, igraph, main = paste0(caption, "\n", clustering))
    }



    if (!is.null(sort_clusters)) {
      # if ((length(which(names(node_clustering) != "")) == 1)
      #     && names(node_clustering)[length(node_clustering)] != "") { # only one color defined, and that for the last cluster
      #   length_except_last <- vapply(node_clustering[1:(length(node_clustering) - 1)], length, 1L)
      #   node_clustering <- node_clustering[c(order(length_except_last, decreasing = TRUE), length(node_clustering))]
      # } else {
      if (typeof(sort_clusters) == "closure" || typeof(sort_clusters) == "builtin") {
        node_clustering <- node_clustering[order(vapply(node_clustering, sort_clusters, 1L), decreasing = TRUE)]
      } else if (sort_clusters == "DDS-SVD") {
        DDS_SVD <- read_data(files = paste0(protein, "_DDS-SVD-1"))
        cluster_weights <- lapply(node_clustering, function(nodes) {
          nodes <- nodes[which((nchar(nodes) >= 3) && (nodes %in% colnames(DDS_SVD)))]
          len <- length(nodes)
          if (len > 0) {
            sum <- sum(DDS_SVD[,nodes])
            return(sum / len)
          } else {
            return(0)
          }
        })
        node_clustering <- node_clustering[order(unlist(cluster_weights), decreasing = TRUE)]
      }
    }


    if (add_cluster_of_conserved_positions && exists("add_clusters")) {
      node_clustering <- c(node_clustering, add_clusters)
    }


    plot_clusters_in_pymol(node_clustering = node_clustering, protein = protein, outpath = outpath,
                           file_separator = file_separator, type_of_clustering = type)
  }
}

remove_singular_clusters <- function(clustering, force = FALSE) {
  clustering_new <- clustering[sapply(clustering,
    function(cluster) {
      length(cluster) > 1
    })]
  if ((length(clustering_new) == 0) && !(force)) {
    return(clustering)
  }
  return(clustering_new)
}

merge_singular_clusters <- function(clustering) {
  non_singular <- sapply(clustering,
                         function(cluster) {
                           length(cluster) > 1
                         })
  others <- unlist(clustering[!non_singular])
  clustering <- clustering[non_singular]
  clustering <- c(clustering, list("#000000" = others))
  # clustering[[length(clustering) + 1]] <- others
  # names(clustering)[[length(clustering)]] <- "#FFFFFF"
  return(clustering)
}

#' Cluster by Causal Effects
cluster_pairwise_effects <- function(results, pairwise_effects, k, cluster_method,
                                     hclust_method, dist_measure, alpha = 0.95,
                                     iterations_pv, protein, outpath, file_separator) {
  results$effects_clustering <- list()
  #k-means
  if (grepl(pattern = "k", cluster_method) && grepl(pattern = "means", cluster_method)) {
    k_m <- kmeans(t(pairwise_effects), k)
    lapply(seq(1:k), function (i) {dim(pairwise_effects[,names(which(k_m$cluster == i))])})

    # apply(pairwise_effects[,names(which(k_m$cluster == 1))], 2, barplot);

    results$pairwise_effects <- pairwise_effects

    clustering_with_duplicates <- k_m$cluster

    cl <- position_clustering_from_clustering_with_duplicates(clustering_with_duplicates = clustering_with_duplicates)

    print(cl)
    results$effects_clustering$k_means <- cl

    type <- "effects-km"
    # names(cl) <- NULL
    plot_clusters_in_pymol(node_clustering = cl, protein = protein, outpath = outpath,
                           file_separator = file_separator, type_of_clustering = type)
  }

  #hierarchical clustering
  # d <- dist(t(pairwise_effects))
  #
  # hc <- hclust(d, method = method)
  #
  # plot.new()
  # plot(hc)
  #
  # hc_cut <- cutree(hc, 5)
  #
  # hc_cut
  #
  # clustering_with_duplicates <- hc_cut
  #
  # cl <- position_clustering_from_clustering_with_duplicates(clustering_with_duplicates = clustering_with_duplicates)
  #
  # print(cl)
  #
  # type <- "effects-hc"
  # # names(cl) <- NULL
  # plot_clusters_in_pymol(node_clustering = cl, protein = protein, outpath = outpath,
  #                        file_separator = file_separator, type_of_clustering = type)
  #
  # pvclust
  # library(pvclust)


  else {
    FUN_pv <- function_set_parameters(pvclust, parameters = list(data = pairwise_effects,
                                                                 method.hclust = hclust_method,
                                                                 method.dist = dist_measure, nboot = iterations_pv))

    effects_pv <- compute_if_not_existent(filename = paste(outpath, "pv", hclust_method,
                                                           substr(dist_measure, 0, 3), iterations_pv, sep="-"),
                                          FUN = FUN_pv,
                                          obj_name = "effects_pv",
                                          fun_loaded_object_ok = function(effects_pv) {return(colnames(pairwise_effects) == effects_pv$hclust$labels)}
    )

    # fit <- pvclust(data = pairwise_effects, method.hclust="ward",
    #                method.dist="euclidean", nboot = 1000)

    # cl <- cutDendrogramAt(effects_pv$hclust, cutat = 15)
    # membershiplist_from_clusterlist(cl)

    plot.new()
    plot(effects_pv) # dendogram with p values

    high_clusterlist <- pvpick(effects_pv, alpha = alpha)$clusters
    alpha_changed = FALSE
    while (is.null(high_clusterlist[[1]])){
      alpha = alpha - 0.1
      high_clusterlist <- pvpick(effects_pv, alpha = alpha)$clusters
      alpha_changed = TRUE
    }
    if (alpha_changed) {
      warning(paste("Alpha was reduced to", alpha, "because no clusters could be found otherwise." ))
    }

    # add rectangles around groups highly supported by the data
    pvrect(effects_pv, alpha = alpha)

    high <- membershiplist_from_clusterlist(high_clusterlist)

    cl_pv <- position_clustering_from_clustering_with_duplicates(clustering_with_duplicates = high)

    print(cl_pv)

    type <- paste("effects-pv", hclust_method, substr(dist_measure, 0, 3), iterations_pv, paste0("iter-alpha=", alpha), sep="-")
    results$effects_clustering$pv[[type]] <- cl_pv

    # names(cl) <- NULL
    plot_clusters_in_pymol(node_clustering = cl_pv, protein = protein, outpath = outpath,
                           file_separator = file_separator, type_of_clustering = type)
  }
  return(results)
}



position_clustering_from_clustering_with_duplicates <- function(clustering_with_duplicates) {
  positions <- unique(sapply(names(clustering_with_duplicates), function(long_name) return(gsub("-.*","",long_name))))
  k = max(clustering_with_duplicates)


  ## averaged_clusters <- sapply(positions, function(position) {
  ##  position_clustering <- clustering_with_duplicates[which(grepl(position, names(clustering_with_duplicates)))]
  ##  mean_cluster <- mean(position_clustering)
  ## names(mean_cluster) <- position
  ##  return(mean_cluster)
  ##})

  ## averaged_clusters <- round(averaged_clusters)
  rb_cols = substr(rainbow(k), 1, 7)
  averaged_clusters <- sapply(positions, function(position) {
    position_clustering <- clustering_with_duplicates[which(grepl(position, names(clustering_with_duplicates)))]
    cols <- rb_cols[position_clustering]

    if(length(cols) > 1) {
      mean_col <- mixcolors(cols)
    } else {
      mean_col = cols[[1]]
    }
    ## weird stuff
    names(mean_col) = NULL
    ##  names(mean_col) <- position
    return(mean_col)
  })

  cl <- clusterlist_from_membershiplist(averaged_clusters)

  return(cl)
  # TODO!
}

