#' Cluster by Causal Structure (Graph Clustering)
#'
#'
#' @return node clustering as an igraph clustering
#' @seealso node_clustering_from_igraph_clustering
protein_graph_clustering <- function(graph, clustering) {
  igraph_cluster <- igraph.from.graphNEL(graph)
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
  return(cl)
}

#' transform an igraph clustering to a node clustering (format?)
#' @param communities_clustering the clustering that is to be transformed
#' @param unname = TRUE Should the names of the resulting object should be removed, to avoid that they are interpreted as colors by subsequent functions?
node_clustering_from_igraph_clustering <- function(communities_clustering, unname = TRUE) {
  node_clustering <- groups(communities_clustering)
  if (unname) {
    node_clustering <- unname(node_clustering)  # otherwise interpreted as colors
  }
  return(node_clustering)
}

#' Plot graph clustering in graph and pymol
#'
#' @param node_clustering the node clustering as a list of lists (?); if the object is named, the names ar interpreted as colors (e.g. "#000000")
#' @param clustering_type a string that decribes the kind of clustering (e.g. the clustering algorithm); used for e.g. the outpath
#' @param graph the graph in which the clustering is to be highlighted; if missing (and communities_clustering missing) or NULL, this step is skipped
#' @param outpath the outpath to which output files should be written
#' @param protein the protein that the data, that the computation relies on, belong to
#' @param mute_all_plots boolean which can be set to avoid that the plot (cluster highlighting in the graph) is shown in R (it is written to the outpath anyway)
#' @param caption caption for the graph plot
#' @param file_separator a character the should be used to speparate directories (e.g. "/" or "\", depending on the operating system).
#' This is used in the path to the pdb-file in the output pymol file.
#' @param remove_singular_clusters = TRUE Should clusters consisting of only one element be removed? (If all clusters are singular, this is not done.)
#' @param merge_singular_clusters = FALSE Should all clusters that consist of only one element be merged into one single cluster? (Only if !remove_singular_clusters)
#' @param sort_clusters = length either a function which can be applied to the clusters and by which result values the clusters are to be sorted in decreasing order
#' or "DDS-SVD" (...)
#' @param ...  parameters passed on to call_plot_igraph for the highlighting in the graph
output_node_clustering <- function(node_clustering, #clustering_type,
                                   graph,
                                   communities_clustering,
                                   outpath, protein, file_separator,
                                   remove_singular_clusters = TRUE, merge_singular_clusters = FALSE,
                                   sort_clusters = length, additional_clusters,
                                   mute_all_plots, ...) {

  if (!missing(communities_clustering)) {
    if (missing(node_clustering)) {
      node_clustering <- node_clustering_from_igraph_clustering(communities_clustering = communities_clustering, unname = TRUE)
    }
    if (missing(graph)) {
      graph <- communities_clustering$graph
    }
    # if (missing(clustering_type)) {
    #   clustering_type <- communities_clustering$algorithm
    # }

  } else {
    if (!missing(graph)) { # otherwise, nnno plotting is done and thus, the communities clustering object is not needed
      communities_clustering <- make_clusters(graph = igraph.from.graphNEL(graph),
                                              membership =  membershipvector_from_clusterlist(node_clustering),
                                              modularity = FALSE)
    }
    # if (missing(clustering_type)) {
    #   warning("No clustering_type given in output_node_clustering! It was set to \"\".")
    # }
  }


  if (remove_singular_clusters) {
    node_clustering <- remove_singular_clusters(node_clustering)
  } else if (merge_singular_clusters) {
    node_clustering <- merge_singular_clusters(node_clustering)
  }

  if (is.null(names(node_clustering))) {
    clustering_colors <- rainbow(length(node_clustering))
    names(node_clustering) <- clustering_colors
  } else {
    clustering_colors <- names(node_clustering)
  }


  ## TODO: save plot, instaed of plotting
  if (!missing(graph) && !is.null(graph)) {
    igraph_plot <- graph
    ## edge.arrow.size determines size of arrows (1 is default), vertex.size determines size of the vertices (15 is default), edge.width determines width of edges (1 is default)
    ## plot(cl, igraph, main = paste0(caption, "\n", clustering), edge.arrow.size=0.2, vertex.size=8, edge.width=0.7)
    # call_plot_igraph(g = igraph_plot, clusters = TRUE, cluster_str = paste0(clustering_type, "-graph"),
    #                  clustering = communities_clustering,
    #                  clustering_colors = sapply(clustering_colors, adjustcolor, alpha.f = 0.1),
    #                  outpath = paste0(outpath,"-", length(node_clustering), "_", "clusters"), ...)
    call_plot_igraph(g = igraph_plot, clusters = TRUE, #cluster_str = clustering_type,
                     clustering = communities_clustering,
                     clustering_colors = sapply(clustering_colors, adjustcolor, alpha.f = 0.1),
                     outpath = function_set_parameters(outpath, list(object = "graph")),
                     mute_all_plots = mute_all_plots, ...)
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


  if (!(missing(additional_clusters))) {
    node_clustering <- c(node_clustering, additional_clusters)
  }


  plot_clusters_in_pymol(node_clustering = node_clustering, protein = protein, outpath = outpath,
                         file_separator = file_separator, type_of_clustering = clustering_type)
}

#' @param monochromatic_removed_cols all removed_cols are given (the same shade of) white as color
#' @param more_levels_of_conservedness (because of too little variacne) removed positions are
#' clored in shades of white with increasing amounts of blue with increasing variance
add_clusters <- function(add_cluster_of_conserved_positions = FALSE, removed_cols,
                         more_levels_of_conservedness = FALSE, monochromatic_removed_cols = TRUE) {
  # is this for coloring intesity by variance (an those below min_pos_var are removed)?
  if (add_cluster_of_conserved_positions) { # should this be add_cluster_of_REMOVED_positions?
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
}

#' Remove all clusters with just one element from the clustering.
#'
#'@param clustering the clustering from which all singular clusters are to be removed.
#'@param force = FALSE Should the clusters be removed, even if the resulting clustering is empty, that is, all clusters are singular?
#'If FALSE, the default, the original clustering is returned in that case.
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

remove_effects_of_isolated_nodes <- function(all_pairwise_effects, isolated_nodes) {
  if (missing(isolated_nodes)) {
    # isolated_rows_names <- sapply(rownames(all_pairwise_effects), function(rowname) {
    #   if (norm(all_pairwise_effects[rowname, ], type="2") == 1) {
    #     return(rowname)
    #   }
    # })
    rows_isolated <- apply(all_pairwise_effects, 1, function(row) {
      if (norm(row, type="2") == 1) {
        return(TRUE)
      }
      return(FALSE)
    })
    nodes_isolated <- sapply(colnames(all_pairwise_effects), function(node_name) {
      rownames_for_node <- rownames(all_pairwise_effects)[grep(paste0("^", node_name, "-|^", node_name, "$"), rownames(all_pairwise_effects))]
      # rownames_for_node <- rownames(all_pairwise_effects)[grep(node_name, rownames(all_pairwise_effects))]
      if (all(rows_isolated[rownames_for_node])) {
        return(TRUE)
      }
      return(FALSE)
    })
  } else {
    nodes_isolated <- colnames(all_pairwise_effects) %in% isolated_nodes
    names(nodes_isolated) <- colnames(all_pairwise_effects)
    rows_isolated <- sapply(rownames(all_pairwise_effects), function(row_name) {
      if (grepl("-", row_name)) {
        node <- str_split(row_name, "-")[[1]][1]
      } else {
        node <- row_name
      }
      return(row_name = unname(nodes_isolated[as.character(node)]))
    }, USE.NAMES = TRUE)
  }
  all_pairwise_effects <- all_pairwise_effects[which(!rows_isolated), which(!nodes_isolated)]

}

#' Cluster by Causal Effects
#' @param pairwise_effects a matrix with one col for each node that is to be clustered. There might be more rows than cols,
#'  e.g. because there are several possible effects (and therefore more than one line) for a node
#' @param pre_fct name of a function that is to be applied to the effects before clustering (e.g. "cor" or "cov); default is "identity"
#' @param remove_isolated_nodes if set, the effects are filtered before the clustering, to remove the lines and column for each position
#'  for which there is no line of effects that effect any other position than the position itself.
#'  (Which means that the node is isolated in the cuasal structure.)
#'  (Actually it is checked whether the lies' 2-norm is 1. This is equivalent because by definition, each node must have the effect 1 on itself,
#'  and thus all other positions must be zero in this case.)
#'
cluster_pairwise_effects <- function(results, pairwise_effects, pre_fct = "identity",
                                     cluster_method, hclust_method,
                                     dist_measure,
                                     number_of_clusters_k,
                                     cut_height_h, alpha = 0.05,
                                     iterations_pv, protein, outpath, file_separator,
                                     remove_isolated_nodes = TRUE,
                                     output_formats, mute_all_plots, ...) {

  if (remove_isolated_nodes) {
    pairwise_effects <- remove_effects_of_isolated_nodes(all_pairwise_effects = pairwise_effects)
    # if ((dim(pairwise_effects)[1] != dim(pairwise_effects)[2]) ||
    #     (sum(sapply(1:92, function(i) {return(pairwise_effects[i,i])})) != dim(pairwise_effects)[1])) { # sum of the main diagonal is other than sum of the rows and cols
    #   warning("Pairwise effects is not a square matrix with value 1 on the main diagonal. Removal of isolated nodes might not work.")
    # }
    # isolated_cols <- apply(pairwise_effects, 1, function(v) (norm(v, type="2") != 0))
    # pairwise_effects <- apply(pairwise_effects, 1, function(v) {if (norm(v, type="2") != 0) {return(v)} else return(NULL)})
    # pairwise_effects <- do.call(rbind, apply(pairwise_effects, 1, function(v) {if (norm(v, type="2") != 0) {return(v)}}))
  }

  pre_cluster_fct <- get(pre_fct)
  pairwise_effects <- pre_cluster_fct(pairwise_effects)
  outpath <- paste0(outpath, "-prefct=", pre_fct)



  results$effects_clustering <- list()
  #k-means
  if (grepl(pattern = "k", cluster_method) && grepl(pattern = "means", cluster_method)) {
    k_m <- kmeans(t(pairwise_effects), number_of_clusters_k)
    lapply(seq(1:number_of_clusters_k), function (i) {dim(pairwise_effects[,names(which(k_m$cluster == i))])})

    # apply(pairwise_effects[,names(which(k_m$cluster == 1))], 2, barplot);

    results$pairwise_effects <- pairwise_effects

    clustering_with_duplicates <- k_m$cluster

    cl <- position_clustering_from_clustering_with_duplicates(clustering_with_duplicates = clustering_with_duplicates)

    print(cl)
    results$effects_clustering$k_means <- cl

    type <- "effects-km"
    # names(cl) <- NULL
    plot_clusters_in_pymol(node_clustering = cl, protein = protein, outpath = outpath,
                           file_separator = file_separator, type_of_clustering = type,
                           length_sort = TRUE)
  } else {
    FUN_pv <- function_set_parameters(pvclust, parameters = list(data = pairwise_effects,
                                                                 method.hclust = hclust_method,
                                                                 method.dist = dist_measure, nboot = iterations_pv))

    effects_pv <- compute_if_not_existent(filename = paste(outpath, "pv", hclust_method,
                                                           substr(dist_measure, 0, 3), iterations_pv, sep="-"),
                                          FUN = FUN_pv,
                                          obj_name = "effects_pv",
                                          fun_loaded_object_ok = function(effects_pv) {return(all(colnames(pairwise_effects) == effects_pv$hclust$labels))}
    )


    # fit <- pvclust(data = pairwise_effects, method.hclust="ward",
    #                method.dist="euclidean", nboot = 1000)

    # cl <- cutDendrogramAt(effects_pv$hclust, cutat = 15)
    # membershipvector_from_clusterlist(cl)


    # TODO: caption
    if ((missing(number_of_clusters_k) || is.null(number_of_clusters_k)) && (missing(cut_height_h) || is.null(cut_height_h))) {
            high_clusterlist <- pvpick(effects_pv, alpha = alpha)$clusters

      plot.new()
      plot(effects_pv, hang = -1)


      alpha_changed = FALSE
      while (is.null(high_clusterlist[[1]])){
        alpha = alpha - 0.1
        high_clusterlist <- pvpick(effects_pv, alpha = alpha)$clusters
        alpha_changed = TRUE
      }
      if (alpha_changed) {
        warning(paste("Alpha was reduced to", alpha, "because no clusters could be found otherwise." ))
      }

      high <- membershipvector_from_clusterlist(high_clusterlist)
      colors <- rainbow(length(unique(high)))
      # add rectangles around groups highly supported by the data
      pvrect(effects_pv, alpha = alpha, border = colors)
      type <- paste("effects-pv", hclust_method, substr(dist_measure, 0, 3), iterations_pv, paste0("iter-alpha=", alpha), sep="-")
      outpath <- paste0(outpath, "-", number_of_clusters_k, pastes("clusters", type, sep = "-"))
      output_dendrogram(cluster_obj = effects_pv, alpha = alpha,
                        clusters = high, colors = colors, caption = caption,
                        outpath = outpath, output_formats = "pdf",
                        mute_all_plots = mute_all_plots)
    } # else {
      # effects_pv <- hclust(d = pairwise_effects, method = hclust_method)
      else if (!(missing(number_of_clusters_k) || is.null(number_of_clusters_k))) {
        high <- cutree(effects_pv$hclust, k = number_of_clusters_k)
        colors <- rainbow(number_of_clusters_k)
        type <- paste("effects-pv", hclust_method, substr(dist_measure, 0, 3), iterations_pv, paste0("iter-k=", number_of_clusters_k), sep="-")
        # rect.hclust(effects_pv$hclust, k = number_of_clusters_k, cluster = high, border = colors)
        outpath <- paste0(outpath, "-", number_of_clusters_k, pastes("clusters", type, sep = "-"))
        output_dendrogram(cluster_obj = effects_pv, k = number_of_clusters_k,
                          clusters = high, colors = colors, caption = caption,
                          outpath = outpath, output_formats = "pdf",
                          mute_all_plots = mute_all_plots)
      } else if (!(missing(cut_height_h) || is.null(cut_height_h))) {
        high <- cutree(effects_pv$hclust, h = cut_height_h)
        colors <- rainbow(length(unique(high)))
        type <- paste("effects-pv", hclust_method, substr(dist_measure, 0, 3), iterations_pv, paste0("iter-h=", cut_height_h), sep="-")
        outpath <- paste0(outpath, "-", number_of_clusters_k, pastes("clusters", type, sep = "-"))
        # abline(h = cut_height_h, lty = 2)
        # rect.hclust(effects_pv$hclust, h = cut_height_h, cluster = high, border = colors)
        output_dendrogram(cluster_obj = effects_pv, h = cut_height_h,
                          clusters = high, colors = colors, caption = caption,
                          outpath = outpath, output_formats = "pdf",
                          mute_all_plots = mute_all_plots)
      }
    # }

    node_clustering <- position_clustering_from_clustering_with_duplicates(clustering_with_duplicates = high)

    print(node_clustering)

    results$effects_clustering$pv[[type]] <- node_clustering

    # names(cl) <- NULL

    # if (remove_singular_clusters) {
    #   node_clustering <- remove_singular_clusters(node_clustering)
    # } else if (merge_singular_clusters) {
    #   node_clustering <- merge_singular_clusters(node_clustering)
    # }


    # plot_clusters_in_pymol(node_clustering = node_clustering, protein = protein, outpath = outpath,
    #                        file_separator = file_separator, type_of_clustering = type,
    #                        length_sort = TRUE)
  }
  output_node_clustering(node_clustering = node_clustering, #clustering_type = type,
                         outpath = outpath,
                         protein = protein, file_separator = file_separator, sort_clusters = length,
                         mute_all_plots = mute_all_plots, ...)
  return(results)
}



position_clustering_from_clustering_with_duplicates <- function(clustering_with_duplicates, cluster_colors) {
  positions <- unique(sapply(names(clustering_with_duplicates), function(long_name) return(gsub("-.*","",long_name))))
  k = max(clustering_with_duplicates)


  ## averaged_clusters <- sapply(positions, function(position) {
  ##  position_clustering <- clustering_with_duplicates[which(grepl(position, names(clustering_with_duplicates)))]
  ##  mean_cluster <- mean(position_clustering)
  ## names(mean_cluster) <- position
  ##  return(mean_cluster)
  ##})

  ## averaged_clusters <- round(averaged_clusters)
  if (missing(cluster_colors)) {
    cluster_colors <- rainbow(k)
  }

  rb_cols = substr(cluster_colors, 1, 7)
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

output_dendrogram <- function(cluster_obj, clusters, colors, caption, alpha, h = NULL,
                              outpath = "", output_formats = "pdf", mute_all_plots = FALSE, ...) {
  if (!mute_all_plots) {
    plot.new()
    plot(cluster_obj, hang = -1)
    if (!missing(h)) {
      abline(h = h, lty = 2)
    }
    if (!missing(alpha)) {
      pvrect(cluster_obj, alpha = alpha, border = colors, ...)
    } else {
      rect.hclust(cluster_obj$hclust, h = h, cluster = clusters, border = colors, ...)
    }
  }

  for(format in output_formats) {
    if (!nchar(outpath) == 0) {
      if (format == "pdf") {
        pdf(paste(outpath, "-dendrogram", ".pdf", sep = ""))
      } else if ((format == "ps") || (format == "postscript")) {
        postscript(paste(outpath, "-dendrogram", ".ps", sep = ""), paper = "special", width = 10, height = 9, fonts=c("serif", "Palatino"))
      } else if(format == "svg") {
        svg(paste(outpath, "-dendrogram", ".svg", sep = ""))
      }

      plot(cluster_obj, hang = -1)
      if (!missing(h)) {
        abline(h = h, lty = 2)
      }
      if (!missing(alpha)) {
        pvrect(cluster_obj, alpha = alpha, border = colors, ...)
      } else {
        rect.hclust(cluster_obj$hclust, h = h, cluster = clusters, border = colors, ...)
      }
      dev.off()
    }
  }
}

