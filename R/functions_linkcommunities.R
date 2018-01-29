library(linkcomm)
library(igraph)
library(colorspace)  # for mixcolor, hex

node_function_for_graph <- function(pie = FALSE, counts, colors = c("#1874CD", "#69A019", "#FFD700", "#CC0000")) {
  if (pie) {
    if (missing(counts)) {
      counts <- read.csv2("Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
      counts <- as.matrix(counts, ncol = 4)
    }
    if (missing(colors)) {
      colors <- colnames(counts)
    }
    
    makeNodeDrawFunction <- function(x) {
      force(x)
      if (sum(x) == 0) {
        return(drawAgNode)
      }
      function(node, ur, attrs, radConv) {
        nc <- getNodeCenter(node)
        name <- name(node)
        # cols <- c("#1874CD", "#69A019", "#FFD700", "#CC0000")
        cols <- colors
        pieGlyph(x[x > 0], 
                 xpos=getX(nc),
                 ypos=getY(nc),
                 labels=name,
                 radius=getNodeRW(node),
                 col=cols[x > 0])
      }
    }
    drawFuns <- apply(counts, 1, makeNodeDrawFunction)
  } else {
    drawFuns <- drawAgNode
  }
  return(drawFuns)
}

# get clusters with meaningful names of members
get_meaningful_names_for_clustering <- function(clusters, edge_list) {
  lapply(clusters, function(cluster) return(edge_list[cluster]))
}

node_colors_pie <- function(clustering, nodes) {
  pie_counts <- sapply(as.list(nodes), function(node) return(
    sapply(clustering, function(edges) return(
      length(edges[grepl(paste("\\<", node, "\\>", sep = ""), edges)])))))
  # for (node in nodes) {
  #   #grepl grepl "16" in "160" :(
  #   lapply(clustering, function(edges) return(length(edges[grepl(node, edges)])))
  # }
  colnames(pie_counts) <- nodes
  pie_counts <- prop.table(t(pie_counts), 1)
  rownames(pie_counts) <- nodes
  return(pie_counts)
}

clusterlist_from_membershiplist <- function(membershiplist, cluster_names) {
  clusters_cut_cuttree <- membershiplist
  if (missing(cluster_names)) {
    cluster_names <- as.list(unique(membershiplist))#list(1,2,3,4)
  }
  clusters <- lapply(cluster_names, function(cluster_name) return(as.numeric(names(clusters_cut_cuttree)[clusters_cut_cuttree == cluster_name])))
  names(clusters) <- cluster_names
  return(clusters)
}

membershiplist_from_clusterlist <- function(clusterlist) {
  elements <- sort(unique(unlist(clusterlist)))
  membershiplist <- sapply(elements, function(element) {
    return(which(sapply(clusterlist, function(cluster){element %in% cluster})))})
  return(membershiplist)
}

# if cut_k, the dendrogram is cut in such a way that k clusters arise (otherwise, k is ignored)
compute_link_communities <- function(graph, k, base_colors, plot_bar_plot = FALSE, 
                                     plot_colored_graph = TRUE, classify_nodes = TRUE, round_categories = 2,
                                     pie_nodes = TRUE, color_edges = TRUE, protein, outpath, 
                                     graph_output_formats = c("ps", "svg", "pdf")) {
  # graph <- results$orig$graph$NEL
  edge_list <- as_edgelist(igraph.from.graphNEL(graph))
  rownames(edge_list) <- paste(edge_list[,1], edge_list[,2], sep = "-")
  link_comm <- getLinkCommunities(edge_list, hcmethod = "ward", plot = TRUE)
  # print(get_meaningful_names_for_clustering(link_comm$clusters))
  # clusters_cut_cutDendrogram <- cutDendrogramAt(link_comm$hclust, lc = link_comm, cutat = 0.9982)
  if (!(missing(k) || is.null(k))) {
    # base_colors <- base_colors[1:k]
    clusters_cut_cuttree <- cutree(link_comm$hclust, k = k)
    clusters <- clusterlist_from_membershiplist(clusters_cut_cuttree)
    # cluster_names <- list(1,2,3,4)
    # clusters <- lapply(cluster_names, function(cluster_name) return(as.numeric(names(clusters_cut_cuttree)[clusters_cut_cuttree == cluster_name])))
  } else {
    clusters <- link_comm$clusters
  }
  
  if ((length(base_colors) == 0) || (!(length(clusters) <= length(base_colors)))) {
    base_colors <- rainbow(length(clusters))
  } else {
    base_colors <- base_colors[1:length(clusters)]
  }
  
  names(clusters) <- base_colors
  
  clustering <- get_meaningful_names_for_clustering(clusters = clusters, edge_list = edgeNames(graph))
  # edgeNames(graph)
  names(clustering) <- base_colors
  print("edge clusters:")
  print(clustering)
  
  nodes <- graph@nodes
  cols <- node_colors_pie(clustering, nodes)
  # cols <- na.omit(cols)
  cols[is.na(cols)] <- 0
  # print(cols)
  
  if (plot_bar_plot) {
    cols_sort <- cols[,order(apply(cols, 2, sum, decreasing = TRUE))]
    # counts_sort <- counts[,order(apply(counts, 2, sum, decreasing = TRUE))]
    # par(mfrow = c(1,2))
    barplot(t(cols_sort))
    # barplot(t(counts_sort))
  }
  
  if (plot_colored_graph) {
    if (classify_nodes) {
      node_clusters <- classify_nodes(cols, round_categories = round_categories, mix = TRUE, base_colors = base_colors)
      print("node clusters:")
      print(node_clusters)
      fillcolor <- colors_for_nodes(node_clusters = node_clusters, colors = names(node_clusters), clustering = TRUE)
    } else {
      fillcolor <- NULL
    }
    
    if (pie_nodes) {
      drawnode <- node_function_for_graph(pie = TRUE, counts = cols, colors = base_colors)
    } else {
      drawnode <- drawAgNode
    }
      # 
      # plot_graph(graph = graph, fillcolor = NULL, drawnode = drawnode, numerical = TRUE, outpath = "")
    if (color_edges) {
      edgecolor <- colors_for_edges(clustering = clusters, colors = base_colors, graph = graph)
    } else {
      edgecolor <- NULL
    }
      # plot(graph, edgeAttrs = eAttrs)
    }
    # plot_graph(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, numerical = TRUE, outpath = "")
    
    subgraphs <- subgraphs_from_node_clusters(node_clusters, graph)
    plot.new()
    plot_graph(graph = graph, fillcolor = fillcolor, edgecolor = edgecolor, drawnode = drawnode, 
               caption = "", numerical = TRUE, outpath = "", plot_as_subgraphs = TRUE, 
               subgraphs = subgraphs, output_formats = graph_output_formats)
  # }
    plot_clusters_in_pymol(node_clustering = node_clusters, protein = protein, outpath = outpath, type_of_clustering = "linkcomm")
  return(cols)
}

# returns a vector with the colnumber, in which the input matrix is maximal for each row
max_cluster_linkcomm <- function(counts) {
  return(apply(counts, 1, function(row) return(which(row == max(row)))))
}

subgraphs_from_node_clusters <- function(node_clusters, graph, protein = "") {
  if (protein == "PDZ" || protein == "pdz") {
    node_clusters[[1]] <- c(node_clusters[[1]], node_clusters[[2]])
    node_clusters[[2]]  <- NULL
  }
  # return(lapply(node_clusters, function(node_vector) {return(list(graph=subGraph(as.character(node_vector), graph)))}))
  return(lapply(node_clusters, function(node_vector) {return(list(graph=subGraph(sapply(node_vector, as.character), graph)))}))
}

# colors must be given when mix
classify_nodes <- function(counts, round_categories = 4, plot = TRUE, base_colors = NULL, mix = TRUE) {
  if (mix && is.null(base_colors)) {
    warning("No colors given, all nodes will be white.")
  }
  # print(counts)
  # rounded_counts <- round(counts * round_categories)
  rounded_counts <- floor((counts * round_categories) + 0.5)  # ceil 0.5, which round does not do!
  colnames(rounded_counts) <- base_colors  # sollte schon passiert sein
  print(rounded_counts)
  
  if (mix) {
    nonzero_positions <- function(vector) {
      ret <- rep(colnames(rounded_counts)[vector != 0], vector[vector != 0])
      if (is.null(ret) || length(ret) == 0) {
        ret <- "#FFFFFF"
      }
      if (length(ret) == 1) {
        return(ret)
      } else {
        return(mixcolors(ret))
      }
    }
    categories <- apply(rounded_counts, 1, nonzero_positions)
    # print(categories)
  } else {  # ever used?
    if (round_categories < 2) {
      round_categories <- 2
    }
    base_round_categories <- rounded_counts * round_categories ^ (col(rounded_counts) - 1)
    # print(base_round_categories)
    categories <- apply(base_round_categories, 1, sum)
    # print(categories)
    # base_round_categories <- cbind(base_round_categories, apply(base_round_categories, 1, sum))
    # print(base_round_categories)
    
    
  }
  clusters <- clusterlist_from_membershiplist(categories)
  # print(clusters)
  return(clusters)
} 

mixcolors <- function(color_vector) {
  color <- mixcolor(0.5, hex2RGB(color_vector[1]), hex2RGB(color_vector[2]))
  n <- length(color_vector)
  alphas <- 1/seq(1:n)
  if (n > 2) {
    for (i in 3:length(color_vector)) {
      color <- mixcolor(alpha = as.numeric(alphas[i]), color, hex2RGB(color_vector[i]))
    }
  }
  return(hex(color))
}

# hex(mixcolor(0.5, hex2RGB("#FFD700"), hex2RGB("#FFD700")))

