library("dagitty")

# results is structured as follows:
# results$orig
# results$sub
# results$anc
# with respective sublists
# results$orig$graph$NEL       (only for orig, sub) 
# results$orig$graph$dagitty   (only for orig, anc)
# results$orig$localTests_results$r
# results$orig$localTests_results$r_interesting --> r_int
# results$orig$localTests_results$r_int_alpha --> r_int_sign
# results$orig$localTests_results$r_int_bad -- r_int_sign_bad
analysis_after_pc <- function(pc, data, outpath, protein, position_numbering, graph_layout = graph_layout, 
                              plot_as_subgraphs = FALSE, plot_only_subgraphs = NULL, 
                              coloring = coloring, colors = colors, stages = c("orig", "anc"), 
                              plot_types = c("localTests", "graphs"), unabbrev_r_to_info, print_r_to_console, 
                              lines_in_abbr_of_r, compute_localTests_anew = FALSE, print = TRUE, plot = TRUE, 
                              caption = "") {
  results <- list()
  results$pc <- pc
  
  results$orig <- list()
  results$orig$graph$NEL <- pc@graph
  
  graph_dagitty <- conv_to_r(pc@graph, type_of_graph = "pdag", nodename_prefix = "P")

  results$orig$graph$dagitty <- graph_dagitty
  
  if (("main" %in% stages) || ("orig" %in% stages)) {
  # if (originalgraph) {
    print("ORIGINAL GRAPH...")
    results$orig$localTests <- evaluate_DAG(data = data, graph = graph_dagitty, results = results, protein = protein, position_numbering = position_numbering, outpath = outpath, compute_localTests_anew = compute_localTests_anew)
  }
  
  if ("sub" %in% stages) {
  # if (subgraph) {
    print("SUBGRAPH...")
    # outpath_subgraph = paste(outpath, "-sub", sep = "") 
    
    subgraph <- subgraph_of_interesting_positions(pc@graph, protein = protein, position_numbering = position_numbering)
    # if (!((plotstages == "sub") && (plot_types == "graph"))) {
    #   garbage <- graphics.off()
    # }
    # plot(subgraph)
    results$sub$graph$NEL <- subgraph
    subgraph_dagitty <- conv_to_r(subgraph, type_of_graph = type_of_graph, nodename_prefix = "P")
    results$sub$localTests <- evaluate_DAG(data = data, graph = subgraph_dagitty, results = results, protein = protein, position_numbering = position_numbering, outpath = outpath, stage = "sub", plot = TRUE, compute_localTests_anew = compute_localTests_anew)
  }
  
  if ("anc" %in% stages) {
  # if (ancestorgraph) {
    print("ANCESTOR GRAPH...")
    # outpath_ancestor = paste(outpath, "-anc", sep = "") 
    
    ancestor_graph_dagitty <- ancestorgraph_of_interesting_positions(graph_dagitty = graph_dagitty, protein = protein, position_numbering = position_numbering, nodename_prefix = "P")
    if (!is.null(ancestor_graph_dagitty)) {
      # if (!((plotstages == "ancestor") && (plot_types == "graph"))) {
      #   garbage <- graphics.off()
      # }
      # plot(dagitty::graphLayout(ancestor_graph_dagitty))
      results$anc$graph$dagitty <- dagitty::graphLayout(ancestor_graph_dagitty)
      results$anc$localTests <- evaluate_DAG(data = data, graph = ancestor_graph_dagitty, results = results, protein = protein, position_numbering = position_numbering, outpath = outpath, stage = "anc", plot = TRUE, compute_localTests_anew = compute_localTests_anew) 
    } else {
      print("Ancestor graph is NULL!")
    }
  }
  
  if (plot) {
    plots(results, stages, plot_types, graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs, coloring = coloring, colors = colors, caption = caption, outpath = outpath)
    plots(results, stages, plot_types, graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs, coloring = coloring, colors = colors, caption = caption, outpath = "")
  }
  if (print) {
    results$r_statistics <- print_evaluation_results_to_info_file(results = results, outpath = outpath, stages, unabbrev_r_to_info = unabbrev_r_to_info, print_r_to_console = print_r_to_console, lines_in_abbr_of_r = lines_in_abbr_of_r)
  }
  return(results)
}


# only for Gauss-distributed data
# bad_alpha_threshold: alpha above which the estimate given by localTests(...) is not regarded
# estimate_margin: margin around zero, in which estimates given by localTests(...) are regarded "bad"
evaluate_DAG <- function(data, graph, results, protein, position_numbering = NULL, outpath, stage = "orig", plot = TRUE, 
                         bad_alpha_threshold = 0.01, estimate_margin = 0.4, compute_localTests_anew) {
  if (!(stage == "" || is.null(stage))) {
    outpath <- paste(outpath, stage, sep = "-")
  }  
  if ((file.exists(paste(outpath, "-localTests.RData", sep = ""))) && !(compute_localTests_anew)) {
    filename <- paste(outpath, "-localTests.RData", sep = "")
    load(filename)
    if (!exists("r")) {
      print("The file did not contain an object of name 'r'!")
    } else {
      print(paste("Result of localTests() loaded from ", filename, ".", sep = ""))
      save(r, file = paste(outpath, "-r.RData", sep = ""))
    }
  } else if ((file.exists(paste(outpath, "-r.RData", sep = ""))) && !(compute_localTests_anew)) {
    filename <- paste(outpath, "-r.RData", sep = "")
    load(filename)
    if (!exists("r")) {
      print("The file did not contain an object of name 'r'!")
    } else {
      print(paste("Result of localTests() loaded from ", filename, ".", sep = ""))
    }
  } else {
    d <- data.frame(data, check.names = FALSE)
    colnames(d) <- paste("P", colnames(d), sep="")
    rownames(d) <- paste("R", rownames(d), sep="")
    
    r <- localTests(graph, d)
    save(r, file = paste(outpath, "-r.RData", sep = ""))
    print("r saved.")
  }
  
  result_position = paste("localTests", stage, "r", sep = "_")
  results[[result_position]] <- r
  
  r$p.value <- p.adjust(r$p.value)
  
  interesting_pos <- interesting_positions(protein, position_numbering, allpositions = colnames(data))
  
  if (length(interesting_pos >= 2)) {
    if ((file.exists(paste(outpath, "-r_int.RData", sep = ""))) && !(compute_localTests_anew)) {
    filename <- paste(outpath, "-r_int.RData", sep = "")
    load(filename)
    if (!exists("r_int")) {
      print("The file did not contain an object of name 'r_int'!")
    } else {
      print(paste("r_int loaded from ", filename, ".", sep = ""))
    }
  } else {
      interesting_pos <- paste("P", interesting_pos, sep = "")
      print("Computing r_int...")
        if (!is.null(interesting_pos)) {
          r_int = c()
          for (i in interesting_pos) {
            for (j in interesting_pos) {
                r_ij <- r[grepl(paste(i, "_|"), rownames(r), fixed=TRUE), ]
                r_ij <- r_ij[grepl(paste("|_", j), rownames(r_ij), fixed=TRUE), ]
                r_int <- rbind(r_int, r_ij)
            }
          }
        }
      print("Computed.")
      r_int <- r_int[order(rownames(r_int)), ] # fkt nicht, wenn r_int empty
      if (dim(r)[1] > 1000) {
        save(r_int, file = paste(outpath, "-r_int.RData", sep = ""))
        print("r_int saved.")
      }
    }
      result_position = paste("localTests", stage, "r_int", sep = "_")
      results[[result_position]] <- r_int
  } else {
    print("No interesting positons known, r_int = r")
    r_int <- r
  }
    
  r_int_signif <- r_int[r_int$p.value < bad_alpha_threshold, ]
  result_position = paste("localTests", stage, "r_int_signif", sep = "_")
  results[[result_position]] <- r_int_signif
  
  r_int_signif_bad_low <- r_int_signif[r_int_signif$estimate > estimate_margin, ]
  r_int_signif_bad_high <- r_int_signif[r_int_signif$estimate < -estimate_margin, ]
  
  r_int_signif_bad <- rbind(r_int_signif_bad_low, r_int_signif_bad_high)
  
  r_int_signif_bad <- r_int_signif_bad[order(rownames(r_int_signif_bad)), ]
  
  result_position = paste("localTests", stage, "r_int_signif_bad", sep = "_")
  results[[result_position]] <- r_int_signif_bad
  
  localTests_results = list(r, r_int, r_int_signif, r_int_signif_bad)
  names(localTests_results) = list("r", "r_int", "r_int_signif", "r_int_signif_bad")
  return(localTests_results)
}

# # noch nicht ausfuehrlich getestet
# # deprecated
# which_pairs_of_postions <- function(r, position = NULL) {
#   pairs <- unique(substr(rownames(r), 0, 14)) # wenn Knotennummern 4-stellig 
#   if (!is.null(position)) {
#     return(grep(paste(position), pairs, value=TRUE))
#   } else {
#     return(pairs)
#   }
# }

pairs_of_pos <- function(r) {
  split_rownames <- strsplit(rownames(r), " ")
  first <- lapply(split_rownames, function(x) {.subset(x, i = 1)})
  second <- lapply(split_rownames, function(x) {.subset(x, i = 3)})
  return(unique(paste(first, second, sep = " _||_ ")))
}

# plottypes: "graphs", "localTests", "both"
# plotstages: "main", "sub", "anc", "all"
plots <- function(results, stages, plot_types = c("localTests", "graphs"), graph_layout, 
                  plot_as_subgraphs = FALSE, plot_only_subgraphs = FALSE, 
                  coloring, colors, outpath = "", caption = "") {
  if (!(outpath == "" || is.null(outpath))) {
    # postscript(paste(outpath, "_analysis.ps",  sep = ""), paper="special", width=10, height=9)
    postscript(paste(outpath, "_analysis.ps",  sep = ""), paper="special", width=10, height=8)
  } else if (!combined_plot) {
    garbage <- graphics.off()
  }
  
  if (!combined_plot) {
    old_par <- par()
    par(mfrow = c(length(stages), length(plot_types)))
  }
  
  for (stage in stages) {
    if ("localTests" %in% plot_types) {
      localTests_results <- results[[stage]]$localTests
      
      r <- localTests_results[[1]]
      r_int <- localTests_results[[2]]
      r_int_signif_bad <- localTests_results[[4]]
      
      #TODO: plotten wenn leer fktiert nicht! bzw. es sein lassen
      
      if (dim(r_int_signif_bad)[1] > 0) {
        plotLocalTestResults(r_int_signif_bad, main = "r_int_signif_bad")
      }  else if (dim(r_int)[1] > 0) {
        plotLocalTestResults(r_int, main = "r_int")
      } else if (dim(r)[1] > 0) {
        plotLocalTestResults(r, main = "r")
      }
    }
    
    if ("graphs" %in% plot_types){
      if (!is.null(results[[stage]]$graph$NEL)) {
        plot_graph(graph = results[[stage]]$graph$NEL, caption = caption, protein = protein, position_numbering = position_numbering, 
                   graph_layout = graph_layout, plot_as_subgraphs = plot_as_subgraphs, plot_only_subgraphs = plot_only_subgraphs, 
                   coloring = coloring, colors = colors, outpath = "")
        
        # colors <- colors_for_graph(protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors)
        # nAttrs <- list()
        # nAttrs$fillcolor <- colors
        # # nAttrs <- colors_for_graph(protein = protein, position_numbering = position_numbering, coloring = coloring, colors = colors)
        # # node shapes (pie)
        # drawNode <- node_function_for_graph(!is.null(coloring) && (coloring == "pie"))
        # pc_graph <- agopen(results[[stage]]$graph$NEL, layoutType = graph_layout, nodeAttrs = nAttrs, name = "pc") # circle produziert cluster
        # plot(pc_graph, nodeAttrs = nAttrs, drawNode = drawNode)  
      } else if (!is.null(results[[stage]]$graph$dagitty)) {
        plot(results[[stage]]$graph$dagitty, main = current_alpha)
      } else {
        print("No graph found.")
      }
    }
  }
  if (!(outpath == "" || is.null(outpath))) {
    dev.off()
  }
  
  if (!combined_plot) {
    par(old_par)
  }
  
  return(results)
}

# if print_r_to_console:
## if (print_exactly_one_r_short): r_int_signif_bad is printed to console, or, if empty, r_int_signif, 
### or otherwise r_int, otherwise r (all abbridged, when longer than lines_in_abbr_of_r lines)
## otherwise: r_int_signif_bad is printed to console (abbridged, when longer than lines_in_abbr_of_r lines)
# otherwise, only the statistics is printed
print_evaluation_results_to_info_file <- function(results, outpath, stages, unabbrev_r_to_info = TRUE, print_r_to_console, lines_in_abbr_of_r, print_exactly_one_r_short = TRUE) {
  print("Writing to info-file.")
  sink(file = paste(outpath, "-info.txt", sep = ""), append = TRUE)
  old_width <- options()$width
  options("width" = 200) # alle columns nebeneinander printen
  print("ANALYSIS")
  # cat("\n")
  
  r_statistics <- list()
  r_print <- list()
  for (stage in stages) {
    print(paste(toupper(stage), ":", sep = ""))
    cat("\n")
    
    localTests_results <- results[[stage]]$localTests
    names(localTests_results) <- c("r", "r_int", "r_int_signif", "r_int_signif_bad")
    
    stage_table <- matrix(c(
      sapply(localTests_results, function(x) {return(dim(x)[1])}),
      sapply(localTests_results, function(x) {return(dim(x[x$estimate < 0,])[1])}),
      sapply(localTests_results, function(x) {return(dim(x[x$estimate > 0,])[1])}),
      sapply(localTests_results, function(x) {return(mean(x$estimate))})
    ), nrow = 4)
    
    colnames(stage_table) <- c("# elements", "# negative elements", "# positive elements", "mean")
    rownames(stage_table) <- c("r", "r_int", "r_int_signif", "r_int_signif_bad")
    
    r_statistics[[stage]] <- stage_table
    
    if (print_exactly_one_r_short) {
      # r_print bestimmen
      if (dim(localTests_results["r_int_signif_bad"][[1]])[1] > 0) {
        r_print[[stage]]$name <- "r_int_signif_bad"
        r_print[[stage]]$value <- localTests_results["r_int_signif_bad"][[1]]
      } else if (dim(localTests_results["r_int_signif"][[1]])[1] > 0) {
        r_print[[stage]]$name <- "r_int_signif"
        r_print[[stage]]$value <- localTests_results["r_int_signif"][[1]]
      } else if (dim(localTests_results["r_int"][[1]])[1] > 0) {
        r_print[[stage]]$name <- "r_int"
        r_print[[stage]]$value <- localTests_results["r_int"][[1]]
      } else if (dim(localTests_results["r"][[1]])[1] > 0) {
        r_print[[stage]]$name <- "r, sorted for estimate"
        r_print[[stage]]$value <- localTests_results["r"][[1]][order(localTests_results["r"][[1]]$estimate, decreasing = TRUE), ]
      }
    }
    
    if (unabbrev_r_to_info) {
      options(max.print = 9999999)
      if (dim(localTests_results["r_int_signif_bad"][[1]])[1] > 0) {
        print(paste(stage, "r_int_signif_bad:", sep = " -- "))
        print(localTests_results["r_int_signif_bad"])
        print("pairs of positions in r_int_signif_bad:")
        cat(paste(pairs_of_pos(localTests_results["r_int_signif_bad"][[1]]), collapse = "\n"))
        cat("\n")
      } else {
        print("r_int_signif_bad: -empty-")
      }
      cat("\n")
      if (dim(localTests_results["r_int"][[1]])[1] > 0) {
        print(paste(stage, "r_int:", sep = " -- "))
        print(localTests_results["r_int"])
        print("pairs of positions in r_int:")
        cat(paste(pairs_of_pos(localTests_results["r_int"][[1]]), collapse = "\n"))
        cat("\n")
      } else {
        print("r_int: -empty-")
      }
      cat("\n")
      if (dim(localTests_results["r"][[1]])[1] > 0) {
        print(paste(stage, print("r:"), sep = " -- "))
        print(localTests_results["r"])
        print("pairs of positions in r:")
        cat(paste(pairs_of_pos(localTests_results["r"][[1]]), collapse = "\n"))
        cat("\n")
      } else {
        print("r: -empty-")
      }
      options(max.print = 99999)
    } else {
      # print_r_short(r_print[[stage]]$value, r_print[[stage]]$name, lines_in_abbr_of_r = lines_in_abbr_of_r)
      print_r_short(localTests_results["r"][[1]], "r", lines_in_abbr_of_r = lines_in_abbr_of_r)
      print_r_short(localTests_results["r_int"][[1]], "r_int", lines_in_abbr_of_r = lines_in_abbr_of_r)
      print_r_short(localTests_results["r_int_signif"][[1]], "r_int_signif", lines_in_abbr_of_r = lines_in_abbr_of_r)
      print_r_short(localTests_results["r_int_signif_bad"][[1]], "r_int_signif_bad", lines_in_abbr_of_r = lines_in_abbr_of_r)
      cat("\n")
    }
    cat("\n")
  }
  print("STATISTICS:", quote = FALSE)
  print(r_statistics)
  options("width" = old_width) # back to default
  sink()
  print("Written.")
  
  # cat("\n")
  print("RESULTS:")
  if (print_r_to_console) {
    for (stage in stages) {
      print(paste(toupper(stage), ":", sep = ""), quote = FALSE)
      print_r_short(r_print[[stage]]$value, r_print[[stage]]$name, lines_in_abbr_of_r = lines_in_abbr_of_r)
      cat("\n")
    }
  }
  
  print("STATISTICS:", quote = FALSE)
  print(r_statistics)
  
  return(r_statistics)
}

print_r_short <- function(r_print, name, lines_in_abbr_of_r) {
  if (dim(r_print)[1] < lines_in_abbr_of_r) {
    print(paste(name, ":", sep =""))
    print(r_print)
  } else {
    print(paste(name, ", first ", lines_in_abbr_of_r, " of ", dim(r_print)[1]," values:", sep = ""))
    print(r_print[1:lines_in_abbr_of_r,])
    cat("...")
    cat("... all pairs of positions:\n")
    print(pairs_of_pos(r_print))
  }
}

# print_r_short <- function(r_print, stage) {
#   if (is.null(r_print[[stage]]$pairs)) {
#     print(paste(r_print[[stage]]$name, ":", sep =""))
#     print(r_print[[stage]]$value)
#   } else {
#     print(paste(r_print[[stage]]$name, ", first ", dim(r_print[[stage]]$value)[1], " of ", r_print[[stage]]$orgi_number_lines," values:", sep = ""))
#     print(r_print[[stage]]$value)
#     cat("...")
#     cat("... all pairs of positions:\n")
#     print(r_print[[stage]]$pairs)
#   }
# }

td_analysis <- function(graph, outpath) {
  if (file.exists(paste0(outpath, "-td.td"))) {
    tree_graphnel <- td_to_graphNEL(g = graph, path = outpath)
    graphNEL_to_tikz(tree_graphnel, outpath, suffix = "-td")
  } else {
    conv_to_DIMACS(g = graph, outpath = outpath)
    print(paste("Max", paste0(outpath, "-graph-DIMACS.txt"), "geben!"))
  }
}