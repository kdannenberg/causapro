# library(dagitty)
# library(pcalg)
# library(graph)

# library(dagitty)
# library(pcalg)
# library(graph)

## gets a pcAlgo object and returns a corresponding dagitty r object
## nodename_prefix" is prepended to each node name
conv_to_r <- function(g, type_of_graph = "pdag", nodename_prefix = "") {
  ## work with graph class
  # g <- pc @ graph
  ## list of nodes
  ln <- names(nodeData(g))
  ## list of edges in format "a|b"
  le <- names(edgeData(g))
  ## construct adjacency matrix
  n <- length(ln)
  am <- matrix(0, nrow = n, ncol = n)
  ## this string contains a specification of a graph in the
  ## dot language
  ## it starts with the graphtype and the description
  ## of the actual graph is then between the curly brackets
  ## this is given as a list of nodes (with a prefix found
  ## in nodename_prefix) followed by an edge list
  ## given as [firstnode edgetype secondnode] where
  ## first- and secondnode are the same as previous nodes
  ## (and thus including the prefix)
  ## the edgetype is either -- for undirected or -> for directed edges
  ## note that it does not matter if the edges/nodes are
  ## seperated by spaces or newlines
  s <- paste(type_of_graph, "{", sep = "")
  ## print nodes
  for(node in ln) {
    s <- paste(s, nodename_prefix)
    s <- paste(s, node, sep = "")
    ## comment in to set variables to exposure
    ## s <- paste(s, "[exposure]")
  }
  ## go through all edges and mark them in adjacency matrix
  for(edge in le) {
    splt <- strsplit(edge, "\\|")
    ## this takes linear time and could be replaced by hashmaps
    from <- which(splt[[1]][1] == ln)[[1]]
    to <- which(splt[[1]][2] == ln)[[1]]
    am[[from, to]] <- 1
  }
  ## go through adjacency matrix and print edges
  ## we distinguish between bidirectional edges and unidirectional ones
  ## note that edges of the first kind will be given to dagitty as edges
  ## whose direction is unknown
  ## this might be changed later on
  for(i in 1:n) {
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        ## edge is bidirectional (second check needed to print edge only once)
        if(am[[j, i]] == 1 && i > j) {
          otp <- paste(nodename_prefix, ln[i], " -- ", nodename_prefix,ln[j], sep = "")
          s <- paste(s, otp)
        } else if(am[[j, i]] != 1) { ## edge is unidirectional
          otp <- paste(nodename_prefix, ln[i], " -> ", nodename_prefix,ln[j], sep = "")
          s <- paste(s, otp)
        }
      }
    }
  }
  s <- paste(s, "}")
  ## if you want to create dagitty object in another session, you could store s
  ## and create object later on with this string
  ## print(s)
  ## create dagitty object and return it
  dg <- dagitty(s)
  return(dg)
}

## gets a pcAlgo object and prints a string in the dagitty web format
conv_to_web <- function(g) {
  ## work with graph class
  # g <- pc @ graph
  ## list of nodes
  ln <- names(nodeData(g))
  ## list of edges in format "a|b"
  le <- names(edgeData(g))
  ## print node names followed by newline
  for(node in ln) {
    cat(paste(node, "1\n"))
  }
  cat("\n")
  ## for each node find outgoing edges
  for(node in ln) {
    ## otp is the string of this nodes edges
    otp <- node
    for(edge in le) {
      splt <- strsplit(edge, "\\|")
      ## if starting node of edge is this node, print end node
      if(splt[[1]][1] == node) {
        otp <- paste(otp, splt[[1]][2])
      }
    }
    cat(paste(otp, "\n"))
  }
}

# library(graph)
# library(dagitty)
## gets a dagitty object g and returns a graphnel object
conv_to_graphnel <- function(g) {
  ## list of nodes which will also be input for graphNEL function
  ln <- names(g)
  ## data structure which is a list of factors or something like that
  le <- dagitty::edges(g)
  ## this list of character vectors stores the edges of the graph
  eList <- list()
  ## go over all nodes
  for (i in 1:length(ln)) {
    node <- ln[i]
    ## edges of this node as character vector
    vc <- c()
    for (j in 1:length(le[[1]])) {
      ## le[[1]] stores the start nodes and le[[2]] the end nodes
      start <- as.character(le[[1]][j])
      end <- as.character(le[[2]][j])
      ## if the node we consider in this iteration is the start of an edge
      ## we append the end node to vc
      if(start == node) {
        ## this copies the vector every time and is therefore possibly slow
        vc <- c(vc, end)
      }
    }
    if(length(vc) > 0) {
      eList[[i]] <- vc
    } else {
      ## if an node has no neighbors put character(0)
      eList[[i]] <- character(0)
    }
    ## name character vector with node name
    names(eList)[i] <- node
  }
  ## create directed graphNEL object with nodes and eList
  ng <- graphNEL(nodes = ln, edgeL = eList, edgemode = "directed")
  ## you can plot the graph directly or return it
  ## plot(ng)
  return(ng)
}




# library(graph)

## source("~/dag/Code/graphNEL_to_tikz.R")

## g is a graphNEL object
conv_to_DIMACS <- function(g, type_of_graph = "undirected", outpath = "") {
  ## store number of nodes
  n <- length(nodes(g))
  el <- edgeL(g)
  ## m will be the number of edges which we have to count
  m <- 0
  am <- matrix(0, nrow = n, ncol = n)
  if(type_of_graph == "undirected") {
    ## filter edges with adjacency matrix
    for(i in 1:n) {
      for(j in el[[i]][[1]]) {
        if(am[[j, i]] == 0) {
          am[[i, j]] = 1
          m = m + 1
        }
      }
    }
  } else {
    for(i in 1:n) {
      for(j in el[[i]][[1]]) {
        am[[i, j]] = 1
        m =  m+1
      }
    }
  }
  ## our output string
  s <- "p tw"
  s <- paste(s, n, m)
  s <- paste(s, "\n", sep = "")
  for(i in 1:n) {
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        s <- paste(s, i, " ", j, sep="")
        s <- paste(s, "\n", sep = "")
      }
    }
  }
  if(outpath == "") {
    return(s)
  } else {
    sink(paste0(outpath, "-graph-DIMACS.txt"))
    cat(s)
    sink()
  }
}

## g is a graphNEL object and path gives the path to a textfile
## in tree decomposition format
td_to_graphNEL <- function(g, path) {
  ## n will be the number of nodes
  n = 0
  node_labels = c()
  g_labels = nodes(g)
  am = matrix(0, nrow = n, ncol = n)
  ## go through file line by line
  con = file(paste0(path, "-td.td"), "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    items = strsplit(line, " ")
    if(items[[1]][1] == "s") {
      n = strtoi(items[[1]][3])
      node_labels = vector(mode="character", length=n)
      am = matrix(0, nrow = n, ncol = n)
    } else if(items[[1]][1] == "b") {
      label = "~"
      for(i in 3:length(items[[1]])) {
        if(length(items[[1]]) < 3) break
        label = paste(label, g_labels[strtoi(items[[1]][i])], "~", sep="")
      }
      node_labels[strtoi(items[[1]][2])] = label
    } else {           
      am[[strtoi(items[[1]][1]), strtoi(items[[1]][2])]] = 1
      am[[strtoi(items[[1]][2]), strtoi(items[[1]][1])]] = 1
    }
  }
  eList <- list()
  for(i in 1:n) {
    vc <- c()
    for(j in 1:n) {
      if(am[[i, j]] == 1) {
        vc <- c(vc, node_labels[j])
      }
    }
    if(length(vc) > 0) {
      eList[[i]] <- vc
    } else {
      ## if an node has no neighbors put character(0)
      eList[[i]] <- character(0)
    }
    ## name character vector with node name
    names(eList)[i] <- node_labels[i]
  }
  ng <- graphNEL(nodes = node_labels, edgeL = eList, edgemode = "undirected")
  ## plot(ng)
  ## toFile(agopen(ng, "test"), filename="test.ps", fileType="ps")
  ## return(ng)
  close(con)
  ## graphNEL_to_tikz(ng, outpath, name)
  return(ng)
}

## g is a graphNEL object and path is the directory to which shall be written which includes a img subdirectory
graphNEL_to_tikz <- function(g, path, suffix = "") {
  dir = path
  for(i in nchar(dir):1) {
    if(substr(dir,i,i) == "/") {
      dir = substr(dir, 1, i)
      break
    }
  }
  sink(paste0(path, suffix, ".tex"))
  ## tex body of document
  cat("\\RequirePackage{luatex85,shellesc}\n")
  cat("\\documentclass[]{article}\n")
  cat("\\nofiles\n")
  cat("\n")
  cat("% LuaLaTeX stuff\n")
  cat("\\usepackage[utf8]{luainputenc}\n")
  cat("\\usepackage{luacode, luatextra, luatexbase}\n")
  cat("\n")
  cat("% TikZ\n")
  cat("\\usepackage{tikz}\n")
  cat("\\usetikzlibrary{graphs, graphdrawing}\n")
  cat("\\usegdlibrary{force, trees, circular}\n")
  cat("\n")
  cat("\\usetikzlibrary{external}\n")
  cat("\\tikzexternalize[prefix=./]\n")
  cat("\n")
  cat("\\tikzset{\n")
  cat("  graph/.style = {semithick},\n")
  cat("  vertex/.style = {draw, circle, inner sep = 0.33ex}\n")
  cat("}")
  cat("\n")
  cat("\\begin{document}\n")
  cat("\n")
  cat("\\tikz\\graph[binary tree layout] {\n")
  ## the binary tree
  for(i in 1:length(nodes(g))) {
    cat(i)
    cat("/ $\\{")
    cat(nodes(g)[i])
    cat("\\}$;\n")
  }
  for(i in 1:length(edgeL(g))) {
    for(j in 1:length(edgeL(g)[[i]][[1]])) {
      cat(i)
      cat(" -- ")
      cat(edgeL(g)[[i]][[1]][j])
      cat(";\n")
    }
  }
  cat("};")
  cat("\n")
  cat("\\end{document}\n")
  ## compiling the document
  system(paste0("cd ", dir, " && ", "lualatex --shell-escape ", path, suffix, ".tex"))
  system(paste0("rm ", path, suffix, ".auxlock"))
  system(paste0("rm ", path, suffix, "-figure0.dpth"))
  system(paste0("rm ", path, suffix, "-figure0.log"))
  system(paste0("rm ", path, suffix, "-figure0.md5"))
  system(paste0("rm ", path, suffix, ".log"))
  system(paste0("rm ", path, suffix, ".pdf"))
  system(paste0("mv ", path, suffix, "-figure0.pdf ", path, suffix, ".pdf"))
  sink()
}