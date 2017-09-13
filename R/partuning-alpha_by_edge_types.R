source("~/.configuration_code.R")

# source("functions_causal_effects.R")
# source("functions_general.R")
# source("functions_conversions.R")
# source("functions_tools.R")

# source("analysis_for_a_set_of_graphs.R")
source("~/.configuration_code.R")
source("functions_analysis_for_a_set_of_graphs.R")

source("compute_DAG_G.R")
source("~/.configuration_code.R")
source("compute_DAG_S.R")

plot_labels_as_rows_and_cols = TRUE
plot_logscale_alpha = FALSE
# if scalar: height of the y-axis = max. possible number of edges (binom(n_nodes, 2)) / ylim_for_plots
# if vector: = ylim_for_plots
# if NULL: automatically set in each plot
ylim_for_plots = c(0,200) #5 # c(0, 200) # or NULL
opt_alpha_double_weight_conflict = FALSE
compute_anew = FALSE

measures <- c("DDS",
  "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")

min_pos_vars <- c(0.0001, 0.001, 0.01)

# alphas <- c(1e-10, 1e-5, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05, 0.1)
alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
# alphas <- c(0.01)

load(file = "edge_types.RData")

if (compute_anew || !exists("edge_types")) {
  edge_types <- list()
  # TODO: DDS hinzufügen
  # for (type_of_data in c("DDG", "DDDG")) {
  #   edge_types[[type_of_data]] <- list()
  #   for (subtype_of_data in c("5", "10", "all")) {
  #     edge_types[[type_of_data]][[subtype_of_data]] <- list()
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    edge_types[[type_of_data]][[subtype_of_data_list]] <- list()
    
    protein_causality_function = get(paste0("protein_causality_", measure))
    
    for (min_pos_var in min_pos_vars) {
      edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]] <- list()
      data <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                 min_pos_var = min_pos_var, graph_computation = FALSE, analysis = FALSE, 
                                                 evaluation = FALSE, data_in_results = TRUE)$data
      edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]]$max_nodes <- dim(data)[2]
      for (alpha in alphas) {
        protein_causality_function <- get(paste0("protein_causality_", measure))
        results <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                            min_pos_var = min_pos_var, alpha = alpha,
                            pc_solve_conflicts = TRUE, pc_u2pd = "relaxed",
                            pc_maj_rule = TRUE, pc_conservative = FALSE,
                            evaluation = FALSE, analysis = FALSE)
        edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]][[as.character(alpha)]] <- conflict_edges(results$pc@graph)
      }
    }
    }
  #   }
  # }
}


if (compute_anew || !exists("numbers_of_nodes")) {
  numbers_of_nodes <- list()
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    numbers_of_nodes[[type_of_data]][[subtype_of_data_list]] <- list()
    
    protein_causality_function = get(paste0("protein_causality_", measure))
    
    for (min_pos_var in min_pos_vars) {
      data <- protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                         min_pos_var = min_pos_var, graph_computation = FALSE, analysis = FALSE, 
                                         evaluation = FALSE, data_in_results = TRUE)$data
      numbers_of_nodes[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]] <- dim(data)[2]
    }
  }
}

# edge_types_by_alpha <- matrix(unlist(edge_types$DDG$`5`$`1e-04`), ncol = 3, byrow = TRUE)
# rownames(edge_types_by_alpha) <- names(edge_types$DDG$`5`$`1e-04`)
# colnames(edge_types_by_alpha) <- names(edge_types$DDG$`5`$`1e-04`[[1]])
# 
# matplot(edge_types_by_alpha, x = rownames(edge_types_by_alpha), type = "l", lty = 1, xlab = "alpha", ylab = "# edges")
# legend('top', legend = colnames(edge_types_by_alpha), col = c('black','green','red'), lwd = 1 )
# 
# # edge_types_by_alpha <- matrix(unlist(edge_types$DDG$`5`$`1e-04`), ncol = length(edge_types$DDG$`5`$`1e-04`))
# # colnames(edge_types_by_alpha) <- names(edge_types$DDG$`5`$`1e-04`)
# # rownames(edge_types_by_alpha) <- names(edge_types$DDG$`5`$`1e-04`[[1]])



analyse_edge_types_by_alpha <- function(edge_types, plot_labels_as_rows_and_cols = TRUE,
                                                      plot_logscale_alpha = FALSE,
                                                      ylim_for_plots = c(0,200), # or NULL
                                                      opt_alpha_double_weight_conflict = FALSE,
                                        print = TRUE, plot = TRUE) {
  best_alphas <- list()
  
  if (print) {
    graphics.off()
    if (plot_labels_as_rows_and_cols) {
      # par(mfrow = c(length(measures) + 1, length(min_pos_vars) + 1))
      par(mfcol = c(length(min_pos_vars) + 1, length(measures) + 1))
    } else {
      # par(mfrow = c(length(measures), length(min_pos_vars)))
      par(mfcol = c(length(min_pos_vars), length(measures)))
    }
    # for (measure_type_sub in c(# "DDS", 
    #   "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")) {
    
    # plot row labels
    if (plot_labels_as_rows_and_cols) {
      plot.new()
      legend('center', legend = c(names(edge_types[[1]][[1]][[1]][[1]]), "sum"), col = c('red','green','orange', 'black'), lty = c(1,1,1,1), lwd = 1 )
      for (min_pos_var in min_pos_vars) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("min_pos_var \n=", min_pos_var), 
             cex = 1.6, col = "black")
      }
    }
  }
  
  for (measure_type_sub in measures) {
    # measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data = "-"
    }
    best_alphas[[measure_type_sub]] <- list()
    if (plot) {
      # plot row labels
      if (plot_labels_as_rows_and_cols) {
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, measure_type_sub, 
             cex = 1.6, col = "black")
      }
    }
    
    for (min_pos_var in min_pos_vars) {
      print(paste(measure_type_sub, " - min_pos_var =", min_pos_var))
      edge_types_by_alpha <- matrix(unlist(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]]), ncol = 3, byrow = TRUE)
     
      rownames(edge_types_by_alpha) <- names(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]])
      colnames(edge_types_by_alpha) <- names(edge_types[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]][[1]])
      
      sum_of_edges <- apply(edge_types_by_alpha, 1, sum)
      edge_types_by_alpha <- cbind(edge_types_by_alpha, sum = sum_of_edges)
      if (opt_alpha_double_weight_conflict) {
        conflict_edge_weight <- 2
      } else {
        conflict_edge_weight <- 1
      }
      quality <- apply(edge_types_by_alpha, 1, quality_of_edge_distribution, weight_of_conflict_edges = conflict_edge_weight)
      best_alpha <- as.numeric(names(quality)[which(quality == max(quality))])
      best_alphas[[measure_type_sub]][[as.character(min_pos_var)]] <- best_alpha
      print(paste("Best alpha, according to quality:", best_alpha))
      
      if (print) {
        print(cbind(edge_types_by_alpha, quality = quality))
        cat("\n")
      }
    
      if (plot) {
        if (plot_logscale_alpha) {
          log <- "x"
        } else { 
          log <- ""
        }
        
        n_nodes <- numbers_of_nodes[[type_of_data]][[subtype_of_data]][[as.character(min_pos_var)]]
        n_edges <- 1/2 * n_nodes * (n_nodes - 1)
        if (length(ylim_for_plots) == 1) {
          current_ylim_for_plots <- c(0, n_edges / ylim_for_plots)
        } else {
          current_ylim_for_plots <- ylim_for_plots
        }
        
        matplot(edge_types_by_alpha, x = rownames(edge_types_by_alpha), xlab = "alpha", ylab = "# edges",
                type = "l", lty = c(1,1,1,1), col = c('red','green','orange', 'black'), ylim = current_ylim_for_plots, log = log)
        abline(h = 15, col = "red", lty = 2)
        # abline(h = n_edges, col = "black", lty = 1)
        abline(v = 0.01, col = "black", lty = 2)
        abline(v = best_alpha, col = "black", lty = 1)
        axis(1, at = best_alpha, labels = best_alpha, mgp = c(10, 2, 0))
        if (!plot_labels_as_rows_and_cols) {
          legend('top', legend = colnames(edge_types_by_alpha), col = c('red','green','orange', 'black'), lty = c(1,1,1,1), lwd = 1 )
          caption <- paste0(measure_type_sub, ", min_pos_var = ", min_pos_var)
          title(caption)
        }
      }
    }
  }
  return(best_alphas)
}

quality_of_edge_distribution <- function(edge_types, weight_of_conflict_edges = 1) {
  if (edge_types["conflict"] > 15) {
    return(0)                 # infeasible
  } else {
    return(edge_types["directed"]/(edge_types["undirected"] + weight_of_conflict_edges * edge_types["conflict"]))
  }
}

# best_alphas_1 <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
#                                              plot_logscale_alpha = FALSE, ylim_for_plots = ylim_for_plots, # or NULL
#                                              opt_alpha_double_weight_conflict = FALSE,
#                                              print = TRUE, plot = TRUE)
best_alphas_2 <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
                                             plot_logscale_alpha = FALSE, ylim_for_plots = ylim_for_plots, # or NULL
                                             opt_alpha_double_weight_conflict = TRUE,
print = TRUE, plot = TRUE)



# pc_function = function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                                                    # mute_all_plots = TRUE))
# analyse_set_of_graphs_fct

# distinct_alphas is a list containing for each measure a list (or vector) with a value alpha for each min_pos_var.
# For all measures in names(distinct_alphas) and all min_pos_vars in names(distinct_alphas$<measure>), the corresponding alpha is chosen, effects computed
# and displayed in grid
# if measures / min_pos_vars are given, only those are used (but need to be existent in distinct_alphas)
effects_for_distinct_alphas <- function(distinct_alphas, with_graphs = FALSE, for_all_alphas = FALSE,
                                        measures, min_pos_vars) {
  effects <- list()
  
  if (with_graphs) {
    if (missing(min_pos_vars)) {
      lines_of_plot <- length(names(distinct_alphas[[1]]))
    } else {
      lines_of_plot <- 2 * length(min_pos_vars)
    }
  } else {
    if (missing(min_pos_vars)) {
      lines_of_plot <- length(names(distinct_alphas[[1]]))
    } else {
      lines_of_plot <- length(min_pos_vars)
    }
  }
  if (for_all_alphas) {
    cols_of_plot = 5  # just sth
  } else {
    if (missing(measures)) {
      cols_of_plot <- length(names(distinct_alphas))
    }  else {
      cols_of_plot <- length(measures)
    }
  }
  
  par(mfcol = c(lines_of_plot, cols_of_plot))
  

  if (missing(measures)) {
    measures <- names(distinct_alphas)
  }
  for (measure_type_sub in measures) {
    measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
    type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
    if (grepl("-", measure_type_sub)) {
      subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
    } else {
      subtype_of_data <- ""
      subtype_of_data_list <- "-"
    }
    
    if (missing(min_pos_vars)) {
      min_pos_vars <- as.numeric(names(distinct_alphas[[measure_type_sub]]))
    }
    for (min_pos_var in min_pos_vars) {
      if (!for_all_alphas && length(alpha) > 1) {
        # alpha <- alpha[1]
        alpha <- alpha[length(alpha)] 
      }
      alphas <- distinct_alphas[[measure_type_sub]][[as.character(min_pos_var)]]
      for (alpha in alphas) {
        
        protein_causality_function = get(paste0("protein_causality_", measure))
    
        pc_function_ <- function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                                                               mute_all_plots = TRUE,
                                                                               alpha = alpha, min_pos_var = min_pos_var))
    
        effects[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]][[as.character(alpha)]] <-
          analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data,
                                direction = "mean", measure = measure, pc_function = pc_function_, alpha = alpha, min_pos_var = min_pos_var,
                                for_combined_plot = TRUE, scale_in_the_end = FALSE, new = FALSE, save = TRUE)
        if (with_graphs) {
          protein_causality_function(type_of_data = type_of_data, subtype_of_data = subtype_of_data, alpha = alpha, min_pos_var = min_pos_var, for_combined_plot = TRUE, pc_maj_rule = TRUE)
        } 
      }
    }
  }
}


# effects_best_alphas <- effects_for_distinct_alphas(best_alphas_1)