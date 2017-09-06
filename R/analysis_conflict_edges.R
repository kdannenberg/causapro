source("~/.configuration_code.R")

# source("functions_causal_effects.R")
# source("functions_general.R")
# source("functions_conversions.R")
# source("functions_tools.R")

# source("functions_analysis_for_a_set_of_graphs.R")

source("functions_analysis_for_a_set_of_graphs.R")

source("compute_DAG_G.R")
source("~/.configuration_code.R")
source("compute_DAG_S.R")

plot_labels_as_rows_and_cols = TRUE
plot_logscale_alpha = FALSE
ylim_for_plots = c(0,200) # or NULL
opt_alpha_double_weight_conflict = FALSE
compute_anew = FALSE

measures <- c("DDS",
  "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")

min_pos_vars <- c(0.0001, 0.001, 0.01)

load(file = "edge_types.RData")

if (compute_anew || !exists("edge_types")) {
  edge_types <- list()

  # alphas <- c(1e-10, 1e-5, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05, 0.1)
  alphas <- c(1e-20, 1e-10, 1e-5, 0.0001, seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01), 0.1, 0.15, 0.2)
  
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
    for (min_pos_var in min_pos_vars) {
      edge_types[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]] <- list()
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
        quality_function <- function(v) {
          if (v["conflict"] > 15) {
            return(0)
          } else {
            return(v["directed"]/(v["undirected"] + 2 * v["conflict"]))
          }
        }
      } else {
        quality_function <- function(v) {
          if (v["conflict"] > 15) {
            return(0)
          } else {
            return(v["directed"]/(v["undirected"] + v["conflict"]))
          }
        }
      }
      quality <- apply(edge_types_by_alpha, 1, quality_function)
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
        
        matplot(edge_types_by_alpha, x = rownames(edge_types_by_alpha), xlab = "alpha", ylab = "# edges",
                type = "l", lty = c(1,1,1,1), col = c('red','green','orange', 'black'), ylim = ylim_for_plots, log = log)
        abline(h = 15, col = "red", lty = 2)
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


# best_alphas_1 <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
#                                              plot_logscale_alpha = FALSE, ylim_for_plots = c(0,200), # or NULL
#                                              opt_alpha_double_weight_conflict = FALSE,
#                                              print = TRUE, plot = TRUE)
# best_alphas_2 <- analyse_edge_types_by_alpha(edge_types = edge_types, plot_labels_as_rows_and_cols = TRUE,
#                                              plot_logscale_alpha = FALSE, ylim_for_plots = c(0,200), # or NULL
#                                              opt_alpha_double_weight_conflict = TRUE,
                                             # print = TRUE, plot = TRUE)

# pc_function = function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                                                    # mute_all_plots = TRUE))
# analyse_set_of_graphs_fct
effects_for_best_alphas <- list()
par(mfcol = c(length(min_pos_vars), length(measures)))

for (measure_type_sub in names(best_alphas_1)) {
  measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
  type_of_data = strsplit(measure_type_sub, "-")[[1]][1]
  if (grepl("-", measure_type_sub)) {
    subtype_of_data_list <- subtype_of_data <- strsplit(measure_type_sub, "-")[[1]][2]
  } else {
    subtype_of_data <- ""
    subtype_of_data_list <- "-"
  }
  
  for (min_pos_var in as.numeric(names(best_alphas_1[[measure_type_sub]]))) {
    alpha <- best_alphas_1[[measure_type_sub]][[as.character(min_pos_var)]]
    
    protein_causality_function = get(paste0("protein_causality_", measure))
    
    pc_function_ <- function_set_parameters(protein_causality_function, parameters = list(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                                                                           mute_all_plots = TRUE,
                                                                           alpha = alpha, min_pos_var = min_pos_var))
    
    effects_for_best_alphas[[type_of_data]][[subtype_of_data_list]][[as.character(min_pos_var)]][[as.character(alpha)]] <- 
      analyse_set_of_graphs(type_of_data = type_of_data, subtype_of_data = subtype_of_data, 
                            direction = "mean", measure = measure, pc_function = pc_function_, alpha = alpha, min_pos_var = min_pos_var,
                            for_combined_plot = TRUE, scale_in_the_end = FALSE, new = FALSE, save = TRUE)
      
    #   analyse_set_of_graphs(
    #   type_of_graph_set = "conflict", protein = "PDZ",
    #   measure = measure,
    #   type_of_data = type_of_data,
    #   subtype_of_data = subtype_of_data,
    #   protein_causality_function = get(paste0("protein_causality_", measure)),
    #   alpha = alpha,
    #   min_pos_var = min_pos_var, 
    #   pc_solve_conflicts = FALSE,
    #   # used when type_of_graph_set = "conflict"
    #   pc_maj_rule_conflict = TRUE,
    #   pc_conservative_conflict = FALSE,
    #   use_scaled_effects_for_each_graph = FALSE,
    #   scale_in_the_end = FALSE,
    #   # results <- protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
    #   #                           pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd,
    #   #                           graph_computation = FALSE, evaluation = FALSE, analysis = FALSE,
    #   #                           data_in_results = TRUE, output_parameters_in_results = TRUE)
    #   ida_percentile = 1 - (11 / dim(data)[2]), # top 11
    #   # weight_effects_on_by = "",          # in der Summe ganz schlecht
    #   # weight_effects_on_by = "var",
    #   # weight_effects_on_by = "mean",
    #   weight_effects_on_by = "median",  # sieht (in der Summe) am besten aus
    #   scale_effects_on_so_372_equals_1 = TRUE,
    #   # if dir == "on": only effects on
    #   # if dir == "of": only effects of
    #   # if dir == c("on", "of") or  dir == "both": both speparately
    #   # if dir == "mean" : mean of effects on and of (first effects on are scaeld such that 372 has value 1)
    #   # for best graphs, dir = "both" and dir = "mean" yield the same
    #   direction = "mean",
    #   # direction = "both",
    #   # which part of the analysis should be plotted?
    #   ## IDA für alle s Graphen berechen und summieren (besser wäre vllt: mitteln, also nochmal durch s teilen. Kannst du gerne machen, Marcel.)
    #   plot = "over_all_graphs",
    #   function_over_all_graphs = "mean",
    #   ## IDA für alle s Graphen brechenen und denjenigen bestimmen, der am besten mit den gewünschten Ergebnissen übereinstimmt
    #   # plot = "best graph"
    #   ## für alle Graphen mit Nummern in <plot> die Abweichung von der Summe (dem zukünftigen Mittelwert) über alle Graphen bestimmen, 
    #   ## also gewissermaßen wie repräsentativ der Graph jeweils für die Menge ist
    #   # plot = 1:25     ## if (is.numeric(plot) && length(plot) > 1) --> deviation from mean for graph(s) nr. <plot>
    #   # pc_function <- function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
    #   #   protein_causality_G(min_pos_var = min_pos_var, alpha = alpha,
    #   #                       pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
    #   #                       evaluation = evaluation, analysis = analysis)
    #   # }
    #   
    #   # pc_function = function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
    #   #   protein_causality_function(min_pos_var = min_pos_var, alpha = alpha,
    #   #                              pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
    #   #                              evaluation = evaluation, analysis = analysis, mute_all_plots = for_combined_plot)
    #   # },
    #   # ida_function = function(results) {
    #   #   causal_effects_ida(data = data, perturbated_position = "372", direction = "both", weight_effects_on_by = weight_effects_on_by,
    #   #                      protein = protein, results = results, coloring = "all", no_colors = FALSE, outpath = outpath,
    #   #                      amplification_exponent = 1, amplification_factor = TRUE, rank_effects = FALSE, effect_to_color_mode = "#FFFFFF",
    #   #                      pymol_bg_color = "grey", caption = caption, show_neg_causation = TRUE, neg_effects = "sep", analysis = TRUE, 
    #   #                      percentile = ida_percentile, mute_all_plots = for_combined_plot)
    #   # },
    #   s = 10,    # sample size
    #   for_combined_plot = TRUE,
    #   caption_as_subcaption = TRUE # for_combined_plot
    # )
  }
}

# all_effects <- list()
# for (measure_type_sub in c("DDS", "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")) {
#   measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
#   subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
#   if (is.na(subtype_of_data)) {
#     subtype_of_data <- ""
#   }
#   all_effects[[measure_type_sub]] <- analyse_graphs(measure = measure, type_of_data = strsplit(measure_type_sub, "-")[[1]][1], 
#                                                     subtype_of_data = subtype_of_data, 
#                                                     alphas = c(0.001, 0.005, 0.01, 0.05, 0.1), min_pos_vars = c(0.0001, 0.001, 0.01))
# }
