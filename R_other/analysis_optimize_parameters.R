setwd("~/Documents/Uni/Viren/R/Code")

source("compute_DAG_numerical.R")
source("general_functions.R")
source("evaluate_DAG.R")

optimize_pc_results <- function(pc_fun, analysis_after_pc_fun, output_dir, filename_w_o_alpha, alpha_start, compute_pc_anew, mute) {
  alpha = alpha_start
  alpha_results <- matrix (ncol = 4, nrow = 0)
  colnames(alpha_results) <- c("alpha", "mean", "# positives", "# negatives")
  
  print(paste("current alpha:", alpha))
  outpath <- paste(output_dir, filename_w_o_alpha, alpha, "/", filename_w_o_alpha, alpha, sep="")
  pc_fun_alpha <- function(outpath) {
    return(pc_fun(alpha = alpha, outpath = outpath))
  }
  pc <- get_pc(pc_fun_alpha, outpath, compute_pc_anew)
  results <- analysis_after_pc_fun(pc = pc, outpath = outpath) 
  alpha_results <- rbind(alpha_results, c(alpha, results$r_statistics$orig["r_int_signif_bad","mean"], results$r_statistics$orig["r_int_signif_bad","# positive elements"], results$r_statistics$orig["r_int_signif_bad","# negative elements"]))
  
  done <- FALSE
  countdown = 0
  factor_below_threshold = 10
  offset_above_threshold = 0.01
  threshold = 0.01
  repeat {
    if (done) {
      countdown <- countdown - 1
      if (countdown == 0) {
        break;
      }
    }
    if ((results$r_statistics$orig["r_int_signif_bad","mean"] > 0) 
        && ((results$r_statistics$orig["r_int_signif_bad","# negative elements"] == 0) || countdown > 0)) {
      # too many edges
      # decrease alpha
      if (alpha <= threshold) {
        alpha <- alpha * factor_below_threshold
      } else {
        alpha <- alpha + offset_above_threshold
      }
    } else if (((results$r_statistics$orig["r_int_signif_bad","mean"] < 0)
                && (results$r_statistics$orig["r_int_signif_bad","# positive elements"] == 0) || countdown > 0)) {
      # edges are missing
      # increase alpha
      if (alpha >= threshold) {
        alpha <- alpha - offset_above_threshold 
      } else {
        alpha <- alpha / factor_below_threshold
      }
    } else {
      done <- TRUE
      countdown = 5
      factor_below_threshold <- factor_below_threshold^2
      offset_above_threshold <- offset_above_threshold*5
      # TODO: überschwingen verhindern!! (wieder kleiner weden)
      # Oder erst groß, und nur nach mean, dann, wenn in den anderen Bereich geschweungen 
      # kleiner werden und so lange weiter machen, wie mean und #pos/#neg konsistent 
      # break;
    }
    
    print(paste("current alpha:", alpha))
    outpath <- paste(output_dir, filename_w_o_alpha, alpha, "/", filename_w_o_alpha, alpha, sep="")
    pc_fun_alpha <- function(outpath) {
      return(pc_fun(alpha = alpha, outpath = outpath))
    }
    pc <- get_pc(pc_fun_alpha, outpath, compute_pc_anew)
    results <- analysis_after_pc_fun(pc = pc, outpath = outpath) 
    alpha_results <- rbind(alpha_results, c(alpha, results$r_statistics$orig["r_int_signif_bad","mean"], results$r_statistics$orig["r_int_signif_bad","# positive elements"], results$r_statistics$orig["r_int_signif_bad","# negative elements"]))
  }
  
  print(alpha_results)
  print("Optimal alpha:")
  print(alpha_results[alpha_results[,"mean"] == min(abs(alpha_results[,"mean"])), "alpha"])
}

linear_optimization_of_pc_results <- function(pc_fun, analysis_after_pc_fun, parameters_for_info_file_fun, output_dir, filename_w_o_alpha, alpha_start, alpha_stop, alpha_step, alphas_vector, compute_pc_anew, mute, logscale = FALSE, plot_all_graphs = FALSE, plot_rows = 3) {
  if (plot_all_graphs) {
    graphics.off()
    par(mfrow = c(ceiling(length(alphas_vector)/plot_rows), plot_rows))
  }
  if (missing(alphas_vector)) {
    alphas_vector <- seq(from = alpha_start, to = alpha_stop, by = alpha_step)
  }
  # alpha = alpha_start
  alpha_results <- matrix (ncol = 4, nrow = 0)
  colnames(alpha_results) <- c("alpha", "mean", "# positives", "# negatives")
  
  fct_per_alpha <- function(alpha) {
    print(paste("current alpha:", alpha))
    outpath <- paste(output_dir, filename_w_o_alpha, alpha, "/", filename_w_o_alpha, alpha, sep="")
    directories <- strsplit(outpath, "/")
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = "/", sep = "/")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print("Directory created.")
    }
    
    pc_fun_alpha <- function(outpath) {
      return(pc_fun(alpha = alpha, outpath = outpath))
    }
    parameters_for_info <- parameters_for_info_file_fun(alpha = alpha, outpath = outpath)
    pc <- get_pc(pc_fun_alpha, outpath, compute_pc_anew, parameters_for_info)
    results <- analysis_after_pc_fun(pc = pc, outpath = outpath, caption = paste("alpha =", alpha)) 
    results_vector_signif <- c(alpha, results$r_statistics$orig["r_int_signif","# negative elements"], results$r_statistics$orig["r_int_signif","# positive elements"], results$r_statistics$orig["r_int_signif","mean"])
    results_vector_signif_bad <- c(alpha, results$r_statistics$orig["r_int_signif_bad","# negative elements"], results$r_statistics$orig["r_int_signif_bad","# positive elements"], results$r_statistics$orig["r_int_signif_bad","mean"])
    print(results_vector_signif)
    print(results_vector_signif_bad)
    
    return(c(results_vector_signif, results_vector_signif_bad[2:4]))
    # alpha_results <- rbind(alpha_results, c(alpha, results$r_statistics$orig["r_int_signif_bad","mean"], results$r_statistics$orig["r_int_signif_bad","# positive elements"], results$r_statistics$orig["r_int_signif_bad","# negative elements"]))
  }
  
  alpha_results <- do.call(rbind, lapply(as.list(alphas_vector), fct_per_alpha))
  colnames(alpha_results) <- c("alpha", "# neg (r_signif)", "# pos (r_signif)", "mean (r_signif)", "# neg (r_signif_bad)", "# pos (r_signif_bad)", "mean (r_signif_bad)")
  # do.call(fct_per_alpha, alphas_vector)
  alpha_results[is.nan(alpha_results)] <- 0
  print(alpha_results)
  if (!plot_all_graphs) {
    par(mfrow = c(1,2))
    if (logscale) {
      plot(x = alpha_results[,"alpha"], y = alpha_results[,"mean (r_signif)"], log = "x", main = "r_signif")
      plot(x = alpha_results[,"alpha"], y = alpha_results[,"mean (r_signif_bad)"], log = "x", main = "r_signif_bad")
    } else {
      plot(x = alpha_results[,"alpha"], y = alpha_results[,"mean (r_signif)"], main = "r_signif")
      plot(x = alpha_results[,"alpha"], y = alpha_results[,"mean (r_signif_bad)"], main = "r_signif_bad")
    }
  }
  print("Optimal alpha r_signif:")
  print(alpha_results[alpha_results[,"mean (r_signif)"] == min(abs(alpha_results[,"mean (r_signif)"])), "alpha"])
  print("Optimal alpha r_signif_bad:")
  print(alpha_results[alpha_results[,"mean (r_signif_bad)"] == min(abs(alpha_results[,"mean (r_signif_bad)"])), "alpha"])
  
  return(alpha_results)
}