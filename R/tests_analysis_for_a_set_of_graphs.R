test_sum_all_effects <- function() {
  scaled_results_1_of <- matrix(c(1,0,0.5,1,1,1,2,2,-0.5), byrow = TRUE, ncol = 3)
  scaled_results_1_on <- matrix(c(1,0,1,1,2,1), byrow = TRUE, ncol = 2)
  rownames(scaled_results_1_of) <- rownames(scaled_results_1_on) <- as.character(1:3)
  # fake_result_1 <- list(ida=list(`372` = list(of = list(scaled_effects = scaled_results_1_of), 
  #                                             on = list(scaled_effects = scaled_results_1_on))))
  fake_result_1 <- list(ida=list(`372` = list(of = list(effects = scaled_results_1_of), 
                                              on = list(effects = scaled_results_1_on))))
  
  scaled_results_2_of <- matrix(c(2,3,1), byrow = TRUE, ncol = 1)
  scaled_results_2_on <- matrix(c(0,0,1,1,2,2), byrow = TRUE, ncol = 2)
  rownames(scaled_results_2_of) <- rownames(scaled_results_2_on) <- as.character(1:3)
  # fake_result_2 <- list(ida=list(`372` = list(of = list(scaled_effects = scaled_results_2_of), 
  # on = list(scaled_effects = scaled_results_2_on))))
  fake_result_2 <- list(ida=list(`372` = list(of = list(effects = scaled_results_2_of), 
                                              on = list(effects = scaled_results_2_on))))
  fake_results <- list(fake_result_1, fake_result_2)
  
  return(compute_over_all_graphs(fake_results, function_over_all_graphs = "sum", weight_effects_on_by = "", 
                                 use_scaled_effects_for_sum = FALSE))
}

desired_sum_of <- c(2.5,4,1.75)
desired_sum_on <- c(0.5,2.0,3.5) 
names(desired_sum_of) <- names(desired_sum_on) <- 1:3
# print(test_sum_all_effects())

test_sum_all_effects()

if (!(test_sum_all_effects()$overAllGraphs_of == desired_sum_of && test_sum_all_effects()$overAllGraphs_on == desired_sum_on)) {
  stop("test_sum_all_effects failed.")
}