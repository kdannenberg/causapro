measure_type_sub = "DDS" # "DDG-10", "DDG-5", "DDG-all", "DDDG-10", "DDDG-5", "DDDG-all")) {
measure = str_sub(strsplit(measure_type_sub, "-")[[1]][1], start = -1)
subtype_of_data = strsplit(measure_type_sub, "-")[[1]][2]
if (is.na(subtype_of_data)) {
  subtype_of_data <- ""
}

alphas <- c(0.000001, 0.00001, 0.005, 0.01, 0.05)
min_pos_vars <- c(0.0001, 0.001, 0.01)

eff <- analyse_graphs(measure = measure, type_of_data = strsplit(measure_type_sub, "-")[[1]][1],
                                                  subtype_of_data = subtype_of_data,
                                                  alphas = alphas, min_pos_vars = min_pos_vars)




# effects <- analyse_set_of_graphs(direction = "mean", measure = "G", alpha = 0.03, min_pos_var = 0.01, new = TRUE)