source("~/.configuration_code.R")

source_all_function_scripts()

# source("compute_DAG_G.R")
# source("compute_DAG_S.R")

source("configuration_data.R")

# best_alphas (example)
# $DDS
# $DDS$`1e-04`
# [1] "0.008" "0.009"
# 
# $DDS$`0.001`
# [1] "1e-20" "1e-10" "1e-05" "1e-04" "0.001" "0.004" "0.005"
# 
# $DDS$`0.01`
# [1] "1e-20" "1e-10" "1e-05" "1e-04" "0.001"
# 
# 
# $`DDG-10`
# $`DDG-10`$`1e-04`
# [1] "0.06" "0.07"
# 
# $`DDG-10`$`0.001`
# [1] "0.06" "0.07"
# 
# $`DDG-10`$`0.01`
# [1] "0.06" "0.07"
# 
# 
# $`DDG-5`
# $`DDG-5`$`1e-04`
# [1] "0.04"
# 
# $`DDG-5`$`0.001`
# [1] "0.04"
# 
# $`DDG-5`$`0.01`
# [1] "0.04"
# 
# 
# $`DDG-all`
# $`DDG-all`$`1e-04`
# [1] "0.003"
# 
# $`DDG-all`$`0.001`
# [1] "0.003"
# 
# $`DDG-all`$`0.01`
# [1] "0.002"
# 
# 
# $`DDDG-10`
# $`DDDG-10`$`1e-04`
# [1] "0.08"
# 
# $`DDDG-10`$`0.001`
# [1] "0.004" "0.07"  "0.08"  "0.09" 
# 
# $`DDDG-10`$`0.01`
# [1] "0.004" "0.04"  "0.05" 
# 
# 
# $`DDDG-5`
# $`DDDG-5`$`1e-04`
# [1] "0.009" "0.01" 
# 
# $`DDDG-5`$`0.001`
# [1] "0.009" "0.01" 
# 
# $`DDDG-5`$`0.01`
# [1] "0.08"
# 
# 
# $`DDDG-all`
# $`DDDG-all`$`1e-04`
# [1] "1e-05"
# 
# $`DDDG-all`$`0.001`
# [1] "1e-05"
# 
# $`DDDG-all`$`0.01`
# [1] "0.001" "0.007" "0.008"

# DDS.0.001.0.01 DDDG-10.0.08.1e-04

graphics.off()

certain_parameters <- list()
# certain_parameters$`DDS`$`0.01` <- 1e-20
certain_parameters$`DDS`$`0.01` <- c(1e-20, 0.001)

certain_parameters$`DDDG-10`$`1e-04` <- 0.08
certain_parameters$`DDDG-10`$`0.01` <- 0.07

effects_for_distinct_alphas(certain_parameters, with_graphs = TRUE, for_all_alphas = FALSE, cols_for_measures = FALSE)
