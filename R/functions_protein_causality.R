f_protein_causality_alpha <- function(measure, ...) {
  if (missing("measure")) {
    return(function(alpha) {protein_causality(alpha = alpha, ...)})
  } else { 
    return(function(alpha) {get(paste0("protein_causality_", measure))(alpha = alpha, ...)})
  }
}

f_protein_causality_pc_parameters_eval_analysis <- function(measure, ...) {
  if (missing("measure")) {
    return(function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
      protein_causality(pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                        evaluation = evaluation, analysis = analysis, ...)})
  } else { 
    return(function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
      get(paste0("protein_causality_", measure))(pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                        evaluation = evaluation, analysis = analysis, ...)})
  }
  
}

f_protein_causality_pc_parameters_eval_analysis_varcutoff_alpha <- function(measure, ...) {
  if (missing("measure")) {
    return(function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
      protein_causality(pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                        evaluation = evaluation, analysis = analysis, alpha, min_pos_var, ...)})
  } else { 
    return(function(pc_solve_conflicts, pc_u2pd, pc_maj_rule, pc_conservative, evaluation, analysis) {
      get(paste0("protein_causality_", measure))(pc_solve_conflicts = pc_solve_conflicts, pc_u2pd = pc_u2pd, pc_maj_rule = pc_maj_rule, pc_conservative = pc_conservative,
                                                 evaluation = evaluation, analysis = analysis, alpha, min_pos_var, ...)})
  }
  
}

# prot_caus <- f_protein_causality_alpha(measure = "G", min_pos_var = 0.00007)
# prot_caus(alpha = 0.123)

# TODO
# # TODO have also the remaining input parameter (alpha) variable and passed to the factory function
# f_protein_causality_alpha <- function(..., parameters = parameters(...)) {
#   # pass the ... as an input argument of the returned function (waht happens is taht it is interpreted as a new ..., 
#   # so whatever parameters as an input to the returned function)
#   return(function(list(...)) {do.call(protein_causality, parameters)})
# }
# 
# pc <- f_protein_causality_alpha("alpha", parameters = list(min_pos_var = 0.0123))
# pc(alpha = 0.05)
