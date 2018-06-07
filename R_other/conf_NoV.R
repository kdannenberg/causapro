## data parameters
## available data:
## "NoV_NMR-Tit_B4S"
## "NoV_NMR-Tit_Fuc"
## and
## "NoV_NMR-Tit_B3S-with-unass"
## type_of_data = "NMR-Tit"
## type_of_data = c("NMR_Tit-Fuc", "NMR_Tit-BTS")
## type_of_data = "Fuc-Tit-only-assigned"
## subtype_of_data = "Fuc-old"
## subtype_of_data = "Fuc"
## TODO: wieder ermöglichen
## subtype_of_data = c("Fuc", "BTS")

type_of_variables = "continuous",
protein = "NoV",
## type_of_data = "NMR-Tit",
type_of_data = "DDS",
subtype_of_data = "",
## subtype_of_data = c("Fuc", "BTS"),
data_set = "",
position_numbering = "",
## analysis parameters
min_pos_var = NULL,
show_variance_cutoff_plot = NULL,
only_cols = NULL,
only_cols_label = "",
alpha = 0.05,
pc_indepTest = NULL,
pc_suffStat = NULL,
ranked = FALSE,
rank_obs_per_pos = FALSE,
cor_cov_FUN = NULL,
pc_solve_conflicts = TRUE,
pc_u2pd = NULL,
pc_conservative = NULL,
pc_maj_rule = NULL,
weight_effects_on_by = NULL, # "var", "mean", ""
## analysis parameters: clustering of effects
effects_cluster_k = NULL,
effects_cluster_method = NULL,
effects_hclust_method = NULL,  #"average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
effects_dist_method = NULL,
effects_pv_nboot = NULL,
effects_cluster_alpha = NULL,
## graphical parameters
graph_output_formats = "pdf",
graph_layout = NULL,
graph_layout_igraph = NULL,
coloring = NULL,
colors = NULL,
plot_as_subgraphs = NULL,
plot_only_subgraphs = NULL, # 1 is another option
plot_ida = NULL,                                  # NEW!
plot_clusters = NULL,                              # NEW!
plot_no_isolated_nodes = NULL,  # TODO: make true possible even for edgeless -> empty graphs #NEW
for_combined_plot = NULL,
mute_all_plots = NULL,
other = NULL, # "cov"
## technical parameters
graph_computation = NULL,
evaluation = NULL,
causal_analysis = NULL,
max_conflict_edges = NULL,
intervention_position = "all",
causal_effects_function = NULL,
ida_direction = NULL,
linkcommunities = NULL,
linkcommunities_k = NULL,
linkcommunities_base_colors = NULL,
stages = NULL, # c("orig", "sub"), "sub"
print_analysis = FALSE,
plot_analysis = TRUE,
plot_types = c("localTests", "graph"),
plot_with_graphviz = FALSE,
pymol_show_int_pos = FALSE,
pymol_sort_connected_components_by_length = NULL, # NEW!
pymol_mix_connected_components = NULL,  # NEW!
print_connected_components = NULL,
compute_pc_anew = NULL,
compute_localTests_anew = NULL,
unabbrev_r_to_info = NULL,
print_r_to_console = NULL,
lines_in_abbr_of_r = NULL,
data_in_results = NULL,
output_parameters_in_results = NULL,
ida_percentile = "11",
file_separator = NULL
