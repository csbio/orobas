# Loads packages
packages <- c("ggplot2", "ggthemes", "devtools", "stringr", "dplyr", "readr")
for (p in packages) {
  library(p, character.only = TRUE)
}

# Prints start time
start_time = Sys.time()
cat(paste("start time:", start_time, "\n"))

# Loads Orobas scoring package
devtools::load_all("orobas-main")

# Sets path to raw read count file
raw_read_count_file = file.path("data","Dataset_all_readcounts.txt") 

# Sets path to screen_replicate_map and condition_control_map tables
screen_replicate_map_file = file.path("data","screen_replicate_map_table.tsv")
condition_control_map_file = file.path("data","condition_control_map_table.tsv")

# Sets path to output folder
parent_folder = file.path("output")

orobas:::run_single_screen_scoring(
parent_folder = parent_folder,
raw_read_count_file = raw_read_count_file,
screen_replicate_map_file = screen_replicate_map_file,
condition_control_map_file = condition_control_map_file,
plot_type = 'png', 
display_numbers = FALSE, 
show_colnames = FALSE, 
show_rownames = FALSE,
filter_names_postfix = 'T0', 
cf1 = 1e6, 
cf2 = 1, 
min_reads = 30, 
max_reads = 10000, 
nonessential_norm = TRUE,
replace_NA = TRUE,
min_guides = 1, 
loess = TRUE, 
ma_transform = FALSE,
control_genes = c("LacZ", "luciferase", "EGFP"),
fdr_method = "BH",
fdr_threshold_positive  = 0.05, 
fdr_threshold_negative = 0.05,
differential_threshold_positive = 0.2, 
differential_threshold_negative = -0.2,
neg_type = "Negative",
pos_type = "Positive", 
label_fdr_threshold = .01, 
save_guide_dlfc = TRUE
)

# Prints start time
end_time = Sys.time()
cat(paste("end time:", end_time, "\n"))
total_time <- as.numeric(end_time - start_time, units = "mins") 
cat(paste("time:", total_time, "mins\n"))
