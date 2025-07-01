# Loads packages
packages <- c("stringr", "dplyr", "readr", "devtools", "optparse", "reticulate", "rappdirs", "sjmisc")
for (p in packages) {
  library(p, character.only = TRUE)
}

# Use the Python version installed in the virtual environment for reticulate
use_condaenv("orobas_env", required = TRUE)
py_config()

#################################################

start_time = Sys.time()
cat(paste("start time:", start_time, "\n"))

drug_list_4x4 <- list(
  "S01_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18","S08_Cisplatin_T18")
  , "S04_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18","S08_Cisplatin_T18")
  , "S08_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18","S08_Cisplatin_T18")
  ,"S05_CPT_T18"=list("S05_CPT_T18","S08_CPT_T18")
  ,"S08_CPT_T18"=list("S05_CPT_T18","S08_CPT_T18")
  ,"S12_HU_T18"=list("S12_HU_T18","S08_HU_T18")
  ,"S08_HU_T18"=list("S12_HU_T18","S08_HU_T18")
)

null_drug_list <- list(
  "S01_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18","S08_Cisplatin_T18")
  , "S04_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18",	"S08_Cisplatin_T18")
  , "S08_Cisplatin_T18"=list("S01_Cisplatin_T18","S04_Cisplatin_T18","S08_Cisplatin_T18")
  , "S05_CPT_T18"=list("S05_CPT_T18","S08_CPT_T18")
  , "S08_CPT_T18"=list("S05_CPT_T18","S08_CPT_T18")  
  , "S12_HU_T18"=list("S12_HU_T18","S08_HU_T18")
  , "S08_HU_T18"=list("S12_HU_T18","S08_HU_T18")
)


devtools::load_all('orobas')

input_folder = file.path("output","run_6_15_25")
output_folder = file.path("output","run_6_15_25")
condition_control_map_file = file.path('data','condition_control_map_table.tsv')
screen_replicate_map_file = file.path("data","screen_replicate_map_table.tsv")
raw_read_count_file <- file.path("data","Dataset_all_readcounts.txt")

orobas:::run_global_normalization(
input_folder = input_folder, 
condition_control_map_file = condition_control_map_file, 
output_folder = output_folder, 
screen_replicate_map_file = screen_replicate_map_file, 
raw_read_count_file = raw_read_count_file, 
filter_names_postfix = 'T0', 
cf1 = 1e6, 
cf2 = 1, 
min_reads = 30, 
max_reads = 10000, 
nonessential_norm = TRUE,
replace_NA = TRUE,
min_guides = 3,
loess = TRUE, 
ma_transform = FALSE,
control_genes = c("LacZ", "luciferase", "EGFP"),  
qc_control_pcc=0.50,
verbose = FALSE, 
black_list=c(),
screen_control_keyword="NT",
remove_screen_list=c(),
drug_list_4x4= drug_list_4x4,
null_drug_list= null_drug_list,
save_intermediate=FALSE,
neg_type = "Negative", 
pos_type = "Positive",
fdr_threshold_positive = 0.05, 
fdr_threshold_negative = 0.05,
differential_threshold_positive = 0.2, 
differential_threshold_negative = -0.2,
plot_type = "png", 
label_fdr_threshold = .001
)

end_time = Sys.time()
cat(paste("end time:", end_time, "\n"))
total_time <- as.numeric(end_time - start_time, units = "mins") 
cat(paste("time:", total_time, "mins\n"))
