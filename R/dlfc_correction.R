######
# Batch correction of differential log fold change scores from numerous screens
######

#' remove screens not to be included in the dLFC score data-frame
#' 
#' The scores matrix has library genes on the row side and compound screens on the column side.
#' A list of screens can be provided to be excluded and those columns will be removed from the scores matrix.
#' 
#' @param scores  differential log fold change scores data-frame (library-genes X screens; each screen column contain dLFC scores)
#' @param exclude_screen_list a list of screens to be removed from the scores data-frame; default is empty list c().  
#' @return  differential log fold change scores data-frame (library-genes X screens) after removing selected screens (columns)
#' 
#' @export
remove_screens<- function(scores, exclude_screen_list=c())
{
  #remove screens listed in exclude_screen_list from dLFC scores matrix 
	return(scores[,!(colnames(scores) %in% exclude_screen_list)])
}

#' Apply standard-deviation scaling to dLFC score data-frame
#' 
#' Calculates per screen standard deviation.
#' Calculate target standard-deviations per screen i.e. sd of scores between 10%-90% percentiles.
#' Normalize target standard-deviations per screen by average of all target standard-deviations.
#' Normalize each screen score with target standard-deviation.
#' 
#' @param scores  differential log fold change scores data-frame (library-genes X screens; each screen column contain dLFC scores)
#' @param sd_table_output_directory  directory path to save the sd_table.tsv (contains per screen standard-deviations pre- and post- sd-scaling)
#' @return  differential log fold change scores data-frame (library-genes X screens) after applying SD scaling per screen (column)
#' 
#' @export
apply_SD_scaling<- function(scores,sd_table_output_directory)
{
	  #calculate standard-deviation of each screen before scaling
	  pre_scaling_sd <- apply(scores, 2, stats::sd, na.rm=TRUE)
	  
	  #calculate target standard deviation per screen(SD of scores between 10%-90% percentiles)
	  lfc_range <- apply(scores, 2, stats::quantile, probs = c(0.1, 0.9), na.rm = TRUE)
	  target_sd <- rep(NA, ncol(scores))
	  for (i in 1:ncol(scores)) {
	    target_sd[i] <- stats::sd(scores[scores[,i] > lfc_range[1,i] & scores[,i] < lfc_range[2,i], i], na.rm = TRUE)
	  }
	  #normalize target standard deviation
	  sd_scale_factor <- mean(target_sd)
	  target_sd <- target_sd / sd_scale_factor
	  
	  #normalize per-screen dlfc scores by target standard deviation 
	  for (i in 1:ncol(scores)) {
	    scores[,i] <- scores[,i] / target_sd[i]
	  }
  
	#sanity check
	#calculate standard-deviation of each screen after scaling and write the scaling values to file
	post_scaling_sd <- apply(scores, 2, stats::sd, na.rm=TRUE)
	sd_table <- data.frame(cbind(pre_scaling_sd, post_scaling_sd))
	colnames(sd_table) <- c("pre_scaling_sd", "post_scaling_sd")
	sd_fname <- file.path(sd_table_output_directory,"sd_table.tsv")
	write.table(sd_table, sd_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


	return(scores)
}

#' remove principal component signal from dLFC scores
#' 
#' apply principal component analysis on the dLFC score matrix.
#' calculate projection of first mentioned number of principal component onto the score matrix.
#' subtract the projection from dLFC score matrix to remove principal component signal.
#' @param scores  differential log fold change scores data-frame (library-genes X screens; each screen column contain dLFC scores)
#' @param pc number of principal components to remove; default 1
#' @return  differential log fold change scores data-frame (library-genes X screens) after removing principal components
#' 
#' @export
remove_principal_component_signal<- function(scores,pc=1)
{	
	scores <- data.matrix(scores) #convert to matrix format
	scores <- scores[complete.cases(scores), ] #only keep rows with no 'NA' values

	pca <- prcomp(scores, center = TRUE, scale. = T) #apply principal component analysis on the score matrix
	rot <- pca$rotation[,pc] 
	projected <- scores %*% rot %*% t(rot) #calculate projection of first #'pc' principal components to the score matrix
	scores <- scores - projected #subtract the the projection matrix to remove signal from principal components
	
	return(data.frame(scores))

}

#' remove signal from control screens from dLFC scores
#' 
#' load differential log fold change scores generated from control screens from mentioned file.
#' apply principal component analysis to the control dLFC scores.
#' create projection of all principal components from control-dLFCs onto the main(condition) dLFC score matrix.
#' subtract the the projection matrix from condition score matrix to remove dLFC signal of control screens 
#' 
#' @param scores  differential log fold change scores data-frame from condition screens (library-genes X screens, each screen column contain dLFC scores)
#' @param control_dlfc_filepath file path to DMSO (control screens) differential log fold change score file (generated by OROBAS pipeline )
#' @return  differential log fold change scores data-frame (library-genes X screens) after removing principal components
#' 
#' @export
remove_control_dlfc_signal<- function(scores,control_dlfc_filepath)
{	
  	#load dLFC score file from control_dlfc_filepath
  	#file format- tab delimited .tsv file with screens as columns and genes as rows. 
	cg_dmso = read.table(control_dlfc_filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
	cg_dmso <- data.matrix(cg_dmso) #convert to matrix format
	cg_dmso <- cg_dmso[complete.cases(cg_dmso), ] #keep rows with no 'NA' values
	
	#get common genes between control and condition score files and get those subset of rows
	#control dlfc and condition dlfc gene sets may not overlap completely
	genes <- intersect(rownames(scores), rownames(cg_dmso))
	scores_ <- scores[genes,]
	cg_dmso <- cg_dmso[genes,]
	
	#apply principal component analysis to control dLFC scores
	pca_cg_dmso_t <- prcomp(as.matrix(t(cg_dmso)), center = TRUE, scale. = T)

	#create projection of all principal components from control-dLFCs onto the main(condition) dLFC score matrix
	cg_dmso_loadings_t <- pca_cg_dmso_t$rotation
	dmso_v_t <- cg_dmso_loadings_t[,1:ncol(cg_dmso_loadings_t)]
	projected_t <- t(scores_) %*% dmso_v_t %*% t(dmso_v_t) 
	
	#subtract the projection matrix from condition score matrix subset to remove dLFC signal of control screens 
	scores_ <- scores_ - data.frame(t(projected_t))
	
	#insert the corrected gene profiles into condition score matrix
	#update the corrected gene profiles, and preserve the other gene profiles (with no control dlfc) instead of removing them	
	scores[match(rownames(scores_), rownames(scores)), ] <- scores_

	return(data.frame(scores))

}

#' wrapper function to call python script screen_batch_correction_LDA.py
#' 
#'  @param scores differential log fold change scores data-frame from condition screens (library-genes X screens, each screen column contain dLFC scores)
#'  @param output_directory path to directory to save output files (roc plot files)
#'  @return differential log fold change scores data-frame (library-genes X screens) after batch correction by LDA
#'  
#'  @export
screen_batch_correction_with_lda<- function(scores, output_directory)
{	
  	#load screen level batch correction LDA python script 
	source_python(system.file("python","screen_batch_correction_LDA.py", package = "orobas"))
	#call function from python script that returns batch corrected dLFC score file
	scores = run_batch_correction(scores, output_directory)	
	return(scores)
}

#'generate wbc (within-between correlation) score of 4x4 screens
#'
#'Score representing how much the replicate screens are correlated compared to non-replicate screens
#'Equation: wbc_per_screen = (average(within_correlation) - average(between_correlation))/sig_b
#'          within_correlation = pearson correlation of screen with other chosen replicates
#'          between_correlation = pearson_correlation of screen with all non-replicate screens
#'          sig_b = sqrt(sum( square(between_correlation - average(between_correlation)) ))/(number of non-replicate screens - 1))
#' <cite paper>         
#' @param data  differential log fold change scores data-frame (library-genes X screens, each screen column contain dLFC scores)
#' @param drug_list_4x4 A named nested list of 4x4 screens; used in the 'within correlation' calculation in wbc score
#' example of one entry: "CHEM015_BORTEZOMIB_T14"=list("CHEM015_BORTEZOMIB_T14","CHEM050_BORTEZOMIB_T14","CHEM051_BORTEZOMIB_T14","CHEM052_BORTEZOMIB_T13") 
#' @param null_drug_list  A named nested list of all screens that are replicates of the 4x4 screens; these should have no effect in the 'between correlation' calculation in wbc score
#' Example of one entry: "CHEM029_NGI1_T16"=list("CHEM031_NGI1_T13","CHEM050_NGI1_T11","CHEM029_NGI1_T16","CHEM050_NGI1_T14","CHEM051_NGI1_T14","CHEM057_NGI1_T13")
#' @return  data-frame with one wbc score for each 4x4 screens
#' 
#' @export
score_wbc <- function(data, drug_list_4x4, null_drug_list) {
 
  #generate pearson correlation matrix across all screens
  sim_net <- cor(data, use = "complete.obs",method="pearson")
  
  all_drugs = colnames(sim_net) #get all screen names
  drug_names = names(drug_list_4x4) #get the 4x4 screen names
  
  #create dataframe with 4x4 screen names and wbc score columns initialized to 0
  drug_wbc = data.frame(drug_names) 
  drug_wbc$wbc = 0
  
  #for each  4x4 screen
  for (drug in drug_names){
    #calculate 'within correlation' - with replicate screens
    w = sim_net[drug,all_drugs %in% drug_list_4x4[[drug]]]
    #calculate 'between correlation' - with non-replicate screens
    b = sim_net[drug,!all_drugs %in% null_drug_list[[drug]]]
    #calculate average of within correlation
    mean_w = mean(w)
    #calculate average of between correlation
    mean_b = mean(b)
    #calculate scaling factor: scaled root-sum-square of mean-centered between correlation
    sig_b = sqrt(sum((b-mean_b)**2)/(length(b)-1))
    #calculate wbc score and add to dataframe
    wbc = (mean_w - mean_b)/sig_b
    drug_wbc$wbc[drug_wbc$drug_names==drug]=wbc
  }
  return(drug_wbc)
}

#' wrapper function to apply several steps of correction to differential LFC scores from numerous screens combined
#' 
#'  @param scores differential log fold change scores data-frame from condition screens (library-genes X screens, each screen column contain dLFC scores).
#'  @param output_folder path to directory to save all output files.
#'  @param control_dlfc_filepath path to control dLFC score file "control_effect_scores.tsv" (in folder 'control') generated by function generate_control_dlfc_scores(). 
#'  The path is hard-coded to reflect the generated file.
#'  @param remove_screen_list list of screens to remove from batch correction.
#' @param drug_list_4x4 A named nested list of 4x4 screens; used in the 'within correlation' calculation in wbc score.
#' example of one entry: "CHEM015_BORTEZOMIB_T14"=list("CHEM015_BORTEZOMIB_T14","CHEM050_BORTEZOMIB_T14","CHEM051_BORTEZOMIB_T14","CHEM052_BORTEZOMIB_T13").
#' @param null_drug_list  A named nested list of all screens that are replicates of the 4x4 screens; these should have no effect in the 'between correlation' calculation in wbc score.
#' Example of one entry: "CHEM029_NGI1_T16"=list("CHEM031_NGI1_T13","CHEM050_NGI1_T11","CHEM029_NGI1_T16","CHEM050_NGI1_T14","CHEM051_NGI1_T14","CHEM057_NGI1_T13").
#'  @param save_intermediate If TRUE, save dLFC score file after each intermediate correction step (default FALSE).
#'  @return differential log fold change scores data-frame (library-genes X screens) after all batch correction steps.
#'  
#'  @export
apply_dlfc_correction<- function(
scores,
output_folder,
control_dlfc_filepath,
remove_screen_list=c(),
drug_list_4x4=c(),
null_drug_list=c(),
save_intermediate=FALSE
)
{
  flag_wbc = TRUE
  if(is_empty(drug_list_4x4))
  {
    flag_wbc = FALSE
  }
  
  if(flag_wbc)
  {
    wbc_4x4 = score_wbc(scores, drug_list_4x4, null_drug_list)
    colnames(wbc_4x4) = c('screens','no_correction')
  }
  
  #### remove listed screens from differencial LFC scores
  if(!is_empty(remove_screen_list))
  {
    scores <- remove_screens(scores, remove_screen_list)
    if(flag_wbc)
    {
      wbc_4x4_temp = score_wbc(scores, drug_list_4x4, null_drug_list)
      wbc_4x4$bad_screens_removed = wbc_4x4_temp$wbc
    }
  }
  
  ####Apply SD scaling to differenctial LFC scores
  scores <- apply_SD_scaling(scores, output_folder)
  if(save_intermediate)
  {
    score_fname <- file.path(output_folder, "dLFC_scores_sd_scaled.tsv")
    write.table(scores, score_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  if(flag_wbc)
  {
    wbc_4x4_temp = score_wbc(scores, drug_list_4x4, null_drug_list)
    wbc_4x4$sd_scaled = wbc_4x4_temp$wbc
  }
  
  ####remove first principal component from the differential LFC scores
  scores <- remove_principal_component_signal(scores,pc=1)
  if(save_intermediate)
  {
    score_fname <- file.path(output_folder, "dLFC_scores_sd_scaled_pc_removed.tsv")
    write.table(scores, score_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  if(flag_wbc)
  {
    wbc_4x4_temp = score_wbc(scores, drug_list_4x4, null_drug_list)
    wbc_4x4$pc_removed = wbc_4x4_temp$wbc
  }
  
  ####remove DMSO signal from the differential LFC scores
  if(file.exists(control_dlfc_filepath))
  {	
    scores <- remove_control_dlfc_signal(scores,control_dlfc_filepath)
    if(save_intermediate)
    {
      score_fname <- file.path(output_folder, "dLFC_scores_sd_scaled_pc_removed_control_removed.tsv")
      write.table(scores, score_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    }
    if(flag_wbc)
    {
      wbc_4x4_temp = score_wbc(scores, drug_list_4x4, null_drug_list)
      wbc_4x4$control_dlfc_removed = wbc_4x4_temp$wbc
    }
  }
  
  ####batch correction using lda
  scores <- screen_batch_correction_with_lda(scores, output_folder)
  score_fname <- file.path(output_folder, "dLFC_scores_sd_scaled_pc_removed_control_removed_batch_corrected.tsv")
  write.table(scores, score_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  if(flag_wbc)
  {
    wbc_4x4_temp = score_wbc(scores, drug_list_4x4, null_drug_list)
    wbc_4x4$lda_batch_corrected = wbc_4x4_temp$wbc
    write.csv(wbc_4x4, file.path(output_folder,'wbc_4x4_all.csv'), quote = FALSE,row.names = FALSE)
  }
  
  return(scores)
}

#' wrapper function to generate differential log fold change from control screens
#' 
#'  @param batch_table_file_path 
#'  @param sample_table_file_path 
#'  @param raw_read_count_data_file_path 
#'  @param output_folder 
#' @param filter_names_postfix Postfix to identify list of screen names to filter based on read counts by. 
#' @param cf1 parameter for \code{normalize_screens}. Scaling factor (default 1e6).
#' @param cf2 parameter for \code{normalize_screens}. Pseudocount (default 1).
#' @param min_reads parameter for \code{normalize_screens}. Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads parameter for \code{normalize_screens}. Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @param nonessential_norm parameter for \code{normalize_screens}. Whether or not to normalize each screen against its
#'   population of core non-essential genes, as defined by Traver et al. 2015 
#'   (default FALSE).
#' @param replace_NA parameter for \code{normalize_screens}.  Whether or not to replace NA and NULL values in non-T0 screens 
#'   with 0's after filtering out T0 guides with too few or NA readcounts 
#'   (default TRUE).
#' @param min_guides parameter for \code{score_controls_batch}.
#' @param loess parameter for \code{score_controls_batch}.
#' @param ma_transform parameter for \code{score_controls_batch}.
#' @param control_genes parameter for \code{score_controls_batch}.
#' @param qc_control_pcc parameter for \code{score_controls_batch}.
#' @param verbose parameter for \code{score_controls_batch}.
#' @param black_list parameter for \code{score_controls_batch}.
#' @param screen_control_keyword parameter for \code{score_controls_batch}.
#' Output: creates file dLFC score file "control_effect_scores.tsv" in folder 'control' 
#'  
#'  @export
generate_control_dlfc_scores <- function(
batch_table_file_path,
sample_table_file_path, 
raw_read_count_data_file_path,
output_folder, 
filter_names_postfix = 'T0', 
cf1 = 1e6, 
cf2 = 1, 
min_reads = 30, 
max_reads = 10000, 
nonessential_norm = TRUE,
replace_NA = TRUE,
min_guides = 3,
loess = TRUE, 
ma_transform = TRUE,
control_genes = c("None", ""),  
qc_control_pcc=0.30,
verbose = FALSE, 
black_list=c('bad_screen_name', ""),
screen_control_keyword="DMSO|MOCK|Mock|WT|NGLY1|BMI1|Control"
)
{	
			   
	batch_table = read.table(batch_table_file_path, sep='\t',  header = T, stringsAsFactors = F)
	sample_table = read.table(sample_table_file_path,sep='\t', header = T, stringsAsFactors = F)

	control_sample_table = sample_table[which(sample_table$Screen %in% unique(batch_table$Control)),]
	sampla_table_NormalizeTo = unique(sample_table$NormalizeTo)
	cur_sample_t0 = sample_table[which(sample_table$Screen %in% sampla_table_NormalizeTo),]
	control_sample_table = unique(rbind(cur_sample_t0,control_sample_table))

	control_table = data.frame('Screen' =unique(batch_table$Control), 'Control' =unique(batch_table$Control), 'Group' = 'control', 'Type' = 'control')
	
	write.table(control_table, file = file.path(output_folder,"control_table.tsv")
            , sep='\t', row.names = F, quote = F)
	write.table(control_sample_table, file = file.path(output_folder,"control_sample_table.tsv")
            , sep='\t', row.names = F, quote = F)		
	
	#read raw-read count file
	raw_reads <- read.csv(raw_read_count_data_file_path, header = TRUE, stringsAsFactors = FALSE, 
               sep = "\t", check.names = FALSE, encoding = "UTF-8")
			   
	# Fix screen names
	colnames(raw_reads) <- format_replicate_names(colnames(raw_reads))
	
	# Process screens
	screens <- add_screens_from_table(file.path(output_folder,"control_sample_table.tsv"))
	
	#Extract gene and relevant screen replicate columns from the raw read-count data
	col_list = c('gene')
	for(item in c(1:length(screens)))
	{
	  reps = screens[[item]][["replicates"]]
	  col_list = c(col_list,reps)
	}
	cols= (colnames(raw_reads))
	raw_reads <- raw_reads[,(cols %in% col_list)]
	
	# normalize raw read-counts to earlier time-point 
	t0_screens <- names(screens)[grepl(paste('_',filter_names_postfix,'$',sep=''), names(screens))] #filter T0 screens while normalizing
	raw_reads <- normalize_screens(raw_reads, 
	screens, 
	filter_names = t0_screens,
	cf1 = cf1, 
	cf2 = cf2, 
    min_reads = min_reads, 
	max_reads = max_reads, 
	nonessential_norm = nonessential_norm,
    replace_NA = replace_NA
	)
	
	# Score dataset
	score_controls_batch(df = raw_reads, 
				screens = screens, 
				batch_file = file.path(output_folder,"control_table.tsv"), 
				output_folder = output_folder,
                control_genes = control_genes,
                min_guides = min_guides, 
				loess = loess, 
				ma_transform = ma_transform, 
                verbose = verbose, 
                qc_control_pcc = qc_control_pcc, 
                black_list = black_list,
                screen_control_keyword=screen_control_keyword)
	
}

#' wrapper function to generate differential log fold change from control screens
#'
#'  @param input_path 
#'  @param batch_table_file_path 
#'  @param sample_table_file_path 
#'  @param raw_read_count_data_file_path 
#'  @param output_folder 
#' @param filter_names_postfix parameter for \code{generate_control_dlfc_scores}.
#' @param cf1 parameter for \code{generate_control_dlfc_scores}.  Scaling factor (default 1e6).
#' @param cf2 parameter for \code{generate_control_dlfc_scores}. Pseudocount (default 1).
#' @param min_reads parameter for \code{generate_control_dlfc_scores}. Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads parameter for \code{generate_control_dlfc_scores}. Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @param nonessential_norm parameter for \code{generate_control_dlfc_scores}. Whether or not to normalize each screen against its
#'   population of core non-essential genes, as defined by Traver et al. 2015 
#'   (default FALSE).
#' @param replace_NA parameter for \code{generate_control_dlfc_scores}. Whether or not to replace NA and NULL values in non-T0 screens 
#'   with 0's after filtering out T0 guides with too few or NA readcounts 
#'   (default TRUE).
#' @param min_guides parameter for \code{generate_control_dlfc_scores}.
#' @param loess parameter for \code{generate_control_dlfc_scores}.
#' @param ma_transform parameter for \code{generate_control_dlfc_scores}.
#' @param control_genes parameter for \code{generate_control_dlfc_scores}.
#' @param qc_control_pcc parameter for \code{generate_control_dlfc_scores}.
#' @param verbose parameter for \code{generate_control_dlfc_scores}.
#' @param black_list parameter for \code{generate_control_dlfc_scores}.
#' @param screen_control_keyword parameter for \code{generate_control_dlfc_scores}.
#' @param remove_screen_list parameter for \code{apply_dlfc_correction}. List of screens to remove from batch correction.
#' @param drug_list_4x4 parameter for \code{apply_dlfc_correction}. A named nested list of 4x4 screens; used in the 'within correlation' calculation in wbc score.
#' example of one entry: "CHEM015_BORTEZOMIB_T14"=list("CHEM015_BORTEZOMIB_T14","CHEM050_BORTEZOMIB_T14","CHEM051_BORTEZOMIB_T14","CHEM052_BORTEZOMIB_T13"). 
#' @param null_drug_list parameter for \code{apply_dlfc_correction}. A named nested list of all screens that are replicates of the 4x4 screens; these should have no effect in the 'between correlation' calculation in wbc score.
#' Example of one entry: "CHEM029_NGI1_T16"=list("CHEM031_NGI1_T13","CHEM050_NGI1_T11","CHEM029_NGI1_T16","CHEM050_NGI1_T14","CHEM051_NGI1_T14","CHEM057_NGI1_T13").
#' @param save_intermediate parameter for \code{apply_dlfc_correction}. If TRUE, save dLFC score file after each intermediate correction step (default FALSE).
#' @param neg_type parameter for \code{plot_drug_response} and \code{call_drug_hits}. Label for significant effects with a negative differential effect
#'   passed to \code{call_drug_hits} (default "Negative").
#' @param pos_type parameter for\code{plot_drug_response} and \code{call_drug_hits}. Label for significant effects with a positive differential effect
#'   passed to \code{call_drug_hits} (default "Positive").
#' @param fdr_threshold_positive parameter for \code{call_drug_hits}. The threshold below which to call gene effects as significant positive hits
#'   (default 0.1).
#' @param fdr_threshold_negative parameter for \code{call_drug_hits}. The threshold below which to call gene effects as significant negative hits
#'   (default 0.1).
#' @param differential_threshold_positive parameter for \code{call_drug_hits}. Threshold on differential effects, 
#'   over which to call gene effects as significant positive hits (default 0).
#' @param differential_threshold_negative parameter for \code{call_drug_hits}. Threshold on differential effects, 
#'   below which gene effects are called as significant negative hits (default 0).
#' @param plot_type parameter for \code{plot_drug_response}. Type of plot to output, one of "png" or "pdf" (default "png").
#' @param label_fdr_threshold parameter for \code{plot_drug_response}. The threshold below which to plot gene labels for significant
#'   hits, or NULL to plot without labels (default NULL).
#'
#' Output: generates corrected dLFC score file and other files 
#'  
#'  @export
correct_dlfc_scores_in_batch <- function(
    input_path, 
	batch_table_file_path, 
	output_folder, 
	sample_table_file_path, 
	raw_read_count_data_file_path,
	filter_names_postfix = 'T0', 
	cf1 = 1e6, 
	cf2 = 1, 
	min_reads = 30, 
	max_reads = 10000, 
	nonessential_norm = TRUE,
	replace_NA = TRUE,
	min_guides = 3,
	loess = TRUE, 
	ma_transform = TRUE,
	control_genes = c("None", ""),  
	qc_control_pcc=0.30,
	verbose = FALSE, 
	black_list=c('bad_screen_name', ""),
	screen_control_keyword="DMSO|MOCK|Mock|WT|NGLY1|BMI1|Control",
    remove_screen_list=c(),
	drug_list_4x4=c(),
	null_drug_list=c(),
	save_intermediate=FALSE,
    neg_type = "Negative", 
	pos_type = "Positive",
    fdr_threshold_positive  = 0.1, 
	fdr_threshold_negative = 0.1,
    differential_threshold_positive = 0, 
	differential_threshold_negative = 0,
    plot_type = "png", 
	label_fdr_threshold = NULL    
    )
{
  # check if batch table exists
  if (!file.exists(batch_table_file_path))
  {	
    stop(paste("ERROR: Could not find file ", batch_table_file_path))
  }
  
  # check if sample table exists
  if (!file.exists(sample_table_file_path))
  {	
    stop(paste("ERROR: Could not find file ", sample_table_file_path))
  }
  
  # check if raw read-count file exists
  if (!file.exists(raw_read_count_data_file_path))
  {	
    stop(paste("ERROR: Could not find file ", raw_read_count_data_file_path))
  }
  
  # read batch file
  batch <- utils::read.csv(batch_table_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  
  #gather scores various screens from different files
  flag = 0
  for (i in 1:nrow(batch)) {
    condition <- batch[i,1] # get current condition screen name
    control <- batch[i,2] # get associated control screen name
    screen_name = strsplit(condition,"_")[[1]][1]
    input_file <- file.path(input_path, screen_name,"condition_gene_calls.tsv")
    
    if (file.exists(input_file) ) {
      curr_screen_score <- read.csv(file.path(input_file),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      col_list <- c('gene',
                    paste0('mean_',control,sep=''),
                    paste0('mean_',condition,sep=''),
                    paste0('differential_',condition,'_vs_',control,sep=''),
                    paste0('fdr_',condition,'_vs_',control,sep=''),
                    paste0('significant_',condition,'_vs_',control,sep=''),
                    paste0('effect_type_',condition,sep=''))
      
      
      if (paste0('mean_',condition,sep='') %in% colnames(curr_screen_score) ) { 
        
        curr_screen_score = curr_screen_score[,col_list]
        if(flag==0){
          all_score <- curr_screen_score
          flag = 1
        }else{
          genes <- intersect(all_score$gene, curr_screen_score$gene)
          all_score <- all_score[all_score$gene %in% genes,]
          curr_screen_score <- curr_screen_score[curr_screen_score$gene %in% genes,]
          all_score <- cbind(all_score, curr_screen_score[,2:ncol(curr_screen_score)])
        }
      }else{
        print(paste('Missing screen data: ',condition))
      }
    }else{
      print(paste('Missing screen file: ',screen_name))
    }
  }
  
  if(length(colnames(all_score))<36){ 
    stop(paste("Not enough data retrieved from screen files. Should be more than 3 screens!"))
  }
  
  #get differential LFC scores
  scores <- all_score[,grepl("gene|differential", colnames(all_score))]
  scores <- abbreviate_names(scores, "differential_", 2:ncol(scores))
  rownames(scores) <- scores$gene #set matrix rownames to gene names ('gene' column (first column) contains gene names at this point)
  scores <- scores[,-1] #remove first column ('gene' column)
  
  #generate control dlfc_score
  generate_control_dlfc_scores(
	batch_table_file_path = batch_table_file_path,
	sample_table_file_path = sample_table_file_path, 
	raw_read_count_data_file_path = raw_read_count_data_file_path,
	output_folder = output_folder,	
	filter_names_postfix = filter_names_postfix, 
	cf1 = cf1, 
	cf2 = cf2, 
	min_reads = min_reads, 
	max_reads = max_reads, 
	nonessential_norm = nonessential_norm,
	replace_NA = replace_NA,
	min_guides = min_guides,
	loess = loess, 
	ma_transform = ma_transform,
	control_genes = control_genes,  
	qc_control_pcc = qc_control_pcc,
	verbose = verbose, 
	black_list = black_list,
	screen_control_keyword = screen_control_keyword
	)
  
  control_dlfc_filepath = file.path(output_folder,'control',"control_effect_scores.tsv")
  
  # apply batch correction to dLFC scores 
  dlfc_score <- apply_dlfc_correction(
    scores = scores,
    output_folder = output_folder,
    control_dlfc_filepath = control_dlfc_filepath,
    remove_screen_list = remove_screen_list,
    drug_list_4x4 = drug_list_4x4,
    null_drug_list = null_drug_list,
    save_intermediate=save_intermediate
    )
  
  dlfc_score$gene <- rownames(dlfc_score)
  genes <- intersect(all_score$gene, dlfc_score$gene)
  all_score <- all_score[all_score$gene %in% genes,]
  dlfc_score <- dlfc_score[dlfc_score$gene %in% genes,]
  
  # Update the dLFC, significance, and effect type columns from all screens based on the corrected dLFC scores
  # Generate LFC scatter plots based on the corrected dLFC scores 
  plots_folder = file.path(output_folder, "plots")
  if (!dir.exists(plots_folder)) { dir.create(plots_folder) }
  
  for (i in 1:nrow(batch)) {
    condition <- batch[i,1] # get current condition screen name
    control <- batch[i,2] # get associated control screen name
    
    if(( paste0("differential_", condition, "_vs_", control) %in% colnames(all_score) ) & (condition %in% colnames(dlfc_score)) )
    {
      col_list <- c('gene',
                    paste0('mean_',control,sep=''),
                    paste0('mean_',condition,sep=''),
                    paste0('differential_',condition,'_vs_',control,sep=''),
                    paste0('fdr_',condition,'_vs_',control,sep=''),
                    paste0('significant_',condition,'_vs_',control,sep=''),
                    paste0('effect_type_',condition,sep=''))
      curr_screen_score = all_score[,col_list]
      curr_screen_score[[paste0("differential_", condition, "_vs_", control)]] <- dlfc_score[[condition]]
      
      # call call_drug_hits() to add information to significant and effect-type columns indicating significant positive and negative interactions that meet the provided cut-offs
      curr_screen_score <- call_drug_hits(scores = curr_screen_score,
                                          control_screen_name = control, 
                                          condition_screen_names = condition,
                                          neg_type = neg_type, 
                                          pos_type = pos_type,
                                          fdr_threshold_positive  = fdr_threshold_positive, 
                                          fdr_threshold_negative = fdr_threshold_negative,
                                          differential_threshold_positive = differential_threshold_positive, 
                                          differential_threshold_negative = differential_threshold_negative
                                          
      )
      
      # plot drug responses 
      plot_drug_response(curr_screen_score, 
                         control_name = control, 
                         condition_name = condition, 
                         output_folder = plots_folder,
                         neg_type = neg_type, 
                         pos_type = pos_type,
                         plot_type = plot_type, 
                         label_fdr_threshold = label_fdr_threshold)
      
      col_list <- c(
        paste0('differential_',condition,'_vs_',control,sep=''),
        paste0('significant_',condition,'_vs_',control,sep=''),
        paste0('effect_type_',condition,sep=''))
      all_score[,col_list] <- curr_screen_score[,col_list]
      
    }
  }   

  # save scores from all screens
  scores_fname <- file.path(output_folder, "scores_all.csv")
  write.table(all_score, scores_fname, sep = ",", row.names = F, col.names = TRUE, quote = FALSE)
  
  #create dataframe with only the FDR columns from all screens
  fdr_scores <- all_score[,grepl("gene|fdr", colnames(all_score))] # get column 'gene' and any columns with name containing string 'fdr'
  fdr_scores <- abbreviate_names(fdr_scores, "fdr_", 2:ncol(scores)) # reformat 'fdr' column names
  rownames(fdr_scores) <- fdr_scores$gene # set matrix rownames to gene names ('gene' column (first column) contains gene names at this point)
  fdr_scores <- fdr_scores[,-1] # remove first column ('gene' column)
  # save FDR score dataframe to file
  scores_fdr_fname <- file.path(output_folder, "fdr_scores_all.tsv")
  write.table(fdr_scores, scores_fdr_fname, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}
