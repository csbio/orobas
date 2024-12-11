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
