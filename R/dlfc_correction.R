######
# Batch correction of differential log fold change scores from numerous screens
######

#' remove screens not to be included in the dLFC score data-frame
#' 
#' The scores matrix has library genes on the row side and compound screens on the column side.
#' The first column is gene names.
#' A list of screens can be provided to be excluded and those columns will be removed from the scores matrix.
#' 
#' @param scores  differential log fold change scores data-frame (library-genes X screens; The first column is gene names, each screen column contain dLFC scores)
#' @param exclude_screen_list a list of screens to be removed from the scores data-frame; default is empty list c().  
#' @return  differential log fold change scores data-frame (library-genes X screens; The first column is gene names) after removing selected screens (columns)
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
#' @param scores  differential log fold change scores data-frame (library-genes X screens; The first column is gene names, each screen column contain dLFC scores)
#' @return  differential log fold change scores data-frame (library-genes X screens; The first column is gene names) after applying SD scaling per screen (column)
#' 
#' @export
apply_SD_scaling<- function(scores)
{
	#calculate standard-deviation of each screen before scaling
	pre_scaling_sd <- apply(scores[,2:ncol(scores)], 2, stats::sd, na.rm=TRUE)

	#calculate target standard deviation per screen(SD of scores between 10%-90% percentiles)
	lfc_range <- apply(scores[,2:ncol(scores)], 2, stats::quantile, probs = c(0.1, 0.9), na.rm = TRUE)
	target_sd <- rep(NA, ncol(scores))
	for (i in 2:ncol(scores)) {
	  target_sd[i] <- stats::sd(scores[scores[,i] > lfc_range[1,i-1] & scores[,i] < lfc_range[2,i-1], i], na.rm = TRUE)
	}
	#normalize target standard deviation
	sd_scale_factor <- mean(target_sd[2:length(target_sd)])
	target_sd <- target_sd / sd_scale_factor

	#normalize per-screen dlfc scores by target standard deviation 
	for (i in 2:ncol(scores)) {
	  scores[,i] <- scores[,i] / target_sd[i]
	}
  
	#sanity check
	#calculate standard-deviation of each screen after scaling and write the scaling values to file
	# post_scaling_sd <- apply(scores[,2:ncol(scores)], 2, stats::sd, na.rm=TRUE)
	# sd_table <- data.frame(rbind(pre_scaling_sd, post_scaling_sd))
	# sd_table$source[1] <- "pre_scaling"
	# sd_table$source[2] <- "post_scaling"
	# sd_fname <- file.path("sd_table.tsv")
	# write.table(sd_table, sd_fname, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	return(scores)
}

#' remove principal component signal from dLFC scores
#' 
#' apply principal component analysis on the dLFC score matrix.
#' calculate projection of first mentioned number of principal component onto the score matrix.
#' subtract the projection from dLFC score matrix to remove principal component signal.
#' @param scores  differential log fold change scores data-frame (library-genes X screens; The first column is gene names, each screen column contain dLFC scores)
#' @param pc number of principal components to remove; default 1
#' @return  differential log fold change scores data-frame (library-genes X screens) after removing principal components
#' 
#' @export
remove_principal_component_signal<- function(scores,pc=1)
{
	rownames(scores) <- scores$gene #set matrix rownames to gene names ('gene' column (first column) contains gene names at this point)
	scores <- scores[,-1] #remove first column ('gene' column)
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
#' apply prncipal component analysis to the control dLFC scores.
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
  
	#apply principal component analysis to control dLFC scores
	pca_cg_dmso_t <- prcomp(as.matrix(t(cg_dmso)), center = TRUE, scale. = T)

	#create projection of all principal components from control-dLFCs onto the main(condition) dLFC score matrix
	cg_dmso_loadings_t <- pca_cg_dmso_t$rotation
	dmso_v_t <- cg_dmso_loadings_t[,1:ncol(cg_dmso_loadings_t)]
	projected_t <- t(scores) %*% dmso_v_t %*% t(dmso_v_t) 
	
	#subtract the the projection matrix from condition score matrix to remove dLFC signal of control screens 
	scores <- scores - t(projected_t)

	return(data.frame(scores))

}
