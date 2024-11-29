######
# SCORING CODE
######

#' @importFrom magrittr "%>%"

# Inner function to scale values between 0 and 1
scale_values <- function(x) {
  val <- (x-min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

#' Finds and removes (set to NA) outlier guides from the conditional residual (differential scores) for each gene
#'
#' @param cond_res	dataframe of condition residuals (differential scores); rows: genes, column: replicate-guide pair 
#' @param threshold	jack-knife threshold to decide if a guide is an outlier (default value 2)
#' @return	dataframe of updated condition residuals (differential scores); rows: genes, column: replicate-guide pair 
jackknife_outliers<-function(cond_res, threshold=2)
{  
	#Find and remove (set to NA) outlier for each gene
  for(i in c(1:nrow(cond_res))){
    x = as.numeric(as.vector(cond_res[i,])) #get differential scores for a gene for all replicate-guide pairs
    var_all = var(x,na.rm = T) #calculate variance of differential scores across all replicate-guide pairs
    
	n <- length(x)
    var_jk <- rep(0, n)
	#For each replicate-guide pair, calculate variance excluding the differential score for that pair
    for(j in 1:n) {
      var_jk[j] <- var(x[ - j],na.rm = T)
    }
	#Find the replicate-guide pairs removing which decreses the variance a certain amount (threshold is applied here) 
    rel = (var_all/var_jk) > threshold 
    x[rel] <- NA #set values for the outliers to NA
    cond_res[i,] <-x #update replicate-guide pair values for gene
  }
  return (cond_res)
}


#' Scores conditions against a single control.
#' Revised one-off score function with individual replicate-level loess fitting (loess smoothing Drug A vs DMSO A, B vs B, C vs C per guide)
#' 
#' Scores guides for any number of drug screens against a control screen
#' (e.g. for directly comparing drug response to DMSO response). After running 
#' this function, pass the resulting dataframe to \code{call_drug_hits} to 
#' call significant effects.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
#' @param sd_scale If TRUE, apply standard-deviation scaling to differential LFC scores.
#'   Only works when test = "moderated-t" (default FALSE).
#' @param return_residuals If FALSE, returns NA instead of residuals dataframe (default TRUE).
#'   This is recommend if scoring large datasets and memory is a limitation.  
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A list containing two dataframes. The first entry, named "scored_data" in the list,
#'   contains scored data with separate columns given by the specified control and condition
#'   names. The second entry, named "residuals" in the list, is a dataframe containing control,
#'   condition and loess-normalized residuals for all guides.
#' @export
score_drugs_vs_control <- function(df, screens, control_screen_name, condition_screen_names, 
                                   control_genes = c("None", ""), min_guides = 3, test = "moderated-t", 
                                   loess = TRUE, ma_transform = TRUE, fdr_method = "BY",
                                   sd_scale = FALSE, return_residuals = TRUE, verbose = FALSE) {
  
  
	if (verbose) {
	cat(paste("Preparing to score...\n"))
	}
	# Get name of control screen and get control screen replicate names (data originally from sample table tsv file)
	control_name <- control_screen_name 
	control_cols <- screens[[control_name]][["replicates"]] 
	
	# Get condition screen names and replicates (There can be multiple condition screens; data originally from sample table tsv file)	
	condition_names <- c()
	condition_cols <- list()
	for (condition in condition_screen_names) {
	condition_names <- c(condition_names, condition)
	condition_cols[[condition]] <- screens[[condition]][["replicates"]]
	}
	
	# Removes control genes from the LFC dataframe
	df <- df[!(df$gene %in% control_genes),]
	
	# Makes output dataframe (add a column with unique genes from LFC dataframe)
	unique_genes <- unique(df$gene)
	n_genes <- length(unique_genes)
	scores <- data.frame(gene = rep(NA, n_genes))
	# Appends additional columns for control screen (n_: # of guides, mean_: average LFC across replicates, variance_: variance of LFC across replicates)
	new_cols <- c(paste0("n_", control_name), 
		paste0("mean_", control_name),
		paste0("variance_", control_name))
	# Appends additional columns for condition screens
	# n_: # of guides, mean_: average LFC across replicates, variance_: variance of LFC across replicates
	# differential_: differential LFC scores calculated against control screen
	# pval_: significance of differential score, fdr_: pval corrected for multiple test, significant_: annotate if DLFC is significant or not
	for (name in condition_names) {
	    new_cols <- c(new_cols, c(
	      paste0("n_", name), 
	      paste0("mean_", name),
	      paste0("variance_", name),
	      paste0("differential_", name, "_vs_", control_name),
	      paste0("pval_", name, "_vs_", control_name),
	      paste0("fdr_", name, "_vs_", control_name),
	      paste0("significant_", name, "_vs_", control_name)
	    ))
	}
	scores[new_cols] <- NA
  
	# Make a dataframe for each condition screen to store differential LFC scores later (if the significance test is moderated-t)
	max_guides <- -1
	condition_residuals <- list()
	if (test == "moderated-t") { # if significance test is moderated-t
		  # Gets maximum number of guides per gene (Gene names in LFC dataframe appear equal to the number of guides for that gene, so taking the frequency per gene will count guide per gene)
		  max_guides <- max(table(df$gene))
		  # Makes residual dataframes with column number equal to the max number of guides
		  for (name in condition_names) { # iterate over condition screens
			  # Extract condition screen replicate names 
			  condition_reps=condition_cols[[name]]
			  # create a dataframe with total rows equal to # of genes and total column equal to (maximum # of guides * # of replicates)
			  residual_df <- data.frame(matrix(nrow = n_genes, ncol = max_guides*length(condition_reps))) # ncol = max_guides*total_replicate
			  # Set column names. "guide_residual_guide#_replicate-name". 
			  # colname distribution example for max_guide 4 and total replicate # 3: 1 1 1 2 2 2 3 3 3 4 4 4. biological replicates (guides) are gathered together
			  colnames(residual_df) <- paste0("guide_residual_", rep(1:max_guides,each=length(condition_reps)),'_', condition_reps)
			  condition_residuals[[name]] <- residual_df
		  }
	}
  
	# Pairwise LOESS smoothing - pair-wise control-condition replicates
	# Compute loess-normalized differential LFC scores if specified and store them in a dataframe (loess_residuals)
	# These values will be later added to the output scores dataframe as differential LFC scores for condition screens
	loess_residuals <- NULL
	if (loess & test == "moderated-t") { # if 'loess' parameter is TRUE  and "moderated-t" test is specified
		  # Create dataframe to store loess normalized differential LFC scores
		  loess_residuals <- data.frame(gene = df$gene) # get gene names from LFC dataframe
		  for (name in condition_names) { # iterate over condition screens
			# Extract condition screen replicate names 
			condition_reps=condition_cols[[name]] 
			for(rep_index in c(1:length(condition_reps))){ # pair-wise control-condition replicates are extracted using index
				#get control replicate LFC scores from LFC dataframe
				control_values <- df[,paste0(control_cols[rep_index])]  
				# get condition replicate LFC scores from LFC dataframe
				rep = paste0(condition_reps[rep_index])
				condition_values <- df[,rep] 
				# Apply loess normalization from qGI_utils.R
				temp <- loess_MA(control_values, condition_values, ma_transform = ma_transform)
				# store loess normalized scores
				loess_residuals[[paste0("loess_residual_", name,'_',rep)]] <- temp[["residual"]]
				loess_residuals[[paste0("loess_predicted_", name,'_',rep)]] <- temp[["predicted"]]
			}
		}
	}
  
	# calculate various scores and update output scores per gene across condition screens
	for (i in 1:n_genes) { # iterate over unique gene indexes
		# Gets gene names and control guide values across replicates and removes NaNs introduced by missing guides
		# 
		gene <- unique_genes[i] # get the current gene
		guide_vals <- df[df$gene == gene,] # get the rows associated with the current gene from the LFC dataframe, each row corresponds to a guide
		scores$gene[i] <- gene # update output scores dataframe to add current gene name
		rep_mean_control <- rowMeans(data.frame(guide_vals[control_cols]), na.rm = TRUE) # mean-collapse control screen LFC scores across replicates - one value per guide
		keep_ind <- !is.nan(rep_mean_control) # get indexes with non-nan mean values
		rep_mean_control <- rep_mean_control[keep_ind] # keep guides with non-nan mean values - in this way guides with all nan values will not be selected
		
		# Skip adding scores for current gene - if total guides that have non-nan control-replicate LFC means are fewer than a cutoff. 
		# If all control replicates for a guide have nan LFC values it will result in a nan control-replicate LFC mean for that guide. 
		# The cut-off is defined by the parameter min_guides (default = 3).
		# The gene will have nan values in the output scores dataframe
		if (length(rep_mean_control) < min_guides) {
			next
		}
		
		# Takes the mean across replicates for all conditions
		for (name in condition_names) { # iterate over condition screens
			# mean-collapse condition screen LFC scores across replicates - one value per guide
			rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]), na.rm = TRUE) 
			rep_mean_condition <- rep_mean_condition[keep_ind] # only keep guides where 'control-replicate LFC means' are non-nan
			rep_mean_condition[is.nan(rep_mean_condition)] <- NA # set NA if condition-replicate LFC mean is nan
			diff <- rep_mean_condition - rep_mean_control # Calculate guide-level differential LFC scores - one value per guide
			
			# Update output scores dataframe to gene-level values
			scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control) # Total number of control guides contributing to scores
			scores[[paste0("n_", name)]][i] <- length(rep_mean_condition) # Total number of condition guides contributing to scores (should be same as Total number of control guides contributing to scores)
			scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control, na.rm = TRUE) # mean-collapse the guide-level control LFC scores
			scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition, na.rm = TRUE) # mean-collapse the guide-level condition LFC scores
			scores[[paste0("variance_", control_name)]][i] <- stats::var(rep_mean_control, na.rm = TRUE) # calculate variance of guide-level control LFC scores
			scores[[paste0("variance_", name)]][i] <- stats::var(rep_mean_condition, na.rm = TRUE)  # calculate variance of guide-level condition LFC scores
			scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff, na.rm = TRUE) # mean-collapse the guide-level differential LFC scores
			# differential LFC score column may be updated later based on different parameter settings
			
			# Perform the specified type of significance testing or store differential LFC scores for later testing
			if (test == "rank-sum") { # If 'rank-sum' test is specified
				# calculate p-values indicating the significance of differential LFC scores using wilcoxon rank-sum test
				# wilcoxon rank-sum test compares the guide-level control and condition mean LFC scores
				# update pval column for condition screen in scores dataframe
				scores[[paste0("pval_", name, "_vs_", control_name)]][i] <- 
				  suppressWarnings(stats::wilcox.test(rep_mean_condition, rep_mean_control))$p.value
			
			} else if (test == "moderated-t") { # If 'moderated-t' test is specified
				#create temporary list to store loess-normalized differential LFC values of current condition screen replicates
				#loess_residual_rep = c() 
				condition_reps = condition_cols[[name]] # get condition replicate names
				for(rep_index in c(1:length(condition_reps))){ # iterate over condition replicates
					
					if (loess) { # if 'loess' parameter is TRUE use the loess-normalized differential LFC scores calculated previously
						# Get the guide-level loess-normalized differential scores for the current gene and the current condition replicate - one score per-guide
						resid <- loess_residuals[[paste0("loess_residual_", name, '_', condition_reps[rep_index])]][loess_residuals$gene == gene]
						predicted <- loess_residuals[[paste0("loess_predicted_", name, '_', condition_reps[rep_index])]][loess_residuals$gene == gene]
						
						# update the condition_residuals matrix (that was created previously to store differential LFC scores) with loess-normalized scores
						if (length(resid) < max_guides) { # Pad temporary vector that stores guide-level dLFC with NAs to make size equal to max-guides
							resid <- c(resid, rep(NA, max_guides - length(resid))) 
						}
						condition_residuals[[name]][i, paste0("guide_residual_", 1:max_guides,'_', rep(condition_reps[rep_index], max_guides))] <- resid
						
						# update temporary list with mean of guide-level loess-normalized differential LFC values of current condition screen replicate 
						#loess_residual_rep <- c(loess_residual_rep, mean(resid, na.rm = TRUE))
					
					}else{ # if 'loess' parameter is FALSE
						# calculate differential LFC by directly subtracting control replicate LFC score from condition replicate LFC score - replicates matched by index
						resid <- guide_vals[condition_reps[rep_index]] - guide_vals[control_cols[rep_index]]
						resid <- resid[[condition_reps[rep_index]]] 
						
						# update the condition_residuals matrix (that was created previously to store differential LFC scores) with dLFC score
						if (length(resid) < max_guides) { # Pad temporary vector that stores guide-level dLFC with NAs to make size equal to max-guides
							resid <- c(resid, rep(NA, max_guides - length(resid))) 
						} 
						condition_residuals[[name]][i, paste0("guide_residual_", 1:max_guides,'_', rep(condition_reps[rep_index], max_guides))] <- resid 
					}
				}
				#if (loess) { # if 'loess' parameter is TRUE use the loess-normalized differential LFC scores calculated previously
					## Update differential LFC values with mean-collapsed Loess-normalized differential LFC values across replicates
					#scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(loess_residual_rep, na.rm = TRUE)
				#}
			}
		}
	}
  
	#remove (set to NA) outlier guides per-gene using jack-knife method for each condition if the significance test is moderated-t 
	if(test == "moderated-t"){
		for (name in condition_names) {
			# remove outlier guides per-gene (set to NA) across replicates for a condition screen
			condition_residuals[[name]] = jackknife_outliers(condition_residuals[[name]])
			# mean-collapse differential LFC scores across guides and replicates for each gene
			# update dLFC column for the condition screen in scores dataframe
			scores[[paste0("differential_", name, "_vs_", control_name)]] <- rowMeans(condition_residuals[[name]],na.rm = TRUE)
		}
	}
  
	  # calculate p-values indicating significance of differential LFC scores using moderated t-test if specified
	if (test == "moderated-t") {
		for (name in condition_names) {
			block <- rep(1:max_guides, each=length(condition_cols[[name]])) #group same guides together - all technical-replicates in each guide block
			dupcor <- limma::duplicateCorrelation(condition_residuals[[name]], block=block) #Estimate the correlation between technical replicates from a series of arrays ( here calculating inter-biological replicate (guides) correlation)
			ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]], design=NULL, block=block, correlation=dupcor$consensus)) # compute moderated t-statistics of dLFC by empirical Bayes moderation of the standard errors towards a common value
			p_val <- ebayes_fit$p.value[,1]
			scores[[paste0("pval_", name, "_vs_", control_name)]] <- p_val # update pval column for condition screen in scores dataframe
		}   
	}

	# Correct p-values for multiple testing and add a new column with FDR values
	# update FDR column for condition screen in scores dataframe
	for (name in condition_names) {
		scores[[paste0("fdr_", name, "_vs_", control_name)]] <- 
		stats::p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = fdr_method)
	}
  
	# Scales moderate effects in top and bottom 10% of data to de-emphasize those if a scaling factor is provided. 
	# The mean to divide SD values by is a pre-computed scalar
	if (test == "moderated-t") {
		if (sd_scale) {
			for (name in condition_names) { # iterate over condition screens
				resid <- condition_residuals[[name]] # get the differential LFC matrix for the current condition screen
				# calculate target standard deviation per screen(SD of scores between 10%-90% percentiles)
				lfc_range <- stats::quantile(resid, probs = c(0.1, 0.9), na.rm = TRUE) 
				target_sd <- stats::sd(resid[resid > lfc_range[1] & resid < lfc_range[2]], na.rm = TRUE)
				# normalize target standard deviation
				sd_scale_factor <- mean(target_sd)
				target_sd <- target_sd / sd_scale_factor
				# normalize per-screen dlfc scores by target standard deviation 
				condition_residuals[[name]] <- resid / target_sd
				# mean-collapse differential LFC scores across guides and replicates for each gene
				# update dLFC column for the condition screen in scores dataframe
				scores[[paste0("differential_", name, "_vs_", control_name)]] <- rowMeans(condition_residuals[[name]],na.rm = TRUE)
			} 
		}
	}
  
  # Explicitly returns scored data
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals & loess) {
    output[["residuals"]] <- loess_residuals
  } else if (return_residuals) {
    cat("WARNING: returning residuals is currently only supported with loess-normalization enabled\n")
    output[["residuals"]] <- NA
  } else {
    output[["residuals"]] <- NA
  }
  return(output)
}

#' Call significant responses for scored chemogenomic data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_drugs_vs_control} or \code{score_drugs_vs_controls}.
#' 
#' @param scores Dataframe returned from \code{score_drugs_vs_control} or 
#'   \code{score_drugs_vs_controls}.
#' @param control_screen_name Name of a control screen to test condition screens against, 
#'   or NULL to score data returned from \code{score_drugs_vs_controls}.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param fdr_threshold_positive Threshold below which to call gene effects as significant positive hits
#'   (default 0.1).
#' @param fdr_threshold_negative Threshold below which to call gene effects as significant negative hits
#'   (default 0.1).
#' @param differential_threshold_positive Threshold on differential effects, 
#'   over which to call gene effects as significant positive hits (default 0).
#' @param differential_threshold_negative Threshold on differential effects, 
#'   below which gene effects are called as significant negative hits (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_drug_hits <- function(scores, control_screen_name = NULL, condition_screen_names = NULL,
                           fdr_threshold_positive  = 0.1, fdr_threshold_negative = 0.1,
			   differential_threshold_positive = 0, differential_threshold_negative = 0,
                           neg_type = "Negative", pos_type = "Positive") {
  
	# Get name of control screen, condition screen, and condition screen replicate names (data originally from sample table tsv file)
	control_name <- control_screen_name
	condition_names <- c()
	for (condition in condition_screen_names) {
		condition_names <- c(condition_names, condition)
	}
  
	# Call significant differences for each condition against the control
	# if the fdr score is less than the fdr_threshold_positive or fdr_threshold_negative parameters 
	for (name in condition_names) {
		scores[[paste0("significant_", name, "_vs_", control_name)]] <- 
		(scores[[paste0("fdr_", name, "_vs_", control_name)]] < fdr_threshold_positive) |
		(scores[[paste0("fdr_", name, "_vs_", control_name)]] < fdr_threshold_negative) 
		}
	
	# Makes thresholded calls for significant negative and positive effects
	for (name in condition_names) {
		response_col <- paste0("effect_type_", name)
		scores[[response_col]] <- "None"
		diffs <- scores[[paste0("differential_", name, "_vs_", control_name)]]
		sig <- scores[[paste0("significant_", name, "_vs_", control_name)]]
		scores[[response_col]][sig & diffs < differential_threshold_negative] <- neg_type
		scores[[response_col]][sig & diffs > differential_threshold_positive] <- pos_type
		}
  
  # Explicitly returns scored data
  return(scores)
}

#' Scores multiple drugs against multiple controls
#' 
#' Takes in an input .tsv file with two columns for "Screen" and "Control" and scores
#' all screens listed in "Screen" against their corresponding screens listed in
#' "Control." Outputs all files and plots in the specified folder. 
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param batch_file Path to .tsv file mapping screens to their controls for scoring, with two 
#'   columns for "Screen" and "Control." Screens to score against dervied null-models with the
#'   combn scoring mode must have their respective control labeled as "combn." 
#' @param output_folder Folder to output scored data and plots to. 
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param sd_scale If TRUE, apply standard-deviation scaling to differential LFC scores.
#'   Only works when test = "moderated-t" (default FALSE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold_positive Threshold below which to call gene effects as significant positive hits
#'   (default 0.1).
#' @param fdr_threshold_negative Threshold below which to call gene effects as significant negative hits
#'   (default 0.1).
#' @param differential_threshold_positive Threshold on differential effects, 
#'   over which to call gene effects as significant positive hits (default 0).
#' @param differential_threshold_negative Threshold on differential effects, 
#'   below which gene effects are called as significant negative hits (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param label_fdr_threshold Threshold below which to plot gene labels for significant
#'   hits, or NULL to plot without labels (default NULL).
#' @param save_residuals If true, saves residuals for each screen to the output folder
#'   (default FALSE).
#' @param plot_residuals If true, plots residual effects for all top hits (default TRUE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @param verbose If true, prints verbose output (default FALSE). 
#' @export
score_drugs_batch <- function(df, screens, batch_file, output_folder, 
                              min_guides = 3, test = "moderated-t", 
                              loess = TRUE, ma_transform = TRUE,
                              control_genes = c("None", ""), sd_scale = FALSE,
                              fdr_method = "BY", 
			      fdr_threshold_positive  = 0.1, fdr_threshold_negative = 0.1,
			      differential_threshold_positive = 0, differential_threshold_negative = 0,
			      neg_type = "Negative", 
                              pos_type = "Positive", label_fdr_threshold = NULL,
                              save_residuals = FALSE, plot_residuals = TRUE, 
                              plot_type = "png", verbose = FALSE) {
  
  # Creates output folder if nonexistent
  if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
  
  # Checks format of batch file
  
  first_file <- utils::read.table(file = batch_file, header = F, nrows = 1, sep = "\t", encoding = "UTF-8")
  batch <- NULL
  if (ncol(first_file) == 2) {
    check_batch_file(batch_file, screens)  
  } else {
    stop(paste("file", batch_file, "must contain exactly 2 columns (see score_drugs_batch documentation)"))
  }
  batch <- utils::read.csv(batch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Scores guides for each batch
	all_scores <- NULL
	for (i in 1:nrow(batch)) {
	
	# Makes output folders if nonexistent
	lfc_folder <- file.path(output_folder, "lfc")
	plot_folder <- file.path(output_folder, "plots")
	if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
	if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
	
	# Scores each screen separately
	condition <- batch[i,1]
	control <- batch[i,2]
	temp <- score_drugs_vs_control(df, screens, control, condition, test = test, 
				     min_guides = min_guides, 
				     loess = loess, 
				     ma_transform = ma_transform, 
				     control_genes = control_genes, 
				     fdr_method = fdr_method, 
				     sd_scale = sd_scale,
				     verbose = verbose)
	scores <- temp[["scored_data"]]
	residuals <- temp[["residuals"]]
	scores <- call_drug_hits(scores, control, condition,
			       neg_type = neg_type, pos_type = pos_type,
			       fdr_threshold_positive = fdr_threshold_positive, 
				 fdr_threshold_negative = fdr_threshold_negative, 
			       differential_threshold_positive = differential_threshold_positive,
				 differential_threshold_negative = differential_threshold_negative
				)
	plot_drug_response(scores, 
			 control_name = control, 
			 condition_name = condition, 
			 output_folder = plot_folder,
			 neg_type = neg_type, 
			 pos_type = pos_type,
			 plot_type = plot_type, 
			 label_fdr_threshold = label_fdr_threshold)
	if (plot_residuals) {
	if (!is.na(residuals)) {
	  plot_drug_residuals(scores, residuals, control, condition, lfc_folder, 
			      neg_type = neg_type, pos_type = pos_type,
			      plot_type = plot_type) 
	} else {
	  cat("WARNING: residuals are set to NA, skipping residual plots\n")
	}
	}
	if (save_residuals) {
	if (!is.na(residuals)) {
	  residuals_file <- paste0(condition, "_vs_", control, "_residuals.tsv")
	  utils::write.table(residuals, file.path(output_folder, residuals_file), sep = "\t",
			     row.names = FALSE, col.names = TRUE, quote = FALSE) 
	} else {
	  cat("WARNING: residuals are set to NA, skipping writing residuals to file\n")
	}
	}
	if (is.null(all_scores)) {
	all_scores <- scores
	} else {
	all_scores <- cbind(all_scores, scores[,3:ncol(scores)])
	}
	}
	if (!is.null(all_scores)) {
	utils::write.table(all_scores, file.path(output_folder, "condition_gene_calls.tsv"), sep = "\t",
			 row.names = FALSE, col.names = TRUE, quote = FALSE) 
	} 
  
}
