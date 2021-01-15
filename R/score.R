######
# SCORING CODE
######

# Inner function to scale values between 0 and 1
scale_values <- function(x) {
  val <- (x-min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

#' Scores conditions against a single control.
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
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
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
                                   loess = TRUE, fdr_method = "BY",
                                   return_residuals = TRUE, verbose = FALSE) {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  control_cols <- screens[[control_name]][["replicates"]]
  condition_names <- c()
  condition_cols <- list()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
    condition_cols[[condition]] <- screens[[condition]][["replicates"]]
  }
  
  # Removes control genes
  df <- df[!(df$gene %in% control_genes),]
  
  # Makes output dataframe
  unique_genes <- unique(df$gene)
  n_genes <- length(unique_genes)
  scores <- data.frame(gene = rep(NA, n_genes))
  
  # Makes residual dataframes if necessary
  max_guides <- -1
  condition_residuals <- list()
  if (test == "moderated-t") {
    
    # Gets max number of guides first
    max_guides <- max(table(df$gene))
    
    # Makes residual dataframes with columns equal to the max number of guides
    for (name in condition_names) {
      residual_df <- data.frame(matrix(nrow = n_genes, ncol = max_guides))
      colnames(residual_df) <- paste0("guide_residual_", 1:max_guides)
      condition_residuals[[name]] <- residual_df
    }
  }
  
  # Makes loess residual dataframe if specified
  loess_residuals <- NULL
  if (loess & test == "moderated-t") {
    loess_residuals <- data.frame(n = rep(0, max_guides*n_genes))
    loess_residuals[[paste0("mean_", control_name)]] <- rep(0, nrow(loess_residuals))
    for (name in condition_names) {
      loess_residuals[[paste0("mean_", name)]] <- rep(0, nrow(loess_residuals))
    }
  }
  
  # Appends additional columns for each condition
  new_cols <- c(paste0("n_", control_name), 
                paste0("mean_", control_name),
                paste0("variance_", control_name))
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
  
  # Scores guides for each condition
  counter <- 1
  for (i in 1:n_genes) {
    
    # Gets gene names and control guide values across replicates and removes 
    # NaNs introduced by missing guides
    gene <- unique_genes[i]
    guide_vals <- df[df$gene == gene,]
    scores$gene[i] <- gene
    rep_mean_control <- rowMeans(data.frame(guide_vals[control_cols]), na.rm = TRUE)
    keep_ind <- !is.nan(rep_mean_control)
    rep_mean_control <- rep_mean_control[keep_ind]
    
    # Skips if too few guides
    if (length(rep_mean_control) < min_guides) {
      next
    }
    
    # Makes loess-normalized residual dataframe if necessary
    ind <- counter:(counter + length(rep_mean_control) - 1)
    if (loess) {
      loess_residuals$n[ind] <- i
      loess_residuals[[paste0("mean_", control_name)]][ind] <- rep_mean_control
    }
    
    # Takes the mean across replicates for all conditions
    for (name in condition_names) {
      
      # Gets residual LFCs across replicates after removing NaNs
      rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]), na.rm = TRUE)
      rep_mean_condition <- rep_mean_condition[keep_ind]
      rep_mean_condition[is.nan(rep_mean_condition)] <- NA
      diff <- rep_mean_condition - rep_mean_control
      
      # Stores gene-level stats
      scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control)
      scores[[paste0("n_", name)]][i] <- length(rep_mean_condition)
      scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control, na.rm = TRUE)
      scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition, na.rm = TRUE)
      scores[[paste0("variance_", control_name)]][i] <- stats::var(rep_mean_control, na.rm = TRUE)
      scores[[paste0("variance_", name)]][i] <- stats::var(rep_mean_condition, na.rm = TRUE)
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff, na.rm = TRUE)
      
      # Appends mean LFCs for loess-normalization if specified
      if (loess) {
        loess_residuals[[paste0("mean_", name)]][ind] <- rep_mean_condition
      }
      
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        scores[[paste0("pval_", name, "_vs_", control_name)]][i] <- 
          suppressWarnings(stats::wilcox.test(rep_mean_condition, rep_mean_control))$p.value
      } else if (test == "moderated-t") {
        if (length(diff) < max_guides) { diff <- c(diff, rep(NA, max_guides - length(diff))) } 
        condition_residuals[[name]][i,1:max_guides] <- diff 
      }
    }
    counter <- counter + length(rep_mean_control)
  }
   
  # Computes loess-normalized residuals if specified
  if (loess & test == "moderated-t") {
    loess_residuals <- loess_residuals[1:counter,]
    control_values <- loess_residuals[[paste0("mean_", control_name)]]
    for (name in condition_names) {
      condition_values <- loess_residuals[[paste0("mean_", name)]]
      temp <- loess_MA(control_values, condition_values)
      loess_residuals[[paste0("loess_residual_", name)]] <- temp[["residual"]]
      loess_residuals[[paste0("loess_predicted_", name)]] <- temp[["predicted"]]
    }
    
    # Replaces residuals with loess-normalized residuals
    for (i in 1:n_genes) {
      for (name in condition_names) {
        resid <- loess_residuals[[paste0("loess_residual_", name)]][loess_residuals$n == i]
        predicted <- loess_residuals[[paste0("loess_predicted_", name)]][loess_residuals$n == i]
        if (length(resid) < max_guides) { 
          resid <- c(resid, rep(NA, max_guides - length(resid))) 
        } 
        condition_residuals[[name]][i,1:max_guides] <- resid
        scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(resid, na.rm = TRUE)
        # scores[[paste0("loess_predicted_", name)]][i] <- mean(predicted)
      }
    }
  } else if (loess) {
    cat("Warning: loess-normalization is only enabled for the test=\"moderated-t\" option\n")
  }
  
  # Scores condition response with moderated t-test
  if (test == "moderated-t") {
    for (name in condition_names) {
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]]))
      p_val <- ebayes_fit$p.value[,1]
      scores[[paste0("pval_", name, "_vs_", control_name)]] <- p_val
    }   
  }
  
  # Computes FDRs
  for (name in condition_names) {
    scores[[paste0("fdr_", name, "_vs_", control_name)]] <- 
      stats::p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = fdr_method)
  }
  
  # Removes genes with too few observations
  scores <- scores[scores[[paste0("n_", control_name)]] >= min_guides,]
  
  # Removes extra zero row from residuals
  loess_residuals <- loess_residuals[1:(nrow(loess_residuals) - 1),]
  
  # Explicitly returns scored data
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals) {
    output[["residuals"]] <- loess_residuals
  } else {
    output[["residuals"]] <- NA
  }
  return(output)
}

#' Call significant responses for scored chemogenomic data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_drugs_vs_control}. 
#' 
#' @param scores Dataframe returned from \code{score_drugs_vs_control}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_drug_hits <- function(scores, control_screen_name, condition_screen_names,
                           fdr_threshold = 0.1, differential_threshold = 0,
                           neg_type = "Negative", pos_type = "Positive",
                           fdr_method = "BH") {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  condition_names <- c()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
  }
  
  # Calls significant differences for each condition against the control
  for (name in condition_names) {
    scores[[paste0("significant_", name, "_vs_", control_name)]] <- 
      scores[[paste0("fdr_", name, "_vs_", control_name)]] < fdr_threshold
  }
  
  # Makes thresholded calls for significant negative and positive effects
  for (name in condition_names) {
    response_col <- paste0("effect_type_", name)
    scores[[response_col]] <- "None"
    diffs <- scores[[paste0("differential_", name, "_vs_", control_name)]]
    sig <- scores[[paste0("significant_", name, "_vs_", control_name)]]
    scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
    scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
  }
  
  # Explicitly returns scored data
  return(scores)
}

# Fits a loess curve to predict y given x
loess_MA <- function(x, y, sp = 0.4, dg = 2, binSize = 100) {
  #this concept is based on pythagoras and cancels out sqrt, square and factor 2
  #it also ignores the factor sqrt(2) as factor between y and x vs distance of x,y from diagonal x = y
  if(all(x == y, na.rm = T)) { #if e.g. wt scored against itself
    gi <- rep(NA, length(x))
  }
  else {
    m <- y - x
    a <- y + x
    A <- (a - stats::median(a, na.rm = T)) / stats::mad(a, na.rm = T) #scale to generate bins along m
    B <- seq(floor(min(A, na.rm = T)), ceiling(max(A, na.rm = T)), .1) #define bins
    b <- c() #indices for model training
    for(i in 1:(length(B) - 1)) {
      temp_b <- which(A > B[i] & A < B[i+1])
      if(length(temp_b) > binSize) { #sample if more events in bin than max bin size
        set.seed(1)
        temp_b <- sample(temp_b, binSize)
      }
      b <- c(b, temp_b)
    }
    I <- is.finite(m[b]) & is.finite(a[b]) #only use finite values
    model <- stats::loess(m[b][I] ~ a[b][I], span = sp, degree = dg) #train model on m ~ a (approx. y ~ x)
    expected <- stats::predict(model, a) #predict expected m ~ a
    gi <- m - expected
  }
  result <- list()
  result[["residual"]] <- gi
  result[["predicted"]] <- expected
  return(result)
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
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param plot_residuals If true, plots residual effects for all top hits (default TRUE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
score_drugs_batch <- function(df, screens, batch_file, output_folder, 
                              min_guides = 3, test = "moderated-t", 
                              loess = TRUE, control_genes = c("None", ""), fdr_method = "BY",
                              fdr_threshold = 0.1, differential_threshold = 0.5, 
                              neg_type = "Negative", pos_type = "Positive",
                              plot_residuals = TRUE, plot_type = "png") {
  
  # Checks batch file and loads it
  check_batch_file(batch_file, screens)
  batch <- utils::read.csv(batch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Makes output folders if nonexistent
  lfc_folder <- file.path(output_folder, "lfc")
  plot_folder <- file.path(output_folder, "plots")
  if (!dir.exists(output_folder)) { dir.create(output_folder) }
  if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
  if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
  
  # Scores guides for each batch
  all_scores <- NULL
  for (i in 1:nrow(batch)) {
    condition <- batch[i,1]
    control <- batch[i,2]
    temp <- score_drugs_vs_control(df, screens, control, condition, test = test, 
                                   min_guides = min_guides, loess = loess, 
                                   control_genes = control_genes, fdr_method = fdr_method)
    scores <- temp[["scored_data"]]
    residuals <- temp[["residuals"]]
    scores <- call_drug_hits(scores, control, condition,
                             neg_type = neg_type, pos_type = pos_type,
                             fdr_threshold = fdr_threshold, 
                             differential_threshold = differential_threshold)
    plot_drug_response(scores, control, condition, plot_folder,
                       neg_type = neg_type, pos_type = pos_type,
                       plot_type = plot_type)
    if (plot_residuals) {
      plot_drug_residuals(scores, residuals, control, condition, lfc_folder, 
                          neg_type = neg_type, pos_type = pos_type,
                          plot_type = plot_type)
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
