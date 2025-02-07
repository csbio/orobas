######
# UTILITY FUNCTIONS FOR CONTROL EFFECT SCORING
######

### Some utility functions directly ported from the qGI scoring pipeline,
### written by Maximilian Billmann and Henry Ward, are located here.
### Updated by Xiang Zhang and Arshia Hassan for computing control effect scores.


#' Correlation check function for control replicates
#' Will pick the screens with poor control replicate correlations given a cutoff
#' 
#' @param df Reads or lfc dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder for output. 
#' @param pcc_cutoff Pearson correlation cutoff below which we define as poorly-correlated controls.
#' @param verbose If true, prints verbose output (default FALSE). 
#' @param control_keyword Keyword (e.g. DMSO) in the control replicate name for identification.
#' @return Two dataframes: first for screens passing the pcc cutoff, second for "bad" screens.
#' @export
batch_control_cor_check <- function(df, screens, output_folder, pcc_cutoff=0.5, 
                                    verbose = TRUE,
                                    control_keyword="DMSO|MOCK|Mock|WT|NGLY1|BMI1|Control") {
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Gets sample groups
  all_cols <- c()
  col_groups <- c()
  i <- 1
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    for (col in screen[["replicates"]]) {
      all_cols <- c(all_cols, col)
      col_groups[i] <- screen_name
      i <- i +1
    }
  }
  
  # Builds dataframe of replicate PCCs
  cor_df <- NULL
  
  # Compares replicates across all screens
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    rep_cols <- screen[["replicates"]]
    if (length(rep_cols) > 1) {
      pairs <- utils::combn(rep_cols, 2)
      for (i in 1:ncol(pairs)) {
        col1 <- pairs[1,i]
        col2 <- pairs[2,i]
        x_label <- paste0(col1, " log fold change")
        y_label <- paste0(col2, " log fold change")
        temp <- plot_samples(df, col1, col2, x_label, y_label, print_cor = verbose)
        p <- temp[[1]]
        pcc <- temp[[2]]
        scc <- temp[[3]]
        #file_name <- paste0(col1, "_vs_", col2, "_replicate_comparison.", plot_type)
        #suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300))
        
        # Stores PCC in dataframe
        if (is.null(cor_df)) {
          cor_df <- data.frame(screen = screen_name, rep1 = col1, rep2 = col2, pcc = pcc, scc = scc,
                               stringsAsFactors = FALSE)
        } else {
          cor_df <- rbind(cor_df, c(screen_name, col1, col2, pcc, scc))
        }
      }   
    }
  }
  
  # Identify control (DMSO) replicates with PCCs less than the cutoff
  cor_df <- cor_df[grepl(control_keyword, cor_df$screen),]
  # Get all unique screen names
  all_unique_screens <- unique(cor_df$screen)
  
  # Store all replicate PCC info
  cor_file <- file.path(output_folder, "replicate_cor.tsv")
  utils::write.table(cor_df, cor_file, quote = FALSE, sep = "\t",
                     row.names = FALSE, col.names = TRUE)
  
  cor_df<- cor_df[cor_df$pcc < pcc_cutoff, ]
  # Get the unique screen names for those using consensus control in the batch mode
  unique_consensus_screens <- unique(cor_df$screen)
  
  # Writes excluded screens based on PCCs to file
  #cor_file <- file.path(output_folder, "excluded_replicate_cor.tsv")
  #utils::write.table(cor_df, cor_file, quote = FALSE, sep = "\t",
  #                   row.names = FALSE, col.names = TRUE)
  
  # Get all other screen names which will be applied in the "match" setting
  unique_match_screens <- unique(all_unique_screens[!(all_unique_screens %in% unique_consensus_screens)])
  
  # Return a nested list consisting the "match" (good) screens and "consensus" (bad) screens
  return(list(unique_match_screens, unique_consensus_screens))
}


#' Computes effect sizes for control screens.
#' 
#' Computes rough estimates of effect sizes for control screens, scored against
#' all control screens passing QC with equal weights given to each.  
#' 
#' @param control_df LFC dataframe for control screens.
#' @param control_names List of screen names for control screens to compute effect
#'   sizes for.
#' @param bad_cor_controls A list of control screen names tested to have bad correlations.
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @return Dataframe of effect sizes for control screens.
#' @export
compute_control_effects <- function(control_df, control_names, bad_cor_controls,
                                    loess = TRUE, ma_transform = TRUE) {
  weights <- compute_control_weights(control_names=control_names, 
                                     condition_names=control_names, 
                                     bad_cor_controls=bad_cor_controls)
  control_scores <- control_df[,colnames(control_df) %in% c("gene", control_names)]
  for (i in 1:length(control_names)) {
    control <- control_names[i]
    temp <- as.matrix(control_df[,2:ncol(control_df)])
    temp <- temp %*% diag(weights[i,])
    temp <- rowSums(temp, na.rm = TRUE) * NA^!rowSums(!is.na(temp))
    temp <- data.frame(gene = control_df$gene, lfc = temp)
    if (!loess) {
      control_scores[[control]] <- control_scores[[control]] - temp$lfc
    } else {
      temp <- loess_MA(temp$lfc,  control_scores[[control]], ma_transform = ma_transform)
      control_scores[[control]] <- temp[["residual"]]
    }
  }
  return(control_scores)
}


#' Computes weights for control screens.
#' 
#' Computes weights for control screens based on matching condition screens
#' for \code{score_controls_vs_controls}.
#' 
#' @param control_names List of screen names for control screens to compute effect
#'   sizes for.
#' @param condition_names List of screen names for control screens to compute effect
#'   sizes for.
#' @param bad_cor_controls A list of control screen names tested to have bad correlations. 
#' @return Dataframe of weights for each condition and control with the number of 
#'   rows equal to length(condition_names) and the number of columns equal to 
#'   length(control_names).
#' @export
compute_control_weights <- function(control_names, condition_names, 
                                    bad_cor_controls) {
  # If there is only 1 control, all weights are equal to 1
  weights <- matrix(1, nrow = length(condition_names), ncol = length(control_names))
  rownames(weights) <- condition_names
  colnames(weights) <- control_names
  if (ncol(weights) > 1) {
    #print(ncol(weights))
    for (i in 1:nrow(weights)) {
      # Gets the current condition, its matched control and all other controls
      condition <- rownames(weights)[i]
      # All good controls are given an equal weight
      weights[i,] <- 1 / (ncol(weights) - length(bad_cor_controls))
      if (length(bad_cor_controls) > 0) {
        weights[i, bad_cor_controls] <- 0
      }
    }
  }
  return(weights)
}


#' Scores conditions (=controls) against consensus controls.
#' 
#' Scores guides for any number of screens against multiple control screens
#' (e.g. for directly comparing control response to consensus control response).
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_names A list of control screen names to test the condition
#'   screens against.
#' @param bad_cor_controls A list of control screen names tested to have bad correlations. 
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A dataframe.
#'   Contains scored control effect data with separate columns for specified control
#'   names. It was aggregated from loess-normalized residuals for all guides.
#' @export
score_controls_vs_controls <- function(df, screens, control_screen_names,
                                       bad_cor_controls,
                                       control_genes = c("None", ""), 
                                       min_guides = 2, loess = TRUE, 
                                       ma_transform = TRUE,
                                       verbose = FALSE) {
  
  # Disables dplyr warnings
  options(dplyr.summarise.inform = FALSE)
  
  # Gets condition names and columns for all specified conditions and controls
  control_names <- c()
  control_cols <- list()
  
  for (control in control_screen_names) {
    control_names <- c(control_names, control)
    control_cols[[control]] <- screens[[control]][["replicates"]]
  }
  #print(control_cols)
  
  # Removes control genes
  df <- df[!(df$gene %in% control_genes),]
  
  # Sorts dataframe by alphabetical order of genes
  unique_genes <- sort(unique(df$gene))
  df <- df[order(match(df$gene, unique_genes)),]
  
  # Takes replicate averages for all controls and conditions
  control_df <- data.frame(gene = df$gene)
  for (control in names(control_cols)) {
    reps <- df[,control_cols[[control]]]
    control_df[[control]] <- rowMeans(reps, na.rm = TRUE)
  }
  
  # Ensures that NaNs are converted to NAs
  control_df[is.na(control_df)] <- NA
  
  # Generate the residual dataframe based on the difference between each condition (control)
  # versus the consensus control
  residual_df <- compute_control_effects(control_df, control_names, bad_cor_controls,
                                         loess = loess, ma_transform = ma_transform)
  
  # Control effect scores
  control_effect_scores <- residual_df %>%
    dplyr::group_by(gene) %>%
    stats::na.omit() %>%
    dplyr::summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  
  return(control_effect_scores)
}


#' Scores multiple control screens against consensus controls
#' 
#' Checks whether the controls' average internal correlation pass the threshold (e.g. pcc=0.50),
#' based on which divide the screens into "good" and "bad".
#' The "good" group will be averaged into a consensus collection,
#' and all screens will be scored against average controls from all good screens.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param batch_file Path to .tsv file mapping screens to their controls for scoring, with four 
#'   columns for "Screen", "Control", "Group", and "Type".
#' @param output_folder Folder to output scored data and plots to. 
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param label_fdr_threshold Threshold below which to plot gene labels for significant
#'   hits, or NULL to plot without labels (default NULL).
#' @param save_residuals If true, saves residuals for each screen to the output folder
#'   (default FALSE).
#' @param verbose If true, prints verbose output (default FALSE). 
#' @param black_list A list of control screens to be excluded beforehand. 
#' @param screen_control_keyword Key words to identify control screen names. 
#' @export
score_controls_batch <- function(df, screens, batch_file, output_folder, 
                                 min_guides = 3, 
                                 loess = TRUE, ma_transform = TRUE,
                                 control_genes = c("None", ""),  
                                 qc_control_pcc=0.30,
                                 verbose = FALSE, 
                                 black_list=c('bad_screen_name', ""),
                                 screen_control_keyword="DMSO|MOCK|Mock|WT|NGLY1|BMI1|Control") {
  
  # Creates output folder if nonexistent
  if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
  
  # Checks format of batch file and loads it
  first_file <- utils::read.table(file = batch_file, header = F, nrows = 1, sep = "\t", encoding = "UTF-8")
  batch <- NULL
  if (ncol(first_file) == 4) {
    check_group_file(batch_file, screens, control_only=TRUE)
  } else {
    stop(paste("file", batch_file, "must contain exactly 4 columns"))
  }
  batch <- utils::read.csv(batch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Check control correlation and split into "good" and "bad" groups
  controls_correlation_list <- batch_control_cor_check(df, screens, output_folder, 
                                                       pcc_cutoff = qc_control_pcc, 
                                                       verbose = verbose,
                                                       control_keyword = screen_control_keyword)
  good_cor_controls <- controls_correlation_list[[1]]
  bad_cor_controls <- controls_correlation_list[[2]]
  if (verbose) {
    cat(paste("\nThese control screens are below the PCC cutoff", qc_control_pcc, "\n"))
    cat(paste("They will not get weights for the consensus control to score against\n"))
    print(bad_cor_controls)
    
    cat(paste("\nThese screens are in the black list and will be excluded if found\n"))
    print(black_list)
  }
  
  # Gets unique groups from batch file and scores each group separately
  groups <- unique(batch$Group)
  for (group in groups) {
    cat(paste("\nScoring group", group, "\n"))
    # Makes output folders for the group if nonexistent
    group_folder <- file.path(output_folder, group)
    if (!dir.exists(group_folder)) { dir.create(group_folder, recursive = TRUE) }
    
    # Get subset of the group batch
    group_ind <- batch$Group %in% group
    batch_group <- batch[group_ind, ]
    batch_group_controls <- batch_group[batch_group$Type %in% "control", ]
    
    # Scores data for the control screens 
    # The variable "control" is equivalent to "conditions" here
    controls <- unique(batch_group_controls$Screen)
    
    # Exclude screens in the black list
    controls <- controls[!(controls %in% black_list)]
    bad_cor_controls_group <- bad_cor_controls[(bad_cor_controls %in% controls)]
    
    control_effect_scores <- score_controls_vs_controls(df, screens, 
                                                        control_screen_names = controls,
                                                        bad_cor_controls = bad_cor_controls_group, 
                                                        min_guides = min_guides, 
                                                        loess = loess, 
                                                        ma_transform = ma_transform, 
                                                        control_genes = control_genes, 
                                                        verbose = verbose)
    utils::write.table(control_effect_scores, file.path(group_folder, "control_effect_scores.tsv"), sep = "\t",
                       row.names = FALSE, col.names = TRUE, quote = FALSE) 
  }
}