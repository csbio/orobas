######
# PRE-PROCESSING CODE
######

#' Normalizes reads for given screens
#' 
#' Log2 and depth-normalizes reads between a given list of columns and a given
#' column of an earlier timepoint (e.g. T0) specified in the "normalize_name"
#' entry of each screen. Screens with NULL for their "normalize_name" entry
#' are log2 and depth-normalized, but not normalized to earlier timepoints.
#' If a screen to normalize against has multiple replicates, those replicates
#' are averaged before normalization. Multiple replicates for a screen being
#' normalized, however, are normalized separately against the provided early
#' timepoint.
#' 
#' @param df Reads dataframe.
#' @param screens List of screens generated with \code{add_screens}. 
#' @param filter_names List of screen names to filter based on read counts by. 
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @param min_reads Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @return Normalized dataframe.
#' @export 
normalize_screens <- function(df, screens, filter_names = NULL, cf1 = 1e6, cf2 = 1, 
                              min_reads = 30, max_reads = 10000) {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Flags guides with too few read counts
  all_names <- names(screens)
  for (name in filter_names) {
    if (!(name %in% all_names)) {
      cat(paste("WARNING: screen", name, "not found, data not filtered by this screen\n"))
    }
  }
  to_remove <- rep(FALSE, nrow(df))
  filter_cols <- sapply(screens[filter_names], "[[", "replicates")
  for (col in filter_cols) {
    to_remove[df[,col] < min_reads] <- TRUE
  }
  sum_low <- sum(to_remove)
  for (col in filter_cols) {
    to_remove[df[,col] > max_reads] <- TRUE
  }
  sum_high <- sum(to_remove) - sum_low
  for (col in filter_cols) {
    to_remove[is.na(df[,col])] <- TRUE
  }
  sum_na <- sum(to_remove) - (sum_high + sum_low)
  removed_guides_ind <- which(to_remove)
  
  # Log2 and depth-normalizes every screen
  for (screen in screens) {
    for (col in screen[["replicates"]]) {
      df[,col] <- normalize_reads(df[,col], cf1, cf2)
    }
  }
  
  # Normalizes specified screens to earlier timepoints
  for (screen in screens) {
    normalize_name <- screen[["normalize_name"]]
    if (!is.null(normalize_name)) {
      if (normalize_name %in% all_names) {
        for (col in screen[["replicates"]]) {
          rep_cols <- screens[[normalize_name]][["replicates"]]
          rep_norm <- df[,rep_cols]
          if (length(rep_cols) > 1) {
            rep_norm <- rowMeans(rep_norm)
          }
          df[,col] <- df[,col] - rep_norm
        }
      } else {
        cat(paste("WARNING: screen", normalize_name, "not found.\n"))
      }
    }
  }
  
  # Removes flagged guides
  df <- df[!to_remove,]
  cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
  cat(paste("Excluded a total of", sum_high, "guides for high t0 representation\n"))
  cat(paste("Excluded a total of", sum_na, "guides with NA t0 values\n"))
  return(df)
}

#' Filters guides with too few read counts.
#' 
#' Filters guides out with too few read counts from a given reads dataframe
#' and a given set of columns (e.g. T0 columns).
#' 
#' @param df Reads dataframe.
#' @param cols Columns to filter by.
#' @param min_reads Minimum number of reads to keep (anything below
#'   this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (anything above
#'   this value will be filtered out).
#' @return Filtered dataframe.
filter_reads <- function(df, cols, min_reads = 30, max_reads = 10000) {
  to_remove <- rep(FALSE, nrow(df))
  for (col in cols) {
    to_remove[df[,col] < min_reads] <- TRUE
  }
  sum_low <- sum(to_remove)
  for (col in cols) {
    to_remove[df[,col] > max_reads] <- TRUE
  }
  sum_high <- sum(to_remove) - sum_low
  removed_guides_ind <- which(to_remove)
  df <- df[!to_remove,]
  cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
  cat(paste("Excluded a total of", sum_high, "guides for low t0 representation\n"))
  return(df)
}

#' Computes essential gene recovery AUC.
#' 
#' Computes area under the curve for ROC curves that measure how well each technical replicate
#' recovers signal for essential-targeting guides and saves results to file.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}. 
#' @param gene_col A column containing gene names for all guides.
#' @param output_folder Folder to which essential gene QC results should be saved.
#' @return All output is saved to the file "essential_PR_QC.txt" in the folder given by
#'   output_folder.
#' @export
essential_lfc_qc <- function(df, screens, gene_col, output_folder) {
  
  # Checks that the given gene_col is in the data
  if (!(gene_col %in% colnames(df))) {
    stop(paste("gene name column", gene_col, "not in df"))
  }
  
  # Loads essential gene standard from internal data
  essentials <- traver_core_essentials
  
  # Gets PR curves for all essential genes and all technical replicates
  output_file <- file.path(output_folder, "essential_PR_QC.txt")
  sink(output_file)
  if (sum(df[[gene_col]] %in% essentials) > 0) {
    for (screen in screens) {
      for (rep in screen[["replicates"]]) {
        temp <- df[!is.na(df[[rep]]),]
        ind <- temp[[gene_col]] %in% essentials
        roc <- PRROC::roc.curve(-temp[[rep]], weights.class0 = as.numeric(ind), curve = TRUE)
        cat(paste(rep, "essential-gene recovery AUC under ROC curve:", roc$auc, "\n")) 
      }
    }
  }
  sink()
}

#' Log-normalizes reads.
#' 
#' Log2-normalizes reads with a given pseudocount and scaling factor, and also
#' depth-normalizes the data. 
#' 
#' @param df List of read counts.
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @return Log- and depth-normalized read counts.
#' @export
normalize_reads <- function(df, cf1 = 1e6, cf2 = 1) {
  log2((df / sum(df, na.rm = TRUE)) * cf1 + cf2)
}

