#' Corrects for GIs shifted in specific chromosomes.
#' 
#' Run this within \code{score_drugs_vs_controls} to implement qGI chromosomal correction 
#' steps.
#' 
#' @param df LFC dataframe passed into \code{score_drugs_vs_controls} with chromosomal location 
#'   specified in a column named "chr", start coordinates specified in "start_loc" and stop 
#'   coordinates specified in "stop_loc." 
#' @param guide_df Dataframe of guide-level qGI scores.
#'   
#' @return Dataframe of chromosomal-corrected, guide-level qGI scores. 
correct_chromosomal_effects <- function(df, guide_df) {
  
  # Gets guide density for each relevant chromosome
  chrom <- names(sort(table(df$chr)))
  chr_features <- array(NA, dim = c(length(chrom), 4),
                        dimnames = list(chrom, c("length", "n_genes", "length_per", "roll_window")))
  for (i in 1:length(chrom)) {
    chr_features[i,"length"] <- max(df$stop_loc[df$chr == chrom[i]]) - min(df$start_loc[df$chr == chrom[i]])
    chr_features[i,"n_genes"] <- length(which(df$chr == chrom[i]))
  }
  chr_features[,3] <- chr_features[,2] / chr_features[,1]
  chr_features[,4] <- ceiling(log2(chr_features[,3] / stats::median(chr_features[,3]) + 2) * 150)
  
  # Ignores chromosomes with too few genes for the rolling window
  to_keep <- rep(TRUE, nrow(chr_features))
  for (i in 1:nrow(chr_features)) {
    if (chr_features[i,"n_genes"] < chr_features[i,"roll_window"]) {
      to_keep[i] <- FALSE
    }
  }
  chr_features <- chr_features[to_keep,]
  
  # Identifies windows where chromosomal shifts occur
  noise_factor <- apply(guide_df[,2:ncol(guide_df)], 2, stats::sd, na.rm = TRUE)
  chr_shifted_genes <- list()
  for (i in 1:nrow(chr_features)) {
    chrom <- row.names(chr_features)[i]
    for (j in 2:ncol(guide_df)) {
      condition <- colnames(guide_df)[j]
      
      # ???
      th1_match <- .09 + noise_factor[condition] * .2
      th2_match <- .03 + noise_factor[condition] * .1
      th3_match <- .04 + noise_factor[condition] * .2
      th4_match <- .03 + noise_factor[condition] * .2
      chrom_guides <- which(df$chr %in% chrom)
      x <- df$start_loc[chrom_guides]
      y <- guide_df[chrom_guides, condition]
      y <- y[order(x)]
      x <- x[order(x)]
      
      # Takes a running mean across genomic coordinates, with windows
      # sized relative to the gene density for the chromosome
      rollwindow <- chr_features[chrom, "roll_window"]
      rmean <- rep(NA, (length(x) - rollwindow))
      for (k in 1:length(rmean)) {
        y_rw <- y[k:(rollwindow + k - 1)]
        y_q <- stats::quantile(y_rw, probs = c(.05, .95), na.rm = TRUE)
        y_rw <- y_rw[y_rw > y_q[[1]] & y_rw < y_q[[2]]]
        rmean[k] <- mean(y_rw, na.rm = TRUE)
      }
      
      # 2. d running mean
      d_rmean <- rmean[1:(length(rmean)-rollwindow)] - rmean[(1+rollwindow):length(rmean)]
      
      # 3.2 maxima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "max", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "max", pi = y)
        }
      }
      x_maxima <- x
      
      # 3.1 minima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "min", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "min", pi = y)
        }
      }
      x_minima <- x
      
      # Identifies all shifted stretches
      chr_shifted_genes <- define_fragments(chr_shifted_genes, 
                                            x_max = x_maxima, x_min = x_minima,
                                            th3 = th3_match, th4 = th4_match,
                                            chromOI = chrom, condition = condition, 
                                            x_ref = rmean, pi = y,
                                            b = rollwindow, chrAnno = names(y))
    }
  }
  
  # Applies chromosomal correction to shifted genes if any exist
  if (length(unlist(chr_shifted_genes)) > 0) {
    for (condition in colnames(guide_df)[2:dim(guide_df)]) {
      for (chrom in names(chr_shifted_genes[[condition]])) {
        if (length(unlist(chrom)) > 0) {
          for (gene in chr_shifted_genes[[condition]][[chrom]]) {
            goi <- chr_shifted_genes[[condition]][[chrom]][[gene]]
            goi <- which(guide_df$gene %in% goi) 
            guide_df[goi, condition] <- df[goi, condition] - mean(df[goi, condition], na.rm = TRUE)
          }
        }
      }
    } 
  }
  
  # Returns corrected scores
  return(guide_df)
}

#' Plots chromosomal QC. 
#' 
#' Plots qGI effects by chromosomal location.
#' 
#' @param df LFC dataframe passed into \code{score_drugs_vs_controls} with chromosomal location 
#'   specified in a column named "chr", start coordinates specified in "start_loc" and stop 
#'   coordinates specified in "stop_loc." 
#' @param guide_df Dataframe of guide-level qGI scores. 
#' @param output_folder Folder to output plots to. 
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
plot_chromosome <- function(df, guide_df, output_folder, plot_type = "png") {
  
  # Gets guide density for each relevant chromosome
  chrom <- names(sort(table(df$chr)))
  chr_features <- data.frame(matrix(nrow = length(chrom), ncol = 4))
  rownames(chr_features) <- chrom
  colnames(chr_features) <- c("length", "n_genes", "length_per", "roll_window")
  for (i in 1:length(chrom)) {
    chr_features[i,"length"] <- max(df$stop_loc[df$chr == chrom[i]]) - min(df$start_loc[df$chr == chrom[i]])
    chr_features[i,"n_genes"] <- length(which(df$chr == chrom[i]))
  }
  chr_features$length_per <- chr_features$n_genes / chr_features$length
  chr_features$roll_window <- ceiling(log2(chr_features$length_per / stats::median(chr_features$length_per) + 2) * 150)
  chr_features$labels <- chrom
  
  # Plots chromosomal features information
  file_name <- paste0("chromosome_genes.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "n_genes")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Number of genes") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  file_name <- paste0("chromosome_adjusted_genes.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "length_per")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Number of guides adjusted by chromosomal length") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  file_name <- paste0("chromosome_rolling_window.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "roll_window")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Width of rolling window") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  
  # Ignores chromosomes with too few genes for the rolling window
  to_keep <- rep(TRUE, nrow(chr_features))
  for (i in 1:nrow(chr_features)) {
    if (chr_features[i,"n_genes"] < chr_features[i,"roll_window"]) {
      to_keep[i] <- FALSE
    }
  }
  chr_features <- chr_features[to_keep,]
  
  # Identifies windows where chromosomal shifts occur
  noise_factor <- apply(guide_df[,2:ncol(guide_df)], 2, stats::sd, na.rm = TRUE)
  chr_shifted_genes <- list()
  for (i in 1:nrow(chr_features)) {
    chrom <- row.names(chr_features)[i]
    for (j in 2:ncol(guide_df)) {
      condition <- colnames(guide_df)[j]
      
      # ???
      th1_match <- .09 + noise_factor[condition] * .2
      th2_match <- .03 + noise_factor[condition] * .1
      th3_match <- .04 + noise_factor[condition] * .2
      th4_match <- .03 + noise_factor[condition] * .2
      chrom_guides <- which(df$chr %in% chrom)
      x <- df$start_loc[chrom_guides]
      y <- guide_df[chrom_guides, condition]
      y <- y[order(x)]
      x <- x[order(x)]
      
      # Takes a running mean across genomic coordinates, with windows
      # sized relative to the gene density for the chromosome
      rollwindow <- chr_features[chrom, "roll_window"]
      rmean <- rep(NA, (length(x) - rollwindow))
      for (k in 1:length(rmean)) {
        y_rw <- y[k:(rollwindow + k - 1)]
        y_q <- stats::quantile(y_rw, probs = c(.05, .95), na.rm = TRUE)
        y_rw <- y_rw[y_rw > y_q[[1]] & y_rw < y_q[[2]]]
        rmean[k] <- mean(y_rw, na.rm = TRUE)
      }
      
      # 2. d running mean
      d_rmean <- rmean[1:(length(rmean)-rollwindow)] - rmean[(1+rollwindow):length(rmean)]
      
      # 3.2 maxima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "max", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "max", pi = y)
        }
      }
      x_maxima <- x
      
      # 3.1 minima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "min", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "min", pi = y)
        }
      }
      x_minima <- x
      
      # Identifies all shifted stretches
      chr_shifted_genes <- define_fragments(chr_shifted_genes, 
                                            x_max = x_maxima, x_min = x_minima,
                                            th3 = th3_match, th4 = th4_match,
                                            chromOI = chrom, condition = condition, 
                                            x_ref = rmean, pi = y,
                                            b = rollwindow, chrAnno = names(y))
      
      # # Plots shifted stretches
      # x <- df$start_loc[chrom_guides]
      # y <- guide_df[chrom_guides, condition]
      # 
      # plot(y[order(x)], 
      #      ylim = range(c(y, rmean, d_rmean), na.rm = TRUE), 
      #      cex = 0.4, pch = 16,
      #      col = ifelse(names(y[order(x)]) %in% unlist(chr_shifted_genes[[condition]][[chrom]]), "#4daf4a", "grey"),
      #      main = condition, bty = "n", 
      #      xlab = paste("Guide RNA on", chrom), 
      #      ylab = "Guide RNA Pi-score")
      # abline(h = c(0, th1_match, -th1_match), lty=c(2,3,3), col = c("black", "blue", "blue"))
      # 
      # for (i in 1:length(chr_shifted_genes[[condition]][[chrom]])) {
      #   abline(h = mean(y[chr_shifted_genes[[condition]][[chrom]][[i]]], na.rm = TRUE), lty = 2, col="orange")
      # }
      # 
      # lines(rmean, lwd = 1)
      # lines(d_rmean, lwd = 1, col = "blue")
      # abline(v = c(x_minima, x_minima+rollwindow, x_maxima, x_maxima+rollwindow), lty = 3, col = "red")
      # 
      # plot(x, y, ylim = range(c(y, rmean, d_rmean), na.rm = TRUE), cex = 0.4, pch = 16,
      #      col = ifelse(names(y) %in% unlist(chr_shifted_genes[[condition]][[chrom]]), "#4daf4a", "grey"),
      #      main = condition, bty = "n", xlab = paste("Location of guide RNA on", chrom), ylab = "Guide RNA Pi-score")
      # abline(h = 0, lty = 2)
      # 
      # for (i in 1:length(chr_shifted_genes[[condition]][[chrom]])) {
      #   abline(h = mean(y[chr_shifted_genes[[condition]][[chrom]][[i]]], na.rm = TRUE), lty = 2, col = "orange")
      # }
      # par(par.shift)
      # dev.off()
    }
  }
  
  # Applies chromosomal correction to shifted genes if any exist
  if (length(unlist(chr_shifted_genes)) > 0) {
    for (condition in colnames(guide_df)[2:dim(guide_df)]) {
      for (chrom in names(chr_shifted_genes[[condition]])) {
        if (length(unlist(chrom)) > 0) {
          for (gene in chr_shifted_genes[[condition]][[chrom]]) {
            goi <- chr_shifted_genes[[condition]][[chrom]][[gene]]
            goi <- which(guide_df$gene %in% goi) 
            guide_df[goi, condition] <- df[goi, condition] - mean(df[goi, condition], na.rm = TRUE)
          }
        }
      }
    } 
  }
  
  # Returns corrected scores
  return(guide_df)
}

