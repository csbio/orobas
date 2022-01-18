######
# qGI UTILITY FUNCTIONS
######

### All utility functions directly ported from the qGI scoring pipeline,
### written by Maximilian Billmann and Henry Ward, are located here.

# Computes control screen effects against all other control screens as if they
# were condition screens, setting the matched control contribution from the 
# closest control screen to 0
compute_control_effects <- function(control_df, control_names, method = "linear",
                                    loess = TRUE, ma_transform = TRUE) {
  weights <- compute_control_weights(control_df, control_df, control_names, control_names,
                                     control_names, method = method, matched_fraction = 0)
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

# Weights control screens for \code{score_drugs_vs_controls} depending on the 
# mode specified by the user.
# matched_controls: if NULL, no matching is performed, otherwise takes a vector
#   of control screen names corresponding to all condition names
# method: "linear" or "exp"
compute_control_weights <- function(control_df, condition_df, control_names, condition_names,
                                    matched_controls, method = "exp", matched_fraction = 0.75) {
  
  # For a 100% matched control weighting, the default type is set to linear
  if (matched_fraction >= 1) {
    cat(paste("Defaulting to linear weighting because matched_fraction >= 1\n"))
  }
  
  # If there is only 1 control, all weights are equal to 1
  weights <- matrix(1, nrow = length(condition_names), ncol = length(control_names))
  rownames(weights) <- condition_names
  colnames(weights) <- control_names
  if (ncol(weights) > 1) {
    for (i in 1:nrow(weights)) {
      
      # Gets the current condition, its matched control and all other controls
      condition <- rownames(weights)[i]
      matched_control <- matched_controls[i]
      other_controls <- control_names[!(control_names == matched_control)]
      
      # Sets weight for matched control
      weights[i, matched_control] <- matched_fraction
      
      # If the method is "linear" all non-matched controls are given an equal weight
      if (method == "linear") {
        weights[i,] <- (1 - matched_fraction) / (ncol(weights) - 1)
        weights[i, matched_control] <- matched_fraction
      } else if (method == "exp") {
        
        # For exponential ranking, we determine weights of non-matched screens depending
        # on the rank of their correlation to the condition according to an exponentially
        # decreasing function
        temp <- cbind(condition_df[,condition], control_df[,other_controls])
        colnames(temp) <- c(condition, other_controls)
        cor <- cor(temp, use = "complete.obs")
        weight_rank <- rank(cor[1,2:ncol(cor)])
        if (length(weight_rank) == 1) {
          names(weight_rank) <- colnames(cor)[2]
        }
        weight_rank <- exp(-weight_rank)
        weight_rank <- weight_rank * ((1 - matched_fraction) / sum(weight_rank))
        weight_rank[[matched_control]] <- matched_fraction
        weight_rank <- weight_rank[colnames(weights)]
        if (!isTRUE(all.equal(sum(weight_rank), 1))) {
          cat(paste("WARNING: weights for screen,", condition, "do not sum to 1\n"))
        }
        weights[i,] <- weight_rank
      }
    }
  }
  return(weights)
}

# Find wildtype (control) weights
#
#' @param x: genetic interaction data array
#' @param y: wildtype (control) foldchange data
#' @param fraction_floor: ?
#
wt_weightID_func <- function(x, y, fraction_floor = 0.05) {
  
  x <- apply((x)^2, 2:3, sum, na.rm = TRUE) #compute square sum of guide-level gi scores for each query-wt pair
  x[x == 0] <- NA #self-gi gets 0 sum, because na.rm = TRUE, set to NA
  y <- apply(y, 2, var, na.rm = TRUE)
  
  for (i in 1:dim(x)[1] ) {
    x[i,] <- x[i,] / y #if wt screens have higher log2FC variance, they tend to also have more gi signal (not linear)
  }
  
  x <- ( x - apply(x, 1, max, na.rm = TRUE) ) * -1 #highest value for smallest sum
  if (fraction_floor > 0) {
    x <- x + apply(x, 1, max, na.rm = TRUE) / (1 / fraction_floor - 1) #add fraction (default 5%) of max weight
  }
  
  x <- x / apply(x, 1, sum, na.rm = TRUE) #scale so sum is 1
  return(x)
}

# Fits a loess curve to predict y given x
loess_MA <- function(x, y, sp = 0.4, dg = 2, binSize = 100, ma_transform = TRUE) {
  #this concept is based on pythagoras and cancels out sqrt, square and factor 2
  #it also ignores the factor sqrt(2) as factor between y and x vs distance of x,y from diagonal x = y
  gi <- NULL
  expected <- NULL
  if(all(x == y, na.rm = T)) { #if e.g. wt scored against itself
    gi <- rep(NA, length(x))
  }
  else {
    if (ma_transform) {
      m <- y - x
      a <- y + x
      A <- (a - stats::median(a, na.rm = T)) / stats::mad(a, na.rm = T)  # scale to generate bins along m
      B <- seq(floor(min(A, na.rm = T)), ceiling(max(A, na.rm = T)), .1) # define bins
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
      model <- stats::loess(m[b][I] ~ a[b][I], span = sp, degree = dg) # train model on m ~ a (approx. y ~ x)
      expected <- stats::predict(model, a) #predict expected m ~ a
      gi <- m - expected
    } else {
      m <- y
      a <- x
      A <- (a - stats::median(a, na.rm = T)) / stats::mad(a, na.rm = T)  # scale to generate bins along m
      B <- seq(floor(min(A, na.rm = T)), ceiling(max(A, na.rm = T)), .1) # define bins
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
      model <- stats::loess(m[b][I] ~ a[b][I], span = sp, degree = dg) # train model on m ~ a (approx. y ~ x)
      expected <- stats::predict(model, a) #predict expected m ~ a
      gi <- m - expected
    }
  }
  result <- list()
  result[["residual"]] <- gi
  result[["predicted"]] <- expected
  return(result)
}

# Chromosomal
fragment_transition <- function(x, x_ref, th1, th2, tb,
                                maxORmin = "max", pi = NA) {
  
  if (mean(pi, na.rm = TRUE) * -2 < max(x) | mean(pi, na.rm = TRUE) * 2 > min(x)) {
    if (missing(th1)) {
      th1 <- max(abs(x_ref))/2 #stringent threshold if not specifically stated
    }
    if (missing(th2)) {
      th2 <- max(abs(x_ref))/4
    }
    
    # While large SD also indicates a shift, strongly fluctuating screens
    # should not be corrected on every chromosome
    #th1 <- th1 + sd(pi, na.rm = TRUE)/5 # ~0.07 per chromosome in median. 
    #th2 <- th2 + sd(pi, na.rm = TRUE)/10
    
    if (maxORmin == "max") {
      if (any(x > th1)) {
        if (missing(tb)) {
          x_max <- x > th1 #detects peaks
          tb <- min(which(x_max), na.rm = TRUE)
          tb <- c(tb, max(which(x_max), na.rm = TRUE))
          if (all(x[tb[1] : tb[2]] > th2, na.rm = TRUE)) {
            localMax <- which(x == max(x[tb[1] : tb[2]], na.rm = TRUE))
            return(localMax)
          } else {
            return(tb) # has twice the length as localMax
          }
        } else if (length(tb) %in% seq(2,50,4)) {
          a <- length(tb) - 1
          b <- length(tb)
          x_max <- x[tb[a] : tb[b]] < th2 #detects breaks between peaks
          tb <- c(tb, min(which(x_max), na.rm = TRUE) + tb[a]) #extend input tb
          tb <- c(tb, max(which(x_max), na.rm = TRUE) + tb[a])
          if (all(x[tb[b+1] : tb[b+2]] < th1, na.rm = TRUE)) {
            tb <- sort(tb)
            localMax <- which(x == max(x[tb[1] : tb[2]], na.rm = TRUE))
            for (i in 2:(length(tb)/2) - 1) {
              localMax <- c(localMax, which(x == max(x[tb[i * 2 + 1] : tb[i * 2 + 2]], na.rm = TRUE)))
            }
            return(localMax)
          } else {
            return(tb) #has twice the length as localMax
          }
        } else if (length(tb) %in% seq(4,48,4)) {
          a <- length(tb) - 1
          b <- length(tb)
          x_max <- x[tb[a] : tb[b]] > th1 #detects peaks
          tb <- c(tb, min(which(x_max), na.rm = TRUE) + tb[a]) #extend input tb
          tb <- c(tb, max(which(x_max), na.rm = TRUE) + tb[a])
          if (all(x[tb[b+1] : tb[b+2]] > th2, na.rm = TRUE)) {
            tb <- sort(tb)
            localMax <- which(x == max(x[tb[1] : tb[2]], na.rm = TRUE))
            for (i in 2:(length(tb)/2) - 1) {
              localMax <- c(localMax, which(x == max(x[tb[i * 2 + 1] : tb[i * 2 + 2]], na.rm = TRUE)))
            }
            return(localMax)
          } else {
            return(tb) #has twice the length as localMax
          }
        }
      }
    } else if (maxORmin == "min") {
      if (any(x < -th1)) {
        if (missing(tb)) {
          x_max <- x < -th1 #detects peaks
          tb <- min(which(x_max), na.rm = TRUE)
          tb <- c(tb, max(which(x_max), na.rm = TRUE))
          if (all(x[tb[1] : tb[2]] < -th2, na.rm = TRUE)) {
            localMax <- which(x == min(x[tb[1] : tb[2]], na.rm = TRUE))
            return(localMax)
          } else {
            return(tb) #has twice the length as localMax
          }
        }
        else if (length(tb) %in% seq(2,50,4)) {
          a <- length(tb) - 1
          b <- length(tb)
          x_max <- x[tb[a] : tb[b]] > -th2 #detects breaks between peaks
          tb <- c(tb, min(which(x_max), na.rm = TRUE) + tb[a]) #extend input tb
          tb <- c(tb, max(which(x_max), na.rm = TRUE) + tb[a])
          if (all(x[tb[b+1] : tb[b+2]] > -th1, na.rm = TRUE)) {
            tb <- sort(tb)
            localMax <- which(x == min(x[tb[1] : tb[2]], na.rm = TRUE))
            for (i in 2:(length(tb)/2) - 1) {
              localMax <- c(localMax, which(x == min(x[tb[i * 2 + 1] : tb[i * 2 + 2]], na.rm = TRUE)))
            }
            return(localMax)
          } else {
            return(tb) #has twice the length as localMax
          }
        } else if (length(tb) %in% seq(4,48,4)) {
          a <- length(tb) - 1
          b <- length(tb)
          x_max <- x[tb[a] : tb[b]] < -th1 #detects peaks
          tb <- c(tb, min(which(x_max), na.rm = TRUE) + tb[a]) #extend input tb
          tb <- c(tb, max(which(x_max), na.rm = TRUE) + tb[a])
          if (all(x[tb[b+1] : tb[b+2]] < -th2, na.rm = TRUE)) {
            tb <- sort(tb)
            localMax <- which(x == min(x[tb[1] : tb[2]], na.rm = TRUE))
            for (i in 2:(length(tb)/2) - 1) {
              localMax <- c(localMax, which(x == min(x[tb[i * 2 + 1] : tb[i * 2 + 2]], na.rm = TRUE)))
            }
            return(localMax)
          } else {
            return(tb) #has twice the length as localMax
          }
        }
      }
    } else { print("what do you want?") }
  }
}

define_fragments <- function(chrShift_genes_temp, x_max = x_maxima, x_min = x_minima,
                             th3, th4, chromOI = chrom, Qoi = qoi, x_ref = rmean,
                             pi = y, b = rollwindow, chrAnno = names(y)) {
  
  sIndex <- 1
  
  # identify all stretches in middle of chromosome
  if (length(x_min) > 0 & length(x_max) > 0) { #run only if shifts detected
    for (i in 1:length(x_min)) {
      a <- x_max - x_min[i]
      if (any(a < 0)) {
        a1 <- x_max[which(a == min(abs(a[a < 0]), na.rm = TRUE) * -1)] #previous max
        a1 <- a1 + 2 * b #add b, because max to min on drmean can only be negative rmean and a max shows were it starts with a shift of length(b)
        # add second b, because indices on pi (not x_ref) are of interest
        if (mean(pi[x_min[i]:a1], na.rm = TRUE) < -th3) {
          if (sIndex == 1) {
            chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
          }
          x <- x_min[i] : a1#save negative shift indices x_min[i]:a1
          chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x]) #get gene names for saved guide indices
          sIndex <- sIndex + 1
        }
      }
      
      if (any(a > 0)) {
        a2 <- x_max[which(a == min(abs(a[a > 0]), na.rm = TRUE))]
        if (mean(pi[x_min[i]:a2 + b], na.rm = TRUE) > th3) {
          if (sIndex == 1) {
            chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
          }
          x <- x_min[i] : a2 + b #save positive shift indices x_min[i]:a2 + b
          chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x]) #get gene names for saved guide indices
          sIndex <- sIndex + 1
        }
      }
    }
  }
  
  x_m <- c(x_min, x_max)
  
  if (length(x_m) > 0) {
    if (min(x_m) %in% x_min) { #define as first stretch
      if (mean(pi[1:(min(x_m)+b)], na.rm = TRUE) < -th3) {
        if (sIndex == 1) {
          chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
        }
        x <- 1:(min(x_m) + b)
        chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x])
        sIndex <- sIndex + 1
      }
    } else if (min(x_m) %in% x_max) {
      if (mean(pi[1:(min(x_m)+b)], na.rm = TRUE) > th3) {
        if (sIndex == 1) {
          chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
        }
        x <- 1:(min(x_m) + b)
        chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x])
        sIndex <- sIndex + 1
      }
    }
    
    if (max(x_m) %in% x_min) { #define as last stretch
      if (mean(pi[b + max(x_m):length(x_ref)], na.rm = TRUE) > th3) {
        if (sIndex == 1) {
          chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
        }
        x <- b + max(x_m):length(x_ref)
        chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x])
        sIndex <- sIndex + 1
      }
    } else if (max(x_m) %in% x_max) {
      if (mean(pi[b + max(x_m):length(x_ref)], na.rm = TRUE) < -th3) {
        if (sIndex == 1) {
          chrShift_genes_temp[[Qoi]][[chromOI]] <- list()
        }
        x <- b + max(x_m):length(x_ref)
        chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno[x])
        sIndex <- sIndex + 1
      }
    }
  }
  if (length(x_m) == 0 & abs(mean(pi, na.rm = TRUE)) > th4) {
    chrShift_genes_temp[[Qoi]][[chromOI]][[sIndex]] <- unique(chrAnno)
  }
  return(chrShift_genes_temp) #return list
}