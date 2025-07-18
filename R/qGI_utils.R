######
# qGI UTILITY FUNCTIONS
######

### All utility functions are ported from the qGI scoring pipeline.
### Billmann, Maximilian, et al. "Quantitative analysis of genetic interactions in human cells from genome-wide CRISPR-Cas9 screens." bioRxiv (2025): 2025-06.


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
