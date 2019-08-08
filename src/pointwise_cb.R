pointwise_cb <- function(dsim, 
                         probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                         y_nam = c('CI_lo', 'Q1', 'Med', 'Q3', 'CI_up'),
                         x_grid = seq(0, 1, 0.01)){
  ## take a data frame with simulated Lorenz curves and construct a point-wise
  ##  confidence band (via linear interpolation on a grid of x-values)
  ## -----------------------------------------------------------------
  ## dsim     data.frame with simulated curves, must contain 'iter', 'crisk', and 'cpinfect'
  ## probs    probability thresholds for the quantiles taken to construct the CIs
  
  get_lorenz <- function(x) {approx(round(x[,'crisk'], 5), x[, 'cpinfect'], xout = x_grid)$y}
  y_all <- matrix(unlist(lapply(split(dsim, dsim$iter), FUN = 'get_lorenz')),
                  nrow = length(x_grid))
  
  dcb <- data.frame(x = x_grid,
                    t(apply(y_all, MAR = 1, FUN = 'quantile', probs = probs)))
  names(dcb)[-1] <- y_nam
  
  return(dcb)
}
