lorenz_boot <- function(data, x_cts, y_bin, wt, 
                        plot_lorenz = TRUE, x_lab = NULL, y_lab = NULL,
                        R
                        ) {
  ## Calc Lorenz curve and Gini coeff,
  ## and derive bootstrap CI for Gini
  ## -------------------------------------------
  ## data       
  ## x_cts      char string, count variable to be used for x-axis
  ## y_bin      char string, binary variable to be used for y-axis
  ## wt         char string, sampling weights
  ## plot_lorenz  logical, should Lorenz curve be plotted?
  ## x_lab      
  ## y_lab
  ## R          integer, nb of bootstrap samples
  ## -------------------------------------------
  ## Requirements: 
  ## - libraries: plyr, boot
  ## - utility functions: lorenz, gini
  ## -------------------------------------------
  ## S. Gsteiger, Mar-2014
  ## -------------------------------------------

  # Calculate the Gini coefficient and plot the Lorenz curve. 
  lng <- lorenz(x_cts = data[, x_cts], y_bin = data[, y_bin], wt = data[, wt])
  
  p0 <- qplot(x_cum, y_cum, data = lng$lorenz, geom = 'line') + geom_abline(aes(intercept = 0, slope = 1), col = 'grey')   
  p0 <- p0 + xlab(x_lab) + ylab(y_lab)
  if(plot_lorenz){
    print(p0)
  } 
  
  # Bootstrap CI.
  boot_gini <- function(data, indices){
    d <- data[indices,]
    lorenz(x_cts = d[, x_cts], y_bin = d[, y_bin], wt = d[, wt])$gini
  }
  gini_boot <- boot(data = data, statistic = boot_gini, R = R)
  gini_BCI <- boot.ci(gini_boot, type = 'perc')
    
  out <- list(gini = round(lng$gini, 3), 
              lorenz = lng$lorenz, 
              plot = p0, 
              gini_boot = gini_boot, 
              gini_BCI = gini_BCI, 
              gini_bCI = round(gini_BCI$percent[, 4:5], 3)
  )
  return(out)
}
