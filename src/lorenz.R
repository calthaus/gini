lorenz <- function(x_cts, y_bin, wt){
  ## Calculate Lorenz curve and Gini coefficient
  ## -------------------------------------------
  ## x_cts      count variable to be used for x-axis
  ## y_bin      binary variable to be used for y-axis
  ## wt         sampling weights
  ## -------------------------------------------
  ## Requirements: 
  ## - libraries: plyr
  ## - utility functions: gini
  ## -------------------------------------------
  ## S. Gsteiger, Mar-2014
  ## -------------------------------------------
  
  # calculate proportions in x-axis groups
  w_cts <- ddply(data.frame(wt, x_cts), .(x_cts), .fun = function(v){c('w_sum' = sum(v$wt, na.rm = TRUE))})
  w_cts <- within(w_cts, {x_prop = w_sum / sum(w_sum)})
  
  # calculate the prevalences within each group 
  y_prev <- ddply(data.frame(x_cts, y_bin, wt), .(x_cts), .fun = function(v) {c('y_prev' = sum(v$wt * v$y_bin, na.rm = TRUE) / sum(v$wt, na.rm = TRUE))})


  # combine and calculate the cummulative proportions to obtain the Lorenz curve values
  dlor <- join(w_cts, y_prev, by = 'x_cts')[, c('x_prop', 'y_prev')]
  dlor <- within(rbind(c(0, 0), dlor), {
    x_cum = cumsum(x_prop)
    y_cts = y_prev * x_prop
    y_cum = cumsum(y_cts) / sum(y_cts)
    x_vals = c(NA, sort(unique(x_cts))) # (added only later on; rm if issues)
  })
  
  out <- list(lorenz = dlor, 
              gini = gini(dlor$x_prop, dlor$y_prev)
              )
  return(out)
}
  

# FCT BELOW IS ORIGINAL VERSION (x_val not given as output)
# 
# lorenz <- function(x_cts, y_bin, wt){
#   ## Calculate Lorenz curve and Gini coefficient
#   ## -------------------------------------------
#   ## x_cts      count variable to be used for x-axis
#   ## y_bin      binary variable to be used for y-axis
#   ## wt         sampling weights
#   ## -------------------------------------------
#   ## Requirements: 
#   ## - libraries: plyr
#   ## - utility functions: gini
#   ## -------------------------------------------
#   ## S. Gsteiger, Mar-2014
#   ## -------------------------------------------
#   
#   # calculate proportions in x-axis groups
#   w_cts <- ddply(data.frame(wt, x_cts), .(x_cts), .fun = function(v){c('w_sum' = sum(v$wt, na.rm = TRUE))})
#   w_cts <- within(w_cts, {x_prop = w_sum / sum(w_sum)})
#   
#   # calculate the prevalences within each group 
#   y_prev <- ddply(data.frame(x_cts, y_bin, wt), .(x_cts), .fun = function(v) {c('y_prev' = sum(v$wt * v$y_bin, na.rm = TRUE) / sum(v$wt, na.rm = TRUE))})
#   
#   
#   # combine and calculate the cummulative proportions to obtain the Lorenz curve values
#   dlor <- join(w_cts, y_prev, by = 'x_cts')[, c('x_prop', 'y_prev')]
#   dlor <- within(rbind(c(0, 0), dlor), {
#     x_cum = cumsum(x_prop)
#     y_cts = y_prev * x_prop
#     y_cum = cumsum(y_cts) / sum(y_cts)
#   })
#   
#   
#   out <- list(lorenz = dlor, 
#               gini = gini(dlor$x_prop, dlor$y_prev)
#   )
#   return(out)
# }

