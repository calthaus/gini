##########
# Gini coefficient
##########
gini <- function(N,p) {
	x<-c(0,N)
	y<-c(0,p)
	area <- 0
	for(i in 2:length(x)) {
		area <- area + x[i]*(cumsum(x[1:i]*y[1:i])[i]/sum(x*y)-cumsum(x[1:i-1]*y[1:i-1])[i-1]/sum(x*y))/2 + x[i]*cumsum(x[1:i-1]*y[1:i-1])[i-1]/sum(x*y)
	}
	return(1 - 2*area)
}

##########
# Gini coefficient
# n: proportion of people in each class
# c: partner change rate in each class
# p: prevalence in each class
# max: maximal number of partners
##########
gini.pois <- function(N,c,p,max=100) {
	pn <- array(0,max+1)
	pp <- array(0,max+1)
	for(i in 0:max) for(j in 1:length(N)) {
		pn[i+1] <- pn[i+1] + N[j]*dpois(i,c[j])
		pp[i+1] <- pp[i+1] + N[j]*dpois(i,c[j])*p[j]
	}
	pp <- pp/pn
	return(gini(pn,pp))
}