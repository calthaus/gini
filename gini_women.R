# Load libraries
library(foreign)
library(reshape)
library(ggplot2)
library(boot)
library(plotrix)
library(grid)
library(gridExtra)
library(deSolve)
library(plyr)
library(knitr)
library(colorspace)

if (!file.exists('out')){
  dir.create('out')
}

for (s in dir('src')){
  source(file.path('src', s))
}
set.seed(3297348)

# Load Natsal-3 data (available at http://data-archive.ac.uk)
dnat3 <- read.dta('data/natsal3.dta')

# Analyse data set
age_lo <- 16
age_up <- 44

durine <- within(subset(dnat3, urintested == 'yes' & dage >= age_lo & dage <= age_up & rsex == 'Female', 
                        select = c('psu', 'dage', 'rsex', 
                                   'ct_posconfirmed', 'mg_posconfirmed',
                                   'HPV_6', 'HPV_11', 'HPV_16', 'HPV_18', 
                                   'urine_wt',
                                   'hetnonew'
                                   )), 
{
  ct = as.integer(ifelse(ct_posconfirmed == 'positive', 1, 0))
  mg = as.integer(ifelse(mg_posconfirmed == 'positive', 1, 0))
  hpv_6 = as.integer(ifelse(HPV_6 == 'Positive', 1, 0))
  hpv_11 = as.integer(ifelse(HPV_11 == 'Positive', 1, 0))
  hpv_16 = as.integer(ifelse(HPV_16 == 'Positive', 1, 0))
  hpv_18 = as.integer(ifelse(HPV_18 == 'Positive', 1, 0))
  uwt_missing = is.na(urine_wt)
  uid = 1:length(psu)
})

dct <- within(subset(durine, ct_posconfirmed %in% c('negative', 'positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  ct_bin = as.integer(ifelse(ct_posconfirmed == 'positive', 1, 0))
})
dhpv_6 <- within(subset(durine, HPV_6 %in% c('Negative', 'Positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  hpv_bin = as.integer(ifelse(HPV_6 == 'Positive', 1, 0))
})
dhpv_11 <- within(subset(durine, HPV_11 %in% c('Negative', 'Positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  hpv_bin = as.integer(ifelse(HPV_11 == 'Positive', 1, 0))
})
dhpv_16 <- within(subset(durine, HPV_16 %in% c('Negative', 'Positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  hpv_bin = as.integer(ifelse(HPV_16 == 'Positive', 1, 0))
})
dhpv_18 <- within(subset(durine, HPV_18 %in% c('Negative', 'Positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  hpv_bin = as.integer(ifelse(HPV_18 == 'Positive', 1, 0))
})
dmg <- within(subset(durine, mg_posconfirmed %in% c('negative', 'positive') & !(hetnonew %in% c(-1, 995, 999)) ), {
  mg_bin = as.integer(ifelse(mg_posconfirmed == 'positive', 1, 0))
})

# Calculate prevalence
dprev <- as.data.frame(matrix(NA,6,4))
names(dprev) <- c("STI","prev","prev_lo","prev_up")
temp <- dct$ct
temp <- binom.test(sum(temp),length(temp))
dprev[1,] <- c("CT", temp$estimate, temp$conf.int[1], temp$conf.int[2])
temp <- dmg$mg
temp <- binom.test(sum(temp),length(temp))
dprev[2,] <- c("MG", temp$estimate, temp$conf.int[1], temp$conf.int[2])
temp <- dhpv_6$hpv_6
temp <- binom.test(sum(temp),length(temp))
dprev[3,] <- c("HPV-6", temp$estimate, temp$conf.int[1], temp$conf.int[2])
temp <- dhpv_11$hpv_11
temp <- binom.test(sum(temp),length(temp))
dprev[4,] <- c("HPV-11", temp$estimate, temp$conf.int[1], temp$conf.int[2])
temp <- dhpv_16$hpv_16
temp <- binom.test(sum(temp),length(temp))
dprev[5,] <- c("HPV-16", temp$estimate, temp$conf.int[1], temp$conf.int[2])
temp <- dhpv_18$hpv_18
temp <- binom.test(sum(temp),length(temp))
dprev[6,] <- c("HPV-18", temp$estimate, temp$conf.int[1], temp$conf.int[2])
save(dprev, file = file.path('out', 'dprev.RData'))

# Calculate Lorenz curves
ct <- lorenz_boot(data = dct, x_cts = 'hetnonew', y_bin = 'ct_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)
hpv_6 <- lorenz_boot(data = dhpv_6, x_cts = 'hetnonew', y_bin = 'hpv_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)
hpv_11 <- lorenz_boot(data = dhpv_11, x_cts = 'hetnonew', y_bin = 'hpv_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)
hpv_16 <- lorenz_boot(data = dhpv_16, x_cts = 'hetnonew', y_bin = 'hpv_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)
hpv_18 <- lorenz_boot(data = dhpv_18, x_cts = 'hetnonew', y_bin = 'hpv_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)
mg <- lorenz_boot(data = dmg, x_cts = 'hetnonew', y_bin = 'mg_bin', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)

# Construct the data sets used for the plotting the Lorenz curve; save the objects, which are used in separate chunks
dfig1 <- rbind(data.frame(STI = 'CT', ct$lorenz),
               data.frame(STI = 'HPV 6', hpv_6$lorenz),
               data.frame(STI = 'HPV 11', hpv_11$lorenz),
               data.frame(STI = 'HPV 16', hpv_16$lorenz),
               data.frame(STI = 'HPV 18', hpv_18$lorenz),
               data.frame(STI = 'MG', mg$lorenz))
dfig1 <- within(dfig1, {
  STI = factor(STI, levels = c('CT', 'MG', 'HPV 6', 'HPV 11', 'HPV 16', 'HPV 18'), ordered = TRUE)
})

save(dfig1, file = file.path('out', 'dfig1.RData'))

dfig1_gini <- rbind(data.frame(STI = 'CT', Gini = ct$gini, CI_lo = ct$gini_bCI[1], CI_up = ct$gini_bCI[2]),
                    data.frame(STI = 'MG', Gini = mg$gini, CI_lo = mg$gini_bCI[1], CI_up = mg$gini_bCI[2]),
                    data.frame(STI = 'HPV 6', Gini = hpv_6$gini, CI_lo = hpv_6$gini_bCI[1], CI_up = hpv_6$gini_bCI[2]),
                    data.frame(STI = 'HPV 11', Gini = hpv_11$gini, CI_lo = hpv_11$gini_bCI[1], CI_up = hpv_11$gini_bCI[2]),
                    data.frame(STI = 'HPV 16', Gini = hpv_16$gini, CI_lo = hpv_16$gini_bCI[1], CI_up = hpv_16$gini_bCI[2]),
                    data.frame(STI = 'HPV 18', Gini = hpv_18$gini, CI_lo = hpv_18$gini_bCI[1], CI_up = hpv_18$gini_bCI[2]))

save(dfig1_gini, file = file.path('out', 'dfig1_gini.RData'))

# For chlamydia, redo derive the bootstrap CI for the Lorenz curve (the function used above only saves the Gini coefs but does not the bootstrapped curves)
dim(dct)
dct_boot <- subset(dct, !uwt_missing)
dim(dct_boot)

n_boot <- 1000
get_wprev <- function(x) {sum(x$urine_wt * x$ct_bin) / sum(x$urine_wt)}
dlall <- matrix(NA, nrow = n_boot * length(unique(dct$hetnonew)), ncol = 6)
last_row <- 0                
for(i in 1:n_boot){
  # cat(i, ' ')
  ids_i <- sample(nrow(dct_boot), replace = TRUE)
  data_i <- dct_boot[ids_i, ]
  
  wh_i <- weighted.hist(x = data_i$hetnonew, w = data_i$urine_wt, breaks = seq(-0.5, 150.5, 1), plot = FALSE)
  risk_i <- data.frame(nb_partners = wh_i$mids,
                       risk = wh_i$counts / sum(wh_i$counts))
  
  prev_i <- rename(ddply(data_i, .(hetnonew), .fun = 'get_wprev'), 
                   replace = c('hetnonew' = 'nb_partners', 'get_wprev' = 'prev'))
  
  drp_i <- merge(prev_i, risk_i, by= 'nb_partners')
  
  dlorenz_i <- within(rbind(c(0, 0), drp_i[, c('risk', 'prev')]), {
    crisk = cumsum(risk)
    pinfect = prev * risk
    cpinfect = cumsum(pinfect) / sum(pinfect)    
  }) 
  dlall[last_row + 1:nrow(dlorenz_i), 1] <- i
  dlall[last_row + 1:nrow(dlorenz_i), 2:6] <- as.matrix(dlorenz_i)
  
  last_row <- last_row + nrow(dlorenz_i)
}
colnames(dlall) <- c('iter', names(dlorenz_i))
dlall <- as.data.frame(dlall)
dcb <- pointwise_cb(dlall)
save(dcb, file = file.path('out', 'dcb.RData'))
save(ct, file = file.path('out', 'ct.RData'))

# Calculate the Lorenz curve and Gini coeffients for chlamydia using Natsal-2
dnat2 <- read.dta('data/natsal2.dta')

dct_nat2 <- subset(dnat2,
                   dage >= age_lo & dage <= age_up & rsex == 1 & !(hetnonew %in% c(-1, 995, 999)) & c_result %in% 0:1, 
                   select = c('dage', 'hetnonew', 'c_result', 'urine_wt'))

ct_nat2 <- lorenz_boot(data = dct_nat2, x_cts = 'hetnonew', y_bin = 'c_result', wt = 'urine_wt', plot_lorenz = FALSE, R = 1000)

dprev_nat2 <- as.data.frame(matrix(NA,1,4))
names(dprev_nat2) <- c("STI","prev","prev_lo","prev_up")
temp <- dct_nat2$c_result
temp <- binom.test(sum(temp),length(temp))
dprev_nat2[1,] <- c("CT", temp$estimate, temp$conf.int[1], temp$conf.int[2])

save(dprev_nat2, file = file.path('out', 'dprev_nat2.RData'))
save(ct_nat2, file = file.path('out', 'ct_nat2.RData'))

load(file = file.path('out', 'dcb.RData'))
load(file = file.path('out', 'ct.RData'))
load(file = file.path('out', 'ct_nat2.RData'))

dct_nat2_nat3 <- rbind(data.frame(Survey = 'Natsal 3', Gini = ct$gini, CI_lo = ct$gini_bCI[1], CI_up = ct$gini_bCI[2]),
                    data.frame(Survey = 'Natsal 2', Gini = ct_nat2$gini, CI_lo = ct_nat2$gini_bCI[1], CI_up = ct_nat2$gini_bCI[2]))
dct_nat2_nat3[, -1] <- round(dct_nat2_nat3[, -1], 2)

par(mfrow=c(3,1))
load(file = file.path('out', 'dfig1.RData'))
STI <- levels(dfig1$STI)
cols <- rainbow_hcl(length(STI))
# Fig. 1A
plot(c(0,1),c(0,1),ty="l",lty=3,xlab="New heterosex. partners last year (cumul. proportion)",ylab="STI infections (cumul. proportion)",frame=FALSE)
for(i in 1:length(STI)) {
    temp <- dfig1[dfig1$STI==STI[i],]
    lines(temp$x_cum,temp$y_cum,lty=i,col=cols[i])
}
legend("topleft",inset=0.05,legend=c("CT","MG","HPV 6","HPV 11","HPV 16","HPV 18"),col=cols,lty=1:length(STI),bty="n")
mtext("A",side=3,adj=0,cex=1,font=2)
# Fig. 1B
plot(c(0,1),c(0,1),ty="l",lty=3,xlab="New heterosex. partners last year (cumul. proportion)",ylab="CT infections (cumul. proportion)",frame=FALSE)
polygon(x = c(dcb$x, rev(dcb$x)), y = c(dcb$CI_lo, rev(dcb$CI_up)), col = rgb(0, 0, 1, alpha=0.2), border = NA)
polygon(x = c(dcb$x, rev(dcb$x)), y = c(dcb$Q1, rev(dcb$Q3)), col = rgb(0, 0, 1, alpha=0.4), border = NA)
lines(ct$lorenz$x_cum,ct$lorenz$y_cum, col = rgb(0, 0, 1))
legend("topleft",inset=0.05,legend="CT",col=rgb(0, 0, 1),lty=1,bty="n")
mtext("B",side=3,adj=0,cex=1,font=2)
# Fig. 1C
cols <- rainbow_hcl(2)
plot(c(0,1),c(0,1),ty="l",lty=3,xlab="New heterosex. partners last year (cumul. proportion)",ylab="CT infections (cumul. proportion)",frame=FALSE)
temp <- data.frame(Survey = 'Natsal 3', ct$lorenz)
lines(temp$x_cum,temp$y_cum, col = cols[1], lty=1)
temp <- data.frame(Survey = 'Natsal 2', ct_nat2$lorenz)
lines(temp$x_cum,temp$y_cum, col = cols[2], lty=2)
legend("topleft",inset=0.05,legend=c("Natsal-3","Natsal-2"),col=cols,lty=1:2,bty="n")
mtext("C",side=3,adj=0,cex=1,font=2)

load(file = file.path('out', 'dfig1_gini.RData'))
dfig1_gini[,2:4] <- round(dfig1_gini[,2:4],2)
levels(dfig1_gini$STI) <- c("Chlamydia trachomatis","Mycoplasma genitalium","HPV 6","HPV 11","HPV 16","HPV 18")
dfig1_gini <- within(dfig1_gini, {
  CI = paste(format(round(CI_lo,2),nsmall=2), ' - ', format(round(CI_up,2),nsmall=2), sep='')
})
names(dfig1_gini) <- c("Infection","Gini coefficient","CI_lo","CI_up","95% confidence interval (CI)")
kable(dfig1_gini[, c("Infection","Gini coefficient","95% confidence interval (CI)")],
    align = "lcc") 

# Transmission model
# Import sexual behaviour data
data <- read.csv("data/sexual_activity.csv",header=TRUE)
data <- na.omit(data)

# Maximal number of partners per year
# max <- max(data$activity)
# Alternative definition of max
max <- length(data$activity) - 1
# Partner change rate
c <- data$activity
# Proportion of people in each host class (weighted or unweighted)
N <- data$w.proportion
# Average partner change rate in the population
cN <- sum(c*N)

# Mixing matrix
mixing <- function(epsilon) {
	f <- matrix(nrow=(max+1),ncol=(max+1))
	for(i in 1:(max+1))
		for(j in 1:(max+1)) {
			if(i == j) f[i,j] <- epsilon
			else f[i,j] <- 0
		}
	rho <- matrix(nrow=(max+1),ncol=(max+1))
	for(i in 1:(max+1)) for(j in 1:(max+1)) rho[i,j] <- f[i,j] + (1.-epsilon)*c[j]*N[j]/cN
	rho
}

# Sex acts matrix
sexacts <- function(c1,c2,c3,epsilon) {
	# The average number of sex acts per partner per year
	s <- (c1 + c2*c^c3)/c/28*365
	s[1] <- 0
	a <- matrix(nrow=(max+1),ncol=(max+1))
	a[1,] <- 0
	a[,1] <- 0
	rho <- mixing(epsilon)
	for(i in (max+1):2) {
		if(i < (max+1)) remaining <- s[i] - sum(rho[i,(i+1):(max+1)]*a[i,(i+1):(max+1)],na.rm=T) else remaining <- s[i]
		sex <- s[2:i]/sum(rho[i,2:i]*s[2:i],na.rm=T)*remaining
		a[i,2:i] <- sex
		a[2:i,i] <- sex
	}
	a
}

# ODE model
model <- function(t, x, parms) {
	with(as.list(parms),{
		dx <- numeric()
		for(i in 1:(max+1)) dx[i] <- mu*sum(x*N) + (1-x[i])*c[i]*sum(b[i,]*rho[i,]*x) - gamma*x[i] - mu*x[i]
		list(dx)
	})
}

# Some conditions
years <- 1e3
times  <- c(0,years)

# Testing data and starting values
# number: total number of tests in each class; pos: number of positive tests in each class
number <- round(data$w.number)
pos <- round(data$w.number*data$w.prevalence)
init <- pos/number
for(i in 1:length(init)) {
	if(is.na(init[i])) init[i] <- 0
}

par(mfrow=c(1,2))
plot(NA,log="x",xlim=c(1e-3,5e-1),ylim=c(0,0.5),xlab="STI prevalence",ylab="Gini coefficient",frame=FALSE,axes=FALSE)
axis(1)
axis(2)
mtext("A",side=3,adj=0,cex=1,font=2)

beta <- c(0.075,0.1,0.15,0.25)
gamma <- 1/c(1,2,4,6)
epsilon <- 0
rho <- mixing(epsilon)
a <- sexacts(0,28/365,1,epsilon)

for(i in beta) {
	p <- numeric()
	g <- numeric()
	for(j in gamma) {
		b <- 1-(1-i)^a
		out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=j,mu=1))
		prev <- out[dim(out)[1],2:(max+2)]
		p <- c(p,sum(N*prev))
		g <- c(g,gini(N,prev))
	}
	lines(p,g,lty=2)
	text(max(p),min(g),paste(100*i,"%"),pos=1,cex=0.75)
}

for(j in gamma) {
	p <- numeric()
	g <- numeric()
	for(i in beta) {
		b <- 1-(1-i)^a
		out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=j,mu=1))
		prev <- out[dim(out)[1],2:(max+2)]
		p <- c(p,sum(N*prev))
		g <- c(g,gini(N,prev))
	}
	lines(p,g,lty=2)	
	text(max(p),max(g),paste(1/j,"years"),pos=4,cex=0.75)
}

# Add the different STIs
load(file = file.path('out', 'dprev.RData'))
load(file = file.path('out', 'dfig1_gini.RData'))
cols <- rainbow_hcl(length(dprev$STI))
for(i in 1:length(dprev$STI)) {
    points(dprev[i,2],dfig1_gini[i,2],pch=c(0,1,2,4,5,6)[i],col=cols[i])
    lines(c(dprev[i,2],dprev[i,2]),c(dfig1_gini[i,3],dfig1_gini[i,4]),col=cols[i])
    lines(c(dprev[i,3],dprev[i,4]),c(dfig1_gini[i,2],dfig1_gini[i,2]),col=cols[i])
    #text(dprev[i,2],dfig1_gini[i,2]+0.01,dprev[i,1],adj=c(0,0),cex=0.75,col=cols[i])
}

legend("topright",inset=0.0,legend=c("CT","MG","HPV 6","HPV 11","HPV 16","HPV 18"),col=cols,pch=c(0,1,2,4,5,6),bty="n")

plot(NA,log="",xlim=c(0.0125,0.03),ylim=c(0.15,0.65),xlab="STI prevalence",ylab="Gini coefficient",frame=FALSE,axes=FALSE)
axis(1)
axis(2)
mtext("B",side=3,adj=0,cex=1,font=2)

# Add simulations
red <- 0.1
beta_start <- 0.282
gamma_start <- 0.955
beta <- beta_start
gamma <- gamma_start
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x0 <- sum(N*prev)
y0 <- gini(N,prev)
beta <- beta_start*(1-red)
gamma <- gamma_start
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x1 <- sum(N*prev)
y1 <- gini(N,prev)
arrows(x0,y0,x1,y1,length=0.1,lty=1,lwd=2,col="darkgray")
beta <- beta_start
gamma <- gamma_start/(1-red)
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x1 <- sum(N*prev)
y1 <- gini(N,prev)
arrows(x0,y0,x1,y1,length=0.1,lty=1,lwd=2,col="darkgray")

beta_start <- 0.225
gamma_start <- 0.665
beta <- beta_start
gamma <- gamma_start
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x0 <- sum(N*prev)
y0 <- gini(N,prev)
beta <- beta_start*(1-red)
gamma <- gamma_start
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x1 <- sum(N*prev)
y1 <- gini(N,prev)
arrows(x0,y0,x1,y1,length=0.1,lty=1,lwd=2,col="darkgray")
beta <- beta_start
gamma <- gamma_start/(1-red)
b <- 1-(1-beta)^a
out <- ode(init, times, model, parms=list(rho=rho,b=b,gamma=gamma,mu=1))
prev <- out[dim(out)[1],2:(max+2)]
x1 <- sum(N*prev)
y1 <- gini(N,prev)
arrows(x0,y0,x1,y1,length=0.1,lty=1,lwd=2,col="darkgray")

# Add CT from Natsal-2 and Natsal-3
cols <- rainbow_hcl(2)
i <- 1
lines(c(dprev[i,2],dprev[i,2]),c(dfig1_gini[i,3],dfig1_gini[i,4]),col=cols[i])
lines(c(dprev[i,3],dprev[i,4]),c(dfig1_gini[i,2],dfig1_gini[i,2]),col=cols[i])
text(dprev[i,2],dfig1_gini[i,2]+0.01," CT Natsal-3",adj=c(0,0),cex=0.75,col=cols[i])
    
load(file = file.path('out', 'dprev_nat2.RData'))
load(file = file.path('out', 'ct_nat2.RData'))

x <- dprev_nat2[1,2] # Prevalence in Natsal-2
xci <- dprev_nat2[1,3:4]
y <- ct_nat2$gini # Gini coefficient in Natsal-2
yci <- ct_nat2$gini_bCI
lines(c(x,x),yci,col=cols[2])
lines(xci,c(y,y),col=cols[2])
text(x,y+0.01," CT Natsal-2",adj=c(0,0),cex=0.75,col=cols[2])

# Add points on top
points(dprev[i,2],dfig1_gini[i,2],pch=16,col=cols[i])
points(x,y,pch=16,col=cols[2])
