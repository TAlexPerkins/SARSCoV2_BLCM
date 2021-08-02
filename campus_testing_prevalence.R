## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(tidyverse,
       ## incidence,
       ## plotrix,
       ## distr,
       ## data.table,
       ## viridis,
       surveillance,
       ## BayesianTools,
       ## grDevices,
       ## plotfunctions
       )

## Load data and organize
covid.data <- fread("./COVID_Surveillance_Impact_20200810-20201204.csv")
dates.onset <- rep(covid.data$CRU_COMPLETED_DATE, covid.data$Positive)
test.type <- rep(covid.data$Test_Type, covid.data$Positive)
incid <- incidence(as.Date(dates.onset), groups = test.type)
dates <- incid$dates
cases <- rowSums(incid$counts)
cases.symp <- incid$counts[,"Symptomatic"]

## parameters for incubation period and testing delay
inc.meanlog <- 1.621
inc.sdlog <- 0.418
delay.mean <- 2
delay.rate <- 1/delay.mean

## Denominator data
## https://www.nd.edu/about/
## https://news.nd.edu/news/we-are-all-nd-reaches-3000-employees/
n.students <- 12681
n.facstaff <- 3000/0.6

## Deconvolve case notifications with incubation period and delay, using surveillance package
dmax <- 50 # truncation of incubation plus pmf
cases <- c(rep(0,dmax),cases)
cases.symp <- c(rep(0,dmax),cases.symp)
names(cases) <- c(seq(min(dates)-dmax,
                      min(dates)-1,length.out=dmax),dates) 
names(cases.symp) <- names(cases)
dates <- as.Date(names(cases))
sts.symp <- new("sts", epoch=1:length(cases.symp),observed=matrix(cases.symp,ncol=1))
## functions to convolve
Inc <- Lnorm(inc.meanlog,inc.sdlog)
Delay <- Pois(delay.mean) 
p.inf2test <- p(Inc+Delay)
## convert pdf to pmf
inc.symp.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))
bpnp.control <- list(k=10,eps=rep(1e-5,2),iter.max=rep(250,2),B=-1)
sts.symp.bp <- backprojNP(sts.symp, incu.pmf=inc.symp.pmf,
                          control=modifyList(bpnp.control,list(eq3a.method="C")))
n.exp <- upperbound(sts.symp.bp)

## Get all cases by multiplying by proportion asymptomatic
n.exp.tot <- n.exp/0.57

## Convert to prevalence
sensitivity <- fread("./grassly_pcr_sensitivity_empirical.csv")
sensitivity[,sensitivity:=sensitivity/100]
prev <- rep(0,length(n.exp.tot))
for (ii in 2:length(prev)){
    prev[ii] <- sum(n.exp[(ii-1):max(1,ii-length(sensitivity$sensitivity))]*sensitivity$sensitivity[1:min(length(sensitivity$sensitivity),ii-1)])
}

## Write to file
fwrite(data.table(Dates=dates,Prevalence=prev/(n.students+n.facstaff)),
       file=paste0("prevalence_with_delay_",delay.mean,".csv"))
