library(BayesianTools)
library(lubridate)



set.seed(46556)

# read in data
load('data.RData')



# extract some basic information about the data set
sum(y$Count) # 853 individuals
sum((y$INCIDENT_TYPE==1)*y$Count) # 41 red passes
sum((y$COVID_PCR_TEST_RESULT!='')*y$Count) # 833 commercial tests
sum((y$COVID_RAPID_TEST_RESULT!='')*y$Count) # 37 antigen tests
sum((y$SALIVA_PCR_RESULT!='')*y$Count) # 846 saliva tests
sum((y$COVID_PCR_TEST_RESULT!='' &
     y$COVID_RAPID_TEST_RESULT=='' &
     y$SALIVA_PCR_RESULT!='')*y$Count) # 799 comm & saliva but not antigen
sum((y$COVID_PCR_TEST_RESULT!='' &
       y$COVID_RAPID_TEST_RESULT!='' &
       y$SALIVA_PCR_RESULT!='')*y$Count) # 27 with all three tests



# make some calculations about positivity

# commercial tests
sum(as.numeric(as.character(y$COVID_PCR_TEST_RESULT)) * y$Count, na.rm=T)
# 12 commercial positives
sum((y$COVID_PCR_TEST_RESULT!='')*y$Count)
# 833 commercial tests
sum(ifelse(
  y$INCIDENT_TYPE==1,
  as.numeric(as.character(y$COVID_PCR_TEST_RESULT)) * y$Count,
  0),na.rm=T)
# 2 commercial positives among red passes
# 10 commercial positives among surveillance
sum(ifelse(
  y$INCIDENT_TYPE==1,
  (y$COVID_PCR_TEST_RESULT!='')*y$Count,
  0),na.rm=T)
# 31 commercial tests among red passes
# 802 commercial tests among surveillance

# rapid antigen tests
sum(as.numeric(as.character(y$COVID_RAPID_TEST_RESULT)) * y$Count, na.rm=T)
# 5 rapid positives
sum((y$COVID_RAPID_TEST_RESULT!='')*y$Count)
# 37 rapid tests
sum(ifelse(
  y$INCIDENT_TYPE==1,
  as.numeric(as.character(y$COVID_RAPID_TEST_RESULT)) * y$Count,
  0),na.rm=T)
# 5 rapid positives among red passes
# 0 rapid positives among surveillance
sum(ifelse(
  y$INCIDENT_TYPE==1,
  (y$COVID_RAPID_TEST_RESULT!='')*y$Count,
  0),na.rm=T)
# 34 rapid tests among red passes
# 3 rapid tests among surveillance

# saliva tests
sum(as.numeric(as.character(y$SALIVA_PCR_RESULT)) * y$Count, na.rm=T)
# 26 saliva positives
sum((y$SALIVA_PCR_RESULT!='')*y$Count)
# 846 saliva tests
sum(ifelse(
  y$INCIDENT_TYPE==1,
  as.numeric(as.character(y$SALIVA_PCR_RESULT)) * y$Count,
  0),na.rm=T)
# 6 saliva positives among red passes
# 20 saliva positives among surveillance
sum(ifelse(
  y$INCIDENT_TYPE==1,
  (y$SALIVA_PCR_RESULT!='')*y$Count,
  0),na.rm=T)
# 41 saliva tests among red passes
# 805 saliva tests among surveillance

# calculate concodance of the three test types
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='0'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='0'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='0'))])
# 99.27%
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='0'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'))])
# 96.30%
sum(y$Count[c(
  which(y$SALIVA_PCR_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$SALIVA_PCR_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'))]) /
sum(y$Count[c(
  which(y$SALIVA_PCR_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$SALIVA_PCR_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='0'),
  which(y$SALIVA_PCR_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$SALIVA_PCR_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'))])
# 97.30%

# calulative sensitivity and specificity in reference to commercial
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='1'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='1'&y$SALIVA_PCR_RESULT=='0'))])
# 0.833
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='0'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='0'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$SALIVA_PCR_RESULT=='1'))])
# 0.995
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='1'),
  which(y$COVID_PCR_TEST_RESULT=='1'&y$COVID_RAPID_TEST_RESULT=='0'))])
# 0.5
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'))]) /
sum(y$Count[c(
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='0'),
  which(y$COVID_PCR_TEST_RESULT=='0'&y$COVID_RAPID_TEST_RESULT=='1'))])
# 1.0



# functions to enable MCMC

# log likelihood function
LL = function(par){
  sens_lc = par[1]; sens_sa = par[2]; sens_an = par[3]
  spec_lc = par[4]; spec_sa = par[5]; spec_an = par[6]
  prev_symp = par[7]; prev_asymp = par[8]
  Pr_lc_pos = ifelse(
    y$COVID_PCR_TEST_RESULT=='',1,
    ifelse(y$COVID_PCR_TEST_RESULT=='1',sens_lc,1-sens_lc))
  Pr_lc_neg = ifelse(
    y$COVID_PCR_TEST_RESULT=='',1,
    ifelse(y$COVID_PCR_TEST_RESULT=='1',1-spec_lc,spec_lc))
  Pr_sa_pos = ifelse(
    y$SALIVA_PCR_RESULT=='',1,
    ifelse(y$SALIVA_PCR_RESULT=='1',sens_sa,1-sens_sa))
  Pr_sa_neg = ifelse(
    y$SALIVA_PCR_RESULT=='',1,
    ifelse(y$SALIVA_PCR_RESULT=='1',1-spec_sa,spec_sa))
  Pr_an_pos = ifelse(
    y$COVID_RAPID_TEST_RESULT=='',1,
    ifelse(y$COVID_RAPID_TEST_RESULT=='1',sens_an,1-sens_an))
  Pr_an_neg = ifelse(
    y$COVID_RAPID_TEST_RESULT=='',1,
    ifelse(y$COVID_RAPID_TEST_RESULT=='1',1-spec_an,spec_an))
  Pr_pos = ifelse(y$INCIDENT_TYPE=='1',prev_symp,prev_asymp)
  Pr_neg = 1 - Pr_pos
  sum(y$Count * log(
    Pr_lc_pos * Pr_sa_pos * Pr_an_pos * Pr_pos +
    Pr_lc_neg * Pr_sa_neg * Pr_an_neg * Pr_neg))
}

# prior density function
prior_density = function(par){
  dbeta(par[1],1,1,log=T) + dbeta(par[2],1,1,log=T) +
  dbeta(par[3],1,1,log=T) + dbeta(par[4],1,1,log=T) +
  dbeta(par[5],1,1,log=T) + dbeta(par[6],1,1,log=T) +
  dbeta(par[7],1,9,log=T) + dbeta(par[8],1,99,log=T)
}

# prior sampler function
prior_sampler = function(n=1){
  cbind(rbeta(n,1,1), rbeta(n,1,1),
        rbeta(n,1,1), rbeta(n,1,1),
        rbeta(n,1,1), rbeta(n,1,1),
        rbeta(n,1,9), rbeta(n,1,99))
}

# set up prior
prior = createPrior(
  prior_density, prior_sampler,
  rep(1e-7,8), rep(1-1e-7,8))

# run and process MCMC
settings = list(iterations = 1e5, nrChains = 3, message = F)
bayesianSetup = createBayesianSetup(
  likelihood = LL, prior = prior)
mcmc = runMCMC(bayesianSetup, settings = settings)
mcmc.coda = getSample(
  mcmc, coda=T, start=1e4, thin=100)
post = getSample(
  mcmc, parametersOnly = T, thin = 100, start = 1e4)
save(list=ls(),file='mcmc.RData')



jpeg('traceplot.jpeg',width=6.5,height=6.5,units='in',res=300)

varnames(mcmc.coda) = c(
  'Se_Comm','Se_Saliva','Se_Antigen',
  'Sp_Comm','Sp_Saliva','Sp_Antigen',
  'Prev_Non-surv','Prev_Surv')
layout(matrix(1:8,4,2,byrow=T))
plot(mcmc.coda,density=F)

dev.off()

gelmanDiagnostics(mcmc)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# par 1       1.00       1.01
# par 2       1.00       1.00
# par 3       1.00       1.00
# par 4       1.01       1.01
# par 5       1.00       1.01
# par 6       1.01       1.01
# par 7       1.00       1.01
# par 8       1.00       1.01
# 
# Multivariate psrf
# 
# 1.01



jpeg('correlationplot.jpeg',width=6.5,height=7.5,units='in',res=300)

correlationPlot(mcmc.coda)

dev.off()



jpeg('params_post.jpeg',width=6.5,height=6.5,units='in',res=300)

layout(matrix(1:4,2,2,byrow=T))
par(mar=c(5,4,2,2))

prev.surv.post = density(post[,8])
prev.surv.prior = list(
  x=c(prev.surv.post$x,max(prev.surv.post$x),0),
  y=c(dbeta(prev.surv.post$x,1,99),0,0))
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Surveillance prevalence',ylab='Density',
     xlim=range(c(prev.surv.post$x)),
     ylim=c(0,1.05)*range(c(prev.surv.prior$y,prev.surv.post$y)))
polygon(prev.surv.prior,col=rgb(0,0,0,0.3))
polygon(prev.surv.post,col=rgb(1,1,0,0.3),border=rgb(1,1,0,1))
mtext('A',side=3,at=min(c(prev.surv.post$x)))
legend('topright',legend=c('Prior','Posterior'),
       fill=c(rgb(0,0,0,0.3),rgb(1,1,0,0.3)),bty='n')

prev.red.post = density(post[,7])
prev.red.prior = list(
  x=c(prev.red.post$x,max(prev.red.post$x),0),
  y=c(dbeta(prev.red.post$x,1,9),0,0))
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Non-surveillance prevalence',ylab='Density',
     xlim=range(c(prev.red.post$x)),
     ylim=c(0,1.05)*range(c(prev.red.prior$y,prev.red.post$y)))
polygon(prev.red.prior,col=rgb(0,0,0,0.3))
polygon(prev.red.post,col=rgb(1,1,0,0.3),border=rgb(1,1,0,1))
mtext('B',side=3,at=min(c(prev.red.post$x)))

d.sens.comm.surv = density(post[,1])
d.sens.salv.surv = density(post[,2])
d.sens.antg.surv = density(post[,3])
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Sensitivity',ylab='Density',
     xlim=range(c(d.sens.comm.surv$x,d.sens.salv.surv$x,d.sens.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.sens.comm.surv$y,d.sens.salv.surv$y,d.sens.antg.surv$y)))
polygon(d.sens.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.sens.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.sens.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('C',side=3,at=min(c(d.sens.comm.surv$x,d.sens.salv.surv$x,d.sens.antg.surv$x)))
legend('topleft',legend=c('Commercial','Saliva','Antigen'),
       fill=c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3)),bty='n')

d.spec.comm.surv = density(post[,4])
d.spec.salv.surv = density(post[,5])
d.spec.antg.surv = density(post[,6])
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Specificity',ylab='Density',
     xlim=c(0.9,max(c(d.spec.comm.surv$x,d.spec.salv.surv$x,d.spec.antg.surv$x))),
     ylim=c(0,1.05)*range(c(d.spec.comm.surv$y,d.spec.salv.surv$y,d.spec.antg.surv$y)))
polygon(d.spec.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.spec.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.spec.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('D',side=3,at=0.9)

dev.off()


# evaluate positive and negative predictive values
ppv.comm.surv = post[,1]*post[,8] / (post[,1]*post[,8] + (1-post[,4])*(1-post[,8]))
ppv.salv.surv = post[,2]*post[,8] / (post[,2]*post[,8] + (1-post[,5])*(1-post[,8]))
ppv.antg.surv = post[,3]*post[,8] / (post[,3]*post[,8] + (1-post[,6])*(1-post[,8]))
npv.comm.surv = post[,4]*(1-post[,8]) / ((1-post[,1])*post[,8] + post[,4]*(1-post[,8]))
npv.salv.surv = post[,5]*(1-post[,8]) / ((1-post[,2])*post[,8] + post[,5]*(1-post[,8]))
npv.antg.surv = post[,6]*(1-post[,8]) / ((1-post[,3])*post[,8] + post[,6]*(1-post[,8]))

jpeg('predvalue_surv.jpeg',width=6.5,height=3.25,units='in',res=300)

layout(matrix(1:2,1,2))
par(mar=c(5,4,2,2))

d.ppv.comm.surv = density(ppv.comm.surv)
d.ppv.salv.surv = density(ppv.salv.surv)
d.ppv.antg.surv = density(ppv.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Positive predictive value',ylab='Density',
     xlim=range(c(d.ppv.comm.surv$x,d.ppv.salv.surv$x,d.ppv.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.ppv.comm.surv$y,d.ppv.salv.surv$y,d.ppv.antg.surv$y)))
polygon(d.ppv.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.ppv.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.ppv.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('A',side=3,at=min(c(d.ppv.comm.surv$x,d.ppv.salv.surv$x,d.ppv.antg.surv$x)))
legend('topleft',legend=c('Commercial','Saliva','Antigen'),
       fill=c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3)),bty='n')

d.npv.comm.surv = density(npv.comm.surv)
d.npv.salv.surv = density(npv.salv.surv)
d.npv.antg.surv = density(npv.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='Negative predictive value',ylab='Density',
     xlim=range(c(d.npv.comm.surv$x,d.npv.salv.surv$x,d.npv.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.npv.comm.surv$y,d.npv.salv.surv$y,d.npv.antg.surv$y)))
polygon(d.npv.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.npv.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.npv.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('B',side=3,at=min(c(d.npv.comm.surv$x,d.npv.salv.surv$x,d.npv.antg.surv$x)))

dev.off()



# predicted outcomes per 1,000 surveillance tests
tp.comm.surv = 1e3 * post[,1]*post[,8]
fp.comm.surv = 1e3 * (1-post[,4])*(1-post[,8])
tn.comm.surv = 1e3 * post[,4]*(1-post[,8])
fn.comm.surv = 1e3 * (1-post[,1])*post[,8]
tp.salv.surv = 1e3 * post[,2]*post[,8]
fp.salv.surv = 1e3 * (1-post[,5])*(1-post[,8])
tn.salv.surv = 1e3 * post[,5]*(1-post[,8])
fn.salv.surv = 1e3 * (1-post[,2])*post[,8]
tp.antg.surv = 1e3 * post[,3]*post[,8]
fp.antg.surv = 1e3 * (1-post[,6])*(1-post[,8])
tn.antg.surv = 1e3 * post[,6]*(1-post[,8])
fn.antg.surv = 1e3 * (1-post[,3])*post[,8]

jpeg('testoutcomes_surv.jpeg',width=6.5,height=6.5,units='in',res=300)

layout(matrix(1:4,2,2,byrow=T))
par(mar=c(5,4,2,2))

d.tp.comm.surv = density(tp.comm.surv)
d.tp.salv.surv = density(tp.salv.surv)
d.tp.antg.surv = density(tp.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='True positives per 1,000 tests',ylab='Density',
     xlim=range(c(d.tp.comm.surv$x,d.tp.salv.surv$x,d.tp.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.tp.comm.surv$y,d.tp.salv.surv$y,d.tp.antg.surv$y)))
polygon(d.tp.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.tp.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.tp.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('A',side=3,at=min(c(d.tp.comm.surv$x,d.tp.salv.surv$x,d.tp.antg.surv$x)))

d.fp.comm.surv = density(fp.comm.surv)
d.fp.salv.surv = density(fp.salv.surv)
d.fp.antg.surv = density(fp.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='False positives per 1,000 tests',ylab='Density',
     xlim=range(c(d.fp.comm.surv$x,d.fp.salv.surv$x,d.fp.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.fp.comm.surv$y,d.fp.salv.surv$y,d.fp.antg.surv$y)))
polygon(d.fp.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.fp.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.fp.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('B',side=3,at=min(c(d.fp.comm.surv$x,d.fp.salv.surv$x,d.fp.antg.surv$x)))
legend('topright',legend=c('Commercial','Saliva','Antigen'),
       fill=c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3)),bty='n')

d.fn.comm.surv = density(fn.comm.surv)
d.fn.salv.surv = density(fn.salv.surv)
d.fn.antg.surv = density(fn.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='False negatives per 1,000 tests',ylab='Density',
     xlim=range(c(d.fn.comm.surv$x,d.fn.salv.surv$x,d.fn.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.fn.comm.surv$y,d.fn.salv.surv$y,d.fn.antg.surv$y)))
polygon(d.fn.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.fn.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.fn.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('C',side=3,at=min(c(d.fn.comm.surv$x,d.fn.salv.surv$x,d.fn.antg.surv$x)))

d.tn.comm.surv = density(tn.comm.surv)
d.tn.salv.surv = density(tn.salv.surv)
d.tn.antg.surv = density(tn.antg.surv)
plot(-100,-100,xaxs='i',yaxs='i',las=1,xlab='True negatives per 1,000 tests',ylab='Density',
     xlim=range(c(d.tn.comm.surv$x,d.tn.salv.surv$x,d.tn.antg.surv$x)),
     ylim=c(0,1.05)*range(c(d.tn.comm.surv$y,d.tn.salv.surv$y,d.tn.antg.surv$y)))
polygon(d.tn.comm.surv,col=rgb(1,0,0,0.3),border=rgb(1,0,0,1))
polygon(d.tn.salv.surv,col=rgb(0,1,0,0.3),border=rgb(0,1,0,1))
polygon(d.tn.antg.surv,col=rgb(0,0,1,0.3),border=rgb(0,0,1,1))
mtext('D',side=3,at=min(c(d.tn.comm.surv$x,d.tn.salv.surv$x,d.tn.antg.surv$x)))

dev.off()




predvalfun = function(prev.in){
  cbind(
    post[,1]*prev.in / (post[,1]*prev.in + (1-post[,4])*(1-prev.in)),
    post[,2]*prev.in / (post[,2]*prev.in + (1-post[,5])*(1-prev.in)),
    post[,3]*prev.in / (post[,3]*prev.in + (1-post[,6])*(1-prev.in)),
    post[,4]*(1-prev.in) / ((1-post[,1])*prev.in + post[,4]*(1-prev.in)),
    post[,5]*(1-prev.in) / ((1-post[,2])*prev.in + post[,5]*(1-prev.in)),
    post[,6]*(1-prev.in) / ((1-post[,3])*prev.in + post[,6]*(1-prev.in)))
}

prev.time = seq(0.01,0.1,length.out=125)
prev.time = read.csv('prevalence_with_delay_2.csv')
prev.time = prev.time[42:167,2]

predvalstime = array(NA,dim=c(length(prev.time),3,6))

for(tt in 1:length(prev.time)){
  predvalstime[tt,,] = apply(
    predvalfun(prev.time[tt]),
    2,
    function(x)quantile(x,c(0.025,0.5,0.975)))
}

jpeg('predvaltime.jpeg',width=6.5,height=3.25,units='in',res=300)

par(mar=c(1.5,1,1.5,1),oma=c(1,4,0.25,4))
layout(matrix(1:6,2,3,byrow=T))

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7,labels=NA)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,las=1,col.ticks=2)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,1],rev(predvalstime[,3,1])),
  border=NA,col=rgb(1,0,0,0.3))
lines(predvalstime[,2,1],col=2)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('A',3,at=1,cex=0.7)
mtext('Commercial',3,cex=0.7,line=0.5)
mtext('Positive predictive value',2,line=3.5,cex=0.7)

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7,labels=NA)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,labels=NA,col.ticks=3)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,2],rev(predvalstime[,3,2])),
  border=NA,col=rgb(0,1,0,0.3))
lines(predvalstime[,2,2],col=3)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('B',3,at=1,cex=0.7)
mtext('Saliva',3,cex=0.7,line=0.5)

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7)
mtext('Prevalence',4,line=3.75,cex=0.7)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,labels=NA,col.ticks=4)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,3],rev(predvalstime[,3,3])),
  border=NA,col=rgb(0,0,1,0.3))
lines(predvalstime[,2,3],col=4)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('C',3,at=1,cex=0.7)
mtext('Antigen',3,cex=0.7,line=0.5)

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7,labels=NA)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0.978,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,las=1,col.ticks=2)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,4],rev(predvalstime[,3,4])),
  border=NA,col=rgb(1,0,0,0.3))
lines(predvalstime[,2,4],col=2)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('D',3,at=1,cex=0.7)
mtext('Negative predictive value',2,line=3.5,cex=0.7)
mtext('Date',1,cex=0.7,line=1.25)

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7,labels=NA)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0.978,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,labels=NA,col.ticks=3)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,5],rev(predvalstime[,3,5])),
  border=NA,col=rgb(0,1,0,0.3))
lines(predvalstime[,2,5],col=3)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('E',3,at=1,cex=0.7)
mtext('Date',1,cex=0.7,line=1.25)

plot(prev.time,xaxt='n',yaxt='n',box='off',type='l',xaxs='i')
axis(4,las=1,cex=0.7)
mtext('Prevalence',4,line=3.75,cex=0.7)
par(new=T)
plot(-100,-100,xlim=c(1,length(prev.time)),ylim=c(0.978,1),
     xaxs='i',las=1,xaxt='n',yaxt='n')
axis(2,labels=NA,col.ticks=4)
polygon(
  c(1:length(prev.time),rev(1:length(prev.time))),
  c(predvalstime[,1,6],rev(predvalstime[,3,6])),
  border=NA,col=rgb(0,0,1,0.3))
lines(predvalstime[,2,6],col=4)
axis(1,at=c(1,1+31,1+31+30,1+31+30+31,1+31+30+31+30),labels=c('','','','',''))
mtext(c('Aug','Sep','Oct','Nov'),1,at=c(1+15,1+15+31,1+15+31+30,1+15+31+30+31),cex=0.7)
mtext('F',3,at=1,cex=0.7)
mtext('Date',1,cex=0.7,line=1.25)

dev.off()


