# set aside data and posterior estimates from the primary analysis of empirical data
y.orig = y
post.orig = post

# replace the posterior samples with independent uniform draws from ranges
for(ii in 1:ncol(post.orig)){
  post.orig[,ii] = runif(
    nrow(post.orig),
    min = quantile(post.orig[,ii],0.025),
    max = quantile(post.orig[,ii],0.975))
}



# number of simulated data sets
n.reps = 100

# allocate storage for summary statistics of analysis of simulated data
post.sim.stats = array(NA,dim=c(n.reps,3,8))

# draw random samples from the posterior distribution to use for simulation
pp.vec = sample(nrow(post.orig),n.reps,replace=F)

# loop across simulated data sets
for(pp in 1:n.reps){
  
  # loop across each type of participant in the empirical study, in terms of
  # reason for testing and types of tests applied
  for(ii in 1:nrow(tested)){
    
    # simulate whether they really have an infection or not
    prev.sim = ifelse(tested[ii,'INCIDENT_TYPE']=='1',post.orig[pp.vec[pp],7],post.orig[pp.vec[pp],8])
    inf.sim = rbinom(tested[ii,'Count'],1,prev.sim)
    
    # simulate the outcome of each type of test, conditional on each individual's true infection status
    which.test = which(tested[ii,]=='1')
    # commercial
    if(2 %in% which.test){
      pos.comm.sim = rbinom(length(inf.sim),1,ifelse(inf.sim,post.orig[pp.vec[pp],1],1-post.orig[pp.vec[pp],4]))
    } else {
      pos.comm.sim = rep('',length(inf.sim))
    }
    # antigen
    if(3 %in% which.test){
      pos.antg.sim = rbinom(length(inf.sim),1,ifelse(inf.sim,post.orig[pp.vec[pp],3],1-post.orig[pp.vec[pp],6]))
    } else {
      pos.antg.sim = rep('',length(inf.sim))
    }
    # saliva
    if(4 %in% which.test){
      pos.salv.sim = rbinom(length(inf.sim),1,ifelse(inf.sim,post.orig[pp.vec[pp],2],1-post.orig[pp.vec[pp],5]))
    } else {
      pos.salv.sim = rep('',length(inf.sim))
    }
    
    # collate data from each type of test
    if(ii == 1){
      # for the first category of participants
      y = data.frame(
        INCIDENT_TYPE = rep(as.character(tested[ii,'INCIDENT_TYPE']),length(inf.sim)),
        COVID_PCR_TEST_RESULT = as.character(pos.comm.sim),
        COVID_RAPID_TEST_RESULT = as.character(pos.antg.sim),
        SALIVA_PCR_RESULT = as.character(pos.salv.sim))
    } else {
      # for all subsequent categories of participants
      y = rbind(
        y,
        data.frame(
          INCIDENT_TYPE = rep(as.character(tested[ii,'INCIDENT_TYPE']),length(inf.sim)),
          COVID_PCR_TEST_RESULT = as.character(pos.comm.sim),
          COVID_RAPID_TEST_RESULT = as.character(pos.antg.sim),
          SALIVA_PCR_RESULT = as.character(pos.salv.sim)))
    }
  }
  # count up the number of individuals with a given combination of reason for testing,
  # types of tests applied, and outcomes of each test
  y = as.data.frame.table(table(
    y$INCIDENT_TYPE,y$COVID_PCR_TEST_RESULT,
    y$COVID_RAPID_TEST_RESULT,y$SALIVA_PCR_RESULT))
  names(y) = c(
    'INCIDENT_TYPE','COVID_PCR_TEST_RESULT',
    'COVID_RAPID_TEST_RESULT','SALIVA_PCR_RESULT','Count')
  y = y[-which(y$Count==0),]
  
  # perform MCMC to obtain posterior parameter estimates, just like in the analysis of empirical data
  mcmc = runMCMC(bayesianSetup, settings = settings)
  mcmc.coda = getSample(
    mcmc, coda=T, start=1e4, thin=100)
  post.sim = getSample(
    mcmc, parametersOnly = T, thin = 100, start = 1e4)
  post.sim.stats[pp,,] = apply(post.sim,2,function(x)quantile(x,c(0.025,0.5,0.975)))
}



# proportion of replicates in which model underestimates
colMeans(post.orig[pp.vec,] > post.sim.stats[,2,])
# par 1 par 2 par 3 par 4 par 5 par 6 par 7 par 8 
# 0.41  0.59  0.58  0.69  0.68  0.69  0.64  0.80
mean(post.orig[pp.vec,] > post.sim.stats[,2,])
# 0.635

# coverage probability
colMeans((post.orig[pp.vec,] >= post.sim.stats[,1,]) & (post.orig[pp.vec,] <= post.sim.stats[,3,]))
# par 1 par 2 par 3 par 4 par 5 par 6 par 7 par 8 
# 0.99  1.00  0.99  0.99  0.99  0.95  0.97  0.98
mean((post.orig[pp.vec,] >= post.sim.stats[,1,]) & (post.orig[pp.vec,] <= post.sim.stats[,3,]))
# 0.9825



# make a figure showing the 95% credible intervals compared to simulated values
jpeg('simulation_validation.jpeg',width=6.5,height=6.5,units='in',res=300)

layout(matrix(1:8))
par(mar=c(0.25,6,0.25,0),oma=c(0.5,0,0,0.25))

plot(
  post.sim.stats[,2,1],ylim=range(post.sim.stats[,,1]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,1],1:n.reps,post.sim.stats[,3,1],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,1],pch=4)
mtext('Commercial',2,line=5,cex=0.7)
mtext('sensitivity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,2],ylim=range(post.sim.stats[,,2]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,2],1:n.reps,post.sim.stats[,3,2],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,2],pch=4)
mtext('Saliva',2,line=5,cex=0.7)
mtext('sensitivity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,3],ylim=range(post.sim.stats[,,3]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,3],1:n.reps,post.sim.stats[,3,3],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,3],pch=4)
mtext('Antigen',2,line=5,cex=0.7)
mtext('sensitivity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,4],ylim=range(post.sim.stats[,,4]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,4],1:n.reps,post.sim.stats[,3,4],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,4],pch=4)
mtext('Commercial',2,line=5,cex=0.7)
mtext('specificity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,5],ylim=range(post.sim.stats[,,5]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,5],1:n.reps,post.sim.stats[,3,5],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,5],pch=4)
mtext('Saliva',2,line=5,cex=0.7)
mtext('specificity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,6],ylim=range(post.sim.stats[,,6]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,6],1:n.reps,post.sim.stats[,3,6],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,6],pch=4)
mtext('Antigen',2,line=5,cex=0.7)
mtext('specificity',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,7],ylim=range(post.sim.stats[,,7]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,7],1:n.reps,post.sim.stats[,3,7],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,7],pch=4)
mtext('Non-surveillance',2,line=5,cex=0.7)
mtext('prevalence',2,line=3.75,cex=0.7)
plot(
  post.sim.stats[,2,8],ylim=range(post.sim.stats[,,8]),
  xlab='',xaxt='n',las=1,ylab='',col=rgb(0,0,0,0.5))
segments(1:n.reps,post.sim.stats[,1,8],1:n.reps,post.sim.stats[,3,8],col=rgb(0,0,0,0.5))
points(post.orig[pp.vec,8],pch=4)
mtext('Surveillance',2,line=5,cex=0.7)
mtext('prevalence',2,line=3.75,cex=0.7)

dev.off()



# make a figure showing the 95% credible intervals compared to simulated values
jpeg('simulation_validation_xy.jpeg',width=6.5,height=3.25,units='in',res=300)

layout(matrix(1:8,2,4,byrow=T))
par(mar=c(1.5,1.75,1.75,1.5),oma=c(2,3,0,0))

plot(
  post.orig[pp.vec,1],
  post.sim.stats[,2,1],ylim=range(post.sim.stats[,,1]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,1],
  post.sim.stats[,1,1],
  post.orig[pp.vec,1],
  post.sim.stats[,3,1],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,1],post.sim.stats[,2,1]) # 0.65
legend('bottomright',legend='r=0.65',bty='n',adj=c(-0.12,1.3))
mtext('Comm. sens.',3,cex=0.7)
mtext('Inferred value',2,cex=0.7,line=3)
mtext('A',3,at=min(post.orig[pp.vec,1]),cex=0.7)
plot(
  post.orig[pp.vec,2],
  post.sim.stats[,2,2],ylim=range(post.sim.stats[,,2]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,2],
  post.sim.stats[,1,2],
  post.orig[pp.vec,2],
  post.sim.stats[,3,2],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,2],post.sim.stats[,2,2]) # 0.40
legend('bottomright',legend='r=0.40',bty='n',adj=c(-0.12,1.3))
mtext('Saliva sens.',3,cex=0.7)
mtext('B',3,at=min(post.orig[pp.vec,2]),cex=0.7)
plot(
  post.orig[pp.vec,3],
  post.sim.stats[,2,3],ylim=range(post.sim.stats[,,3]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,3],
  post.sim.stats[,1,3],
  post.orig[pp.vec,3],
  post.sim.stats[,3,3],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,3],post.sim.stats[,2,3]) # 0.57
legend('bottomright',legend='r=0.57',bty='n',adj=c(-0.12,1.3))
mtext('Antigen sens.',3,cex=0.7)
mtext('C',3,at=min(post.orig[pp.vec,3]),cex=0.7)
plot(
  post.orig[pp.vec,4],
  post.sim.stats[,2,4],ylim=range(post.sim.stats[,,4]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,4],
  post.sim.stats[,1,4],
  post.orig[pp.vec,4],
  post.sim.stats[,3,4],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,4],post.sim.stats[,2,4]) # 0.58
legend('bottomright',legend='r=0.58',bty='n',adj=c(-0.12,1.3))
mtext('Comm. spec.',3,cex=0.7)
mtext('D',3,at=min(post.orig[pp.vec,4]),cex=0.7)
plot(
  post.orig[pp.vec,5],
  post.sim.stats[,2,5],ylim=range(post.sim.stats[,,5]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,5],
  post.sim.stats[,1,5],
  post.orig[pp.vec,5],
  post.sim.stats[,3,5],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,5],post.sim.stats[,2,5]) # 0.76
legend('bottomright',legend='r=0.76',bty='n',adj=c(-0.12,1.3))
mtext('Saliva spec.',3,cex=0.7)
mtext('Simulated value',1,cex=0.7,line=2.25)
mtext('Inferred value',2,cex=0.7,line=3)
mtext('E',3,at=min(post.orig[pp.vec,5]),cex=0.7)
plot(
  post.orig[pp.vec,6],
  post.sim.stats[,2,6],ylim=range(post.sim.stats[,,6]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,6],
  post.sim.stats[,1,6],
  post.orig[pp.vec,6],
  post.sim.stats[,3,6],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,6],post.sim.stats[,2,6]) # 0.54
legend('bottomright',legend='r=0.54',bty='n',adj=c(-0.12,1.3))
mtext('Antigen spec.',3,cex=0.7)
mtext('Simulated value',1,cex=0.7,line=2.25)
mtext('F',3,at=min(post.orig[pp.vec,6]),cex=0.7)
plot(
  post.orig[pp.vec,7],
  post.sim.stats[,2,7],ylim=range(post.sim.stats[,,7]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,7],
  post.sim.stats[,1,7],
  post.orig[pp.vec,7],
  post.sim.stats[,3,7],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,7],post.sim.stats[,2,7]) # 0.66
legend('bottomright',legend='r=0.66',bty='n',adj=c(-0.12,1.3))
mtext('Non-surv. prev.',3,cex=0.7)
mtext('Simulated value',1,cex=0.7,line=2.25)
mtext('G',3,at=min(post.orig[pp.vec,7]),cex=0.7)
plot(
  post.orig[pp.vec,8],
  post.sim.stats[,2,8],ylim=range(post.sim.stats[,,8]),
  las=1,col=rgb(0,0,0,0.5),
  xlab='',ylab='')
segments(
  post.orig[pp.vec,8],
  post.sim.stats[,1,8],
  post.orig[pp.vec,8],
  post.sim.stats[,3,8],col=rgb(0,0,0,0.5))
abline(0,1,lty=2)
cor(post.orig[pp.vec,8],post.sim.stats[,2,8]) # 0.73
legend('bottomright',legend='r=0.73',bty='n',adj=c(-0.12,1.3))
mtext('Surv. prev.',3,cex=0.7)
mtext('Simulated value',1,cex=0.7,line=2.25)
mtext('H',3,at=min(post.orig[pp.vec,8]),cex=0.7)

dev.off()



# save outputs
save(list=ls(),file='mcmc_simulated.RData')
