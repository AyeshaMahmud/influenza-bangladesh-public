######################################
### Computation of center of gravity with the 'circular' package
### Juliette Paireau, 2015
######################################
library(circular)
# Returns CG mean, CIlower, CIupper
center.gravity<-function(data,var,period, reps = 1000){
  
  # var is the name of the cases column
  # period is 12 for monthly data, 52 for weekly data
  
  stats<-as.data.frame(matrix(NA,ncol=3,nrow=1))
  names(stats)<-c('mean','CIlower','CIupper')
  
  # compute monthly sum of cases over the whole time period (check that you have only FULL years, or it will be biased) 
  data <- aggregate(data[,which(names(data)==var)],by=list(data$month),FUN=sum,na.rm=T)
  names(data) <- c('month','cases')   
  
  # repeat each month times the number of cases 
  x<-rep(data$month,times=data$cases)
  
  # create circular object in radians
  x<-circular(x*2*pi/period,modulo="2pi")
  
  # compute center of gravity in months
  stats$mean[1]<-mean.circular(x)/(2*pi)*period
  
  # compute confidence interval in radians
  ci<-mle.vonmises.bootstrap.ci(x, bias = T, alpha = 0.05, reps = reps)
  # fix problem with CI
  if (ci$mu.ci[1] > ci$mu.ci[2]) {
    ci$mu.ci[1]<-ci$mu.ci[1]-2*pi
  }
  # convert from radian to month 
  stats$CIlower[1]<-ci$mu.ci[1]/(2*pi)*period
  stats$CIupper[1]<-ci$mu.ci[2]/(2*pi)* period
  
  return(stats)
}

