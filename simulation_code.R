source('~/Desktop/Papers/functional_clustering/functions_temp.R', chdir = TRUE)



################   Simulation 1  ################
nvar=15
n=99
Lt<-list()
tt=c(1:nvar)
for(j in 1:n){  Lt[[j]]=tt }

u_sigma=0.1
num_sim=100
yj = mixed=func =funHD =matrix(nrow= num_sim, ncol=2)
for(simm in 1:num_sim){
  
  print(paste(simm,'th simulation'))
  nvar=15
  n=99
  Lt<-list()
  tt=c(1:nvar)
  for(j in 1:n){  Lt[[j]]=tt }
  Ly<-list()
  type='pois1'
  Ly[[1]]=data_generation(n, u_sigma=0.1, type, tt)
  
  type='pois2'
  Ly[[2]]=data_generation(n, u_sigma=0.1, type, tt)
  
  type='binom2'
  Ly[[3]]=data_generation(n, u_sigma=0.05, type, tt)
  
  dis=c('poisson(log)','poisson(log)', 'binomial(logit)')
  
  num_clus=3
  cls=c(rep(c(1: num_clus),each=(99/num_clus)) ) 
  true=permn(num_clus)  
  
  ################ My method  ################
  result=FPCA_y(Ly, Lt, optnsSW=list(), K=num_clus, dis) 
  
  yj[simm, 1]=perc_correct(true, result)
  yj[simm, 2]=adjustedRandIndex(cls, result)  

  ################ MixAK  ################
  output1 <- as.vector(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = F))
  output2 <- as.vector(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = F))
  output3 <- as.vector(matrix(unlist(Ly[[3]]), ncol = nvar, byrow = F))
  tt=as.vector(matrix(unlist(Lt), ncol = nvar, byrow = F))
  mod <- GLMM_MCMC(y = cbind(output1, output2, output3),  dist = dis,  id = rep(c(1:99),each=nvar), 
                   x = list(  "empty", "empty", "empty"), z = list(tt, tt, tt) ,  
                   prior.b = list(Kmax = num_clus), nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),  parallel = FALSE)
  
  fit_mix <- apply(mod[[1]]$poster.comp.prob_u, 1, which.max)
  
  
  mixed[simm, 1]= perc_correct(true, fit_mix)
  mixed[simm, 2]=adjustedRandIndex(cls, fit_mix) 
  
  
  
  
  ################ funclust  ################
   CWtime<- 1:nvar
  output1 <- t(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd1 <- smooth.basisPar( CWtime, output1,  lambda=1e-2)$fd
  output2 <- t(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd2 <- smooth.basisPar( CWtime, output2,  lambda=1e-2)$fd
  output3 <- t(matrix(unlist(Ly[[3]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd3 <- smooth.basisPar( CWtime, output3,  lambda=1e-2)$fd

  CWfd=list(CWfd1,CWfd2,CWfd3)
  
  res=funclust(CWfd,K=num_clus)
  fit_fun=res$cls
  func[simm, 1]= perc_correct(true, fit_fun)
  func[simm, 2]=adjustedRandIndex(cls, fit_fun) 
  
  ################ funHDDC ################
  
    res.multi<-funHDDC(CWfd,K=num_clus,model="AkjBQkDk",init="kmeans",threshold=0.2)
   funHDDC_result=res.multi$class
   funHD[simm, 1]= perc_correct(true, funHDDC_result)
   funHD[simm, 2]=adjustedRandIndex(cls, funHDDC_result) 
  
}



yj_out=round(apply(yj,2, mean, na.rm=T),2)
mixed_out=round(apply(mixed,2, mean, na.rm=T),2)
func_out=round(apply(func,2, mean, na.rm=T),2)
HD_out=round(apply(funHD,2, mean, na.rm=T),2)


yj_sd=round(apply(yj,2, sd, na.rm=T),2)
mixed_sd=round(apply(mixed,2, sd, na.rm=T),2)
func_sd=round(apply(func,2, sd, na.rm=T),2)
HD_sd=round(apply(funHD,2, sd, na.rm=T),2)






############################  Simulation 2 (binomal/ poisson-bivariate)  ############################

nvar=15
n=99
Lt<-list()
tt=c(1:nvar)
for(j in 1:n){  Lt[[j]]=tt }


num_sim=100
yj = mixed=func =funHD =matrix(nrow= num_sim, ncol=2)
for(simm in 1:num_sim){
  
  nvar=15
  n=99
  Lt<-list()
  tt=c(1:nvar)
  for(j in 1:n){  Lt[[j]]=tt }
  
  print(paste(simm,'th simulation'))
  
  Ly<-list()
  type='pois1'
  Ly[[1]]=data_generation(n, u_sigma=0.1, type, tt)
  type='binom1'
  Ly[[2]]=data_generation(n, u_sigma=0.05, type, tt)


  dis=c('poisson(log)','binomial(logit)')
  
  num_clus=3
  cls=c(rep(c(1: num_clus),each=(99/num_clus)) ) 
  true=permn(num_clus)  
  
  ################ My method  ################
  result=FPCA_y(Ly, Lt, optnsSW=list(), K=num_clus, dis) 
  
  yj[simm, 1]=perc_correct(true, result)
  yj[simm, 2]=adjustedRandIndex(cls, result)  
  
  
  
  ################ MixAK  ################
  output1 <- as.vector(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = F))
  output2 <- as.vector(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = F))
  tt=as.vector(matrix(unlist(Lt), ncol = nvar, byrow = F))
  mod <- GLMM_MCMC(y = cbind(output1, output2),  dist = dis,  id = rep(c(1:99),each=nvar), 
                   x = list(  "empty", "empty"), z = list(tt, tt) ,  
                   prior.b = list(Kmax = num_clus), nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),  parallel = FALSE)
  
  fit_mix <- apply(mod[[1]]$poster.comp.prob_u, 1, which.max)
  
  
  mixed[simm, 1]= perc_correct(true, fit_mix)
  mixed[simm, 2]=adjustedRandIndex(cls, fit_mix) 
  
  
  
  
  ################ funclust  ################
  CWtime<- 1:nvar
  output1 <- t(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd1 <- smooth.basisPar( CWtime, output1,  lambda=1e-2)$fd
  output2 <- t(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd2 <- smooth.basisPar( CWtime, output2,  lambda=1e-2)$fd
  
  CWfd=list(CWfd1,CWfd2)
  
  res=funclust(CWfd,K=num_clus)
  fit_fun=res$cls
  func[simm, 1]= perc_correct(true, fit_fun)
  func[simm, 2]=adjustedRandIndex(cls, fit_fun) 
  
  ################ funHDDC ################
  
  res.multi<-funHDDC(CWfd,K=num_clus,model="AkjBQkDk",init="kmeans",threshold=0.2)
  funHDDC_result=res.multi$class
  funHD[simm, 1]= perc_correct(true, funHDDC_result)
  funHD[simm, 2]=adjustedRandIndex(cls, funHDDC_result) 
  
}



yj_out=round(apply(yj,2, mean, na.rm=T),2)
mixed_out=round(apply(mixed,2, mean, na.rm=T),2)
func_out=round(apply(func,2, mean, na.rm=T),2)
HD_out=round(apply(funHD,2, mean, na.rm=T),2)


yj_sd=round(apply(yj,2, sd, na.rm=T),2)
mixed_sd=round(apply(mixed,2, sd, na.rm=T),2)
func_sd=round(apply(func,2, sd, na.rm=T),2)
HD_sd=round(apply(funHD,2, sd, na.rm=T),2)




############################  Simulation 3 (4 variate)  ############################



num_sim=100
yj = mixed=func =funHD =matrix(nrow= num_sim, ncol=2)
for(simm in 15:num_sim){
  
  print(paste(simm,'th simulation'))
  
  
  nvar=15
  n=99
  Lt<-list()
  tt=c(1:nvar)
  for(j in 1:n){  Lt[[j]]=tt }
  Ly<-list()
  type='pois1'
  Ly[[1]]=data_generation(n, u_sigma=0.1, type, tt)
  
  type='pois2'
  Ly[[2]]=data_generation(n, u_sigma=0.1, type, tt)
  
  type='binom1'
  Ly[[3]]=data_generation(n, u_sigma=0.05, type, tt)
  
  type='binom2'
  Ly[[4]]=data_generation(n, u_sigma=0.05, type, tt)

  dis=c('poisson(log)', 'poisson(log)','binomial(logit)','binomial(logit)')
  
  num_clus=3
  cls=c(rep(c(1: num_clus),each=(99/num_clus)) ) 
  true=permn(num_clus)  
  
  ################ My method  ################
  result=FPCA_y(Ly, Lt, optnsSW=list(), K=num_clus, dis) 
  
  yj[simm, 1]=perc_correct(true, result)
  yj[simm, 2]=adjustedRandIndex(cls, result)  
  
  
  
  ################ MixAK  ################
  output1 <- as.vector(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = F))
  output2 <- as.vector(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = F))
  output3 <- as.vector(matrix(unlist(Ly[[3]]), ncol = nvar, byrow = F))
  output4 <- as.vector(matrix(unlist(Ly[[4]]), ncol = nvar, byrow = F))
  tt=as.vector(matrix(unlist(Lt), ncol = nvar, byrow = F))
  mod <- GLMM_MCMC(y = cbind(output1, output2,output3, output4),  dist = dis,  id = rep(c(1:99),each=nvar), 
                   x = list(  "empty", "empty", "empty", "empty"), z = list(tt, tt,tt,tt) ,  
                   prior.b = list(Kmax = num_clus), nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),  parallel = FALSE)
  
  fit_mix <- apply(mod[[1]]$poster.comp.prob_u, 1, which.max)
  
  
  mixed[simm, 1]= perc_correct(true, fit_mix)
  mixed[simm, 2]= adjustedRandIndex(cls, fit_mix) 
  
  
  
  
  ################ funclust  ################
  CWtime<- 1:nvar
  output1 <- t(matrix(unlist(Ly[[1]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd1 <- smooth.basisPar( CWtime, output1,  lambda=1e-2)$fd
  output2 <- t(matrix(unlist(Ly[[2]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd2 <- smooth.basisPar( CWtime, output2,  lambda=1e-2)$fd
  output3 <- t(matrix(unlist(Ly[[3]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd3 <- smooth.basisPar( CWtime, output3,  lambda=1e-2)$fd
  output4 <- t(matrix(unlist(Ly[[4]]), ncol = nvar, byrow = TRUE))  ## nrow=15, ncol=curve
  CWfd4 <- smooth.basisPar( CWtime, output4,  lambda=1e-2)$fd
  
  CWfd=list(CWfd1,CWfd2,CWfd3,CWfd4)
  
  res=funclust(CWfd,K=num_clus)
  fit_fun=res$cls
  func[simm, 1]= perc_correct(true, fit_fun)
  func[simm, 2]=adjustedRandIndex(cls, fit_fun) 
  
  ################ funHDDC ################
  
  res.multi<-funHDDC(CWfd,K=num_clus,model="AkjBQkDk",init="kmeans",threshold=0.2)
  funHDDC_result=res.multi$class
  funHD[simm, 1]= perc_correct(true, funHDDC_result)
  funHD[simm, 2]=adjustedRandIndex(cls, funHDDC_result) 
  
}



yj_out=round(apply(yj,2, mean, na.rm=T),2)
mixed_out=round(apply(mixed,2, mean, na.rm=T),2)
func_out=round(apply(func,2, mean, na.rm=T),2)
HD_out=round(apply(funHD,2, mean, na.rm=T),2)


yj_sd=round(apply(yj,2, sd, na.rm=T),2)
mixed_sd=round(apply(mixed,2, sd, na.rm=T),2)
func_sd=round(apply(func,2, sd, na.rm=T),2)
HD_sd=round(apply(funHD,2, sd, na.rm=T),2)



