##########################################################################################################
######## 			Generate simulation data 'N' times and perform proposed method and other methods		  			########
########   		(desginate N at code)																																########
######## 	   output :  mean and standard deviation of correct classification rate (CCR) and aRand values over N simulations. 	########
##########################################################################################################


library(cluster)	
library(fields)
library(mclust)
library(EbayesThresh)
library(waveslim)
library(curvclust)
library(pracma)
library(QuantPsyc)
library(combinat)
library(funcy)
library(funFEM)
library(gamlss)
library(extraDistr)
library(pscl)
library(gtools)
library(mixtools)


## load source functions
setwd("/Users/yaejilim/Desktop/Papers/JY_clustering/homepage/sending_code/")
source('functions.R')

#############################################################################
#########   Sinusoidal function   (different amount )   										  #########
#############################################################################

N <- 5  #number of times of simulation
bin3= bin2 = bin1 = km= kmh= funff=curv= kmeansEn = curvEn = funFEM= funFEMEn <- matrix(nrow = N, ncol = 2)

for(simm in 1:N) {
	
	
		  print(paste(simm,'th simulation'))
		
		### data generation  ###
		  num_clus=4
		  nn=1024	
		   cls = true=c(rep(c(1: num_clus),each=(200/num_clus)) )
		   
		  data <-matrix(nrow=nn ,ncol=200)
		  
		  t=c(1:nn)/(nn/5)
		  y=sin(t)
		  for(i in 1:200){
		    sigma=rep(c(0.1, 0.4, 0.7, 1), each=50)
		    data[,i]= abs( y+rnorm(nn, 0, sigma[i]))
		  }
		

		data = t(data)
		
		
		  ########## proposed method  - bin1 ########################
		  ### using starting point 1 for all bins
		  ####################################################################
		  print ('bin1')
		  res <- matrix(nrow = 6, ncol = 200)
		  ind <- 1
		  for(B in c(1024, 512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin(t(data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		
		  bin1[simm, 1]=ccr(result, true, 4)
		  bin1[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		
		
		
		
		  ########## proposed method  - bin2 ########################
		  ###  2 additional random starting points
		  ####################################################################
		  print ('bin2')
		  res <- matrix(nrow = 16, ncol = 200)
		  res[1,] <- kmeansBin(t(data), 1024, startPoint = 1, K = num_clus)
		  ind <- 2
		
		  for(B in c(512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin(t(data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		    sp_ind <- sample(B, 2)
		    res[ind,] <- kmeansBin(t(data), B, startPoint = sp_ind[1], K = num_clus)
		    ind <- ind + 1
		    res[ind,] <- kmeansBin(t(data), B, startPoint = sp_ind[2], K = num_clus)
		    ind <- ind + 1
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		
		  bin2[simm, 1]=ccr(result, true, 4)
		  bin2[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		
		  ########## proposed method  - bin3 ########################
		  res <- matrix(nrow = 9, ncol = 200)
		  ind <- 1
		  for(B in c(1024, 512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin(t(data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		  }
		  for(B in c(4, 2, 1)){
		    basis <- create.fourier.basis(c(0, 1024/B), nbasis=11)
		    fdobj <- smooth.basis(1:(1024/B), 
		                          (bin(t(data), B, startPoint = 1)),
		                          basis)$fd
		    result = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
		    res[ind,] <- result$cls
		    ind <- ind + 1
		  }
		  
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  bin3[simm, 1]=ccr(result, true, 4)
		  bin3[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		  ########## KMeans ########################
		  result=kmeans((data),num_clus )$cluster
		  km[simm, 1]=ccr(result, true, num_clus)
		  km[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		  		  ########## KMeans h=500 ########################
		  result=kmeansBin(t(data), 500, startPoint = 1, K = num_clus)
		  kmh[simm, 1]=ccr(result, true, num_clus)
		  kmh[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		    ######## K-means Ensemble #####################
		  print('K-means Ensemble')
		  #data1
		  res <- matrix(nrow = 15, ncol = 200)
		  ind <- 1
		  for(ind in 1:15) {
		    res[ind,] <- kmeans(data, centers = num_clus, nstart = 1)$cluster
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  kmeansEn[simm, 1]=ccr(result, true, 4)
		  kmeansEn[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		
		  ############# funFEM C. Bouveyron, E. Côme and J. Jacques,  ###########
		  basis <- create.fourier.basis(c(0, 1024), nbasis=11)
		  fdobj <- smooth.basis(1: 1024,t(data),basis)$fd
		  res = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
		  pam=res$cls
		  funFEM[simm, 1]=ccr(pam, true, num_clus)
		  funFEM[simm, 2]=round(adjustedRandIndex(true, pam)  ,4)
		  
		    ### funFEM Ensemble #####
		  print('funFEM Ensemble')
		  
		  res <- matrix(nrow = 12, ncol = 200)
		  modelList <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB",
		                 "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
		  for(ind in 1:12) {
		    basis <- create.fourier.basis(c(0, 1024), nbasis=11)
		    fdobj <- smooth.basis(1: 1024,t(data),basis)$fd
		    result = funFEM(fdobj,K= num_clus,model = modelList[ind],
		                    init="kmeans",lambda=0,disp=TRUE)
		    res[ind,] <- result$cls
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  funFEMEn[simm, 1]=ccr(result, true, num_clus) #0.83
		  funFEMEn[simm, 2]=round(adjustedRandIndex(true, result)  ,4) #0.669
		  
		
		#############  curvclust  ######################
		
			 print('cuvclust')
			fdat= list()
			K= num_clus
			 for (j in 1:ncol(t(data))) fdat[[j]] =t(data)[,j]
			  CCD = new("CClustData",Y=fdat,filter.number=1)
		 	 CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
		 	 CCO = new("CClustO")
		 	 #CCO["nbclust"] = K
		 	 #CCO["Gamma2.structure"] = "none"
		 	CCO@nbclust = K
		 	CCO@Gamma2.structure = "none"
		 	 CCR = getFCM(CCDred,CCO)
		
			cluster = apply(CCR@Tau,1,which.max)
			curv[simm, 1]=ccr(true, cluster, num_clus)
			curv[simm, 2]=( round( adjustedRandIndex(cls, cluster)  ,4))
		
		  
		  #############  curvclust Ensmeble ######################
		  
		  print('cuvclust Ensemble')
		  
		  #data1
		  fdat= list()
		  K= num_clus
		  for (j in 1:ncol(t(data))) fdat[[j]] =t(data)[,j]
		  
		  res <- matrix(nrow = 10, ncol = 200)
		  ind <- 1
		  for(ind in 1:10) {
		    CCD = new("CClustData",Y=fdat,filter.number=ind)
		    CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
		    CCO = new("CClustO")
		    #CCO["nbclust"] = K
		    #CCO["Gamma2.structure"] = "none"
		    CCO@nbclust = K
		    CCO@Gamma2.structure = "none"
		    CCR = getFCM(CCDred,CCO)
		    res[ind,] <- apply(CCR@Tau,1,which.max)
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  curvEn[simm, 1]=ccr(result, true, 4)
		  curvEn[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		

  
}  #end of simulation

output <- list(Kmeans = km, curvclust = curv, funFEM = funFEM,
                    bin1 = bin1, bin2 = bin2, bin3 = bin3, kmh=kmh
                    kmeansEn = kmeansEn, curvEn = curvEn, funFEMEn = funFEMEn)
                    
for(k in 1:2){
print(paste(
round(  apply(curvEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(curvEn,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(curv,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(curv,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(funFEMEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(funFEMEn,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(funFEM,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(funFEM,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(kmeansEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(kmeansEn,2, sd,na.rm=TRUE) ,3)[k]  , ') & '   , 
round(  apply(km,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(km,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(kmh,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(kmh,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(bin1,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin1 ,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(bin2,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin2 ,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(bin3, 2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin3 ,2, sd,na.rm=TRUE) ,3)[k]  , ') '  
)
)}



#############################################################################
#########   Block function   (different pattern )   										  #########
#############################################################################

N <- 50  #number of times of simulation
kmh= bin3= bin2 = bin1 = km=funff=curv= kmeansEn = curvEn = funFEM= funFEMEn <- matrix(nrow = N, ncol = 2)

for(simm in 1:N) {
		  print(paste(simm,'th simulation'))
		  
		  ####### data generation  ##########
		  num_clus=3
		  n=1024
		  sigma=3
		  data <-matrix(nrow=n ,ncol=150)
		  cls = true=c(rep(c(1: num_clus),each=(150/num_clus)) )
		  
		  x <- seq(1, n)/n
		  
		  for(j in 1:50){  
		    t <-(sample(   seq(1, n)/(3*n) , 5))
		    h1 <- sort(runif(5, 0,20))[ c( 3,2, 4 ,1,5)] * c(-1,1,-1,-1,1)
		    
		    blocks <- rep(0, n)
		    for (i in seq(1, length(h1))) {
		      blocks <- blocks + (h1[i] * (1 + sign(x - t[i])))/2
		    }
		    data[,j] = abs(blocks+rnorm(n, 0, sigma))
		  }
		  
		  for(j in 51:100){  
		    t <-(sample(   seq( round(n/3), round(2*n/3 ) ) /n  , 5))
		    h1 <- sort(runif(5, 0,20))[ c( 3,2, 4 ,1,5)] * c(-1,1,-1,-1,1)
		    
		    blocks <- rep(0, n)
		    for (i in seq(1, length(h1))) {
		      blocks <- blocks + (h1[i] * (1 + sign(x - t[i])))/2
		    }
		    data[,j] = abs(blocks+rnorm(n, 0, sigma))
		  }
		  
		  for(j in 101:150){  
		    t <-(sample(   seq( round(2*n/3), n ) /n  , 5))
		    h1 <- sort(runif(5, 0,20))[ c( 3,2, 4 ,1,5)] * c(-1,1,-1,-1,1)
		    
		    blocks <- rep(0, n)
		    for (i in seq(1, length(h1))) {
		      blocks <- blocks + (h1[i] * (1 + sign(x - t[i])))/2
		    }
		    data[,j] = abs(blocks+rnorm(n, 0, sigma))
		  }
		
		
		  ########## proposed method  - bin1 ########################
		  print ('bin1')
		  res <- matrix(nrow = 6, ncol = 150)
		  ind <- 1
		  for(B in c(1024, 512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin((data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		
		  bin1[simm, 1]=ccr(result, true, 4)
		  bin1[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		
		
		
		
		  ########## proposed method  - bin2 ########################
		  print ('bin2')
		  res <- matrix(nrow = 16, ncol = 150)
		  res[1,] <- kmeansBin((data), 1024, startPoint = 1, K = num_clus)
		  ind <- 2
		
		  for(B in c(512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin((data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		    sp_ind <- sample(B, 2)
		    res[ind,] <- kmeansBin((data), B, startPoint = sp_ind[1], K = num_clus)
		    ind <- ind + 1
		    res[ind,] <- kmeansBin((data), B, startPoint = sp_ind[2], K = num_clus)
		    ind <- ind + 1
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		
		  bin2[simm, 1]=ccr(result, true, 4)
		  bin2[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		
		  ########## proposed method  - bin3 ########################
		  res <- matrix(nrow = 9, ncol = 150)
		  ind <- 1
		  for(B in c(1024, 512, 256, 128, 64, 32)) {
		    res[ind,] <- kmeansBin((data), B, startPoint = 1, K = num_clus)
		    ind <- ind + 1
		  }
		  for(B in c(4, 2, 1)){
		    basis <- create.fourier.basis(c(0, 1024/B), nbasis=11)
		    fdobj <- smooth.basis(1:(1024/B), 
		                          (bin((data), B, startPoint = 1)),
		                          basis)$fd
		    result = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
		    res[ind,] <- result$cls
		    ind <- ind + 1
		  }
		  
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  bin3[simm, 1]=ccr(result, true, 4)
		  bin3[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		  ########## KMeans ########################
		  result=kmeans(t(data),num_clus )$cluster
		  km[simm, 1]=ccr(result, true, num_clus)
		  km[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		  		  		  ########## KMeans h=300 ########################
		  result=kmeansBin((data), 500, startPoint = 1, K = num_clus)
		  		  kmh[simm, 1]=ccr(result, true, num_clus)
		  kmh[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		    ######## K-means Ensemble #####################
		  print('K-means Ensemble')
		  #data1
		  res <- matrix(nrow = 15, ncol = 150)
		  ind <- 1
		  for(ind in 1:15) {
		    res[ind,] <- kmeans(t(data), centers = num_clus, nstart = 1)$cluster
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  kmeansEn[simm, 1]=ccr(result, true, 4)
		  kmeansEn[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))
		  
		
		  ############# funFEM C. Bouveyron, E. Côme and J. Jacques,#######
		  
		  basis <- create.fourier.basis(c(0, 1024), nbasis=11)
		  fdobj <- smooth.basis(1: 1024,(data),basis)$fd
		  res = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
		  pam=res$cls
		  funFEM[simm, 1]=ccr(pam, true, num_clus)
		  funFEM[simm, 2]=round(adjustedRandIndex(true, pam)  ,4)
		  
		    ### funFEM Ensemble #####
		  print('funFEM Ensemble')
		  
		  res <- matrix(nrow = 12, ncol = 150)
		  modelList <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB",
		                 "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
		  for(ind in 1:12) {
		    basis <- create.fourier.basis(c(0, 1024), nbasis=11)
		    fdobj <- smooth.basis(1: 1024,(data),basis)$fd
		    result = funFEM(fdobj,K= num_clus,model = modelList[ind],
		                    init="kmeans",lambda=0,disp=TRUE)
		    res[ind,] <- result$cls
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  funFEMEn[simm, 1]=ccr(result, true, num_clus) #0.83
		  funFEMEn[simm, 2]=round(adjustedRandIndex(true, result)  ,4) #0.669
		  
		
		#############  curvclust  ######################
		
			 print('cuvclust')
			fdat= list()
			K= num_clus
			 for (j in 1:ncol((data))) fdat[[j]] =(data)[,j]
			  CCD = new("CClustData",Y=fdat,filter.number=1)
		 	 CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
		 	 CCO = new("CClustO")
		 	 #CCO["nbclust"] = K
		 	 #CCO["Gamma2.structure"] = "none"
		 	CCO@nbclust = K
		 	CCO@Gamma2.structure = "none"
		 	 CCR = getFCM(CCDred,CCO)
		
			cluster = apply(CCR@Tau,1,which.max)
			curv[simm, 1]=ccr(true, cluster, num_clus)
			curv[simm, 2]=( round( adjustedRandIndex(cls, cluster)  ,4))
		
		  
		  #############  curvclust Ensmeble ######################
		  
		  print('cuvclust Ensemble')
		  
		  #data1
		  fdat= list()
		  K= num_clus
		  for (j in 1:ncol((data))) fdat[[j]] =(data)[,j]
		  
		  res <- matrix(nrow = 10, ncol = 150)
		  ind <- 1
		  for(ind in 1:10) {
		    CCD = new("CClustData",Y=fdat,filter.number=ind)
		    CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
		    CCO = new("CClustO")
		    #CCO["nbclust"] = K
		    #CCO["Gamma2.structure"] = "none"
		    CCO@nbclust = K
		    CCO@Gamma2.structure = "none"
		    CCR = getFCM(CCDred,CCO)
		    res[ind,] <- apply(CCR@Tau,1,which.max)
		  }
		  voteRes <- voting(res, num_clus)
		  result <- voteToCluster(voteRes)
		  
		  curvEn[simm, 1]=ccr(result, true, 4)
		  curvEn[simm, 2]=(round( adjustedRandIndex(cls, result) ,4))


  
}  #end of simulation

output <- list(Kmeans = km, curvclust = curv, funFEM = funFEM,
                    bin1 = bin1, bin2 = bin2, bin3 = bin3,
                    kmeansEn = kmeansEn, curvEn = curvEn, funFEMEn = funFEMEn)
                    






#############################################################################
#########   Count data simulation      									  #########
#############################################################################
 load('step_full.RData')
newStep= stepAll
#newStep <- read.csv("step.csv") #real data
fitPar <- matrix(ncol = 2, nrow = 1440)
colnames(fitPar) <- c("mu", "p")


for(i in 1:1440) {
  fitPar[i, ] <- fitZIP(as.numeric(newStep[i,]))
}



mixRes <- normalmixEM(fitPar[,1])
orderedMu <- sort(fitPar[,1])
mixResAmount1 <- normalmixEM(orderedMu[1:(1440*0.4)])
mixResAmount2 <- normalmixEM(orderedMu[(1440*0.2+1):(1440*0.6)])
mixResAmount3 <- normalmixEM(orderedMu[(1440*0.4+1):(1440*0.8)])
mixResAmount4 <- normalmixEM(orderedMu[(1440*0.6+1):(1440*1)])




kmh=kmeansEn = curvEn = funFEMEn =bin3 = bin1 = bin2 =Kmeans = curv = funFEM = matrix(nrow = 50, ncol = 6)

for(simm in 1:50) {
	
				      tryCatch({	
		  print(paste(simm,'th simulation'))
		  #############################################################################
		  #########   Count Time Series- Pattern(K=4) : #새벽활발(00-06) #아침활발(6-12) #오후활발(13-18) #저녁활발(18-24)							#########
		  ##############################################################################

		  N <- 1440
		  mu_p=0.99

		  components <- sample(1:2, prob = mixRes$lambda, size = N, replace = TRUE)
		  lambda <- rnorm(n = N, mean = mixRes$mu[components], sd = mixRes$sigma[components])
		  lambda <- ifelse(lambda < 0, 0.01, lambda)
		  
		  num_clus = 4
		  n = 1440
		  data <- matrix(nrow = n, ncol = 200)
		  true=c(rep(c(1: num_clus),each=(200/num_clus)) )
		  for(j in 1:50) {   #active during 1:480
		    p <- rnorm(n = N, mean = mu_p, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixRes$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixRes$mu[components]*2, sd = mixRes$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- c(seq(sample(1:116, size = 1), by = 1, length.out = 5),
		                seq(sample(121:236, size = 1), by = 1, length.out = 5),
		                seq(sample(241:356, size = 1), by = 1, length.out = 5),
		                seq(sample(360:476, size = 1), by = 1, length.out = 5))
		    data[signal, j] <- 50 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                min = quantile(lambda, 0.75),
		                                                                max = quantile(lambda, 0.95)))
		  }
		  for(j in 51:100) { #active during 320:800
		    p <- rnorm(n = N, mean = mu_p, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixRes$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixRes$mu[components]*2, sd = mixRes$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- c(seq(sample(320:476, size = 1), by = 1, length.out = 5),
		                seq(sample(480:596, size = 1), by = 1, length.out = 5),
		                seq(sample(600:716, size = 1), by = 1, length.out = 5),
		                seq(sample(720:796, size = 1), by = 1, length.out = 5))
		    data[signal, j] <- 50 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		    
		  }
		  for(j in 101:150) { #active during 640:1120
		    p <- rnorm(n = N, mean = mu_p, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixRes$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixRes$mu[components]*2, sd = mixRes$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- c(seq(sample(640:756, size = 1), by = 1, length.out = 5),
		                seq(sample(760:876, size = 1), by = 1, length.out = 5),
		                seq(sample(880:996, size = 1), by = 1, length.out = 5),
		                seq(sample(1000:1116, size = 1), by = 1, length.out = 5))
		    data[signal, j] <- 50 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		  }
		  for(j in 151:200) { #active during 960:1440
		    p <- rnorm(n = N, mean = mu_p, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixRes$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixRes$mu[components]*2, sd = mixRes$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- c(seq(sample(960:1076, size = 1), by = 1, length.out = 5),
		                seq(sample(1080:1196, size = 1), by = 1, length.out = 5),
		                seq(sample(1200:1316, size = 1), by = 1, length.out = 5),
		                seq(sample(1320:1436, size = 1), by = 1, length.out = 5))
		    data[signal, j] <- 50 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		  }
		  data <- t(data)
		  
		  #############################################################################
		  #########   Count Time Series  - Amount(K=4): Count mean 차이로 4개			      #########
		  ##############################################################################
		
		  num_clus = 4
		  n = 1440
		  mu_p=.99
		  data2 <- matrix(nrow = n, ncol = 200)
		  cls = true=c(rep(c(1: num_clus),each=(200/num_clus)) )
		  for(j in 1:50) {   #lowest activity
		    p <- rnorm(n = N, mean = mu_p, sd = 0.05)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount1$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount1$mu[components], 
		                    sd = mixResAmount1$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data2[,j] <- rzip(n = n, lambda = lambda, pi = p)
		  }
		  for(j in 51:100) { #2nd lowest
		    p <- rnorm(n = N, mean = mu_p, sd = 0.05)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount2$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount2$mu[components]*1, sd = mixResAmount2$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data2[,j] <- rzip(n = n, lambda = lambda, pi = p)
		  }
		  for(j in 101:150) { #3rd lowest
		    p <- rnorm(n = N, mean = mu_p, sd = 0.05)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount3$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount3$mu[components]*1, sd = mixResAmount3$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data2[,j] <- rzip(n = n, lambda = lambda, pi = p)
		  }
		  for(j in 151:200) { #most active
		    p <- rnorm(n = N, mean = mu_p, sd = 0.05)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount4$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount4$mu[components]*1, sd = mixResAmount4$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data2[,j] <- rzip(n = n, lambda = lambda, pi = p)
		  }
		
		  data2 <- t(data2)
		  
		  #############################################################################
		  #########   Count Time Series  - Pattern + Amount (K=4) #저활동군 #아침활발 #저녁활발 #전체활발#######
		  ##############################################################################
		
		  num_clus = 4
		  n = 1440
		  data3 <- matrix(nrow = n, ncol = 200)
		  true=c(rep(c(1: num_clus),each=(200/num_clus)) )
		  for(j in 1:50) {   #lowest activity
		    p <- rnorm(n = N, mean = 0.99, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount1$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount1$mu[components]*1, 
		                    sd = mixResAmount1$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data3[,j] <- rzip(n = n, lambda = lambda, pi = p)
		  }
		  for(j in 51:100) { #Active Morning 8:14 (481:840)
		    p <- rnorm(n = N, mean = 0.99, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount2$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount2$mu[components]*2, sd = mixResAmount2$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data3[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- unique(c(seq(sample(481:540, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(541:600, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(601:660, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(661:720, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(720:780, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(781:840, size = 1), by = 1, length.out = sample(3:10, size = 1))))
		    data3[signal, j] <- 100 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		  }
		  for(j in 101:150) { #Active 18:24 (1080:1440)
		    p <- rnorm(n = N, mean = 0.99, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount3$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount3$mu[components]*4, sd = mixResAmount3$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data3[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- unique(c(seq(sample(1081:1140, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(1141:1200, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(1201:1260, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(1261:1320, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(1321:1380, size = 1), by = 1, length.out = sample(3:10, size = 1)),
		                       seq(sample(1381:1431, size = 1), by = 1, length.out = sample(3:10, size = 1))))
		    data3[signal, j] <- 100 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		    
		  }
		  for(j in 151:200) { #most active (480:1440)
		    p <- rnorm(n = N, mean = 0.99, sd = 0.1)
		    p <- ifelse(p >= 1, 0.99, p); 		   p <- ifelse(p <= 0, 0, p)
		    components <- sample(1:2, prob = mixResAmount4$lambda, size = N, replace = TRUE)
		    lambda <- rnorm(n = N, mean = mixResAmount4$mu[components]*5, sd = mixResAmount4$sigma[components])
		    lambda <- ifelse(lambda < 0, 0.01, lambda)
		    data3[,j] <- rzip(n = n, lambda = lambda, pi = p)
		    signal <- unique(c(seq(sample(481:540, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(541:600, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(601:660, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(661:720, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(720:780, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(781:840, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(841:900, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(901:960, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(961:1020, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1021:1080, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1081:1140, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1141:1200, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1201:1260, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1261:1320, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1321:1380, size = 1), by = 1, length.out = sample(1:10, size = 1)),
		                       seq(sample(1381:1431, size = 1), by = 1, length.out = sample(1:10, size = 1))))
		    data3[signal, j] <- 100 + rpois(n = length(signal), lambda = runif(length(signal),
		                                                                     min = quantile(lambda, 0.75),
		                                                                     max = quantile(lambda, 0.95)))
		  }
		  data3 <- t(data3)
		  
		full_data=list(data=data, data2=data2, data3=data3)
		

		
		  ####### Binning 1 ########
		  
		  for(mm in 1:3){
			  	use_data=full_data[[mm]]
				 print ('bin1')
			  res <- matrix(nrow = 13, ncol = 200)
			  ind <- 1
			  for(B in c(1440, 720, 360, 240, 180, 144, 120, 90, 72, 60, 40, 24, 20)) {
			    res[ind,] <- kmeansBin(t(use_data), B, startPoint = 1, K = num_clus)
			    ind <- ind + 1
			  }
			  voteRes <- voting(res, num_clus)
			  result <- voteToCluster(voteRes)
			
			  bin1[simm, (2*(mm-1)+1 )]=ccr(result, true, 4)
			  bin1[simm, (2*(mm-1)+2 )]=(round( adjustedRandIndex(cls, result) ,4))
		
			}
		
		
		  ####### Binning 2 ########
		  for(mm in 1:3){
			  		use_data=full_data[[mm]]
			 	 print ('bin2')
			  res <- matrix(nrow = 40, ncol = 200)
			  res[1,] <- kmeansBin(t(use_data), 1440, startPoint = 1, K = num_clus)
			  ind <- 2
			
			  for(B in c(720, 360, 240, 180, 144, 120, 90, 72, 60, 40, 24, 20)) {
			    res[ind,] <- kmeansBin(t(use_data), B, startPoint = 1, K = num_clus)
			    ind <- ind + 1
			    sp_ind <- sample(B, 2)
			    res[ind,] <- kmeansBin(t(use_data), B, startPoint = sp_ind[1], K = num_clus)
			    ind <- ind + 1
			    res[ind,] <- kmeansBin(t(use_data), B, startPoint = sp_ind[2], K = num_clus)
			    ind <- ind + 1
			  }
			  voteRes <- voting(res, num_clus)
			  result <- voteToCluster(voteRes)
			
			  bin2[simm,  (2*(mm-1)+1 )]=ccr(result, true, 4)
			  bin2[simm,  (2*(mm-1)+2 )]=(round( adjustedRandIndex(cls, result) ,4))
		
			}
		  ####### Binning 3 ########
		  
		  for(mm in 1:3){
		  	use_data=full_data[[mm]]
			  res <- matrix(nrow = 16, ncol = 200)
			  ind <- 1
			  for(B in c(1440, 720, 360, 240, 180, 144, 120, 90, 72, 60, 40, 24, 20)) {
			    res[ind,] <- kmeansBin(t(use_data), B, startPoint = 1, K = num_clus)
			    ind <- ind + 1
			  }
			  for(B in c(3, 2, 1)){
			    basis <- create.fourier.basis(c(0, 1440/B), nbasis=11)
			    fdobj <- smooth.basis(1:(1440/B),
			                          (bin(t(use_data), B, startPoint = 1)),
			                          basis)$fd
			    result = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
			    res[ind,] <- result$cls
			    ind <- ind + 1
			  }
			
			  voteRes <- voting(res, num_clus)
			  result <- voteToCluster(voteRes)
			
			  bin3[simm, (2*(mm-1)+1 ) ]=ccr(result, true, 4)
			  bin3[simm, (2*(mm-1)+2 )]=(round( adjustedRandIndex(cls, result) ,4))
		
		}
		
		
		
		  ########## K-means #####################
		  print('K-means')
		    for(mm in 1:3){
			  	use_data=full_data[[mm]]
			  	cluster = kmeans(use_data, centers = num_clus, nstart = 10)$cluster
			  	Kmeans[simm, (2*(mm-1)+1 )]=ccr(cluster, true, num_clus)
			 	 Kmeans[simm, (2*(mm-1)+2 )]=round( adjustedRandIndex(true, cluster)  ,4)
			 	 
			}
			
			
	 ########## KMeans h=300 ########################
	 
	    for(mm in 1:3){
			  	use_data=full_data[[mm]]
		 		 result=kmeansBin(t(use_data), 300, startPoint = 1, K = num_clus)
		  		kmh[simm, (2*(mm-1)+1 )]=ccr(result, true, num_clus)
			 	 kmh[simm, (2*(mm-1)+2 )]=round( adjustedRandIndex(true, result)  ,4)
			 	 print(ccr(result, true, num_clus) )
		  }
		  
		  #############  curvclust  ######################
		
		  print('cuvclust')
		   for(mm in 1:3){
		   	  use_data=full_data[[mm]]
			  fdat= list()
			  K= num_clus
			  for (j in 1:ncol(t(use_data))) fdat[[j]] =t(use_data)[,j]
			  CCD = new("CClustData",Y=fdat,filter.number=4)
			  CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
			  CCO = new("CClustO")
			  CCO@nbclust = K
			  CCO@Gamma2.structure = "none"
			  CCR = getFCM(CCDred,CCO)
			
			  cluster = apply(CCR@Tau,1,which.max)
			  curv[simm, (2*(mm-1)+1 )]=ccr(cluster, true, num_clus)
			  curv[simm, (2*(mm-1)+2 )]=round( adjustedRandIndex(true, cluster)  ,4)
		
			}
		
		  ############ funFEM C. Bouveyron, E. Côme and J. Jacques,
		   for(mm in 1:3){
		   	   	  use_data=full_data[[mm]]
				  basis <- create.fourier.basis(c(0, 1440), nbasis=11)
				  fdobj <- smooth.basis(1: 1440,t(use_data),basis)$fd
				  res = funFEM(fdobj,K= num_clus,model="AkjBk",init="kmeans",lambda=0,disp=TRUE)
				  pam=res$cls
				  funFEM[simm, (2*(mm-1)+1 )]=ccr(pam, true, num_clus) #0.83
				  funFEM[simm, (2*(mm-1)+2 )]=round(adjustedRandIndex(true, pam)  ,4) #0.669
			}
		  
		  ######### K-means Ensemble #####################
		  print('K-means Ensemble')
		   for(mm in 1:3){
		   	   	  use_data=full_data[[mm]]
				  res <- matrix(nrow = 15, ncol = 200)
				  ind <- 1
				  for(ind in 1:15) {
				    res[ind,] <- kmeans(use_data, centers = num_clus, nstart = 1)$cluster
				  }
				  voteRes <- voting(res, num_clus)
				  result <- voteToCluster(voteRes)
				
				  kmeansEn[simm, (2*(mm-1)+1 )]=ccr(result, true, 4)
				  kmeansEn[simm, (2*(mm-1)+2 )]=(round( adjustedRandIndex(cls, result) ,4))
		
		 	}
		
		
		  #############  curvclust Ensmeble ######################
		
		  print('cuvclust Ensemble')
		   for(mm in 1:3){
		   	   	   	use_data=full_data[[mm]]
				  fdat= list()
				  K= num_clus
				  for (j in 1:ncol(t(use_data))) fdat[[j]] =t(use_data)[,j]
				
				  res <- matrix(nrow = 10, ncol = 200)
				  ind <- 1
				  for(ind in 1:10) {
				    CCD = new("CClustData",Y=fdat,filter.number=ind)
				    CCDred = getUnionCoef(CCD)  #Dimension reduction is then performed
				    CCO = new("CClustO")
				    #CCO["nbclust"] = K
				    #CCO["Gamma2.structure"] = "none"
				    CCO@nbclust = K
				    CCO@Gamma2.structure = "none"
				    CCR = getFCM(CCDred,CCO)
				    res[ind,] <- apply(CCR@Tau,1,which.max)
				  }
				  voteRes <- voting(res, num_clus)
				  result <- voteToCluster(voteRes)
				
				  curvEn[simm,  (2*(mm-1)+1 )]=ccr(result, true, 4)
				  curvEn[simm,  (2*(mm-1)+2 )]=(round( adjustedRandIndex(cls, result) ,4))
		}  
		  ### funFEM Ensemble #####
		 	 print('funFEM Ensemble')
		   for(mm in 1:3){
		   	   	   	use_data=full_data[[mm]]
				  res <- matrix(nrow = 12, ncol = 200)
				  ind <- 1
				  modelList <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB",
				                 "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
				  for(ind in 1:12) {
				    basis <- create.fourier.basis(c(0, 1440), nbasis=11)
				    fdobj <- smooth.basis(1: 1440,t(use_data),basis)$fd
				    result = funFEM(fdobj,K= num_clus,model = modelList[ind],
				                    init="kmeans",lambda=0,disp=TRUE)
				    res[ind,] <- result$cls
				  }
				  voteRes <- voting(res, num_clus)
				  result <- voteToCluster(voteRes)
				
				  funFEMEn[simm, (2*(mm-1)+1 )]=ccr(result, true, num_clus)
				  funFEMEn[simm, (2*(mm-1)+2 )]=round(adjustedRandIndex(true, result)  ,4)
			}
			
			
		        
 }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}  #end of simulation


for(k in 5:6){
print(paste(
round(  apply(curvEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(curvEn,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(curv,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(curv,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(funFEMEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(funFEMEn,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(funFEM,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(funFEM,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(kmeansEn,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(kmeansEn,2, sd,na.rm=TRUE) ,3)[k]  , ') & '   , 
round(  apply(Kmeans,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(Kmeans,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(kmh,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(kmh,2, sd,na.rm=TRUE) ,3)[k]  , ') &'  , 
round(  apply(bin1,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin1 ,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(bin2,2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin2 ,2, sd,na.rm=TRUE) ,3)[k]  , ') & '  , 
round(  apply(bin3, 2, mean,na.rm=TRUE) ,3)[k]  , '(', round(  apply(bin3 ,2, sd,na.rm=TRUE) ,3)[k]  , ') '  
)
)}






