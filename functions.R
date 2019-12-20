library(pacman)
library(MASS)
library(fdapace)
library(cluster)	
library(fields)
library(modelcf)
library(mclust)
library(EbayesThresh)
library(waveslim)
library(curvclust)
library(pracma)
library(QuantPsyc)
library(combinat)
library(funcy)
library(funFEM)
library(optCluster)
library(mixAK)
library(funHDDC)

data_generation<-function(n, u_sigma, type, tt){
  if(type=='pois1'){
  Ly<-list()
  u<-rnorm(n, 0,u_sigma)
  for(j in 1:33){  Ly[[j]]=rpois(nvar, exp(u[j]) ) }
  for(j in 34:66){ 
    Ly[[j]]<-vector() 
    for(t in 1:length(tt)){
      Ly[[j]][t]=rpois(1, exp(0.08*tt[t] +  u[j]) ) 
    }
  }
  
  for(j in 67:99){ 
    Ly[[j]]<-vector() 
    for(t in 1:length(tt)){
      Ly[[j]][t]=rpois(1, exp(-0.2*tt[t]+  u[j]) ) 
    }
  }
  
}

  if(type=='pois2'){   
    # group1 
    Ly<-list()
    u<-rnorm(n, 0,u_sigma)
    for(j in 1:33){  
      for(t in 1:length(tt)){
        Ly[[j]]=rpois(15, exp(1- tt[t]/15 + u[j]) )		
      }}
    
    
    for(j in 34:66){ 
      Ly[[j]]<-vector() 
      for(t in 1:length(tt)){
        Ly[[j]][t]=rpois(1, exp(1-0.2*tt[t] +  u[j]) ) 
      }
    }
    
    
    for(j in 67:99){ 
      Ly[[j]]<-vector() 
      for(t in 1:length(tt)){
        Ly[[j]][t]=rpois(1, exp(1 -0.2*tt[t]+0.25* ( tt[t] -7) *I(tt[t] >= 7) +   u[j]) ) 
      }
    }
    
    
  }  


  if(type=='binom1'){  
    # group1 
    Ly<-list()
    u<-rnorm(n, 0,u_sigma)
    for(j in 1:33){  
      Ly[[j]]<-vector()
      Ly[[j]]=rbinom(nvar, 1, exp(tt/15  + u[j])/nvar )		#rbinom(1, nvar, exp(tt[k]/15 + u[j])/nvar )		
      
    }
    
    
    for(j in 34:66){ 
      Ly[[j]]<-vector() 
      Ly[[j]]=rbinom(nvar, 1,exp(-tt/15 +  u[j]) /nvar )   #rbinom(1, nvar,exp(-tt/15 +  u[j])/nvar ) 
      
    }
    
    
    for(j in 67:99){ 
      Ly[[j]]<-vector() 
      
      Ly[[j]]=  rbinom(nvar, 1, exp(-4 *tt/15 +u[j])/nvar)  # rbinom(1, nvar,exp(-4*tt/15 +   u[j])/nvar ) 
      
    }
  }

  
  if(type=='binom2'){  
    # group1 
    Ly<-list()
    u<-rnorm(n, 0,u_sigma)
    for(j in 1:33){  
      Ly[[j]]<-vector()
      max.value=max(exp(2+ u[j])   )+1
      Ly[[j]]=rbinom(nvar, 1, exp(2  + u[j])/nvar )		#rbinom(1, nvar, exp(tt[k]/15 + u[j])/nvar )		
      
    }
    
    
    for(j in 34:66){ 
      Ly[[j]]<-vector() 
      # max.value=max(exp(2-0.23* (tt-5)* I(tt>=5)  +  u[j])   )+1
      Ly[[j]]=rbinom(nvar, 1,exp(2-0.23* (tt-5)* I(tt>=5)  +  u[j]) /nvar )   #rbinom(1, nvar,exp(-tt/15 +  u[j])/nvar ) 
      
    }
    
    
    for(j in 67:99){ 
      Ly[[j]]<-vector() 
      
      Ly[[j]]=  rbinom(nvar, 1, exp(2.1 -0.1*tt  +0.13 *(tt-11)*I(tt>=11))/nvar)  # rbinom(1, nvar,exp(-4*tt/15 +   u[j])/nvar ) 
      
    }
  }
  
  
  return(Ly)
}



source('~/Desktop/Papers/acc_clustering(JRSS_C)/R-code(for sending)/functions.R', chdir = TRUE)

pos= function(xx){
  y<-vector()
  for(i in 1: length(xx)){
    x=xx[i]
  if( x >=0){
    y[i]=x}
  if (x<0){
    y[i]=0
  }  
  }
  return(y)
}
hidden <- setdiff(p_funs("fdapace", TRUE), p_funs("fdapace"))

invisible(lapply(hidden, function(x) {

    a <- strtrim(x, 1) == "%" 
    b <- substring(x, nchar(x)) == "%"

    if (a && b) {
        x2 <- paste0("`", x, "`")
    } else {
        x2 <- x
    }

    assign(x, eval(parse(text=paste0("fdapace", ":::", x2))), 
        envir = .GlobalEnv)
}))






########################################################################### 
###################### functions   ####################################
############################################################################

  
 FPCA_y<-function (Ly, Lt, optnsSW, K=3, dis) {
   
   	maxIter=20
    #optnsSW = list(methodMuCovEst = "smooth",  FVEthreshold = 0.9, methodBwCov = "GCV", methodBwMu = "GCV", verbose=FALSE)

    #initialClustering <- kmeans(fpcaObjY$xiEst, centers = K,  algorithm = "MacQueen")
   	
   #	for(k in 1:length(Ly)){
   	#  Ly[[k]]=Ly[[k]][which(is.na(Ly[[k]])==FALSE)  ]  
   	#  Lt[[k]]=Lt[[k]][which(is.na(Ly[[k]])==FALSE)  ]  
   	#}

   	output1 <- as.vector(unlist(Ly[[1]]))
   	tt=as.vector(unlist(Lt))
   	mod <- GLMM_MCMC(y = (output1),  dist = dis[1],  id = rep(c(1:length(Ly[[1]])), times=unlist(lapply(Ly[[1]],function(x) length(x)))  ), 
   	                 x = list(  "empty"), z = list(tt) ,  
   	                 prior.b = list(Kmax = K), nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),  parallel = FALSE)
   	
   	fit_mix <- apply(mod[[1]]$poster.comp.prob_u, 1, which.max)
   # mod<-kCFC(Ly[[1]], Lt, k= K)
   	clustConf0<-  fit_mix#mod$cluster
    ClustIds1<-ClustIds0 <-  lapply(unique(clustConf0), function(u) which(clustConf0 ==  u)  )
 	  
    clustlist<-list()
    clustlist[[1]]= clustConf0
    j=1
    if (  min(table(clustConf0))  <  3 | length(unique(clustlist[[j]]) ) < K   ) { 
      print('Stop at j=1 & unbalabced result')
      } 	
    if (   min(table(clustConf0))  >=  3  && length(unique(clustlist[[j]]) ) >= K  ) { 
   
      for (j in 2:(maxIter)) {
  			print(paste(j,'th iteration'))
        ref_hatyc<-list()
        for(u in 1:K){
          uu= ClustIds1[[u]]
          yy<-list()
          for(p in 1:length(Ly)){
             yy[[p]]= (lapply ( uu,  function(u) Ly[[p]][[u]]))
          }
          tt= (lapply ( uu,  function(u) Lt[[u]]))
          ref_hatyc[[u]]= FPCA_yj(yy,tt ,optnsSW=list(maxK=2), dis)  ## return fitted X
        }
        
  		  for(m in 1:length(Ly) ){## index for curve
         # print(m)
  		    ith_yy<-ith_tt<-list()
  		    for(p in 1:length(Ly)){
		    	 ith_yy[[p]]=Ly[[p]][[m]]
		  	   ith_tt[[p]]=Lt[[m]]
  		    }

		  	   hatyc<-list()
		     
				  for(u in 1:K){	## index for cluster group
				  	
				    	uu= ClustIds1[[u]]  ## index of curves in cluster u
				    	
					   if( sum(m == uu)== 1  ){ ## If m-th curve is in the group u
					     
					     yy<-list()
					     for(p in 1:length(Ly)){
					       yy[[p]]= (lapply ( uu,  function(u) Ly[[p]][[u]]))
					       yy[[p]][[ which(m==uu)  ]]=NULL
					     }
					     tt=(lapply ( uu,  function(u) Lt[[u]]))
					     tt[[ which(m==uu)  ]]=NULL
					     
							  temp_r= FPCA_yj(yy,tt ,optnsSW=list(maxK=2), dis)  ## return fitted X
							  hatyc[[u]]= predicted_yhat(yy ,ith_yy,  temp_r$mux,tt ,tt[[ which(m==uu)  ]],  temp_r$pi,temp_r$eigvalc, temp_r$muxtemp, temp_r$muytemp, temp_r$gamma, dis, temp_r$workGrid, temp_r$n_M, Lt)
		  
						}
					
					  if( sum(m == uu)== 0  ){## If m-th curve is not in the group u
					    yy<-list()
					    for(p in 1:length(Ly)){
					      yy[[p]]= (lapply ( uu,  function(u) Ly[[p]][[u]]))
					    }
					    tt= (lapply ( uu,  function(u) Lt[[u]]))
               hatyc[[u]]= predicted_yhat(yy,ith_yy,  ref_hatyc[[u]]$mux,tt ,Lt[[ m ]]  ,ref_hatyc[[u]]$pi,  ref_hatyc[[u]]$eigvalc,  ref_hatyc[[u]]$muxtemp, ref_hatyc[[u]]$muytemp, 
                                          ref_hatyc[[u]]$gamma, dis, workGrid=ref_hatyc[[u]]$workGrid, ref_hatyc[[u]]$n_M, Lt)
					   }
				    
				   } 
				    seCost<-vector()
				    for(u in 1:K){
				      y<-list()
				      for(p in 1:length(Ly)){
					     y[[p]]= (glink(hatyc[[u]][[p]][,1], dis[p] )-  (  ith_yy[[p]] )  )^2
				      }
					    seCost[u] <- sum(Reduce("+",  y) ) #   trapzRcpp(X = fpcaObjY$workGrid, Y = y)
				    }
				     clustConf0[m]=which.min(seCost)
				     ClustIds1<- lapply(unique(clustConf0), function(u) which(clustConf0 ==  u)  )	
				   #  if (   min(summary(ClustIds1))   <  3  ) { break} 				      
		  }			       
		      clustlist[[j]]= clustConf0
		       ##break information
		       
		 if( length(unique(clustlist[[j]]) ) < K  ){ 
		 	print("num of group is less than K . Error!")
		 	break
		 	}
		      
		if ( (j > 1) && any(sapply(clustlist[1:(j - 1)], function(u) all(u == clustlist[[j]])))  ) {
		   print("We made it!")
       break
      }
    if (   min(table(clustlist[[j]]))  <  3  ) { 
        print('Less than 3 curve in a group. Error')
        break} 				
      }	
    }
  return(clustlist[[j]])
    
}


glink<-function(alpha,family){
 
    if(family=='binomial(logit)'){
        ginv = exp(alpha)/ (1+exp(alpha))
        }
    if(family=='poisson(log)'){
        ginv = exp(alpha)
    }
  if(family=='gaussian'){
    ginv = (alpha)
  }
	return(ginv)

}



glink_inv<-function(alpha,family){
 
    if(family=='binomial(logit)'){
        alpha[alpha <=0] = exp(-5)
        alpha[alpha>=1] = 1-exp(-5)
        ginv = log(alpha/(1-alpha))
        }
    if(family=='poisson(log)'){
        alpha[alpha <=0] = exp(-5)
        ginv = log(alpha)
    }
  
  if(family=='gaussian'){
    ginv = (alpha)
  }
	return(ginv)
}



glink_der<-function(x,family){
 
    if(family=='binomial(logit)'){
        gder = exp(x)/((1+exp(x))^2)
        }
    if(family=='poisson(log)'){
        gder = exp(x)
    }
  
  if(family=='gaussian'){
    gder = rep(1, length(x))
  }
  
	return(gder)
}





FPCA_yj<-function(yy,tt, optnsSW, dis){
	
  mux_c=pi_c=eigvalc_c=  muxtemp_c= muytemp_c= error0_c=workGrid_c= gamma_c=n_M=list()
  
  ########################################
	########  compute mean & var of Y   ########
	########################################
	for(j in 1:length(yy)){
   res <- FPCA(yy[[j]], tt, optns=optnsSW)
								            
  	muy=res$mu           
	  covy =res$smoothedCov
	
  ########################################
	########  compute mean & var of X   ########
	########################################
	
  	mux=glink_inv( muy, dis[j] )
  	covx = diag(1/glink_der(mux,dis[j]))%*%covy%*%diag(1/glink_der(mux,dis[j]))
  	covx = (covx+t(covx))/2
	
	 ########################################
	########  PCA on X              ########
	########################################
	
	   Cmat=covx
     result <- eigen(Cmat)
     eigvalc <- result$values
 
     eigvecc <- as.matrix(result$vectors)
     sumvecc <- apply(eigvecc, 2, sum)
     eigvecc[, sumvecc < 0] <- -eigvecc[, sumvecc < 0]
    
    
        phi <- apply(eigvecc, 2, function(x) {
     				  	 x <- x/sqrt(trapzRcpp( res $workGrid, x^2))
			   return(x)
        				 })
    
    
 ## est    xiEst 
 
     muxtemp= mux
     muytemp=muy
     
     
     M<-gamma2<-vector()
     M[1]=2
     gamma2[1]= gamma_cv(yy[[j]], tt,  mux,  muy  , phi, eigvalc, M=M[1], dis[j], res$workGrid)
     
     l=2
     repeat{
     	M[l]=no_eigen(yy[[j]],tt, mux,  muy  , phi, eigvalc , gamma=gamma2[l-1] ,  maxk=5, dis[j], res$workGrid)
    	 gamma2[l]= gamma_cv(yy[[j]],tt,   mux,  muy  , phi, eigvalc, M[l],dis[j],res$workGrid)
    	 
    	 if(M[l]==M[l-1]  | l==10){break}	 
    	 l=l+1
	}      
	 gamma2= gamma2[l]#;print(gamma2)
	 M=M[l]
	 print(M)
   eigvalc= eigvalc[1:M]
  	pi= phi[,1:M]
	 n_M[[j]]=M
     error0 = gamma2*diag(muytemp)
     
     mux_c[[j]]=mux
     pi_c[[j]]=pi
     eigvalc_c[[j]]=eigvalc
     muxtemp_c[[j]]=muxtemp
     muytemp_c[[j]]=muytemp
     error0_c[[j]]=error0
     workGrid_c[[j]]=res$workGrid
     gamma_c[[j]]=gamma2

	}
  
     
	(list(mux=mux_c, pi=pi_c, eigvalc=eigvalc_c, muxtemp=muxtemp_c,muytemp=muytemp_c, error0=error0_c, workGrid=workGrid_c, gamma=gamma_c,n_M=n_M))					
}	




predicted_yhat<-function(yy, ith_yy, mux, t,target_t, pi, eigvalc, muxtemp,muytemp, gamma, dis, workGrid,n_M, Lt){ ## y do not have ith
  ## yy : curves in group c
  ## i th yy : curve that I want to predict
  ## mux, muxtemp, ,... : information from yy (i.e. curves in group c)
 
  
#  mux=temp_r$mux;t=tt ;target_t=tt[[ which(m==uu)  ]];  pi=temp_r$pi;eigvalc=temp_r$eigvalc; muxtemp=temp_r$muxtemp; muytemp=temp_r$muytemp; 
  #gamma=temp_r$gamma; workGrid=temp_r$workGrid;n_M= temp_r$n_M 
  
  
  xiEst<-fittedX<-xiEst_t2<-list()
 for(j in 1:length(yy)){
   muxtemp=spline(workGrid[[j]], mux[[j]], xout=target_t)$y
   muytemp2=spline(workGrid[[j]],muytemp[[j]], xout=target_t)$y
  
   pi_temp<-matrix(nrow=length(muytemp2), ncol=ncol(pi[[j]]) )
  for(l in 1:ncol(pi[[j]])){
     pi_temp[,l]=spline(workGrid[[j]] ,pi[[j]][, l], xout=target_t)$y
   }
  
  error0 = gamma[[j]]*diag(muytemp2, nrow=length(muytemp2))
  
  dev_data<-list()
  for(k in 1:length(yy[[j]])){
    muytemp22=spline(workGrid[[j]],muytemp[[j]], xout=t[[k]])$y
    dev_data[[k]]=(yy[[j]][[k]] - (muytemp22))   ## y_i - g(mu)
  }
 
  output= matrix(unlist(Ly[[1]]), ncol = length(target_t), byrow = F)
   xiEst_t=  diag(eigvalc[[j]])%*%t(pi_temp)%*%diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) ) %*%ginv(diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) )%*%pi_temp%*%diag(eigvalc[[j]])%*%t(pi_temp)
            %*%diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) )+error0)  %*%t(output ) ## (N-1)*n_M  (for all curves in group c)

  
  xiEst_t2[[j]]=  diag(eigvalc[[j]])%*%t(pi_temp)%*%diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) ) %*%ginv(diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) )%*%pi_temp%*%diag(eigvalc[[j]])%*%t(pi_temp)
               %*%diag(glink_der(muxtemp,dis[j]), nrow=length(muxtemp) )+error0) %*% (ith_yy[[j]] -muytemp2)  ## xi_i (for i th patient) - length n_M
  xiEst[[j]]=t(xiEst_t)
  
 }
  
  full_xiEst<-vector()
  for(p in 1:length(yy)){
   full_xiEst=cbind(full_xiEst, xiEst[[p]])
  }
  
  Z= t(full_xiEst) %*% full_xiEst / (nrow(full_xiEst)-1)
  

  index=0
  for(j in 1: length(yy)   ){
    muxtemp=spline(workGrid[[1]], mux[[j]], xout=Lt[[j]])$y
    index=seq(from=(max(index)+1 ) , length.out=n_M[[j]])
   
     pi_temp<-matrix(nrow=length(muytemp2), ncol=ncol(pi[[j]]) )
    for(l in 1:ncol(pi[[j]])){
      pi_temp[,l]=spline(workGrid[[j]] ,pi[[j]][, l], xout=target_t)$y
    }
   
    multi_eigen=(svd(Z)$u)[index,index] %*% t(pi_temp)  ## time * n_M
    multi_score=t(xiEst_t2[[j]])%*%(svd(Z)$u)[index,index] ## 1 * n_M
   fittedX[[j]] = muxtemp + (t(multi_eigen) %*% t(multi_score))  ## time*1
  }
  return(fittedX)
}



gamma_cv<- function(yy,tt,   mux,  muy  , phi, eigvalc, M, dis,workGrid){
	
	
	pe<-vector()
	ga_index=seq(1, 10, length.out=5)
	
	eigvalc= eigvalc[1:M]
	phi= phi[,1:M]
	
	for(k in 1:length(ga_index)){
		
		yprederr<-matrix(nrow=length(yy), ncol=length(workGrid))
		gamma1=ga_index[k]
		
		for(i in 1:length(yy)){
			
			y=yy[[i]]
			t=tt[[i]]
		for(j in 1:length(y)){
		  
		 
			muxtemp=spline(workGrid, mux, xout=t)$y
			muytemp=spline(workGrid,muy, xout=t)$y
			
			phi_temp<-matrix(nrow=length(muytemp), ncol=M)
			for(l in 1:M){
			  phi_temp[,l]=spline(workGrid ,phi[, l], xout=t)$y
			}
			error0 = gamma1*diag(muytemp[-j], nrow=length(muytemp[-j]))
			
			if(length(muxtemp)==2){
			xi_predt =  diag(eigvalc)%*%(phi_temp[-j,])%*%diag(glink_der(muxtemp[-j],dis) , nrow=length(muxtemp[-j])) %*%ginv(diag(glink_der(muxtemp[-j],dis), nrow=length(muxtemp[-j]))%*%phi_temp[-j,]%*%diag(eigvalc)%*%(phi_temp[-j,])%*%diag(glink_der(muxtemp[-j],dis), nrow=length(muxtemp[-j]))+error0)  %*%(y[-j] -muytemp[-j])
			}
			if(length(muxtemp)>2){
			xi_predt =  diag(eigvalc)%*%t(phi_temp[-j,])%*%diag(glink_der(muxtemp[-j],dis) , nrow=length(muxtemp[-j])) %*%ginv(diag(glink_der(muxtemp[-j],dis), nrow=length(muxtemp[-j]))%*%phi_temp[-j,]%*%diag(eigvalc)%*%t(phi_temp[-j,])%*%diag(glink_der(muxtemp[-j],dis))+error0)  %*%(y[-j] -muytemp[-j])
			}
	        temp2 = muxtemp[j]+t(xi_predt)%*%phi_temp[j,]
	         ypred = glink(temp2, family=dis)
	        yprederr[i,j]= ( ypred- y[j])^2/ ypred
	     }    
		}
		
		
	pe[k]=mean(yprederr, na.rm=TRUE)

		}
		
		opt_gamma= ga_index[which.min(pe)]
		return(opt_gamma)
		
}
	
	

no_eigen<-function(yy, tt,mux,  muy  , phi, eigvalc , gamma, maxk=5, dis,workGrid){
	
	ql<-vector()
	 for( k  in  2:maxk){
	 	
	 		eigvalca= eigvalc[1:k]
			phia= phi[,1:k]
			
   temp<-vector()	
      for(i in 1:length(yy)){
			
			y=yy[[i]]
			t=tt[[i]]
			
		muxtemp=spline(workGrid, mux, xout=t)$y
			muytemp=spline(workGrid,muy, xout=t)$y
			
			phi_temp<-matrix(nrow=length(muytemp), ncol=k)
			for(l in 1:k){
			  phi_temp[,l]=spline(workGrid ,phia[, l], xout=t)$y
			}
			
			 error0 = gamma*diag(muytemp, nrow=length(muytemp))
			xi_predt =  diag(eigvalca)%*%t(phi_temp)%*%diag(glink_der(muxtemp,dis), nrow=length(muxtemp)  ) %*% ginv(diag(glink_der(muxtemp,dis), nrow=length(muxtemp)  )%*%phi_temp%*%diag(eigvalca)%*%t(phi_temp)%*%diag(glink_der(muxtemp,dis), nrow=length(muxtemp)  )+error0)  %*%(y -muytemp)
			
			
	        temp2 = muxtemp+phi_temp%*% xi_predt
	         ypred = glink(temp2, family=dis)
	         
	         if(is.na(ypred[1])==TRUE) {temp[i]= NA}
	         if(is.na(ypred[1])==FALSE){ temp[i]= QL_param(ypred, y, family=dis, sigma=gamma)}
	     }    
	
       		 
            ql[k]= sum(temp, na.rm=TRUE)
       }
       
	cret_val = -2*ql+2*c(1:maxk)
	
	M=which.min(cret_val)
	return(M)
}



QL_param<-function(x,y,family,sigma){

    
    if(family=='binomial(logit)'){
        gder =sum(y*log(x)+log((1-x)^(1-y)))/sigma
        }
    if(family=='poisson(log)'){
    	
    	myfac <-function(y){
    			output<-vector()
    			for(k in 1:length(y)){
    				x=y[k]
					 fact<-0
					 i<-x
					while(i>1){
						fact<-fact+log(i)
						 i<-i-1
					 }
					 output[k]=fact
					 }
		 return (output)
		 }

        gder = sum(y*log(x)-x-myfac(y)   ) /sigma
    }
  
  if(family=='gaussian'){
    gder = sum(  -(y-x)^2/2 )/sigma
  }
  
	return(gder)

}



















