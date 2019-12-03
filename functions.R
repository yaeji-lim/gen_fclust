### functions ###
library(permute)


bin_borrowing_next_day <- function(data, binSize, startPoint = 1) {
  n <- nrow(data)
  p <- ncol(data)
  if(binSize > n) {warning("binSize cannot be greater than size of data vector")}
  if(startPoint > binSize | startPoint < 1) {warning("startPoint not between 1 and B")}
  L <- ceiling(n/binSize)
  if (L == 1) return(colSums(data))
  binData <- matrix(nrow = L, ncol = p) 
  if (startPoint == 1) {
    for(j in 1:p) {
      for(i in 1:(L-1)) {
          binData[i,j] <- sum(data[((binSize*(i-1)+1):(binSize*i)), j])
      }
      binData[L, j] <- sum(data[(binSize*(L-1)+1):n, j])
    }
  } else { #starting point between 2 ~ B
    for(j in 1:(p-1)) {
      if(L > 2) {
        for(i in 1:(L-2)) {
          binData[i,j] <- sum(data[((binSize*(i-1)+ startPoint):((binSize*i) + startPoint - 1)), j])
        }
      }
      if((startPoint-1 > (n%%binSize)) & (n%%binSize != 0)) {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):n  , j]) +      sum(data[(1:(startPoint-(n%%B)-1) ) ,j+1])  
        binData[L, j] <- sum(data[(startPoint-(n%%B)):(startPoint-1), j+1])
      } else {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):(binSize*(L-1)+startPoint-1), j])
        binData[L, j] <- sum(data[((binSize*(L-1)+startPoint):n),j]) + sum(data[(1:(startPoint-1)),j+1])   # + sum(data[(1:(startPoint-1)),j])
      }
    }
        for(j in p) {
      if(L > 2) {
        for(i in 1:(L-2)) {
          binData[i,j] <- sum(data[((binSize*(i-1)+ startPoint):((binSize*i) + startPoint - 1)), j])
        }
      }
      if((startPoint-1 > (n%%binSize)) & (n%%binSize != 0)) {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):n  , j]) +      sum(data[(1:(startPoint-(n%%B)-1) ) ,1])  
        binData[L, j] <- sum(data[(startPoint-(n%%B)):(startPoint-1), 1])
      } else {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):(binSize*(L-1)+startPoint-1), j])
        binData[L, j] <- sum(data[((binSize*(L-1)+startPoint):n),j]) + sum(data[(1:(startPoint-1)),1])   # + sum(data[(1:(startPoint-1)),j])
      }
    }
  }
  return(binData)
}

kmeansBin_next_day <- function(data, binSize, startPoint, K) {
  binnedData <- t(bin_borrowing_next_day(data, binSize, startPoint))
  if(nrow(binnedData) == 1) {
    return(kmeans(t(binnedData), centers = K, nstart = 10)$cluster)
  } else {
    return(kmeans(t(bin_borrowing_next_day(data, binSize, startPoint)), centers = K, nstart = 10)$cluster)
  }
} 


best_match_factor=function(x,y){
  library(plyr)
  x=factor(x)
  A=NA
  for(k in 1:length(permn(length(table(x)))) ){
    new_b= as.factor(y)
    newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[k][[1]] ) )
    A=c(A,  length(which(newy==x) ) )
  }
  A=A[-1]
  newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[which.max(A)][[1]] ) )
  return(as.numeric(as.character(newy)) )
}


fitZIP <- function(vec) {
  res <- zeroinfl(vec~1, dist = "poisson") #prob of zero-inflation
  mu <- exp(res$coefficients$count) #mean of count distn
  zero <- inv.logit(res$coefficients$zero)
  return(c(mu, zero))
}

###     Permute clustering result by a given rule
###			Input : clst - length n vector : result of classification
###						  rule - length K vector : class i -> class rule[i] 
###     Ex: permute(c(1,3,4,2,2), c(4,3,2,1)) : 4 2 1 3 3
permute <- function(clst, rule) {
  result = vector(length = length(clst))
  for(i in 1:length(rule)) {
    result[clst == rule[i]] <- i
  }
  return(result)
}
#permute(c(1,3,4,2,2), c(4,3,2,1))

###     Clustering Result to Membership Matrix
###			Input : clst - length n vector : result of classification
###             K    - number of clusters
###     Output: n*K matrix(n: number of data) a_ij = 1 
###             if ith object’s cluster is j
###     
clstToMembership <- function(clst, K) {
  member <- matrix(0, nrow = length(clst), ncol = K)
  for(i in 1:length(clst)) {
    member[i, clst[i]] <- 1
  }
  return(member)
}

#tr(A)
tr <- function(A) {
  return(ifelse(nrow(A) == ncol(A), sum(diag(A)), NA))
}

###     Merging M clustering Results by the algorithm
###			Input : clstRes - M by n matrix of clusting results
###             K       - number of clusters
###     Output: n*K matrix(n: number of data) 
###             a_ij:likelihood that ith object’s cluster is j
###     

voting <- function(clstRes, K) {
  M <- nrow(clstRes)
  N <- ncol(clstRes)
  perm <- rbind(1:K, allPerms(K))
  P <- clstToMembership(clstRes[1,], K)
  for(m in 2:M) {
    Um <- clstToMembership(clstRes[m,], K)
    trcMax <- tr(t(P)%*%Um)
    indMax <- 1
    for(i in 1:nrow(perm)) {
      trc <- tr(t(P)%*%Um[,perm[i,]])
      if(trc > trcMax) {
        trcMax <- trc
        indMax <- i
      }
    }
    P <- ((m-1)/m)*P + Um[,perm[indMax,]]/m
  }
  return(P)
}


###   Membership matrix to clustering result
###   Input:  votemat - n*K membership matrix returned from voting
###   Output: length n vector, return the cluster with maxium likelihood
voteToCluster <- function(voteMat) {
  return(apply(voteMat, 1, which.max))
}

###   Membership matrix to clustering result
###   Input:  votemat - n*K membership matrix returned from voting
###   Output: length n vector, return the maxium likelihood
voteToSureness <- function(voteMat) {
  return(apply(voteMat, 1, max))
}

###   Membership matrix to clustering result
###   Input:  clst1 - clustering result (기준)
###           clst2 - clustering result 
###           K     - number of clusters  
###   Output: Best permuation for cluster2
###           (i.e. permute(c2, findBestRule(c1, c2, 4))하면 c1과 일치)

# findBestRule(c(1,3,4,2,1,3,4,2), c(2,4,1,3,2,4,1,3), 4)
# permute(c(2,4,1,3,2,4,1,3), findBestRule(c(1,3,4,2,1,3,4,2), c(2,4,1,3,2,4,1,3), 4))
# permute(c(1,2,3,4), c(4,3,2,1))
# permute(c(1,3,2,4), c(4,3,2,1))
findBestRule <-function(clst1, clst2, K) {
  perm <- rbind(1:K, allPerms(K))
  P <- clstToMembership(clst1, K)
  Um <- clstToMembership(clst2, K)
  trcMax <- tr(t(P)%*%Um)
  indMax <- 0
  for(i in 1:nrow(perm)) {
    trc <- tr(t(P)%*%Um[,perm[i,]])
    if(trc > trcMax) {
      trcMax <- trc
      indMax <- i
    }
  }
  if(indMax > 0) {
    return(perm[indMax, ])
  } else {
    return(c(1:K))
  }
}
# m <- matrix(c(0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1), nrow = 4, ncol = 4, byrow = T)
# m[,c(1,4,3,2)]
# m[,c(4,3,2,1)]


###------------------------------------------------
###  Binnig a given data given a binSize and startPoint
###  Input:   data  -  ncol = number of subjects(p) , nrow = length of data(512, 1440,... etc) (n)
###           binSize - B, size of the bin
###           startPoint - number between 1,..,B (default = 1)
###  Output:  binned data (ncol = p, nrow = ceiling(n/B))
bin <- function(data, binSize, startPoint = 1) {
  n <- nrow(data)
  p <- ncol(data)
  if(binSize > n) {warning("binSize cannot be greater than size of data vector")}
  if(startPoint > binSize | startPoint < 1) {warning("startPoint not between 1 and B")}
  L <- ceiling(n/binSize)
  if (L == 1) return(colSums(data))
  binData <- matrix(nrow = L, ncol = p) 
  if (startPoint == 1) {
    for(j in 1:p) {
      for(i in 1:(L-1)) {
          binData[i,j] <- sum(data[((binSize*(i-1)+1):(binSize*i)), j])
      }
      binData[L, j] <- sum(data[(binSize*(L-1)+1):n, j])
    }
  } else { #starting point between 2 ~ B
    for(j in 1:p) {
      if(L > 2) {
        for(i in 1:(L-2)) {
          binData[i,j] <- sum(data[((binSize*(i-1)+ startPoint):((binSize*i) + startPoint - 1)), j])
        }
      }
      if((startPoint-1 > (n%%binSize)) & (n%%binSize != 0)) {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):n, j]) +
                             sum(data[(1:(startPoint-(n%%B)-1)),j])
        binData[L, j] <- sum(data[(startPoint-(n%%B)):(startPoint-1), j])
      } else {
        binData[(L-1), j] <- sum(data[(binSize*(L-2)+startPoint):(binSize*(L-1)+startPoint-1), j])
        binData[L, j] <- sum(data[((binSize*(L-1)+startPoint):n),j]) + sum(data[(1:(startPoint-1)),j])
      }
    }
  }
  return(binData)
}

#bin(t(data), binSize = 256, startPoint = 4)

# fitZIP <- function(vec) {
#   res <- zeroinfl(vec~1, dist = "poisson")
#   mu <- exp(res$coefficients$count)
#   zero <- inv.logit(res$coefficients$zero)
#   return(c(mu, zero))
# }
# 
# fitPois <- function(vec) {
#   res <- glm(vec~1, family = poisson(link = log))
#   mu <- exp(res$coefficients)
#   return(mu)
# }
# 
# fitPoisson <- vector(length = 1440)
# for(i in 1:1440) {
#   fitPoisson[i] <- fitPois(newStep[i,])
# }


### Doppler signal 
### t: time, t0 : phase
doppler<-function (t,t0) {
  -0.025+0.6*sqrt(t * (1 - t)) * sin((2 * pi * 1.05)/(t -t0))
}

### K-means binned data given bin size and starting point
### Input:  data  -  ncol = number of subjects(p) , nrow = length of data(512, 1440,... etc) (n)
###        

kmeansBin <- function(data, binSize, startPoint, K) {
  binnedData <- t(bin(data, binSize, startPoint))
  if(nrow(binnedData) == 1) {
    return(kmeans(t(binnedData), centers = K, nstart = 10)$cluster)
  } else {
    return(kmeans(t(bin(data, binSize, startPoint)), centers = K, nstart = 10)$cluster)
  }
} 

#kmeans(t(bin(t(data), binSize = 512, startPoint = 1)), centers = 4, nstart = 10)




### CCR after matching

ccr <- function(result, true, K) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == permute(result,findBestRule(true, result, K)) )/length(result))
}




###  Function for FClust
"fitfclust" <-
  function(x=NULL,curve=NULL,timeindex=NULL,data=NULL, q = 5, h = 1, p = 5,
           K = 2, tol = 0.001, maxit = 20, pert =  
             0.01, grid = seq(0, 1, length = 100), hard = F, plot= F,trace=F)
  {
    # This is the main function to implement the FClust procedure.
    library(splines)
    if (is.null(data))
      data <- list(x=x,curve=curve,timeindex=timeindex)
    initfit <- fclustinit(data = data, pert = pert, grid = grid, h = h,
                          p = p, q = q, K = K)
    # Initialize the parameters
    parameters <- initfit$parameters
    vars <- initfit$vars
    S <- initfit$S
    FullS <- initfit$FullS
    sigma.old <- 0
    sigma.new <- 1
    ind <- 1
    # Main loop. Iterates between M and E steps and stops when
    # sigma  has converged.
    while(abs(sigma.old - sigma.new)/sigma.new > tol & (ind <= maxit)) {
      parameters <- fclustMstep(parameters, data, vars, S, tol, p, hard)
      vars <- fclustEstep(parameters, data, vars, S, hard)
      sigma.old <- sigma.new
      sigma.new <- parameters$sigma[1]
      if (trace)
        print(paste("Iteration", ind,": Sigma = ",sigma.new))
      # Plot cluster mean curves.  
      if(plot){
        cluster.mean <- matrix(0,K,dim(FullS)[1])
        for(k in 1:K)
          cluster.mean[k,] <- FullS %*% (parameters$lambda.zero + parameters$Lambda %*% parameters$alpha[k,  ])
        plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
        for(k in 1:K)
          lines(grid,cluster.mean[k,], col = 4, lwd = 2)
      }
      ind <- ind + 1
    }
    # Enforce parameter constraint (7)
    cfit <- fclustconst(data, parameters, vars, FullS)
    list(data=data,parameters = cfit$parameters, vars = cfit$vars, FullS = FullS,grid=grid)
  }

"fclustinit" <-
  function(data, pert = 0, grid = seq(0.01, 1, length = 100), h = 1, p = 2, q = 5,K = K){
    S <- FullS <- NULL
    # This function initializes all the parameters.
    # Produce spline basis matrix
    FullS <- cbind(1, ns(grid, df = (q - 1)))
    FullS <- svd(FullS)$u
    S <- FullS[data$timeindex,  ]
    N <- length(table(data$curve))
    # Compute initial estimate of basis coefficients.
    points <- matrix(0,N,sum(q))
    for (i in 1:N){
      Si <- S[data$curve==i,]
      xi <- data$x[data$curve==i]
      points[i,] <- solve(t(Si) %*% Si + pert * diag(q)) %*% t(Si) %*%xi}
    # Use k-means to get initial cluster memberships from points.
    if(K > 1)
      class <- kmeans(points, K, 10)$cluster
    else class <- rep(1, N)
    # Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, N, K)
    piigivej[col(piigivej) == class] <- 1
    # Calculate coefficeints for cluster means.
    classmean <- matrix(0,K,q)
    for (k in 1:K)
      classmean[k,] <- apply(points[class==k,],2,mean)
    # Initialize lambda.zero, Lambda and alpha as defined in paper.
    lambda.zero <- apply(classmean, 2, mean)
    Lambda <- as.matrix(svd(scale(classmean, scale = F))$v[, 1:h])
    alpha <- scale(classmean, scale = F)%*%Lambda
    # Calculate estimates for gamma and gprod.
    gamma <- t(t(points) - lambda.zero - (Lambda %*% t(alpha[class,  ])))
    gprod <- NULL
    for(i in 1:N)
      gprod <- cbind(gprod, (gamma[i,  ]) %*% t(gamma[i,  ]))
    gamma <- array(gamma[, rep(1:sum(q), rep(K, sum(q)))], c(N, K, sum(q)))
    gcov <- matrix(0, sum(q), N * sum(q))
    list(S = S, FullS = FullS, parameters = list(lambda.zero = lambda.zero,
                                                 Lambda = Lambda, alpha = alpha), vars = list(gamma = gamma,
                                                                                              gprod = gprod, gcov = gcov, piigivej = piigivej))
  }

"fclustMstep" <-
  function(parameters, data, vars, S, tol, p, hard)
  {
    # This function implements the M step of the EM algorithm.
    K <- dim(parameters$alpha)[1]
    alpha <- parameters$alpha
    Lambda <- parameters$Lambda
    gamma <- vars$gamma
    gcov <- vars$gcov
    curve <- data$curve
    piigivej <- vars$piigivej
    N <- dim(gamma)[1]
    K <- dim(alpha)[1]
    h <- dim(alpha)[2]
    n <- length(curve)
    q <- dim(S)[2]
    sigma.old <- 2
    sigma <- 1
    # Compute pi.
    if(hard)
      parameters$pi <- rep(1/K, K)
    else parameters$pi <- apply(vars$piigivej, 2, mean)
    # Compute rank p estimate of Gamma
    ind <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    gsvd <- svd(vars$gprod %*% ind/N)
    gsvd$d[ - (1:p)] <- 0
    parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
    # This loop iteratively calculates alpha and then Lambda and stops
    # when they have converged.
    loopnumber <- 1
    while((abs(sigma.old[1] - sigma[1])/sigma[1] > tol) & (loopnumber <10)){
      sigma.old <- sigma
      # Calculate lambda.zero.
      gamma.pi <- diag(S %*% t(apply(gamma * as.vector(piigivej),
                                     c(1, 3), sum)[curve,  ]))
      alpha.pi <- t(matrix(apply(alpha, 2, function(x, piigivej,K)
      {rep(1, K) %*% (piigivej * x)}
      , t(piigivej), K), N, h)[curve,  ])
      lambda.zero <- solve(t(S) %*% S) %*% t(S) %*% (data$x - diag(
        S %*% Lambda %*% alpha.pi) - gamma.pi)
      x <- data$x - S %*% lambda.zero
      # Calculate alpha.
      for(i in 1.:K) {
        S.Lam <- S %*% Lambda
        S.Lam.pi <- S.Lam * piigivej[curve, i]
        if(sum(piigivej[, i]) > 10^(-4))
          alpha[i,  ] <- solve(t(S.Lam.pi) %*% S.Lam) %*%
          t(S.Lam.pi) %*% (x - diag(S %*% t(gamma[curve, i,  ])))
        else print("Warning: empty cluster")
      }
      # Calculate Lambda given alpha. This is done by iterating
      # through each column of Lambda holding the other columns fixed.
      for(m in 1:h) {
        pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)[curve]
        pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)[curve]
        S.Lambda <- t(S %*% Lambda)
        if(h != 1) {
          temp <- NULL
          for(i in 1:K) {
            temp <- cbind(temp, as.vector(rep(1, h - 1) %*% matrix((S.Lambda *
                                                                      alpha[i,  ])[ - m,  ], h - 1,dim(S)[1])) * alpha[i, m])
          }
          otherLambda <- (temp * piigivej[curve,  ])%*%rep(1, K)
        }
        else otherLambda <- 0
        gamma.pi.alpha <- apply(gamma * as.vector(piigivej) *
                                  rep(alpha[, m], rep(N, K)), c(1, 3), sum)[curve,  ]
        Lambda[, m] <- solve(t(S * pi.alphasq) %*% S) %*% t(S) %*%
          (x * pi.alpha - otherLambda - (S *gamma.pi.alpha) %*% rep(1, sum(q)))
      }
      # Calculate sigma 
      ind <- matrix(rep(c(rep(F, q), rep(T, N * q)),N)
                    [1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, table(curve)),]
      mat1 <- matrix(rep(S, N), n, N * q)
      mat3 <- t(mat1)
      mat3[t(ind)] <- 0
      ind2 <- matrix(rep(c(rep(F, q), rep(T, N * q)),
                         N)[1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, rep(q, N)),]
      mat2 <- matrix(rep(t(gcov), N), N * q, N * q,byrow = T)
      mat2[ind2] <- 0
      econtrib <- 0
      for(i2 in 1:K) {
        vect1 <- x - S %*% Lambda %*% alpha[
          i2,  ] - (S * gamma[curve, i2,  ]) %*% rep(1, q)
        econtrib <- econtrib + t(piigivej[curve,i2] * vect1) %*% vect1
      }
      sigma <- as.vector((econtrib + sum(diag(mat1 %*% mat2 %*% mat3)))/n)
      loopnumber <- loopnumber + 1
    }
    parameters$lambda.zero <- as.vector(lambda.zero)
    parameters$alpha <- alpha
    parameters$Lambda <- Lambda
    parameters$sigma <- sigma
    parameters
  }

"fclustEstep" <-
  function(parameters, data, vars, S, hard)
  {
    # This function performs the E step of the EM algorithm by
    # calculating the expected values of gamma and gamma %*% t(gamma)
    # given the current parameter estimates.
    par <- parameters
    N <- dim(vars$gamma)[1]
    K <- dim(vars$gamma)[2]
    q <- dim(vars$gamma)[3]
    Gamma <- par$Gamma
    Lambda.alpha <- par$lambda.zero + par$Lambda %*% t(par$alpha)
    vars$gprod <- vars$gcov <- NULL
    for(j in 1:N) {
      # Calculate expected value of gamma.
      Sj <- S[data$curve == j,  ]
      nj <- sum(data$curve == j)
      invvar <- diag(1/rep(par$sigma, nj))
      Cgamma <- Gamma - Gamma %*% t(Sj) %*% solve(diag(nj) + Sj %*%
                                                    Gamma %*% t(Sj) %*% invvar) %*% invvar %*% Sj %*% Gamma
      centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
      vars$gamma[j,  ,  ] <- t(Cgamma %*% t(Sj) %*% invvar %*% centx)
      # Calculate pi i given j.
      covx <- Sj %*% par$Gamma %*% t(Sj) + solve(invvar)
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * par$pi
      vars$piigivej[j,  ] <- d/sum(d)
      if(hard) {
        m <- order( - d)[1]
        vars$piigivej[j,  ] <- 0
        vars$piigivej[j, m] <- 1
      }
      # Calculate expected value of gamma %*% t(gamma).
      vars$gprod <- cbind(vars$gprod, t(matrix(vars$gamma[j,  ,  ],
                                               K, q)) %*% (matrix(vars$gamma[j,  ,  ], K, q) * vars$
                                                             piigivej[j,  ]) + Cgamma)
      vars$gcov <- cbind(vars$gcov, Cgamma)
    }
    vars
  }


"fclustconst" <-
  function(data, parameters, vars, S)
  {
    # This function enforces the constraint (7) from the paper on the
    # parameters. This means that the alphas can be interpreted as the
    # number of standard deviations that the groups are apart etc.
    par <- parameters
    A <- t(S) %*% solve(par$sigma * diag(dim(S)[1]) + S %*% par$Gamma %*%
                          t(S)) %*% S
    svdA <- svd(A)
    sqrtA <- diag(sqrt(svdA$d)) %*% t(svdA$u)
    negsqrtA <- svdA$u %*% diag(1/sqrt(svdA$d))
    finalsvd <- svd(sqrtA %*% par$Lambda)
    par$Lambda <- negsqrtA %*% finalsvd$u
    if(dim(par$Lambda)[2] > 1)
      par$alpha <- t(diag(finalsvd$d) %*% t(finalsvd$v) %*% t(par$alpha))
    else par$alpha <- t(finalsvd$d * t(finalsvd$v) %*% t(par$alpha))
    meanalpha <- apply(par$alpha, 2, mean)
    par$alpha <- t(t(par$alpha) - meanalpha)
    par$lambda.zero <- par$lambda.zero + par$Lambda %*% meanalpha
    list(parameters = par, vars = vars)
  }

"nummax" <-
  function(X)
  {
    ind <- rep(1, dim(X)[1])
    m <- X[, 1]
    if(dim(X)[2] > 1)
      for(i in 2:dim(X)[2]) {
        test <- X[, i] > m
        ind[test] <- i
        m[test] <- X[test, i]
      }
    list(ind = ind, max = m)
  }

"fclust.pred" <-
  function(fit,data=NULL,reweight=F)
  {
    # This function produces the alpha hats used to provide a low
    # dimensional pictorial respresentation of each curve. It also
    # produces a class prediction for each curve. It takes as
    # input the fit from fldafit (for predictions on the original data)
    # or the fit and a new data set (for predictions on new data).
    if (is.null(data))
      data <- fit$data
    FullS <- fit$FullS
    par <- fit$parameters
    curve <- data$curve
    time <- data$time
    N <- length(table(curve))
    h <- dim(par$alpha)[2]
    alpha.hat <- matrix(0, N, h)
    K <- dim(fit$par$alpha)[1]
    distance <- matrix(0, N, K)
    Calpha <- array(0, c(N, h, h))
    for(i in 1:N) {
      Sij <- FullS[time[curve == i],  ]
      xij <- data$x[curve == i]
      n <- length(xij)
      Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
      # Calculate covariance for each alpha hat.
      InvCalpha <- t(par$Lambda) %*% t(Sij) %*% solve(Sigma) %*% Sij %*%
        par$Lambda
      Calpha[i,  ,  ] <- solve(InvCalpha)
      # Calculate each alpha hat.
      alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
        Sij) %*% solve(Sigma) %*% (xij - Sij %*% par$lambda.zero)
      # Calculate the matrix of distances, relative to the
      # appropriate metric of each curve from each class centroid. 
      for (k in 1:K){
        y <- as.vector(alpha.hat[i,])-fit$par$alpha[k,]
        distance[i,k] <- t(y)%*%InvCalpha %*%y}}
    # Calculate final class predictions for each curve.
    class.pred <- rep(1, N)
    log.pi <- log(fit$par$pi)
    if (!reweight)
      log.pi <- rep(0,K)
    probs <- t(exp(log.pi-t(distance)/2))
    probs <- probs/apply(probs,1,sum)
    m <- probs[,1]
    if(K != 1)
      for(k in 2:K) {
        test <- (probs[, k] > m)
        class.pred[test] <- k
        m[test] <- probs[test, k]
      }
    list(Calpha = Calpha, alpha.hat = alpha.hat, class.pred = class.pred,
         distance = distance, m = m,probs=probs)
  }

"fclust.curvepred" <-
  function(fit, data=NULL, index=NULL, tau = 0.95, tau1 = 0.975)
  {
    if (is.null(data))
      data <- fit$data
    if (is.null(index))
      index <- 1:length(table(data$curve))
    tau2 <- tau/tau1
    sigma <- fit$par$sigma
    Gamma <- fit$par$Gamma
    Lambda <- fit$par$Lambda
    alpha <- fit$par$alpha
    lambda.zero <- as.vector(fit$par$lambda.zero)
    S <- fit$FullS
    N <- length(index)
    upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,N,nrow(S))
    etapred <- matrix(0,N,ncol(S))
    ind <- 1
    Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
    for (i in index){
      y <- data$x[data$curve == i]
      Si <- S[data$time[data$curve == i],  ]
      ni <- dim(Si)[1]
      invvar <- diag(1/rep(sigma, ni))
      covx <- Si %*% Gamma %*% t(Si) + solve(invvar)
      centx <- data$x[data$curve == i] - Si %*% Lambda.alpha
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * fit$par$pi
      pi <- d/sum(d)
      K <- length(pi)
      mu <- lambda.zero + Lambda %*% t(alpha * pi) %*% rep(1, K)
      cov <- (Gamma - Gamma %*% t(Si) %*% solve(diag(ni) + Si %*% Gamma %*%
                                                  t(Si)/sigma) %*% Si %*% Gamma/sigma)/sigma
      etapred[ind,] <- mu + cov %*% t(Si) %*% (y - Si %*% mu)
      ord <- order( - pi)
      numb <- sum(cumsum(pi[ord]) <= tau1) + 1
      v <- diag(S %*% (cov * sigma) %*% t(S))
      pse <- sqrt(v + sigma)
      se <- sqrt(v)
      lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(S), numb)
      for(j in 1:numb) {
        mu <- lambda.zero + Lambda %*% alpha[ord[j],  ]
        mean <- S %*% (mu + cov %*% t(Si) %*% (y - Si %*% mu))
        upper.p[, j] <- mean + qnorm(tau2) * pse
        lower.p[, j] <- mean - qnorm(tau2) * pse
        upper.c[, j] <- mean + qnorm(tau2) * se
        lower.c[, j] <- mean - qnorm(tau2) * se
      }
      upci[ind,] <- nummax(upper.c)$max
      lowci[ind,] <-  - nummax( - lower.c)$max
      uppi[ind,] <- nummax(upper.p)$max
      lowpi[ind,] <-  - nummax( - lower.p)$max
      gpred[ind,] <- as.vector(S %*%etapred[ind,])
      ind <- ind+1
    }
    meancurves <- S%*%Lambda.alpha
    list(etapred = etapred, gpred = gpred,  upci = upci,lowci = lowci,  uppi = uppi, lowpi = lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
  }

"fclust.discrim" <-
  function(fit,absvalue=F){
    S <- fit$FullS
    sigma <- fit$par$sigma
    n <- nrow(S)
    Gamma <- fit$par$Gamma
    Sigma <- S%*%Gamma%*%t(S)+sigma*diag(n)
    Lambda <- fit$par$Lambda
    discrim <- solve(Sigma)%*%S%*%Lambda
    if (absvalue)
      discrim <- abs(discrim)
    n <- ncol(discrim)
    nrows <- ceiling(sqrt(n))
    par(mfrow=c(nrows,nrows))
    for (i in 1:n){
      plot(fit$grid,discrim[,i],ylab=paste("Discriminant Function ",i),xlab="Time",type='n')
      lines(fit$grid,discrim[,i],lwd=3)
      abline(0,0)}}

"fclust.plotcurves" <-
  function(object=NULL,fit=NULL,index=NULL,ci=T,pi=T,clustermean=F){
    if (is.null(object))
      object <- fclust.curvepred(fit)
    if (is.null(index))
      index <- 1:length(table(object$data$curve))
    r <- ceiling(sqrt(length(index)))
    par(mfrow=c(r,r))
    for (i in index){
      grid <- object$grid
      upci <- object$upci[i,]
      uppi <- object$uppi[i,]
      lowci <- object$lowci[i,]
      lowpi <- object$lowpi[i,]
      gpred <- object$gpred[i,]
      meancurves <- (object$mean)
      if (clustermean)
        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
      else
        yrange <- c(min(lowpi),max(uppi))
      plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
           main=paste("Curve ",i))
      if (clustermean)
        for (k  in 1:ncol(meancurves))
          lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
      if (ci){
        lines(grid,upci,col=3)
        lines(grid,lowci,col=3)}
      if (pi){
        lines(grid,uppi,col=4)
        lines(grid,lowpi,col=4)}
      lines(grid,gpred,col=2,lwd=2)
      lines(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],lwd=2)
      points(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],pch=19,cex=1.5)
    }
  }

####plot curve
plotCur <- function(dat, yl = 190) {
  plot(dat, type = "l", ylim = c(0, yl))
}

binIter <- function(voteRes, result, data, numRes) { #numRes = nrow(res)
  #data: data의 개수 * data의 dimension (N * d)
  cent <- list()
  for(i in unique(result)) { 
    if(sum(result == i) > 1) {
      cent[[i]] <- apply(data[result == i,], 2, mean)
    } else { #cluster에 원소가 한개 있을 때
      cent[[i]] <- data[result == i, ]
    }
  }
  dist_mat <- matrix(nrow = length(result), ncol = 4) #N * K matrix
  for(j in 1:length(result)) {
    for(i in unique(result)) {
      dist_mat[j, i] <- dist(rbind(cent[[i]], data[j, ]))
    }
  }
  min_clst <- apply(dist_mat, 1, which.min)
  
  for(i in 1:length(result)) {
    if(min_clst[i] != result[i]) {
      if(abs(voteRes[i, min_clst[i]] - voteRes[i,result[i]]) <= 1/numRes) { 
        result[i] = min_clst[i]
      }
    }
  }
  return(result)
}
