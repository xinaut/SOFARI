#return the normal vector of u1 u2 u3 and d1 d2 d3
normal_u <- function(n, p, dU0, D0, u1, u2, u3, d1, d2, d3, varu1, varu2, varu3, vard1, vard2, vard3){
  u01 <- dU0[,1];  u02 <- dU0[,2];  u03 <- dU0[,3]
  d01 <- D0[1,1]; d02 <- D0[2,2]; d03 <- D0[3,3]
  nu1 = nu2 = nu3 = rep(0,p)
  nd1 = nd2 = nd3 = 0
  u1 = as.vector(u1); u2 = as.vector(u2); u3 = as.vector(u3)
  for(i in 1:p){
    nu1[i] <- sqrt(n)*(u1[i] - u01[i])/sqrt(varu1[i])
    nu2[i] <- sqrt(n)*(u2[i] - u02[i])/sqrt(varu2[i])
    nu3[i] <- sqrt(n)*(u3[i] - u03[i])/sqrt(varu3[i])
    nd1 <- sqrt(n)*(d1 - d01^2)/sqrt(vard1)
    nd2 <- sqrt(n)*(d2 - d02^2)/sqrt(vard2)
    nd3 <- sqrt(n)*(d3 - d03^2)/sqrt(vard3)
  }
  return(list(nu1 = nu1, nu2 = nu2, nu3 = nu3, nd1 = nd1, nd2 = nd2, nd3 = nd3))
}




#generate t(U)%*%cov(X)%*%U data X Y
porth.sim <- function(U,D,V,n,snr,Xsigma,rho=0){
  
  ## finding basis along more number of columns of data vector 
  basis.vec =function(x){
    # require(Matrix)
    if(diff(dim(x))<0) x <- t(x)
    qd <- qr(x)
    k <- qr.Q(qd) %*% qr.R(qd)[,1:qd$rank]
    k[abs(k)<1e-6] <- 0
    b.ind <- vector()
    for(i in 1:qd$rank)
      b.ind[i] <- which(apply(x,2,function(x,y)sum(abs(x-y)),k[,i])<1e-6)[1]
    return(list(ind=b.ind,vec = x[,b.ind]))
  }
  
  p <- nrow(U);q <- nrow(V);nrank <- ncol(U)
  
  U.t <- diag(max(dim(U)))
  U.t <- U.t[,-basis.vec(U)$ind]
  P <- cbind(U,U.t)
  UtXsUt <- t(U.t)%*%Xsigma%*%U.t
  UtXsU <- t(U.t)%*%Xsigma%*%U
  UXsU <- t(U)%*%Xsigma%*%U
  UXsUinv <- solve(UXsU)
  ##sigma.X2 <- t(U.t)%*%Xsigma%*%U.t - t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%t(U)%*%Xsigma%*%U.t
  sigma.X2 <- UtXsUt - UtXsU%*%UXsUinv%*%t(UtXsU)
  sigma.X2 <- (sigma.X2+t(sigma.X2))/2
  
  
  X1 <- matrix(nrow=nrank,ncol=n,rnorm(n*nrank))
  ##X1 <- t(mvrnorm(n,rep(0,ncol(U)),diag(ncol(U)) ))
  mean.X2 <- UtXsU%*%UXsUinv%*%X1
  ##mean.X2 <- t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%X1
  X2 <- mean.X2 + t(MASS::mvrnorm(ncol(mean.X2),rep(0,nrow(mean.X2)),sigma.X2))
  X <- t(solve(t(P))%*%rbind(X1,X2))#/sqrt(n)
  # crossprod(X%*%U)
  
  sige0 <- rho^abs(outer(1:q, 1:q,FUN="-"))
  svdrr <- eigen(sige0)
  svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
  
  UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
  
  # UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))
  C <- U%*%D%*%t(V)
  Y3 <- X%*%U[,nrank]%*%t(V[,nrank])*D[nrank,nrank]
  sigma <- sqrt(sum(Y3^2)/sum(UU^2))/snr    ## recheck 
  UU <- UU*sigma
  Y <- X%*%C + UU      ## data prepration ends
  sigmae <- sigma^2*sige0
  return(list(Y=Y,X=X, sigmae = sigmae))
}

##X is I
porth.simide <- function(U,D,V,n,snr,Xsigma,rho=0){
  
  ## finding basis along more number of columns of data vector 
  basis.vec =function(x){
    # require(Matrix)
    if(diff(dim(x))<0) x <- t(x)
    qd <- qr(x)
    k <- qr.Q(qd) %*% qr.R(qd)[,1:qd$rank]
    k[abs(k)<1e-6] <- 0
    b.ind <- vector()
    for(i in 1:qd$rank)
      b.ind[i] <- which(apply(x,2,function(x,y)sum(abs(x-y)),k[,i])<1e-6)[1]
    return(list(ind=b.ind,vec = x[,b.ind]))
  }
  
  p <- nrow(U);q <- nrow(V);nrank <- ncol(U)
  
  X <- diag(p)
  
  sige0 <- rho^abs(outer(1:q, 1:q,FUN="-"))
  svdrr <- eigen(sige0)
  svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
  
  UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
  
  # UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))
  C <- U%*%D%*%t(V)
  Y3 <- X%*%U[,nrank]%*%t(V[,nrank])*D[nrank,nrank]
  sigma <- sqrt(sum(Y3^2)/sum(UU^2))/snr    ## recheck 
  UU <- UU*sigma
  Y <- X%*%C + UU      ## data prepration ends
  sigmae <- sigma^2*sige0
  return(list(Y=Y,X=X, sigmae = sigmae))
}

#the covariance matrix of error E
sigmaee <- function(q){
  rho = 0.3
  siee = rho^abs(outer(1:q, 1:q, FUN = "-"))
  sigma = 1
  return(sigma*siee)
}


# make.data: function generating data matrix

make.data <- function(Sigma.half, n, p){
  ##  set.seed(seed)  
  
  X = matrix(rnorm(n*p),n,p)%*%Sigma.half
  return(X)
}
###########################################################################

# blk.mat: function generating block diagonal Omega

blk.mat <- function(a=0.5, p, permu){
  Omega.blk = matrix(c(1,a,a,1),2,2)
  Sigma.blk = Omega.blk
  Sigma.half.blk =  chol(Sigma.blk)
  Omega = Omega.blk
  Sigma = Sigma.blk
  Sigma.half = Sigma.half.blk
  
  for(j in 1:((p/2)-1)){
    Omega=bdiag(Omega, Omega.blk)
    Sigma = bdiag(Sigma, Sigma.blk)
    Sigma.half = bdiag(Sigma.half, Sigma.half.blk)
  }
  
  Sigma = Sigma[permu,permu]
  Omega = Omega[permu,permu]
  Sigma.half = Sigma.half[permu,permu]
  data = list(Sigma = Sigma, Omega=Omega, Sigma.half = Sigma.half)
  return(data)
}
###########################################################################

# ar1.mat: function generating AR(1) Omega

ar1.mat <- function(a, p, permu=c(1:p)){
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  Omega<- a^H
  Sigma = Omega 
  Sigma = Sigma*(abs(Sigma)>1e-4)
  Sigma.half = chol(Sigma)
  Sigma.half = Sigma.half*(abs(Sigma.half)>1e-4)
  Sigma = Sigma[permu,permu]
  Omega = Omega[permu,permu]
  Sigma.half = Sigma.half[permu,permu]
  obj = list(Sigma=Sigma, Omega = Omega, Sigma.half = Sigma.half)
}
###########################################################################

# band.mat: function generating band Omega

band.mat <- function(a, p, K=1, permu=c(1:p)){
  ones = rep(1,p)
  Omega0 = a*ones%*%t(ones)
  diag(Omega0) = rep(1,p)
  Omega = 1*band(Omega0,-K,K)
  Sigma = Omega
  Sigma = Sigma*(abs(Sigma)>1e-4)
  Sigma.half=chol(Sigma)
  Sigma.half = Sigma.half*(abs(Sigma.half)>1e-4)
  Sigma = Sigma[permu,permu]
  Omega = Omega[permu,permu]
  Sigma.half = Sigma.half[permu,permu]
  obj = list(Sigma=Sigma, Omega = Omega, Sigma.half = Sigma.half)
}


# generating function

rrr.simnew <-
  function(n = 50,
           p = 25,
           q = 25,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0,
           Xsigma,
           D = diag(c(100, 50, 20))){
    # Sigma <- function(p, rho)
    # {
    #   Sigma <- matrix(nrow = p, ncol = p, NA)
    #   for (i in seq_len(p)) {
    #     for (j in seq_len(p)) {
    #       Sigma[i, j] <- rho ^ (abs(i - j))
    #     }
    #   }
    #   Sigma
    # }
    
    U <- matrix(ncol = nrank, nrow = p)
    V <- matrix(ncol = nrank, nrow = q)
    
    U[, 1] <- c(sample(c(-1, 1), 5, replace = T), rep(0, p - 5))
    U[, 2] <-
      c(rep(0, 3), -U[4, 1], U[5, 1], sample(c(-1, 1), 3, replace = T), rep(0, p -
                                                                              8))
    U[, 3] <- c(rep(0, 8), sample(c(-1, 1), 2, replace = T), rep(0, p - 10))
    U[, 1] <- U[, 1] / lnorm(U[, 1], 2)
    U[, 2] <- U[, 2] / lnorm(U[, 2], 2)
    U[, 3] <- U[, 3] / lnorm(U[, 3], 2)
    
    V[, 1] <- c(sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1), rep(0, q -
                                                                           5))
    V[, 2] <-
      c(rep(0, 5),
        sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1),
        rep(0, q - 10))
    V[, 3] <-
      c(rep(0, 10),
        sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1),
        rep(0, q - 15))
    V[, 1] <- V[, 1] / lnorm(V[, 1], 2)
    V[, 2] <- V[, 2] / lnorm(V[, 2], 2)
    V[, 3] <- V[, 3] / lnorm(V[, 3], 2)
    
    #D <- diag(c(100, 50, 20))
    
    
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    ## X <- diag(p)
    
    # UU <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
    # if (rho_E != 0) {
    #   for (t in 1:n)
    #     UU[t, ] <- arima.sim(list(order = c(1, 0, 0), ar = rho_E), n = q)
    # }
    
    sige0 <- rho_E^abs(outer(1:q, 1:q,FUN="-"))
    svdrr <- eigen(sige0)
    svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
    
    UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
    
    C <- U %*% D %*% t(V)
    C3 <- U[, 3] %*% t(V[, 3]) * D[3, 3]
    
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    sigmae <- sigma^2*sige0
    
    list(
      Y = Y,
      X = X,
      C = C,
      U = U,
      V = V,
      D = D,
      Xsigma = Xsigma,
      sigmae = sigmae
    )
    
  }

rrr.simnew2 <-
  function(n = 50,
           p = 25,
           q = 25,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0,
           Xsigma,
           D = diag(c(100, 50, 20))){
    # Sigma <- function(p, rho)
    # {
    #   Sigma <- matrix(nrow = p, ncol = p, NA)
    #   for (i in seq_len(p)) {
    #     for (j in seq_len(p)) {
    #       Sigma[i, j] <- rho ^ (abs(i - j))
    #     }
    #   }
    #   Sigma
    # }
    
    U <- matrix(ncol = nrank, nrow = p)
    V <- matrix(ncol = nrank, nrow = q)
    
    ld1 = 5
    ld2 = 5
    ld3 = 5
    
    U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
    # U[,1]<-c(sample(c(1,-1),ld1,replace=TRUE),rep(0,p-ld1))
    # #U[,2]<-c(rep(0,ld1-2), -U[4, 1], U[5, 1], sample(c(1,-1),ld2-2,replace=TRUE),rep(0,p-ld1-ld2+2))
    # U[,2]<-c(rep(0,ld1), sample(c(1,-1),ld2,replace=TRUE),rep(0,p-ld1-ld2))
    # U[,3]<-c(rep(0,ld1+ld2),sample(c(1,-1),ld3,replace=TRUE),rep(0,p-ld1-ld2-ld3))
    
    U[,1]<-c(c(1,  1,  1, -1,  1),rep(0,p-ld1))
    #U[,2]<-c(rep(0,ld1-2), -U[4, 1], U[5, 1], sample(c(1,-1),ld2-2,replace=TRUE),rep(0,p-ld1-ld2+2))
    U[,2]<-c(rep(0,ld1), c(-1,  -1,  -1, 1,  1),rep(0,p-ld1-ld2))
    U[,3]<-c(rep(0,ld1+ld2),c(1,  -1,  1, 1,  1),rep(0,p-ld1-ld2-ld3))
    
    
    # U[, 1] <- c(sample(c(-1, 1), 5, replace = T), rep(0, p - 5))
    # U[, 2] <-
    #   c(rep(0, 3), -U[4, 1], U[5, 1], sample(c(-1, 1), 3, replace = T), rep(0, p -
    #                                                                           8))
    # U[, 3] <- c(rep(0, 8), sample(c(-1, 1), 2, replace = T), rep(0, p - 10))
    U[, 1] <- U[, 1] / lnorm(U[, 1], 2)
    U[, 2] <- U[, 2] / lnorm(U[, 2], 2)
    U[, 3] <- U[, 3] / lnorm(U[, 3], 2)
    
    V[,1]<-c(sample(c(1,-1),ld1,replace=TRUE)*runif(ld1,0.3,1), rep(0,q-ld1))
    V[,2]<-c(rep(0,ld1), sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.3,1), rep(0,q-ld1-ld2))
    V[,3]<-c(rep(0,ld1+ld2),  sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.3,1), rep(0,q-ld1-ld2-ld3))
    
    # V[,1]<-c(c(0.9741170, -0.9316093, -0.7834937,  0.8568272, -0.3172296), rep(0,q-ld1))
    # V[,2]<-c(rep(0,ld1), c(0.3999600,  0.5901824,  0.5896070 , 0.5581918, -0.4067113), rep(0,q-ld1-ld2))
    # V[,3]<-c(rep(0,ld1+ld2), c( 0.3320818,  0.6095401, -0.8592474,  0.3853295,  0.6926636), rep(0,q-ld1-ld2-ld3))
    # 
    U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
    V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
    
    
    #D <- diag(c(100, 50, 20))
    
    
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    ## X <- diag(p)
    
    # UU <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
    # if (rho_E != 0) {
    #   for (t in 1:n)
    #     UU[t, ] <- arima.sim(list(order = c(1, 0, 0), ar = rho_E), n = q)
    # }
    
    sige0 <- rho_E^abs(outer(1:q, 1:q,FUN="-"))
    svdrr <- eigen(sige0)
    svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
    
    UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
    
    C <- U %*% D %*% t(V)
    C3 <- U[, 3] %*% t(V[, 3]) * D[3, 3]
    
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    sigmae <- sigma^2*sige0
    
    list(
      Y = Y,
      X = X,
      C = C,
      U = U,
      V = V,
      D = D,
      Xsigma = Xsigma,
      sigmae = sigmae
    )
    
  }

lnorm <- function(a, p = 1)
{
  (sum(abs(a) ^ p)) ^ (1. / p)
}



rrr.simnew3 <-
  function(n = 50,
           p = 25,
           q = 25,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0,
           Xsigma,
           D = diag(c(100, 50, 20))){
    
    U <- matrix(ncol = nrank, nrow = p)
    V <- matrix(ncol = nrank, nrow = q)
    
    ld1 = 8
    ld2 = 12
    ru1 <- sample(c(1, -1), ld1, replace = TRUE)
    ru2 <- sample(c(1, -1), ld1, replace = TRUE)
    ru3 <- sample(c(1, -1), ld1, replace = TRUE)
    
    ruu1 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    ruu2 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    ruu3 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    
    U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
    
    U[,1]<-c(ru1, rep(0,p-ld1-ld2), ruu1)
    #U[,2]<-c(rep(0,ld1-2), -U[4, 1], U[5, 1], sample(c(1,-1),ld2-2,replace=TRUE),rep(0,p-ld1-ld2+2))
    U[,2]<-c(rep(0,ld1/2), ruu2, ru2, rep(0,p-ld1-ld2-ld1/2))
    U[,3]<-c(rep(0,ld1+ld2), ruu3, ru3, rep(0,p-2*ld1-2*ld2))
    
    
    # U[, 1] <- c(sample(c(-1, 1), 5, replace = T), rep(0, p - 5))
    # U[, 2] <-
    #   c(rep(0, 3), -U[4, 1], U[5, 1], sample(c(-1, 1), 3, replace = T), rep(0, p -
    #                                                                           8))
    # U[, 3] <- c(rep(0, 8), sample(c(-1, 1), 2, replace = T), rep(0, p - 10))
    U[, 1] <- U[, 1] / lnorm(U[, 1], 2)
    U[, 2] <- U[, 2] / lnorm(U[, 2], 2)
    U[, 3] <- U[, 3] / lnorm(U[, 3], 2)
    
    ld3 <- 8
    rv1 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    rv2 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    rv3 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    
    ld4 <- 3
    rvv1 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    rvv2 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    rvv3 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    
    
    V[,1]<-c(rv1, rep(0,q-ld3-ld4), rvv1)
    V[,2]<-c(rep(0,ld3-2), rvv2, rv2,  rep(0,q-2*ld3-ld4+2))
    V[,3]<-c(rep(0,2*ld3+ld4), rv3, rvv3, rep(0,q-3*ld3-2*ld4))
    
    # V[,1]<-c(c(0.9741170, -0.9316093, -0.7834937,  0.8568272, -0.3172296), rep(0,q-ld1))
    # V[,2]<-c(rep(0,ld1), c(0.3999600,  0.5901824,  0.5896070 , 0.5581918, -0.4067113), rep(0,q-ld1-ld2))
    # V[,3]<-c(rep(0,ld1+ld2), c( 0.3320818,  0.6095401, -0.8592474,  0.3853295,  0.6926636), rep(0,q-ld1-ld2-ld3))
    # 
    U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
    V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
    
    
    #D <- diag(c(100, 50, 20))
    
    
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    ## X <- diag(p)
    
    # UU <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
    # if (rho_E != 0) {
    #   for (t in 1:n)
    #     UU[t, ] <- arima.sim(list(order = c(1, 0, 0), ar = rho_E), n = q)
    # }
    
    sige0 <- rho_E^abs(outer(1:q, 1:q,FUN="-"))
    svdrr <- eigen(sige0)
    svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
    
    UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
    
    C <- U %*% D %*% t(V)
    C3 <- U[, 3] %*% t(V[, 3]) * D[3, 3]
    
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    sigmae <- sigma^2*sige0
    
    list(
      Y = Y,
      X = X,
      C = C,
      U = U,
      V = V,
      D = D,
      Xsigma = Xsigma,
      sigmae = sigmae
    )
    
  }

rrr.simnew4 <-
  function(n = 50,
           p = 25,
           q = 25,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0,
           Xsigma,
           D = diag(c(100, 50, 20))){
    
    U <- matrix(ncol = nrank, nrow = p)
    V <- matrix(ncol = nrank, nrow = q)
    
    ld1 = 10
    ld2 = 20
    ru1 <- sample(c(1, -1), ld1, replace = TRUE)
    ru2 <- sample(c(1, -1), ld1, replace = TRUE)
    ru3 <- sample(c(1, -1), ld1, replace = TRUE)
    
    ruu1 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    ruu2 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    ruu3 <- sample(c(1,-1),ld2,replace=TRUE)*runif(ld2,0.01,0.1)
    
    U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
    
    U[,1]<-c(ru1, rep(0,p-ld1-ld2), ruu1)
    #U[,2]<-c(rep(0,ld1-2), -U[4, 1], U[5, 1], sample(c(1,-1),ld2-2,replace=TRUE),rep(0,p-ld1-ld2+2))
    U[,2]<-c( ruu2, ru2,  rep(0,p-ld1-ld2))
    U[,3]<-c(rep(0,p-ld1-ld2), ruu3, ru3)
    
    
    # U[, 1] <- c(sample(c(-1, 1), 5, replace = T), rep(0, p - 5))
    # U[, 2] <-
    #   c(rep(0, 3), -U[4, 1], U[5, 1], sample(c(-1, 1), 3, replace = T), rep(0, p -
    #                                                                           8))
    # U[, 3] <- c(rep(0, 8), sample(c(-1, 1), 2, replace = T), rep(0, p - 10))
    U[, 1] <- U[, 1] / lnorm(U[, 1], 2)
    U[, 2] <- U[, 2] / lnorm(U[, 2], 2)
    U[, 3] <- U[, 3] / lnorm(U[, 3], 2)
    
    ld3 <- 10
    rv1 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    rv2 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    rv3 <- sample(c(1,-1),ld3,replace=TRUE)*runif(ld3,0.6,1)
    
    ld4 <- 10
    rvv1 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    rvv2 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    rvv3 <- sample(c(1,-1),ld4,replace=TRUE)*runif(ld4,0.01,0.1)
    
    
    V[,1]<-c(rv1, rep(0,q-ld3-ld4), rvv1)
    V[,2]<-c(rvv2, rv2,  rep(0,q-ld3-ld4))
    V[,3]<-c(rep(0,q-ld3-ld4), rvv3, rv3)
    
    # V[,1]<-c(c(0.9741170, -0.9316093, -0.7834937,  0.8568272, -0.3172296), rep(0,q-ld1))
    # V[,2]<-c(rep(0,ld1), c(0.3999600,  0.5901824,  0.5896070 , 0.5581918, -0.4067113), rep(0,q-ld1-ld2))
    # V[,3]<-c(rep(0,ld1+ld2), c( 0.3320818,  0.6095401, -0.8592474,  0.3853295,  0.6926636), rep(0,q-ld1-ld2-ld3))
    # 
    U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
    V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
    
    
    #D <- diag(c(100, 50, 20))
    
    
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    ## X <- diag(p)
    
    # UU <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
    # if (rho_E != 0) {
    #   for (t in 1:n)
    #     UU[t, ] <- arima.sim(list(order = c(1, 0, 0), ar = rho_E), n = q)
    # }
    
    sige0 <- rho_E^abs(outer(1:q, 1:q,FUN="-"))
    svdrr <- eigen(sige0)
    svdrrinv <- svdrr$vectors%*%diag(svdrr$values^0.5,nrow=q)%*%t(svdrr$vectors)
    
    UU <- matrix(nrow=n,ncol=q,rnorm(n*q,0,1))%*%svdrrinv
    
    C <- U %*% D %*% t(V)
    C3 <- U[, 3] %*% t(V[, 3]) * D[3, 3]
    
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    sigmae <- sigma^2*sige0
    
    list(
      Y = Y,
      X = X,
      C = C,
      U = U,
      V = V,
      D = D,
      Xsigma = Xsigma,
      sigmae = sigmae
    )
    
  }



bhq <- function(u, alpha = 0.05) {
  n = length(u)
  r = rank(u, ties.method = "max")
  bh = max(c(r[u <= (r/n) * alpha], 0), na.rm = T)
  su <- sort(u)
  jj <- which(u == 1)
  if (length(jj) != 0) 
    pi0 <- 1
  else pi0 <- min((-1/n) * sum(log(1 - u)), 1)
  if (bh == 0) {
    FDR_BH <- 0
  }
  else {
    FDR_BH <- round((pi0 * su[bh])/(ecdf(u)(su[bh])), 
                    4)
  }
  ad.p = numeric(n)
  ad.p[n] <- sort(u)[n]
  for (i in (n - 1):1) {
    ad.p[i] <- min(sort(u)[i] * (n/i), ad.p[i + 1])
  }
  names(ad.p) = names(sort(u))
  return(c(list(Rejections = bh, FDR = min(FDR_BH, 1), 
                Adjusted.pvalues = sort(ad.p))))
}

##P-VALUE
bhq_adj_p <- function(uk, varuk, q = 0.1, n, p){
  u1_std = rep(0,p)
  pval = rep(0,p)
  for(j in 1:p){
    sej <- sqrt(varuk[j])/sqrt(n)
    u1_std[j] = uk[j]/sej
    pval[j] <- 2 * pnorm( abs(u1_std[j]), lower.tail = FALSE)
  }
  names(pval) = c(1:p)
  res = bhq(pval, q)
  num_rej = res$Rejections
  num_rej
  ad.pva = res$Adjusted.pvalues
  
  sele = ad.pva[1:num_rej]
  num = rep(0, num_rej)
  for(j in 1:num_rej){
    num[j] = which(names(pval) == names(sele)[j])
  }
  S1 = num
  St = sort(S1)
  return(list(St = St, num = length(St), ad.pva = ad.pva, selection = sele))
  
}




sum_simu_result_para <- function(res){
  #this function summarizes the simulation results of parallel computing
  ##24 mean results for all
  result1 <- res[,1]
  ##the normal variables for u1, u2, u3 and d1, d2, d3
  re_u1 <- res[,2]; re_u2 <- res[,3]; re_u3 <- res[,4]
  re_d <- res[,5]
  ##individual CI and length for u1j, u2j, u3j
  ci_u1j <- res[,6]; len_u1j <- res[,7]
  ci_u2j <- res[,8]; len_u2j <- res[,9]
  ci_u3j <- res[,10]; len_u3j <- res[,11]
  #########first index
  iter = length(result1)
  print(iter)
  result1.matrix <- matrix(0, iter, 24)
  for (i in 1:iter) {
    result1.matrix[i,] = result1[[i]]
  }
  colnames(result1.matrix)<- c("U1covs0","U1covs1","U1covall","U1len_s0", "U1len_s1","U1lenall",
                               "U2covs0","U2covs1","U2covall","U2len_s0", "U2len_s1","U2lenall",
                               "U3covs0","U3covs1","U3covall","U3len_s0", "U3len_s1","U3lenall",
                               "d1covs0","d1len_s0",
                               "d2covs0","d2len_s0",
                               "d3covs0","d3len_s0")
  #cat(round(colMeans(result1.matrix),3), "\n")
  print(round(colMeans(result1.matrix),3))
  #######normal plot
  pp = length(re_u1[[1]])
  re_u1.m <- matrix(0, iter, pp); re_u2.m <- matrix(0, iter, pp); re_u3.m <- matrix(0, iter, pp)
  for (i in 1:iter) {
    re_u1.m[i,] = re_u1[[i]]
    re_u2.m[i,] = re_u2[[i]]
    re_u3.m[i,] = re_u3[[i]]
  }
  # hist(re_u2.m[,4], freq = F)
  # hist(re_u1.m[,1], xlim = c(-3,3),  freq = F)
  # hist(re_u2.m[,1], xlim = c(-3,3), breaks =  20, freq = F)
  # hist(re_u3.m[,1], xlim = c(-3,3), breaks =  10, freq = F)
  
  re_d.m <- matrix(0, iter, 3)
  for (i in 1:iter) {
    re_d.m[i,] = re_d[[i]]
  }
  # hist(re_d.m[,1],  breaks =  20, freq = F)
  # hist(re_d.m[,2], xlim = c(-3,3), breaks =  20, freq = F)
  # hist(re_d.m[,3], xlim = c(-3,3), breaks =  15, freq = F)
  
  
  ###CI and length
  ci_u1j.matrix <- matrix(0, iter, pp); ci_u2j.matrix <- matrix(0, iter, pp);  ci_u3j.matrix <- matrix(0, iter, pp)
  len_u1j.matrix <- matrix(0, iter, pp); len_u2j.matrix <- matrix(0, iter, pp);  len_u3j.matrix <- matrix(0, iter, pp)
  for (i in 1:iter) {
    ci_u1j.matrix[i,] = ci_u1j[[i]]
    ci_u2j.matrix[i,] = ci_u2j[[i]]
    ci_u3j.matrix[i,] = ci_u3j[[i]]
    len_u1j.matrix[i,] = len_u1j[[i]]
    len_u2j.matrix[i,] = len_u2j[[i]]
    len_u3j.matrix[i,] = len_u3j[[i]]
  }
  
  cat("Coverage of u1:\n")
  print(colMeans(ci_u1j.matrix))
  cat("Coverage of u2:\n")
  print(colMeans(ci_u2j.matrix))
  cat("Coverage of u3:\n")
  print(colMeans(ci_u3j.matrix))
  
  cat("Length of u1:\n")
  print(round(colMeans(len_u1j.matrix),3))
  cat("Length of u2:\n")
  print(round(colMeans(len_u2j.matrix),3))
  cat("Length of u3:\n")
  print(round(colMeans(len_u3j.matrix),3))
  # print(colMeans(ci_u1j.matrix))
  # print(colMeans(ci_u2j.matrix))
  # print(colMeans(ci_u3j.matrix))
  # print(round(colMeans(len_u1j.matrix),3))
  # print(round(colMeans(len_u2j.matrix),3))
  # print(round(colMeans(len_u3j.matrix),3))
  
  
}

