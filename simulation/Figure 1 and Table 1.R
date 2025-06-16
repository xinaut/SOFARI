###example 1 
### Set working directory to the `SOFARI_code` folder of the reproducibility_materials
### Please replace the path below with your own path containing the `SOFARI_code` folder, remove the `#`, and then run the following line
# setwd("~/SOFARI_code")  ### <-- Replace this with your own path

source("functions/func_strongly.R")
source("functions/func_weakly.R")
source("functions/adaptive_thresholding.R")
source("functions/STRS.R")
source("functions/nodewise_lasso.R")
source("functions/functions_helper.R")
library(rrpack)
library(secure)
library(MASS)
library(glmnet)
library(foreach)
library(doParallel) 
library(flare)
library(expm)

##setting of simulation example 1
simu_setup_1 <- function(sss){
  n <- 200; p <- 25; q <- 15; nrank <- 3
  
  ld1 = 5
  ld2 = 5
  ld3 = 5
  
  
  U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
  U[,1]<-c(sample(c(1,-1),ld1,replace=TRUE),rep(0,p-ld1))
  U[,2]<-c(rep(0,ld1), sample(c(1,-1),ld2,replace=TRUE),rep(0,p-ld1-ld2))
  U[,3]<-c(rep(0,ld1+ld2),sample(c(1,-1),ld3,replace=TRUE),rep(0,p-ld1-ld2-ld3))
  
  V[,1]<-c(sample(c(1,-1),3,replace=TRUE)*runif(3,0.3,1), rep(0,q-3))
  V[,2]<-c(0, 0, 0, sample(c(1,-1),3,replace=TRUE)*runif(3,0.3,1), rep(0,q-6))
  V[,3]<-c(rep(0,6), sample(c(1,-1),3,replace=TRUE)*runif(3,0.3,1), rep(0,q-9))
  
  
  U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
  V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
  D <- diag(c(100,15,5))
  C <- U%*%D%*%t(V)
  
  U0 = U; D0 = D; V0 = V; dU0 = U0%*%D0
  
  #generate X and Y
  xrho = 0.3
  #obj_Xsigma_real <- blk.mat(0.3, p, 1:p)
  #obj_Xsigma_real <- band.mat(0.3, p, 1:p)
  obj_Xsigma_real <- ar1.mat(0.3, p, 1:p)
  Xsigma_real <- obj_Xsigma_real$Sigma
  #Xsigma_real <- diag(p)
  #Xsigma <- Xsigma_real
  #Xsigma_real <- xrho^abs(outer(1:p, 1:p,FUN="-"))
  sim.sample <- porth.sim(U0,D0,V0,n,snr = 1,Xsigma_real,rho=0.3)
  #sim.sample <- porth.simide(U0,D0,V0,n,snr = 1,Xsigma_real,rho=0.3)
  Y <- sim.sample$Y; X <- sim.sample$X; Sigmae_real <- sim.sample$sigmae
  est_rank <- STRS(Y, X)
  # est_rank
  # fit1 <- sofar(Y, X, ic.type = "BIC", nrank =  est_rank, 
  #               control = list(methodA = "adlasso", methodB = "adlasso",
  #              nlam = 100, lam.max.factor = 2000, lam.min.factor = 1e-6,epsilon = 1e-6)); summary(fit1)
  #generate X and Y   #initial estimate
  fit1 <- sofar(Y, X, ic.type = "BIC", nrank = est_rank,
                control = list(methodA = "adlasso", methodB = "adlasso",
                               nlam = 100, lam.max.factor = 2000, lam.min.factor = 1e-6,epsilon = 1e-6))
  #summary(fit1)
  Ur = fit1$U; Dr = fit1$D; V = fit1$V;
  Ur = transuv(Ur,U0)
  V = transuv(V,V0)
  dU = Ur%*%diag(Dr)
  r = fit1$rank
  
  C = Ur%*%diag(Dr)%*%t(V)
  E_est = Y - X%*%C
  obj_sigmae <-  ada.thresholding(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_real)
  Sigmae <- obj_sigmae$cov_est
  
  
  Xsigma <- t(X)%*%X/n
  #estimate the precision matrix
  #Theta <- InverseLinfty(Xsigma, n, resol=1.5, maxiter=50, threshold=1e-6, verbose = F)
  #Theta <- solve(Xsigma)
  # ghat=cv.glasso1(X,0)
  # Theta=ghat$wi
  obj<- score.nodewiselasso(X,  wantTheta = TRUE,
                            verbose = FALSE,
                            lambdaseq = "quantile",
                            parallel = FALSE,
                            ncores = 10,
                            oldschool = FALSE,
                            lambdatuningfactor = 1,
                            cv.verbose = FALSE,
                            do.ZnZ = TRUE)
  Theta = obj$out
  
  
  ##############start debias procedure
  #de-biased 
  
  #the first layer
  k = 1
  u_wh = list()
  ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
  #ukdk1_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
  uk1 <- ukdk1_weak$u
  dk1 <- ukdk1_weak$d
  
  
  uk0 <- dU0[,k]
  dk0 <- D0[k,k]
  
  varobj1 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae)
  varuk1 <- varobj1$varu
  vardk1 <- varobj1$vard
  objavg1 <- avgcov(n, uk1, varuk1, uk0)
  objd1 <- avgd(n, dk1, vardk1, dk0)
  objavg1$cov_s0; objavg1$cov_s1; objavg1$len_s0; objavg1$len_s1
  objd1$num; objd1$len
  
  U1_cs0 <- objavg1$cov_s0; U1_cs1 <- objavg1$cov_s1; U1_ls0 <- objavg1$len_s0
  U1_ls1 <- objavg1$len_s1; U1_call <- objavg1$cov_all;   U1_lall <- objavg1$len_all
  d1_cs0 = objd1$num; d1_ls0 = objd1$len
  U1_num = objavg1$num; U1_len = objavg1$len
  U1_num; U1_len
  
  
  
  #the second layer
  k = 2
  u_wh <- list(dU[ ,1])
  ukdk2_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
  #ukdk2_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
  uk2 <- ukdk2_weak$u
  dk2 <- ukdk2_weak$d
  
  
  uk0 <- dU0[,k]
  dk0 <- D0[k,k]
  
  varobj2 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae)
  #varobj2 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
  varuk2 <- varobj2$varu
  vardk2 <- varobj2$vard
  objavg2 <- avgcov(n, uk2, varuk2, uk0)
  objd2 <- avgd(n, dk2, vardk2, dk0)
  objavg2$cov_s0; objavg2$cov_s1; objavg2$len_s0; objavg2$len_s1
  objd2$num; objd2$len
  
  U2_cs0 <- objavg2$cov_s0; U2_cs1 <- objavg2$cov_s1; U2_ls0 <- objavg2$len_s0
  U2_ls1 <- objavg2$len_s1; U2_call <- objavg2$cov_all; U2_lall <- objavg2$len_all
  d2_cs0 = objd2$num; d2_ls0 = objd2$len
  U2_num = objavg2$num; U2_len = objavg2$len
  U2_num; U2_len
  
  
  
  #the third layer
  k = 3
  u1_wh <- dU[,1]
  u2_wh <- dU[,2]
  u_wh <- list(u1_wh, u2_wh)
  ukdk3_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
  #ukdk3_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
  uk3 <- ukdk3_weak$u
  dk3 <- ukdk3_weak$d
  
  
  
  uk0 <- dU0[,k]
  dk0 <- D0[k,k]
  
  varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae)
  #varobj3 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
  varuk3 <- varobj3$varu
  vardk3 <- varobj3$vard
  objavg3 <- avgcov(n, uk3, varuk3, uk0)
  objd3 <- avgd(n, dk3, vardk3, dk0)
  objavg3$cov_s0; objavg3$cov_s1; objavg3$len_s0; objavg3$len_s1
  objd3$num; objd3$len
  
  U3_cs0 <- objavg3$cov_s0; U3_cs1 <- objavg3$cov_s1; U3_ls0 <- objavg3$len_s0
  U3_ls1 <- objavg3$len_s1; U3_call <- objavg3$cov_all; U3_lall <- objavg3$len_all
  d3_cs0 = objd3$num; d3_ls0 = objd3$len
  U3_num = objavg3$num; U3_len = objavg3$len
  U3_num; U3_len
  
  
  reud <- normal_u(n, p, dU0, D0, uk1, uk2, uk3, dk1, dk2, dk3, varuk1, varuk2, varuk3, vardk1, vardk2, vardk3)
  u1_nor <- reud$nu1; u2_nor <- reud$nu2; u3_nor <- reud$nu3
  d1_nor <- reud$nd1; d2_nor <- reud$nd2; d3_nor <- reud$nd3
  
  result1 =c(round(mean(U1_cs0),3),round(mean(U1_cs1),3),round(mean(U1_call),3),round(mean(U1_ls0),3), round(mean(U1_ls1),3), round(mean(U1_lall),3),
             round(mean(U2_cs0),3),round(mean(U2_cs1),3),round(mean(U2_call),3),round(mean(U2_ls0),3), round(mean(U2_ls1),3), round(mean(U2_lall),3),
             round(mean(U3_cs0),3),round(mean(U3_cs1),3),round(mean(U3_call),3),round(mean(U3_ls0),3), round(mean(U3_ls1),3), round(mean(U3_lall),3),
             round(mean(d1_cs0),3),round(mean(d1_ls0),3),
             round(mean(d2_cs0),3),round(mean(d2_ls0),3),
             round(mean(d3_cs0),3),round(mean(d3_ls0),3)
  )
  result_d = c(d1_nor, d2_nor, d3_nor)
  result = list(result1, u1_nor, u2_nor, u3_nor, result_d, 
                U1_num, U1_len, U2_num, U2_len, U3_num, U3_len)
  
  
  
}




###parallel
start <- (proc.time())[3][[1]]
no_cores <- detectCores(no_cores)-4
cl <- makeCluster(no_cores)
registerDoParallel(cl)
set.seed(12345)
iter = 1000
result <- foreach( ii = 1:iter, .combine = 'rbind', .errorhandling = "remove",
                .packages = c("rrpack", "MASS","glmnet","parallel"))%dopar% simu_setup_1(ii)
stopCluster(cl)
end <- (proc.time())[3][[1]]
print(paste(' time = ', (end-start), 's', sep=''))




###Table 1
temp_res <- sum_simu_result_para(result)



###Figure 1
##obtain the normal variables for u1, u2, u3 and d1, d2, d3
res <- result
re_u1 <- res[,2]; re_u2 <- res[,3]; re_u3 <- res[,4]
re_d <- res[,5]
pp = length(re_u1[[1]])
re_u1.m <- matrix(0, iter, pp); re_u2.m <- matrix(0, iter, pp); re_u3.m <- matrix(0, iter, pp)
re_d.m <- matrix(0, iter, 3)
for (i in 1:iter) {
  re_u1.m[i,] = re_u1[[i]]
  re_u2.m[i,] = re_u2[[i]]
  re_u3.m[i,] = re_u3[[i]]
  re_d.m[i,] = re_d[[i]]
}

library(latex2exp)

par(mfrow=c(3, 3)) #par(mfrow=c(1, 1))
plot(density(re_u1.m[,1]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(1, ',', 1)])),
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_u1.m[,25]), xlim = c(-3,3), ylim = c(0,0.4), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(1, ',', 25)])), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_d.m[,1]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T)[d[1]]), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_u2.m[,4]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(2, ',', 4)])), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_u2.m[,25]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(2, ',', 25)])), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_d.m[,2]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T)[d[2]]), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_u3.m[,7]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(3, ',', 7)])), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_u3.m[,25]), xlim = c(-3,3), ylim = c(0,0.4), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T[paste(3, ',', 25)])), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)

plot(density(re_d.m[,3]), xlim = c(-3,3), col='blue', lwd = 2, axes=F,
     xlab = expression(symbol(T)[d[3]]), 
     ylab = "", main = "",cex.lab=1.8)
axis(1,at=(-3):3,  lwd = 1.5, cex.lab=1, line = -0.8)
title( ylab = "Density", cex.lab=1.3)
axis(2,at=seq(0,0.4,0.1),  lwd = 1.5)
curve(dnorm, from=-3.2, to=3.2, add=T, col="red", lwd=2)




