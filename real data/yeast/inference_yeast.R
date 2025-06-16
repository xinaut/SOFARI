####estimation of all data
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
library(glmnet)
library(parallel)


load("real data/yeast/yeast_preprocess_data.RData")
n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]


##sofar estimation
est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 300, rank_tol = 1e-3, C = 2.01); est_rank

r = est_rank 
fit1 <- sofar(Y, X, ic.type = "GIC", nrank = r,
              control = list(methodA = "adlasso", methodB = "adlasso", lam.AB.factor = 0.1,
                             nlam = 300, lam.max.factor = 5000, lam.min.factor = 1e-12,epsilon = 1e-8), screening = TRUE); summary(fit1)
U <- fit1$U
V <- fit1$V
length(which(U!=0)) 
length(which(V!=0)) 
fit1$D
##rank is 3


#####inference from the all data#####
Ur = fit1$U; Dr = fit1$D; V = fit1$V;
dU = Ur%*%diag(Dr)
r = fit1$rank

C = Ur%*%diag(Dr)%*%t(V)
E_est = Y - X%*%C

Sigmae_null = matrix(0, q, q) #no meaning
obj_sigmae <-  ada.thresholding(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
Sigmae <- obj_sigmae$cov_est


Xsigma <- t(X)%*%X/n
#estimate the precision matrix
# obj<- score.nodewiselasso(X, wantTheta = TRUE,
#                           verbose = FALSE,
#                           lambdaseq = "quantile",
#                           parallel = FALSE,
#                           ncores = 12,
#                           oldschool = FALSE,
#                           lambdatuningfactor = 1,
#                           cv.verbose = FALSE,
#                           do.ZnZ = TRUE)
# Theta = obj$out
#write.csv(Theta, file = "finalXY_theta.csv", row.names = F)
Theta <- InverseLinfty(Xsigma, n,  resol=1.3, maxiter=100, threshold=1e-6)
#asd = cv.glasso_clean(X,0)
#Theta = asd$wi
# asd = CVglasso(X,  cores = 5)

#norm(diag(p)- Xsigma%*%Theta, "M")


n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
r = 3
Xsigma <- t(X)%*%X/n

##############start debias procedure


Ur = fit1$U; Dr = fit1$D; V = fit1$V;
dU = Ur%*%diag(Dr)
r = fit1$rank

C = Ur%*%diag(Dr)%*%t(V)
E_est = Y - X%*%C

Sigmae_null = matrix(0, q, q) #no meaning
obj_sigmae <-  ada.thresholding(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
Sigmae <- obj_sigmae$cov_est


#############read the output data to quickly show inference results#############
load("real data/yeast/yeast_inference.RData")

uk0 = rep(0,p)
dk0 = 0
#de-biased 
#the first layer
k = 1
u_wh = list()
ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
#ukdk1_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
uk1 <- ukdk1_weak$u
dk1 <- ukdk1_weak$d
dk1

varobj1 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae,  S_u)

# e = 1
# a1 = rep(0, 2-1)
# a2 = rep(0, p-2)
# b = c(a1,e,a2)
# varobj1 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
varuk1 <- varobj1$varu
mean(varuk1)
vardk1 <- varobj1$vard

length(which(Ur[,1]!=0))
length(which(Ur[,2]!=0))
length(which(Ur[,3]!=0))
uk = uk1; varuk = varuk1
index = NULL
len2 = 0
zz = 2.81



dk = dk1; vardk = vardk1
low_ci <- dk - sqrt(vardk)*zz/sqrt(n)
up_ci <- dk + sqrt(vardk)*zz/sqrt(n)
len_d1 <- abs(up_ci - low_ci)


low <- rep(0,p); up <- rep(0,p)
for(j in 1:p){
  low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n) #[which(dU[,1]==0)]
  low[j] <- low_ci
  up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
  up[j] <- up_ci
  len2 = len2 + abs(up_ci - low_ci)
  if( up_ci < 0 || low_ci >0) {index = c(index, j)}
}
length(index)
len2/p
index
len_uk1 = len2/p
nonzero_uk1 = index



#the second layer
k = 2
length(which(Ur[,2]!=0))
uk1_t <- uk1
uk1_t[abs(uk1) < log(n)/sqrt(n)] = 0 #no use
S_u1 = which(uk1_t != 0)  # S_u is no used parameter 
S_u1 # S_u is no used parameter
S_u <- list(S_u1) # S_u is no used parameter
u1_wh <- uk1_t
u1_wh <- as.vector(dU[,1]) 
u_wh <- list(u1_wh)
ukdk2_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)  # u_wh is no used parameter
#ukdk2_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
uk2 <- ukdk2_weak$u
dk2 <- ukdk2_weak$d
dk2


varobj2 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u) # S_u is no used parameter
# e = 1
# a1 = rep(0, 2-1)
# a2 = rep(0, p-2)
# b = c(a1,e,a2)
# varobj2 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
varuk2 <- varobj2$varu
vardk2 <- varobj2$vard

uk = uk2; varuk = varuk2
index = NULL
len2 = 0
for(j in 1:p){
  low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n)
  up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
  len2 = len2 + abs(up_ci - low_ci)
  if( up_ci < 0 || low_ci > 0) {index = c(index, j)}
}
length(index)
len2/p
index
nonzero_uk2 = index
len_uk2 = len2/p

dk = dk2; vardk = vardk2
low_ci <- dk - sqrt(vardk)*zz/sqrt(n)
up_ci <- dk + sqrt(vardk)*zz/sqrt(n)
len_d2 <- abs(up_ci - low_ci)


####the third layer
k = 3
length(which(Ur[,3]!=0))
uk1_t <- uk1
uk1_t[abs(uk1) < log(n)/sqrt(n)] = 0 
S_u1 = which(uk1_t != 0) # S_u is no used parameter 

u1_wh <- dU[,1]
uk2_t <- uk2 # u_wh is no used parameter
uk2_t[abs(uk2) < log(n)/sqrt(n)] = 0 # u_wh is no used parameter
S_u2 = which(uk2_t != 0)
S_u <- list(S_u1, S_u2)
u2_wh <- dU[,2]
u_wh <- list(u1_wh, u2_wh)
ukdk3_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh) # u_wh is no used parameter
uk3 <- ukdk3_weak$u
dk3 <- ukdk3_weak$d


varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u) # S_u is no used parameter
varuk3 <- varobj3$varu
vardk3 <- varobj3$vard

uk = uk3; varuk = varuk3
index = NULL
len2 = 0
for(j in 1:p){
  low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n)
  up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
  len2 = len2 + abs(up_ci - low_ci)
  if( up_ci < 0 || low_ci > 0) {index = c(index, j)}
}
length(index)
len2/p
index
nonzero_uk3 = index
len_uk3 = len2/p

dk = dk3; vardk = vardk3
low_ci <- dk - sqrt(vardk)*zz/sqrt(n)
up_ci <- dk + sqrt(vardk)*zz/sqrt(n)
len_d3 <- abs(up_ci - low_ci)


qq = 0.05
res1 <- bhq_adj_p(uk1, varuk1, qq , n, p)
res1$St
res1$num
# res1$ad.pva
# res1$selection

res2 <- bhq_adj_p(uk2, varuk2, qq , n, p)
res2$St
res2$num
# res2$ad.pva
# res2$selection
which(Ur[,2]!=0)


res3 <- bhq_adj_p(uk3, varuk3, qq , n, p)
res3$St
res3$num
# res3$ad.pva
# res3$selection

St12 = union( res1$St, res2$St)
S123 = union(St12,  res3$St)



#############output results##########

################## Table 8 #################
cat("singular values:\n",dk1,dk2,dk3,"\n")
cat("length of singular values:\n",len_d1,len_d2,len_d3,"\n")

###number of features
cat("number of features in three layers:", res1$num, res2$num, res3$num,"\n")
cat("all nonzeors:", res1$num + res2$num + res3$num,"\n")
cat("distinct nonzeros:", length(S123))


save.image("yeast_inference.RData")
