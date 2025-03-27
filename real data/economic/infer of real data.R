load(file = "realdata1.RData")
#estimate
X <- read.csv("final_XX.csv")
Y <- read.csv("final_YY.csv")
X <- as.matrix(X); Y <- as.matrix(Y)
dim(X)
dim(Y)
 

#####estimate and inference from the all data#####
#X <- X0; Y<-Y0
est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 200, rank_tol = 1e-2, C = 2.01); est_rank

n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
r = est_rank #r = 3
fit1 <- sofar(Y, X, ic.type = "BIC", nrank = r,
              control = list(methodA = "adlasso", methodB = "adlasso",
                             nlam = 200, lam.max.factor = 1500, lam.min.factor = 1e-12,
                             penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
              ), screening = FALSE);  summary(fit1)
fit1$lam.id

Ur = fit1$U; Dr = fit1$D; V = fit1$V;
dU = Ur%*%diag(Dr)
r = fit1$rank

C = Ur%*%diag(Dr)%*%t(V)
E_est = Y - X%*%C

Sigmae_null = matrix(0, q, q) #no meaning
obj_sigmae <-  sim(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
Sigmae <- obj_sigmae$cov_est
write.csv(Sigmae, file = "final_Sigmae.csv", row.names = F)
#Sigmae <- t(E_est)%*%E_est/q

Xsigma <- t(X)%*%X/n
#estimate the precision matrix
obj<- score.nodewiselasso(X, wantTheta = TRUE,
                          verbose = FALSE,
                          lambdaseq = "quantile",
                          parallel = FALSE,
                          ncores = 6,
                          oldschool = FALSE,
                          lambdatuningfactor = 1,
                          cv.verbose = FALSE,
                          do.ZnZ = TRUE)
Theta = obj$out
#write.csv(Theta, file = "finalXY_theta.csv", row.names = F)
#Theta <- InverseLinfty(Xsigma, n,  resol=1.5, maxiter=50, threshold=1e-2)
#norm(diag(p)- Xsigma%*%Theta, "M")
write.csv(Theta, file = "final_theta.csv", row.names = F)

#############read data#############

n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
r = 3
Xsigma <- t(X)%*%X/n
Sigmae <-  read.csv("final_Sigmae.csv")
Sigmae <- as.matrix(Sigmae)

Theta <- read.csv("final_theta.csv")
Theta <- as.matrix(Theta)

##read the initital SOFAR estimate
dU = read.csv("dU.csv")
V = read.csv("V.csv")
dU = as.matrix(dU)
V = as.matrix(V)

##############start debias procedure
uk0 = rep(0,p)
dk0 = 0
#de-biased 
#the first layer
k = 1
u_wh = list()
ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
uk1 <- ukdk1_weak$u
dk1 <- ukdk1_weak$d
dk1

varobj1 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae)
varuk1 <- varobj1$varu
mean(varuk1)
vardk1 <- varobj1$vard


uk = uk1; varuk = varuk1
index = NULL
len2 = 0
zz = 2.81


dk = dk1; vardk = vardk1
low_ci <- dk - sqrt(vardk)*zz/sqrt(n)
up_ci <- dk + sqrt(vardk)*zz/sqrt(n)
abs(up_ci - low_ci)


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
u1_wh <- as.vector(dU[,1]) 
u_wh <- list(u1_wh)
ukdk2_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)  
#ukdk2_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
uk2 <- ukdk2_weak$u
dk2 <- ukdk2_weak$d
dk2


varobj2 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae) 
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
abs(up_ci - low_ci)


####the third layer
k = 3

u1_wh <- dU[,1]
u2_wh <- dU[,2]
u_wh <- list(u1_wh, u2_wh)
ukdk3_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh) 
#ukdk3_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
uk3 <- ukdk3_weak$u
dk3 <- ukdk3_weak$d
dk3


varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae) 
#varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)
#varobj3 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
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
abs(up_ci - low_ci)

qq = 0.05
res1 <- bhq_adj_p(uk1, varuk1, qq , n, p)
res1$St
res1$num
res1$ad.pva
res1$selection

res2 <- bhq_adj_p(uk2, varuk2, qq , n, p)
res2$St
res2$num
res2$ad.pva
res2$selection

res3 <- bhq_adj_p(uk3, varuk3, qq , n, p)
res3$St
res3$num
res3$ad.pva
res3$selection


St12 = union( res1$St, res2$St)
S123 = union(St12,  res3$St)
length(S123)

save.image(file = "realdata1.RData")



len_uk1
len_uk2
len_uk3

nonzero_uk1 
nonzero_uk2 
nonzero_uk3 
length(nonzero_uk1)
length(nonzero_uk2)
length(nonzero_uk3)

##save data
ukK = cbind(uk1,uk2,uk3)
write.csv(ukK, file = "final_Uk.csv", row.names = F)
write.csv(dU, file = "dU.csv", row.names = F)
write.csv(V, file = "V.csv", row.names = F)

nonzero_uk1 = res1$St
nonzero_uk2 = res2$St
nonzero_uk3 = res3$St


datau <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
datau[1, nonzero_uk1] = 1
datau[2, nonzero_uk2] = 1
datau[3, nonzero_uk3] = 1
write.csv(datau, file = "plotX.csv")
uk1_f <- uk1; uk1_f[-nonzero_uk1] = 0
uk2_f <- uk2; uk2_f[-nonzero_uk2] = 0
uk3_f <- uk3; uk3_f[-nonzero_uk3] = 0
datauu <- cbind(uk1_f,uk2_f,uk3_f)

#write.csv(datauu, file = "plotUU.csv")

save.image(file = "realdata2.RData")