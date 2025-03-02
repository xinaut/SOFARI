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
    #pval[j] <- 2 * pnorm( abs(u1_std[j]), lower.tail = FALSE)
    pval[j] <- 2 *( 1 - pnorm( abs(u1_std[j])))
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


load(file = "real20230912.RData")
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
                          ncores = 12,
                          oldschool = FALSE,
                          lambdatuningfactor = 1,
                          cv.verbose = FALSE,
                          do.ZnZ = TRUE)
Theta = obj$out
#write.csv(Theta, file = "finalXY_theta.csv", row.names = F)
Theta <- InverseLinfty(Xsigma, n,  resol=1.3, maxiter=100, threshold=1e-6)
asd = cv.glasso_clean(X,0)
Theta = asd$wi
# asd = CVglasso(X,  cores = 5)
# Theta= asd$Omega 
# regfactor = "log"
# npermu = 1
# sis.use = 0
# bia.cor = 0
# obj = isee(X, regfactor, 5, 0, 0) 
# if (bia.cor == 1){
#   # ISEE with bias correction
#   Omega = obj$Omega.isee.c
# } else {
#   #  # ISEE with no bias correction
#   Omega = obj$Omega.isee
# }

norm(diag(p)- Xsigma%*%Theta, "M")
write.csv(Theta, file = "final_theta.csv", row.names = F)
Theta_weig_vard = Theta
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


Ur = fit1$U; Dr = fit1$D; V = fit1$V;
dU = Ur%*%diag(Dr)
r = fit1$rank

C = Ur%*%diag(Dr)%*%t(V)
E_est = Y - X%*%C

Sigmae_null = matrix(0, q, q) #no meaning
obj_sigmae <-  sim(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
Sigmae <- obj_sigmae$cov_est

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
zz = 1.96
# zz= 2.33

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
abs(up_ci - low_ci)


####the third layer
k = 3
length(which(Ur[,3]!=0))
uk1_t <- uk1
uk1_t[abs(uk1) < log(n)/sqrt(n)] = 0 
S_u1 = which(uk1_t != 0) # S_u is no used parameter 
# S_u1 
# uk1[S_u1]

#u1_wh <- uk1_t
u1_wh <- dU[,1]
#(Dr[1]/Dr[2])*(log(n)/sqrt(n))

uk2_t <- uk2 # u_wh is no used parameter
uk2_t[abs(uk2) < log(n)/sqrt(n)] = 0 # u_wh is no used parameter
S_u2 = which(uk2_t != 0)
#S_u2
#uk2[S_u2]
S_u <- list(S_u1, S_u2)
#u2_wh <- uk2_t
u2_wh <- dU[,2]
u_wh <- list(u1_wh, u2_wh)
ukdk3_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh) # u_wh is no used parameter
#ukdk3_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
uk3 <- ukdk3_weak$u
dk3 <- ukdk3_weak$d
dk3
#varobj3 <- var_weakly_apply(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)

varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u) # S_u is no used parameter
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

uk = uk2
varuk = varuk2
u1_std = rep(0,p)
p_values = rep(0,p)
for(j in 1:p){
  sej <- sqrt(varuk[j])/sqrt(n)
  u1_std[j] = uk[j]/sej
  #pval[j] <- 2 * pnorm( abs(u1_std[j]), lower.tail = FALSE)
  p_values[j] <- 2 *( 1 - pnorm( abs(u1_std[j])))
}
bhq1 <- p.adjust(p_values, method = "BH")
# 输出调整后的p值
print(sort(bhq1))
# 找出那些调整后的p值小于等于q的假设，表示可以拒绝原假设
significant_hypotheses <- bhq1 <= 0.05
length(which(significant_hypotheses == TRUE))

St12 = union( res1$St, res2$St)
S123 = union(St12,  res3$St)
length(S123)

save.image(file = "real20230912.RData")



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
# uk1_f[nonzero_uk1] <- standarize_data(as.matrix(uk1_f[nonzero_uk1]))
# uk2_f[nonzero_uk2] <- standarize_data(as.matrix(uk2_f[nonzero_uk2]))
# uk3_f[nonzero_uk3] <- standarize_data(as.matrix(uk3_f[nonzero_uk3]))
datauu <- cbind(uk1_f,uk2_f,uk3_f)

#write.csv(datauu, file = "plotUU.csv")

save.image(file = "realdata20231006.RData")