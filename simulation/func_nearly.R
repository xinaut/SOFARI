##1. strongly case
#M
Mk_nearly <- function(k, Xsigma, dU, V){
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  C_k <- dU[, -k]%*%t(V[, -k])
  M <- -Xsigma%*%C_k/z_kk 
  return(M)
}

#W
wk_nearly <- function(k, r, p, Theta, dU, Xsigma){
  bTheta = Theta
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  A = diag(r-1) - t(dU[, -k])%*%Xsigma%*%dU[, -k]/z_kk
  a1 = (Xsigma%*%dU[, -k])/z_kk 
  w0 = diag(p) + a1%*%solve(A)%*%t(dU[, -k])
  w = bTheta%*%w0 
  return(w)
}

#modified score function
psik_tilde <-function(k, X, Y, dU, V){
  n = nrow(X)
  Xsigma <- t(X)%*%X/n
  XYrho <- t(X)%*%Y/n
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  derL_uk <- Xsigma%*%dU[ ,k] - XYrho%*%V[ ,k]
  derL_vk <- V[ ,k]*z_kk - t(XYrho)%*%dU[ ,k]
  M = Mk_nearly(k, Xsigma, dU, V) #get M
  psi <- derL_uk - M%*%derL_vk
  return(psi)
}


#the de-biased estimate of u_k and d_k
uk_nearly <- function(k, r, p, Theta, dU, V, Xsigma, X, Y){
  uk_tilde <- dU[ ,k]
  W = wk_nearly(k, r, p, Theta, dU, Xsigma)
  psi_tilde = psik_tilde(k, X, Y, dU, V)
  uk_hat <- uk_tilde - W%*%psi_tilde
  dk_hat <- t(uk_tilde)%*%uk_tilde - 2*t(uk_tilde)%*%W%*%psi_tilde
  return(list(u = uk_hat, d = dk_hat))
}


var_nearly <- function(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b){
  M = Mk_nearly(k, Xsigma, dU, V)
  W = wk_nearly(k, r, p, Theta, dU, Xsigma)
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  v_kk = as.numeric(t(V[, k])%*%Sigmae%*%V[, k])
  
  #var of uk
  varufinal = rep(0,p)
  for(i in 1:p){
    e = 1
    a1 = rep(0, i-1)
    a2 = rep(0, p-i)
    a = c(a1,e,a2)
  
  aw = t(a)%*%W
  awm = t(a)%*%W%*%M
  var1 = z_kk*awm%*%Sigmae%*%t(awm)
  var2 = v_kk*aw%*%Xsigma%*%t(aw)
  var3 = -2*aw%*%Xsigma%*%dU[,k]%*%t(V[,k])%*%Sigmae%*%t(awm)
  
  varuk = var1 + var2 + var3
  varufinal[i] = varuk
  }
  #var of given b
  a = b
  aw = t(a)%*%W
  awm = t(a)%*%W%*%M
  var1 = z_kk*awm%*%Sigmae%*%t(awm)
  var2 = v_kk*aw%*%Xsigma%*%t(aw)
  var3 = -2*aw%*%Xsigma%*%dU[,k]%*%t(V[,k])%*%Sigmae%*%t(awm)
  
  varukb = var1 + var2 + var3
  
  #var of dk
  aw = 2*t(dU[,k])%*%W
  awm = 2*t(dU[,k])%*%W%*%M
  var1 = z_kk*awm%*%Sigmae%*%t(awm)
  var2 = v_kk*aw%*%Xsigma%*%t(aw)
  var3 = -2*aw%*%Xsigma%*%dU[,k]%*%t(V[,k])%*%Sigmae%*%t(awm)
  
  vardk = var1 + var2 + var3
  return(list(varu = varufinal, varub = varukb, vard = vardk))
}

avgcov <- function(n, uk, varuk, uk0){
  p0 <- length(uk0)
  length0 <- sum(uk0!=0)
  length1 <- p0 - length0
  index1 <- which(uk0!=0)
  index2 <- which(uk0==0)
  len1 = 0
  len2 = 0
  len = rep(0,p0)
  num1 = 0
  num2 = 0
  num = rep(0,p0)
  for(j in index1){
    #if(uk[j]* uk0[j] < 0){ uk[j] = -uk[j] }
    low_ci <- uk[j] - sqrt(varuk[j])*1.96/sqrt(n)
    up_ci <- uk[j] + sqrt(varuk[j])*1.96/sqrt(n)
    len1 = len1 + abs(up_ci - low_ci)
    len[j] =  abs(up_ci - low_ci)
    if(uk0[j] <= up_ci && uk0[j] >=low_ci) {num[j] = 1; num1 = num1 + 1 }
  }
  for(j in index2){
    low_ci <-  uk[j]-sqrt(varuk[j])*1.96/sqrt(n)
    up_ci <-  uk[j] + sqrt(varuk[j])*1.96/sqrt(n)
    len2 = len2 + abs(up_ci - low_ci)
    len[j] =  abs(up_ci - low_ci)
    if(uk0[j] <= up_ci && uk0[j] >=low_ci) {num[j] = 1; num2 = num2 + 1}
  }
  cov_s0 <- num1/length0 
  cov_s1 <- num2/length1
  len_s0 <- len1/length0
  len_s1 <- len2/length1
  cov_all <- (num1 + num2)/p0
  len_all <- (len1 + len2)/p0
  return(list(cov_s0 = cov_s0, cov_s1 = cov_s1, cov_all = cov_all, len_s0 = len_s0, len_s1 = len_s1, len_all = len_all, num = num, len = len))
}

avgd <- function(n, dk, vardk, dk0){
  dk = as.numeric(dk)
  vardk = as.numeric(vardk)
  num = 0
  len = 0
  low_ci <- dk - sqrt(vardk)*1.96/sqrt(n)
  up_ci <- dk + sqrt(vardk)*1.96/sqrt(n)
  len =  abs(up_ci - low_ci)
  if( dk0^2 <= up_ci &&  dk0^2 >= low_ci) {num =  1}
  
  return(list(num = num, len = len))
}


# transuv <- function(A, A0){#rank 3 transformation
#   index1 <- which(A0[,1]!=0)
#   index2 <- which(A0[,2]!=0)
#   index3 <- which(A0[,3]!=0)
#   ind <- list(index1, index2, index3)
#   for(i in 1:3){
#     for(j in ind[[i]]){
#       if(A[j,i]* A0[j,i] < 0){ A[,i] = -A[,i] }
#     }
#   }
#   return(A)
# }

transuv <- function(A, A0){#angle transformation
  cola <- NCOL(A)
  cola0 <- NCOL(A0)
  if(cola < cola0){
    cat("rank is underestimate as", cola, "\n")
  }
  for(i in 1:cola0){
    index <- which(A0[,i]!=0)
    for(j in index){
      if(A[j,i]* A0[j,i] < 0){ A[,i] = -A[,i] }
    }
  }
  return(A)
}





