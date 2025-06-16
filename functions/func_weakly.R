##functions for the weakly general form

#M
Mk_weakly <- function(k, r, Xsigma, dU, V){
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  if(k < r){
  C_k <- dU[, -(1:k)  ]%*%t(V[, -(1:k)]) 
  M <- -Xsigma%*%C_k/z_kk 
  }else{
    #C_k <- dU[, k  ]%*%t(V[, k])
    M <- matrix(0, nrow = dim(Xsigma)[1], ncol = dim(V)[1])
  }
  return(M)
}

#W
wk_weakly <- function(k, r, p, Theta, dU, Xsigma){
  if(k == 1){
    w = wk_nearly(k, r, p, Theta, dU, Xsigma)
  }else if(k > 1 && k < r){
    bTheta = Theta
    z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
    dU2 = dU[, (k+1):r]
    A = diag(r-k) - t(dU2)%*%Xsigma%*%dU2/z_kk
    a1 = (Xsigma%*%dU2)/z_kk 
    w0 = diag(p) + a1%*%solve(A)%*%t(dU2)
    w = bTheta%*%w0 
  }else if(k == r){
    w = Theta
    # w1r = Theta
    # z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
    # a0r = Xsigma%*%dU[, k]%*%t(dU[, k])
    # w2r = diag(p) - a0r/(2*z_kk)
    # w = w1r %*% w2r
  }
  return(w)
}

#DML score function 
psik_weakly <-function(k, r, X, Y, dU, V, u_wh){
  if(k == 1){
    psi = psik_tilde(k, X, Y, dU, V)
  }else  if(k > 1){
    n = nrow(X)
    Xsigma <- t(X)%*%X/n
    XYrho <- t(X)%*%Y/n
    z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
    C_wh = matrix(0, nrow = dim(X)[2], ncol = dim(Y)[2])
    for(i in 1:(k-1)){
      C_wh <- C_wh + as.vector(u_wh[[i]])%*%t(V[ ,i])
    }
    #if(k == 2){ C_wh = u1_wh%*%t(V[ ,1]) }else if(k == 3){C_wh = u1_wh%*%t(V[ ,1]) + u2_wh%*%t(V[ ,2]) }
    temp1 = Xsigma%*%C_wh%*%V[ ,k]
    temp2 = t(C_wh)%*%Xsigma%*%dU[ ,k]
    derL_uk <- Xsigma%*%dU[ ,k] - XYrho%*%V[ ,k] + temp1
    derL_vk <- V[ ,k]*z_kk - t(XYrho)%*%dU[ ,k] + temp2
    M = Mk_weakly(k, r, Xsigma, dU, V) #get M
    psi <- derL_uk - M%*%derL_vk
    
  }
  return(psi)
}


#the de-biased estimate of u_k and d_k
uk_weakly <- function(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh){
  uk_tilde <- dU[ ,k]
  W = wk_weakly(k, r, p, Theta, dU, Xsigma)
  psi_tilde = psik_weakly(k, r, X, Y, dU, V, u_wh)
  uk_hat <- uk_tilde - W%*%psi_tilde
  dk_hat <- t(uk_tilde)%*%uk_tilde - 2*t(uk_tilde)%*%W%*%psi_tilde
  return(list(u = uk_hat, d = dk_hat))
}




var_weakly <- function(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae){
  M_k = Mk_weakly(k, r, Xsigma, dU, V)
  W_k = wk_weakly(k, r, p, Theta, dU, Xsigma)
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  v_kk = as.numeric(t(V[, k])%*%Sigmae%*%V[, k])
  
  #var of uk
  varufinal = rep(0,p)
  for(ii in 1:p){
    e = 1
    a1 = rep(0, ii-1)
    a2 = rep(0, p-ii)
    a = c(a1,e,a2)
    
    aw = t(a)%*%W_k
    awm = t(a)%*%W_k%*%M_k
    
    
    var1 = z_kk*awm%*%Sigmae%*%t(awm)
    var2 = v_kk*aw%*%Xsigma%*%t(aw)
    var3 = -2*aw%*%Xsigma%*%dU[,k]%*%t(V[,k])%*%Sigmae%*%t(awm)
    cov_final1  = var1 + var2 + var3
    var1; var2;  var3
    cov_final1
    
    cov_final = cov_final1
    
    
    varufinal[ii] = cov_final
  }
  
  
  ###var of dk
  M_k = Mk_weakly(k, r, Xsigma, dU, V)
  W_k = wk_weakly(k, r, p, Theta, dU, Xsigma)
  z_kk = as.numeric(t(dU[, k])%*%Xsigma%*%dU[, k])
  v_kk = as.numeric(t(V[, k])%*%Sigmae%*%V[, k])
  a = 2*dU[,k]
  aw = t(a)%*%W_k
  awm = t(a)%*%W_k%*%M_k
  
  var1 = z_kk*awm%*%Sigmae%*%t(awm)
  var2 = v_kk*aw%*%Xsigma%*%t(aw)
  var3 = -2*aw%*%Xsigma%*%dU[,k]%*%t(V[,k])%*%Sigmae%*%t(awm)
  cov_final1  = var1 + var2 + var3
  
  cov_final = cov_final1
  
  vardk = cov_final
  
  return(list(varu = varufinal, vard = vardk))
}

