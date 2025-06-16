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

#load the preprocessed data
load("real data/yeast/yeast_preprocess_data.RData")
n <- nrow(Y)
q <- ncol(Y)
p <- ncol(X)

n_splits <- 100
result_mat <- matrix(NA, nrow = n_splits, ncol = 2)
colnames(result_mat) <- c("sofar_rv2", "sofar_ori")


train_ratio <- 0.8  
n_train <- round(train_ratio * n)
n_test <- n - n_train



for (i in 1:n_splits) {
  train_idx <- sample(seq_len(n), size = n_train)
  test_idx <- setdiff(seq_len(n), train_idx)
  
  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  Y_test  <- Y[test_idx, , drop = FALSE]
  
  set_all <- sofari(X_train, Y_train)
  p_index <- set_all$iU_set
  q_index <- set_all$V_set
  C_sofar <- set_all$C_sofar
  
  X_train_sele <- X_train[, p_index, drop = FALSE]
  
  fit1 <- sofar(Y_train, X_train_sele, ic.type = "BIC", nrank = r,
                control = list(methodA = "adlasso", methodB = "adlasso",
                               nlam = 200, lam.max.factor = 1500, lam.min.factor = 1e-12,
                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
                ), screening = FALSE);  summary(fit1)
  
  
  Ur = fit1$U; Dr = fit1$D; V = fit1$V;
  C2 = Ur%*%diag(Dr)%*%t(V)
  coef_sofar2 <- matrix(0, nrow = p, ncol = q)
  coef_sofar2[p_index, ] <- C2
  
  loss_sofar2 <- (norm(Y_test-X_test%*%coef_sofar2,"F"))^2/(n_test*q) 
  
  if(length(q_index)!= q){
    
    Y_train_sele <- Y_train[, q_index, drop = FALSE]
    est_rank <- STRS(Y_train_sele, X_train_sele, type = "STRS-MC", rep_MC = 200, rank_tol = 1e-2, C = 2.01)
    
   
    fit1 <- sofar(Y_train_sele, X_train_sele, ic.type = "BIC", nrank = r,
                  control = list(methodA = "adlasso", methodB = "adlasso",
                                 nlam = 200, lam.max.factor = 1500, lam.min.factor = 1e-12,
                                 penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
                  ), screening = FALSE);  summary(fit1)
    
    
    Ur = fit1$U; Dr = fit1$D; V = fit1$V;
    dU = Ur%*%diag(Dr)
    r = fit1$rank
    
    C0 = Ur%*%diag(Dr)%*%t(V)
    coef_sofar <- matrix(0, nrow = p, ncol = q)
    coef_sofar[p_index, q_index] <- C0
    
    loss_sofar <- (norm(Y_test-X_test%*%coef_sofar,"F"))^2/(n_test*q)
    
  }else{
    loss_sofar <- loss_sofar2
  }
  
  loss_sofar_est <- (norm(Y_test - X_test %*% C_sofar, "F"))^2 / (n_test * q)
  
  result_mat[i, ] <- c(loss_sofar, loss_sofar_est)
  
  avg_result <- colMeans(result_mat[1:i, , drop = FALSE], na.rm = TRUE)
  cat(sprintf(
    "Iteration %d:\n  sofar_r2  = %.4f  (avg = %.4f)\n  sofar_ori = %.4f  (avg = %.4f)\n\n",
    i,
    loss_sofar,    avg_result[1],
    loss_sofar_est, avg_result[2]
  ))
}

se_values <- apply(result_mat, 2, sd, na.rm = TRUE)/ sqrt(nrow(result_mat))
se_values

sofari <- function(X, Y){
  est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 200, rank_tol = 1e-2, C = 2.01); est_rank
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
  r = est_rank 
  fit1 <- sofar(Y, X, ic.type = "BIC", nrank = r,
                control = list(methodA = "adlasso", methodB = "adlasso",
                               nlam = 200, lam.max.factor = 1500, lam.min.factor = 1e-12,
                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
                ), screening = FALSE);
  
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
  #Theta <- InverseLinfty(Xsigma, n,  resol=1.3, maxiter=100, threshold=1e-6)
  # asd = cv.glasso_clean(X,0)
  # Theta = asd$wi
  # asd = CVglasso(X,  cores = 5)
  # Theta= asd$Omega 
  
  sets_list <- list()
  qq = 0.05
  for(k in 1:r){
    
    if(k == 1){
      u_wh = list()
      ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
      #ukdk1_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
      uk1 <- ukdk1_weak$u
      dk1 <- ukdk1_weak$d
      
      varobj1 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae,  S_u)
      varuk1 <- varobj1$varu
      vardk1 <- varobj1$vard
      
      res1 <- bhq_adj_p(uk1, varuk1, qq , n, p)
      set1 <- res1$St
      sets_list[[1]] <- set1
    }
    else{
      
      u_wh <- list()
      for(kk in 1:k){
        u_wh_temp <- as.vector(dU[,kk]) 
        u_wh[[kk]] <- u_wh_temp
      }
      #u_wh <- lapply(seq_len(k), function(kk) as.vector(dU[, kk]))
      
      ukdk2_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)  # u_wh is no used parameter
      uk2 <- ukdk2_weak$u
      
      
      varobj2 <- var_weakly_nosu(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae) # S_u is no used parameter
      varuk2 <- varobj2$varu
      
      
      res2 <- bhq_adj_p(uk2, varuk2, qq , n, p)
      set2 <- res2$St
      sets_list[[k]] <- set2
    }
    
  }#end for(k in 1:r)
  
  sofari_selected <- Reduce(union, sets_list)
  sofari_selected <- sort(sofari_selected)
  
  
  nonzero_list <- apply(dU, 2, function(col) {
    which(col != 0)  
  })
  combined_nonzero_values <- unlist(nonzero_list)
  sofar_selected <- unique(combined_nonzero_values)
  
  nonzero_list_v <- apply(V, 2, function(col) {
    which(col != 0)  
  })
  combined_nonzero_values2 <- unlist(nonzero_list_v)
  v_selected <- unique(combined_nonzero_values2)
  v_selected <- sort(v_selected)
  
  return(list(iU_set = sofari_selected, C_sofar = C, V_set = v_selected) )
}
