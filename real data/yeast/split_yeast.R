
n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
# Number of splits
n_splits <- 100  

# 创建一个 n_loop x 5 的矩阵，存放 5 个数在每次循环中的值
result_mat <- matrix(NA, nrow = n_splits, ncol = 5)
colnames(result_mat) <- c("rrr", "sofar_r", "rrr_v2", "sofar_rv2", "sofar_ori")

# Define the training ratio (here it's 70% for training, 30% for testing)
train_ratio <- 0.8  
n_train <- round(train_ratio * n)
n_test <- n - n_train
# Perform 50 random splits
set.seed(2024)
for (i in 1:n_splits) {
  # 1. Randomly sample indices for the training set
  train_idx <- sample(seq_len(n), size = n_train)
  # 2. The remaining indices form the test set
  test_idx <- setdiff(seq_len(n), train_idx)
  
  # 3. Subset X and Y accordingly
  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  Y_test  <- Y[test_idx, , drop = FALSE]
  
  ##apply sofari and sofar
  set_all <- sofari(X_train, Y_train)
  p_index <- set_all$iU_set
  q_index <- set_all$V_set
  C_sofar <- set_all$C_sofar
  
  X_train_sele <- X_train[, p_index, drop = FALSE]
  
  ## refitting
  #2.1 RRR no screening q 
  fit_rrr <-  rrr.fit(Y_train, X_train_sele, nrank = est_rank)
  coef_rrr <- matrix(0, nrow = p, ncol = q)
  coef_rrr[p_index, ] <- fit_rrr$coef
  
  loss_rrr2 <- (norm(Y_test-X_test%*%coef_rrr,"F"))^2/(n_test*q)
  
  #2.2 sofar no screening q 
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
  
  #1.1 RRR screening q 
  fit_rrr <-  rrr.fit(Y_train_sele,X_train_sele,nrank = est_rank)
  coef_rrr <- matrix(0, nrow = p, ncol = q)
  coef_rrr[p_index, q_index] <- fit_rrr$coef
  
  loss_rrr <- (norm(Y_test-X_test%*%coef_rrr,"F"))^2/(n_test*q)
  
  #1.2 sofar screening q 
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
    loss_rrr <- loss_rrr2
    loss_sofar <- loss_sofar2
  }
  
  
  #3. sofar 
  loss_sofar_est <- (norm(Y_test-X_test%*%C_sofar,"F"))^2/(n_test*q) 
  
  result_mat[i, ] <- c(loss_rrr2, loss_sofar2, loss_rrr, loss_sofar, loss_sofar_est)
  
  avg_result <- colMeans(result_mat[1:i, , drop = FALSE])
  cat(sprintf(
    "Iteration %d:\n  rrr       = %.4f  (avg = %.4f)\n  sofar_r   = %.4f  (avg = %.4f)\n  rrr2      = %.4f  (avg = %.4f)\n  sofar_r2  = %.4f  (avg = %.4f)\n  sofar_ori = %.4f  (avg = %.4f)\n\n",
    i,
    loss_rrr2,       avg_result[1],
    loss_sofar2,   avg_result[2],
    loss_rrr,      avg_result[3],
    loss_sofar,  avg_result[4],
    loss_sofar_est, avg_result[5]
  ))
}
std_values <- apply(result_mat, 2, sd)
std_values
se_values <- std_values / sqrt(nrow(result_mat)  )
se_values
sofari <- function(X, Y){
  est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 200, rank_tol = 1e-2, C = 2.01); est_rank
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
  r = est_rank #r = 3
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
  obj_sigmae <-  sim(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
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

