
sim <- function(p = p, n = n, delta = 2, thresholding = 'hard', data, cov_true) {
  # compute theta
  get.theta <- function(i = i, j = j) {  
    X_bar_i <- mean(data[, i]);
    X_bar_j <- mean(data[, j]);
    
    sigma <- cov(data)[i,j];
    
    out <- sum(
      unlist(
        lapply(
          X = 1:n, # X = 1:n instead of 1:p
          FUN = function(z) { 
            ((data[z, i] - X_bar_i)*(data[z, j] - X_bar_j) - sigma)^2;
          }
        )
      )
    );
    
    
    return(out/n); # return out/n instead of out/p
    
  }
  
  # compute the estimated covariance matrix
  get.est.cov <- function(delta = delta, thresholding = thresholding) {
    
    datagrid <- expand.grid(i = 1:p, j = 1:p);
    res <- apply(datagrid, 
                 MARGIN = 1, 
                 FUN = function(z){
                   get.theta(i = z["i"], j = z["j"])
                 }
    );
    
    theta <- matrix(res, nrow = p, ncol = p);
    lambda <- delta*(sqrt(theta*log(p)/n));
    
    
    # hard thresholding
    if('hard' == thresholding) {
      cov_est <- sample_cov*(abs(sample_cov) > lambda)*1;
    }
    
    # adaptive lasso thresholding
    if('al' == thresholding) {
      d <- 1- (abs(lambda/sample_cov))^4;
      d <- d*(d > 0);
      cov_est <- sample_cov*d;
    }
    
    
    return(cov_est);
    
  }
  
  #if('m1' == model) {cov_true <- m1.data(p = p, n = n); }else{ cov_true <- m2.data(p = p, n = n); }
  
  # data <- rmvnorm(n = n, mean = rep(0, q), sigma = cov_true);
  sample_cov <- cov(data);
  
  cov_hard <- get.est.cov(delta = delta, thresholding = thresholding);
  
  return( list(cov_est = cov_hard, cov_true = cov_true, 
          loss2 = norm(     x = cov_hard - cov_true,   type = '2' ),
     lossF = norm( x = cov_hard - cov_true,    type = 'F')
    ) )
  
}


sim.find.delta <- function(delta = delta, 
                           fold.id = fold.id, 
                           p = p, 
                           n = n, 
                           thresholding = 'hard', 
                           model = 'm1',
                           measure = 'f',
                           data, cov_true
) {
  
  # split data into 5 folds
  #set.seed(12345);
  sample <- sample(1:n, size = n, replace = FALSE);
  fold <- rep(1:5, length = n);
  
  train_id <- list();
  test_id <- list();
  
  for(i in 1:5) {
    train_id[[i]] <- sample[fold != i]; 
    test_id[[i]] <- sample[fold == i]; 
  }
  
  #if('m1' == model) {cov_true <- m1.data(p = p, n = n); }else{ cov_true <- m2.data(p = p, n = n); }
  #data <- rmvnorm(n = n, mean = rep(0, p), sigma = cov_true);
  sample_cov <- cov(data);
  
  train <- train_id[[fold.id]];
  data_train <- data[train, ];
  
  # compute theta
  get.theta <- function(i = i, j = j) {
    
    X_bar_i <- mean(data_train[, i]);
    X_bar_j <- mean(data_train[, j]);
    sigma <- cov(data_train)[i,j];
    
    out <- sum(
      unlist(
        lapply(
          X = 1:length(train),
          FUN = function(z) { 
            jj <- train[z]; 
            ((data[jj, i] - X_bar_i)*(data[jj, j] - X_bar_j) - sigma)^2;
          }
        )
      )
    );
    
    
    return(out/length(train));
    
  }
  
  # find the best delta
  datagrid <- expand.grid(i = 1:p, j = 1:p);
  res <- mcmapply(get.theta, 
                  datagrid$i, 
                  datagrid$j,
                  mc.cores = 1 
  );
  
  theta <- matrix(res, nrow = p, ncol = p);
  lambda <- delta*sqrt(theta*log(p)/length(train));
  
  # sample variance from the training samples
  sample_sigma_train <- cov(data_train);
  
  # sample variance from the testing samples
  
  sample_sigma_test <- cov(data[test_id[[fold.id]], ]);
  
  # hard thresholding
  if('hard' == thresholding) {
    if('f' == measure) {
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = 'F' #Frobenius norm
      );
    }else{
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = '2'
      );
    }
    
  }
  
  # adaptive lasso thresholding
  if('al' == thresholding) {
    d <- 1- (abs(lambda/sample_sigma_train))^4;
    d <- d*(d > 0);
    
    if('f' == measure) {
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = 'F'
      );
    }else{
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = '2'
      );
    }
    
  }
  
  return(norm_loss)
  
}

  
sim.cv <- function(n = n, p = p, data, cov_true){
  measure <- '2';
  # N is the fixed integer, refer paper for details
  N <- 10;
  
  deltagrid <- expand.grid(fold.id = 1:5, 
                           delta = seq(0, 4, 1/N), 
                           thresholding = c("hard","al"),
                           p = p,
                           n = n,
                           model = 'm1',
                           measure = measure,
                           data = data,
                           cov_true = cov_true
  );
  find_delta <- deltagrid;
  
  
  find_delta$loss <- mcmapply(
    sim.find.delta, 
    delta = deltagrid$delta, 
    fold.id = deltagrid$fold.id, 
    thresholding = deltagrid$thresholding,
    n = deltagrid$n,
    p = deltagrid$p,
    model = deltagrid$model,
    measure = deltagrid$measure,
    data = deltagrid$data,
    cov_true = deltagrid$cov_true,
    mc.cores = 1
  );
  
  
  # al thresholding 
  find_delta_al <- aggregate( 
    loss ~ delta, 
    data = find_delta["al" == find_delta$thresholding, ],
    FUN = mean 
  );
  
  delta_opt_al <- find_delta_al[which.min(find_delta_al$loss), 'delta'];
  
  # hard thresholding
  find_delta_hard <- aggregate( 
    loss ~ delta, 
    data = find_delta["hard" == find_delta$thresholding, ],
    FUN = mean 
  );
  
  delta_opt_hard <- find_delta_hard[which.min(find_delta_hard$loss), 'delta'];
  
  # # print optimal delta
  # message('the best delta in hard thresholding is ', delta_opt_hard, ' ','under norm', ' ',measure);
  # message('the best delta in al thresholding is ', delta_opt_al,' ','under norm', ' ', measure);
  return(list(delta_opt_al = delta_opt_al, delta_opt_hard = delta_opt_hard ))
}
