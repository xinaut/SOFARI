
calculate.Z <- function(x,
                        parallel,
                        ncores,
                        verbose,
                        Z,
                        do.ZnZ = FALSE,
                        debug.verbose = FALSE)
{
  if(is.null(Z)){
    message("Nodewise regressions will be computed as no argument Z was provided.")
    message("You can store Z to avoid the majority of the computation next time around.")
    message("Z only depends on the design matrix x.")
    nodewiselasso.out <- score.nodewiselasso(x = x,
                                             wantTheta = FALSE,
                                             parallel = parallel,
                                             ncores = ncores,
                                             cv.verbose = verbose || debug.verbose,
                                             do.ZnZ = do.ZnZ,
                                             verbose = debug.verbose)
    Z <- nodewiselasso.out$out$Z
    scaleZ <- nodewiselasso.out$out$scaleZ
  }else{
    scaleZ <- rep(1,ncol(Z))
    
    ## Check if normalization is fulfilled
    if(!isTRUE(all.equal(rep(1, ncol(x)), colSums(Z * x) / nrow(x), tolerance = 10^-8))){
      ##no need to print stuff to the user, this is only an internal detail
      rescale.out <- score.rescale(Z = Z, x = x)
      Z <- rescale.out$Z
      scaleZ <- rescale.out$scaleZ
    }
  }
  list(Z = Z,
       scaleZ = scaleZ)
}

score.nodewiselasso <- function(x,
                                wantTheta = FALSE,
                                verbose = FALSE,
                                lambdaseq = "quantile",
                                parallel = FALSE,
                                ncores = 8,
                                oldschool = FALSE,
                                lambdatuningfactor = 1,
                                cv.verbose = FALSE,
                                do.ZnZ = TRUE)
{
  ## Purpose:
  ## This function calculates the score vectors Z OR the matrix of nodewise
  ## regression coefficients Thetahat, depending on the argument wantTheta.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  ## First, a sequence of lambda values over all nodewise regressions is created

  lambdas <-
    switch(lambdaseq,
           "quantile" = nodewise.getlambdasequence(x), ## this is preferred
           "linear"   = nodewise.getlambdasequence.old(x,verbose),
           ## otherwise
           stop("invalid 'lambdaseq': ", lambdaseq))
  
  if(verbose){
    cat("Using the following lambda values:", lambdas, "\n")
  }

  ## 10-fold cv is done over all nodewise regressions to calculate the error
  ## for the different lambda
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x,
                                      parallel = parallel, ncores = ncores,
                                      oldschool = oldschool,
                                      verbose = cv.verbose)
  if(verbose){
    cat(paste("lambda.min is", cvlambdas$lambda.min), "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), "\n")
  }

  if(do.ZnZ)
  {
    bestlambda <- improve.lambda.pick(x = x,
                                      parallel = parallel,
                                      ncores = ncores,
                                      lambdas = lambdas,
                                      bestlambda = cvlambdas$lambda.min,
                                      verbose = verbose)
    if(verbose)
    {
      cat("Doing Z&Z technique for picking lambda\n")
      cat("The new lambda is",bestlambda,"\n")

      cat("In comparison to the cross validation lambda, lambda = c * lambda_cv\n")
      cat("c=",bestlambda/cvlambdas$lambda.min,"\n")## Some info on how much smaller lambda truly is
      }
  }else{
    if(lambdatuningfactor == "lambda.1se"){
      if(verbose)
        cat("lambda.1se used for nodewise tuning\n")
      ## We use lambda.1se for bestlambda now!!!
      bestlambda <- cvlambdas$lambda.1se
    }else{
      if(verbose)
        cat("lambdatuningfactor used is", lambdatuningfactor, "\n")
      
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  }
  
  if(verbose){
    cat("Picked the best lambda:", bestlambda, "\n")
    ##print("with the error ")
    ##print(min(err))
  }

  ## Having picked the 'best lambda', we now generate the final Z or Thetahat
  if(wantTheta){
    out <- score.getThetaforlambda(x = x,
                                   lambda = bestlambda,
                                   parallel = parallel,
                                   ncores = ncores,
                                   oldschool = TRUE,
                                   verbose = verbose)
  }else{
    Z <- score.getZforlambda(x = x, lambda = bestlambda, parallel = parallel,
                             ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out,
                     bestlambda = bestlambda)
  return(return.out)
}


improve.lambda.pick <- function(x,
                                parallel,
                                ncores,
                                lambdas,
                                bestlambda,
                                verbose)
{
  ## let's improve the lambda choice by using the Z&Z procedure
  ## (But picking the one new lambda value not one for each j)
  
  lambdas <- sort(lambdas, decreasing = TRUE) ## Such that we have a natural ordering of decreasing lambdas --> ordering by increasing variance

  ## Compute the mean theoretical variance for every lambda value
  M <- calcM(x = x,
             lambdas = lambdas,
             parallel = parallel,
             ncores = ncores)
  
  Mcv <- M[which(lambdas == bestlambda)] ## Get the current M alue
  
  if(length(which(M < 1.25*Mcv))>0){
    ## There exists a 25% inflated M in this range of lambda choices
    lambdapick <- min(lambdas[which(M < 1.25*Mcv)])
  }else{
    if(verbose)
      cat("no better lambdapick found\n")
     lambdapick <- bestlambda
   }

  if(max(which(M < 1.25*Mcv)) < length(lambdas))
  {
    ## There is an interval of potential lambda values that give close to the exact 25% inflation
    ## Let's refine the interval to try to find a more accurate solution
    if(verbose)
      cat("doing a second step of discretisation of the lambda space to improve the lambda pick\n")

    lambda.number <- max(which(M < 1.25*Mcv))##assumes the lambdas are sorted from big to small
    newlambdas <- seq(lambdas[lambda.number],lambdas[lambda.number+1], (lambdas[lambda.number+1]-lambdas[lambda.number])/100)
    newlambdas <- sort(newlambdas,decreasing=TRUE)##convention
    M2 <- calcM(x=x,
                lambdas=newlambdas,
                parallel=parallel,
                ncores=ncores)
    
    if(length(which(M2 < 1.25*Mcv))>0){
      evenbetterlambdapick <- min(newlambdas[which(M2 < 1.25*Mcv)])
    }else{
      if(verbose)
        cat("no -even- better lambdapick found\n")
      evenbetterlambdapick <- lambdapick
    }
    
    if(is.infinite(evenbetterlambdapick))
    {##it is possible that when redoing the fits the M2 value is all of a sudden higher for all lambdas in our newlambdas range

      ##just leave lambdapick as it is
      if(verbose)
      {
        cat("hmmmm the better lambda pick after the second step of discretisation is Inf\n")
        cat("M2 is\n")
        cat(M2,"\n")
        cat("Mcv is\n")
        cat(Mcv,"\n")
        cat("and which(M2 < 1.25* Mcv) is\n")
        cat(which(M2 < 1.25*Mcv),"\n")
      }
    }else{
      lambdapick <- evenbetterlambdapick
    }
  }
  
  return(lambdapick)
}

calcM <- function(x,
                  lambdas,
                  parallel,
                  ncores)
{
  ## Compute the theoretical variances for each j
  ## for each choice for lambda
  ## Return the means over the j for each distinct lambda value
  
  if(parallel)
    {
      M <- mcmapply(FUN=calcMforcolumn,
                    x=list(x=x),
                    j=1:ncol(x),
                    lambdas=list(lambdas=lambdas),
                    mc.cores=ncores)
    }else{
      M <- mapply(FUN=calcMforcolumn,
                  j=1:ncol(x),
                  x=list(x=x),
                  lambdas=list(lambdas=lambdas))
    }
  ##every row of M contains for a particular lambda choice the value of M for a particular nodewise regression
  ##M should be of dimension nxp (100x500)
  M <- apply(M,1,mean)
  return(M)
}

calcMforcolumn <- function(x,
                           j,
                           lambdas)
{
  ## do a nodewise regression and calc the l2 norm ^2 of the residuals
  ## and the inner product of Zj and Xj
  glmnetfit <- glmnet(x[,-j],x[,j],
                      lambda=lambdas)
  predictions <- predict(glmnetfit,x[,-j],s=lambdas)
  Zj <- x[,j] - predictions##the columns are the residuals for the different lambda

  Znorms <- apply(Zj^2,2,sum)
  Zxjnorms <- as.vector(crossprod(Zj,x[,j])^2)
  return(Znorms/Zxjnorms)
}


score.getThetaforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                    oldschool = FALSE, verbose = FALSE,
                                    oldtausq = TRUE)
{
  ## Purpose:
  ## This function is for calculating Thetahat once the desired tuning
  ## parameter is known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  message("Calculating Thetahat by doing nodewise regressions and dropping the unpenalized intercept")
  n <- nrow(x)
  p <- ncol(x)
  C <- diag(rep(1,p))
  T2 <- numeric(p)
  
  if(oldschool){
    message("doing getThetaforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i], x[,i])
      coeffs <- as.vector(predict(glmnetfit,x[,-i], type = "coefficients",
                                  s = lambda))[-1]
      ## we just leave out the intercept
      
      C[-i,i] <- -as.vector(coeffs)
      if(oldtausq){
        ## possibly quite incorrect,it ignores the intercept, but intercept is
        ## small anyways. Initial simulations show little difference.
        T2[i] <- as.numeric(crossprod(x[,i]) / n - x[,i] %*% (x[,-i] %*%
                                                              coeffs) / n)
        }else{   
          ##print("now doing the new way of calculating tau^2")
          T2[i] <- as.numeric((x[,i] %*%
                               (x[,i] - predict(glmnetfit,x[,-i],s =lambda)))/n)
        }
    }
  }else{
    stop("not implemented yet!")
  }
  thetahat <- C %*% solve(diag(T2))
  if(verbose){
    cat("1/tau_j^2:", solve(diag(T2)), "\n")
  }
  ##this is thetahat ^ T!!
  thetahat <- t(thetahat)
  
  if(all(thetahat[lower.tri(thetahat)] == 0,
         thetahat[upper.tri(thetahat)] == 0) && verbose)
    cat("Thetahat is a diagonal matrix!\n")
  
  return(thetahat)
}

score.getZforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                oldschool = FALSE)
{
  ## Purpose:
  ## This function is for calculating Z once the desired tuning parameter is
  ## known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version)
    
  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n*p),n)
  
  if(oldschool){
    message("doing getZforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i],x[,i])
      prediction <- predict(glmnetfit,x[,-i],s=lambda)
      Z[,i] <- x[,i] - prediction
    }
  }else{
    ## REPLACING THE FOR LOOP
    if(parallel){
      Z <- mcmapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                    lambda = lambda, mc.cores = ncores)
      
    }else{
      Z <- mapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                  lambda = lambda)
    }
  }
  ## rescale Z such that t(Zj) Xj/n = 1 \-/ j
  Z <- score.rescale(Z,x)
  return(Z)
}

score.getZforlambda.unitfunction <- function(i, x, lambda)
{
  ## Purpose:
  ## Calculate the residuals of a nodewise regression of one column vs the
  ## other columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  glmnetfit  <- glmnet(x[,-i],x[,i])
  prediction <- predict(glmnetfit,x[,-i],s=lambda)
  return(x[,i] - prediction)
}

score.rescale <- function(Z, x)
{
  ## Purpose:
  ## Rescale the Z appropriately such that such that t(Zj) Xj/n = 1 for all j
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  scaleZ <- diag(crossprod(Z,x)) / nrow(x)
  Z      <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z=Z,scaleZ=scaleZ))
}

nodewise.getlambdasequence <- function(x)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## Equidistant quantiles of the complete set of lambda values are returned.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  
  nlambda <- 100 ## use the same value as the glmnet default
  p <- ncol(x)
  
  lambdas <- c()
  for(c in 1:p){
    lambdas <- c(lambdas,glmnet(x[,-c],x[,c])$lambda)
  }
  
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

cv.nodewise.err.unitfunction <- function(c, K, dataselects, x, lambdas,
                                         verbose, p) {
  ## Purpose:
  ## this method returns the K-fold cv error made by the nodewise regression
  ## of the single column c of x vs the other columns for all values of the
  ## tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  
  if(verbose){ ##print some information out about the progress
      ##report every 25%
    interesting.points <- round(c(1/4,2/4,3/4,4/4)*p)
    names(interesting.points) <- c("25%","50%","75%","100%")
    if(c %in% interesting.points){
      message("The expensive computation is now ",
              names(interesting.points)[c == interesting.points],
              " done")
    }
  }
  
  ## return  'totalerr' 
  cv.nodewise.totalerr(c = c,
                       K = K,
                       dataselects = dataselects,
                       x = x,
                       lambdas = lambdas)
}

## gets the standard error for one particular lambda
cv.nodewise.stderr <- function(K, x, dataselects, lambda, parallel, ncores)
{
  ## Purpose:
  ## this method returns the standard error 
  ## of the average K-fold cv error made by the nodewise regression
  ## of each column vs the other columns for a single lambda value.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  p <- ncol(x)
  if(parallel){
    totalerr <- mcmapply(cv.nodewise.totalerr,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = lambda,
                         mc.cores = ncores)
  }else{
    totalerr <- mapply(cv.nodewise.totalerr,
                       c = 1:p,
                       K = K,
                       dataselects = list(dataselects = dataselects),
                       x = list(x = x),
                       lambdas = lambda)
  }
  ## get the mean over the variables
  totalerr.varmean <- rowMeans(totalerr)
  
  ## get the stderror over the K;  return stderr.forlambda
  sd(totalerr.varmean) / sqrt(K)
}

cv.nodewise.totalerr <- function(c, K, dataselects, x, lambdas)
{
  ## Purpose:
  ## this method returns the error made for each fold of a K-fold cv
  ## of the nodewise regression of the single column c of x vs the other
  ## columns for all values of the tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  
  for(i in 1:K){ ## loop over the test sets
    whichj <- dataselects == i ##the test part of the data
    
    glmnetfit <- glmnet(x = x[!whichj,-c, drop = FALSE],
                        y = x[!whichj, c, drop = FALSE],
                        lambda = lambdas)
    predictions  <- predict(glmnetfit, newx = x[whichj, -c, drop = FALSE],
                            s = lambdas)
    totalerr[, i] <- apply((x[whichj, c] - predictions)^2, 2, mean)
  }

  totalerr
}


cv.nodewise.bestlambda <- function(lambdas, x, K = 10, parallel = FALSE,
                                   ncores = 8, oldschool = FALSE,
                                   verbose = FALSE)
{
  ## Purpose:
  ## this function finds the optimal tuning parameter value for minimizing
  ## the K-fold cv error of the nodewise regressions.
  ## A second value of the tuning parameter, always bigger or equal to the
  ## former, is returned which is calculated by allowing the cv error to
  ## increase by the amount of
  ## 1 standard error (a similar concept as to what is done in cv.glmnet).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)
  
  ## Based on code from cv.glmnet for sampling the data
  dataselects <- sample(rep(1:K, length = n))
  
  if(oldschool){
    message("doing cv.nodewise.error oldschool")
    totalerr <- numeric(l)
    for(c in 1:p){ ## loop over the nodewise regressions
      for(i in 1:K){ ## loop over the test sets
        whichj <- dataselects == i ## the test part of the data
        
        glmnetfit <- glmnet(x[!whichj,-c,drop=FALSE], x[!whichj,c,drop=FALSE],
                            lambda=lambdas)
        predictions <- predict(glmnetfit,x[whichj, -c, drop = FALSE],
                               s = lambdas)
        totalerr <- totalerr + apply((x[whichj,c]-predictions)^2, 2, mean)
      }
    }
    totalerr <- totalerr / (K * p)
    stop("lambda.1se not implemented for oldschool cv.nodewise.bestlamba")
  }else{
    ## REPLACING THE FOR LOOP

    ##totalerr <- matrix(nrow = l, ncol = p)
    
    if(parallel){
      totalerr <- mcmapply(cv.nodewise.err.unitfunction,
                           c = 1:p,
                           K = K,
                           dataselects = list(dataselects = dataselects),
                           x = list(x = x),
                           lambdas = list(lambdas = lambdas),
                           mc.cores = ncores,
                           SIMPLIFY = FALSE,
                           verbose = verbose,
                           p=p)
    }else{
      totalerr <- mapply(cv.nodewise.err.unitfunction,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = list(lambdas = lambdas),
                         SIMPLIFY = FALSE,
                         verbose = verbose,
                         p = p)
    }
    ## Convert into suitable array (lambda, cv-fold, predictor)
    err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
    err.mean   <- apply(err.array, 1, mean) ## 1 mean for each lambda

    ## calculate mean for every lambda x fold combination (= average over p)
    ## for every lambda then get the standard errors (over folds)
    err.se     <- apply(apply(err.array, c(1, 2), mean), 1, sd) / sqrt(K)
    ##totalerr <- apply(totalerr, 1, mean)
  }

  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]

  stderr.lambda.min <- err.se[pos.min]
##-   stderr.lambda.min <- cv.nodewise.stderr(K = K,
##-                                           x = x,
##-                                           dataselects = dataselects,
##-                                           lambda = lambda.min,
##-                                           parallel = parallel,
##-                                           ncores = ncores)
  list(lambda.min = lambda.min,
       lambda.1se = max(lambdas[err.mean < (min(err.mean) + stderr.lambda.min)]))
}

##DEPRECATED, not really the best choice
nodewise.getlambdasequence.old <- function(x,verbose=FALSE)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## It returns a __linear__ interpolation of lambda values between the max and
  ## min lambda value found.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  nlambda <- 100#take the same value as glmnet does automatically
  p <- ncol(x)
  maxlambda <- 0
  minlambda <- 100
  
  for(c in 1:p){
    lambdas <- glmnet(x[,-c],x[,c])$lambda
    
      ##DEBUG
    if(verbose || is.nan(max(lambdas))){
      cat(paste("c:", c, "\n"))
      cat("lambdas:", lambdas, "\n")
      cat("max(lambdas) max(lambdas,na.rm=TRUE) maxlambda: ",
          max(lambdas), max(lambdas,na.rm=TRUE), maxlambda, "\n")
    }
    if(max(lambdas,na.rm=TRUE) > maxlambda){
      maxlambda <- max(lambdas,na.rm=TRUE)
    }
    if(min(lambdas,na.rm=TRUE) < minlambda){
      minlambda <- min(lambdas, na.rm = TRUE)
    }
  }
  
  lambdas <- seq(minlambda, maxlambda, by = (maxlambda-minlambda)/nlambda)
  ## return
  sort(lambdas, decreasing=TRUE)
}


###############
#
# Authors: Hamid Javadi, Adel Javanmard, Andrea Montanari, Sven Schmit
# Date: April 24th, 2014
#
# Confidence intervals for high-dimensional regression using the method of
# A. Javanmard, A. Montanari
# "Confidence intervals and hypothesis testing for high-dimensional regression"
# 2013, arXiv:1306.3171
#
library(Matrix);
library(glmnet);
library(expm);
library(flare);

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

NoiseSd <- function( yh, A, n ){
  ynorm <- sqrt(n)*(yh/sqrt(diag(A)));
  sd.hat0 <- mad(ynorm);
  
  zeros <- (abs(ynorm)<3*sd.hat0);
  y2norm <- sum(yh[zeros]^2);
  Atrace <- sum(diag(A)[zeros]);
  sd.hat1 <- sqrt(n*y2norm/Atrace);
  
  ratio <- sd.hat0/sd.hat1;
  if (max(ratio,1/ratio)>2)
    print("Warning: Noise estimate problematic");
  
  s0 <- sum(zeros==FALSE);
  return (list( "sd" = sd.hat1, "nz" = s0));
}


Lasso <- function( X, y, lambda = NULL, intercept = TRUE){
  #
  # Compute the Lasso estimator:
  # - If lambda is given, use glmnet and standard Lasso
  # - If lambda is not given, use square root Lasso
  #
  p <- ncol(X);
  n <- nrow(X);
  
  if  (is.null(lambda)){
    lambda <- sqrt(qnorm(1-(0.1/p))/n);
    outLas <- slim(X,y,lambda=c(lambda),method="lq",q=2,verbose=FALSE);
    # Objective : sqrt(RSS/n) +lambda *penalty
    if (intercept==TRUE) {
      return (c(as.vector(outLas$intercept),as.vector(outLas$beta)))
    }  else {
      return (as.vector(outLas$beta));
    }
  } else {
    outLas <- glmnet(X, y, family = c("gaussian"), alpha =1, intercept = intercept );
    # Objective :1/2 RSS/n +lambda *penalty
    if (intercept==TRUE){
      return (as.vector(coef(Las,s=lambda)));
    } else {
      return (as.vector(coef(Las,s=lambda))[2:(p+1)]);
    }
  }
}

SSLasso <- function (X, y, alpha=0.05, lambda = NULL, mu = NULL, intercept = TRUE, 
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
  #
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X     :  design matrix
  #   y     :  response
  #   alpha :  significance level
  #   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
  #   mu    :  Linfty constraint on M (if null, searches)
  #   resol :  step parameter for the function that computes M
  #   maxiter: iteration parameter for computing M
  #   threshold : tolerance criterion for computing M
  #   verbose : verbose?
  #
  # Returns:
  #   noise.sd: Estimate of the noise standard deviation
  #   norm0   : Estimate of the number of 'significant' coefficients
  #   coef    : Lasso estimated coefficients
  #   unb.coef: Unbiased coefficient estimates
  #   low.lim : Lower limits of confidence intervals
  #   up.lim  : upper limit of confidence intervals
  #   pvals   : p-values for the coefficients						 
  #
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);
  
  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);
  
  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  
  if ((n>=2*p)&&(tmp>=1e-4)){
    M <- solve(sigma.hat)
  }else{
    M <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  }
  
  unbiased.Lasso <- as.numeric(htheta + (M%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  
  A <- M %*% sigma.hat %*% t(M);
  noise <- NoiseSd( unbiased.Lasso, A, n );
  s.hat <- noise$sd;
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));
  
  if  (is.null(lambda)){
    lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  }
  
  addlength <- rep(0,pp);
  MM <- M%*%sigma.hat - diag(pp);
  for (i in 1:pp){
    effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
    effectivemuvec <- effectivemuvec[0:(noise$nz-1)];
    addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  }  
  
  htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm;
  interval.sizes <- interval.sizes*col.norm;
  addlength <- addlength*col.norm;
  
  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
  p.vals <- 2*(1-pnorm(sqrt(n)*abs(unbiased.Lasso)/(s.hat*col.norm[(pp-p+1):pp]*sqrt(diag(A[(pp-p+1):pp,(pp-p+1):pp])))))
  
  returnList <- list("noise.sd" = s.hat,
                     "norm0" = noise$nz,
                     "coef" = htheta,
                     "unb.coef" = unbiased.Lasso,
                     "low.lim" = unbiased.Lasso - interval.sizes - addlength,
                     "up.lim" = unbiased.Lasso + interval.sizes + addlength,
                     "pvals" = p.vals
  )
  return(returnList)
}


