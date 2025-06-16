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


####### do prediction
X0 <- read.csv("real data/economic/final_XX.csv")
Y0 <- read.csv("real data/economic/final_YY.csv")
X0 <- as.matrix(X0); Y0 <- as.matrix(Y0)
q <- dim(Y0)[2]
#trian the data
sam <- 474
train_sam <- 1:sam
Y<- Y0[train_sam, ]; X <- X0[train_sam, ]

est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
r = est_rank
fit10 <- sofar(Y, X, ic.type = "GIC", nrank = r, 
               control = list(methodA = "adlasso", methodB = "adlasso",
                              nlam = 1000, lam.max.factor = 2000, lam.min.factor = 1e-10, 
                              penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
               ), screening = FALSE);  summary(fit10)
fit30 <- sofar(Y, X, ic.type = "GIC", nrank = r, 
               control = list(methodA = "adglasso", methodB = "adglasso",
                              nlam = 1000, lam.max.factor = 3000, lam.min.factor = 1e-12, 
                              penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
               ), screening = FALSE);  summary(fit30)


fit2 <- rssvd(Y,X, nrank = r); summary(fit2)
fit3 <- srrr(Y,X, method = "adglasso", ic.type = "GIC", nrank = r); summary(fit3)
fit4 <- rrs.fit(Y, X,  nrank = r); summary(fit4)




fit_all <- list()
fit_all[[1]] <- fit10
fit_all[[2]] <- fit30
fit_all[[3]] <- fit2
fit_all[[4]] <- fit3
fit_all[[5]] <- fit4
loss <- rep(0,5)
for(i in 1:5){
  fit1 <- fit_all[[i]]
  if(i==4){
    Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
    C = Ur%*%(Dr)%*%t(V)
    loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
    cat(i, loss, "\n")
  }
  if(i <=4 && i !=4){
    Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
    if(fit1$rank <2){C = Ur%*%t(V)*Dr
    }else{
      C = Ur%*%diag(Dr)%*%t(V)
    }
    loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
    cat(i, loss, "\n")
  }
  
  if(i ==5){
    C = fit1$coef
    loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
    cat(i, loss, "\n")
  }
}

#######Table 7
cat("adlasso", loss[1], "\n",
    "adglasso", loss[2], "\n",
    "rssvd", loss[3], "\n",
    "srrr", loss[4], "\n",
    "rrr", loss[5], "\n" )

