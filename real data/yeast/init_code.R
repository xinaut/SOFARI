setwd(".../yeast")
source("func_nearly.R")
source("func_weakly.R")
source("sim.R")
source("lasso_inference.r")
source("STRS.R")
source("helpers.nodewise.R")
source("functions.R")

library(rrpack)
library(secure)
library(MASS)
library(foreach)
library(doParallel)
library(parallel)
library(glmnet)



load(".../yeast_preprocess_data.RData")
dim(X)
dim(Y)
n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
#n = 112  p = 605 q = 54

##sofar estimation
est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 200, rank_tol = 1e-2, C = 2.01); est_rank


r = est_rank #r = 3
fit1 <- sofar(Y, X, ic.type = "GIC", nrank = 3,
              control = list(methodA = "adlasso", methodB = "adlasso",  lam.AB.factor = 0.32,
      nlam = 300, lam.max.factor = 5000, lam.min.factor = 1e-12,epsilon = 1e-8), screening = TRUE); summary(fit1)
U <- fit1$U
V <- fit1$V
length(which(U!=0)) 
length(which(V!=0)) 
fit1$D


# fit1 <- sofar(Y, X, ic.type = "GIC", nrank = 3,
#               control = list(methodA = "adlasso", methodB = "adlasso", lam.AB.factor = 0.1,
#                              nlam = 300, lam.max.factor = 5000, lam.min.factor = 1e-12,epsilon = 1e-8), screening = TRUE); summary(fit1)
# U2 <- fit2$U
# V2 <- fit2$V
# length(which(U2!=0)) #140
# length(which(V2!=0)) #40
