rm(list = ls())
setwd("C:/Users/Xin Zhou/SOFARI")
source("func_nearly.R")
source("func_weakly.R")
source("sim.R")
source("lasso_inference.r")
source("STRS.R")
source("helpers.nodewise.R")

library(rrpack)
library(secure)
library(MASS)
library(foreach)
library(doParallel)
library(parallel)


DIM <- function(x) {
  if( is.null( dim(x) ) )
    return( length(x) )
  dim(x)
}

##find the id of the name
change_name_id <- function(res_name, cname){
  id_re <- NULL
  for(i in 1:length(res_name)){
    id_c <- which(cname == res_name[i])
    id_re <- c(id_re, id_c)
  }
  id_re <- id_re[order(id_re)]
  return(id_re)
}

##standarized X and Y
standarize_data <- function(mdata){
  n0 <- dim(mdata)[1]; p0 <- dim(mdata)[2]
  x <- mdata - matrix(1,NROW(mdata),1)%*%colMeans(mdata)
  mX <- matrix(0,n0,p0)
  for (i in 1:ncol(mdata)){
    mX[,i] <- x[,i]/sd(x[,i])
  }
  return(mX)
}

#simple person test for the corrlation of X
test.corr <- function(X, ratio){
  n <- dim(X)[1]; p <- dim(X)[2]
  corx = matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      corx[i,j] = as.numeric(cor.test(X[,i],X[,j],method="pearson")$estimate) 
    }
  }
  diag(corx) = 0
  #high_cor <- matrix(0,p,p)
  ini_index <- seq(1,p,1)
  high_cor <- list()
  for(i in 1:p){
    high_cor[[i]] <- which(corx[,i] > ratio)
    #temp <-  which(corx[,i] > 0.9)
    # if(length(high_cor[[i]]) > 0){
    #   temp = high_cor[[i]]
    # }
  }
  return(list(corx = corx, fcor = high_cor))
}

#return the final X replaced with residuals
regress_on_X <- function(yX, X0, id){
  #yX = group_X; X0 <- fix_X; id <- fg_id
  p0 = dim(yX)[2]
  for(i in 1:p0){
    y_temp <- yX[,i]
    id_temp <- which(id == i)
    for(j in id_temp){
      X_temp <- X0[,j]
      fit <- lm(y_temp ~ X_temp)
      X0[,j] <- fit$residuals
    }
  }
  f_X <- cbind(yX, X0)
  return(f_X)
}

#add lag to X and Y
AR.lag <-function(lag, X, Y){
  data_organized <- X
  response <- Y
  # construct the design
  AR <- list()
  len_res <- dim(response)[1]
  for (t in 1:lag) {
    AR[[t]] <- data_organized[t:(len_res - lag - 1 + t),]
    AR[[t]] <- cbind(response[t:(len_res - lag - 1 + t), ],  AR[[t]])
  }
  
  
  t <- lag + 1
  response <- response[t:(len_res - lag - 1 + t), ]
  # record the dates of responses; from 7/1/1960 to 12/1/2014
  #date <- date[t:(length(date) - lag - 1 + t)]
  # date[475] # 1/1/2000; so we have 474 observations for trainning the model
  X <- NULL
  for(t in 1:lag){
    X <- cbind(X, AR[[t]])
  }
  return(list(X = X, Y = response))
  
}


file <- 'C:/Users/Xin Zhou/OneDrive/Rcode/SOFARI/real data/rd_code/2015-04.csv'
data_ <- read.csv(file)
# DIM(data_)
# the first row is the category index and the last row is the column numbering I added to the file.

# record the categories of each column
categories <- data_[1,-1]
cname <- colnames(data_[,-1])
# use data from 1960:1 to 2014:12; in the meantime the first row, which records the cotegories, is deleted
data_ <- data_[1:673,]
data_ <- data_[14:DIM(data_)[1],]

# the first column is the date column. Delete it.
date <- data_[, 1]
data_ <- data_[, -1]

# ---- data transformation ----
data_organized <- array(0, c(dim(data_)[1] - 2, dim(data_)[2]))
colnames(data_organized) <- colnames(data_)
for (j in 1:length(categories)) {
  # if (j == 113) {
  #   # if (j == 106) {
  #   # 113 is the CPI for all items
  #   # perchange changes of CPI: ((x_{t} - x_{t-1}) / x_{t-1}) * 100
  #   x <- data_[, j]
  #   
  #   x <- ((x[2:length(x)] / x[1:(length(x) - 1)]) - 1) * 100
  #   # 
  #   data_organized[, j] <- x[2:length(x)]
  #   
  # } else {
  if (categories[j] == 7) {
    x <- data_[, j]
    x <- (x[2:length(x)] / x[1:(length(x) - 1)]) - 1
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    data_organized[, j] <- x
  }
  if (categories[j] == 6) {
    x <- log(data_[, j])
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    data_organized[, j] <- x
  }
  if (categories[j] == 2) {
    x <- data_[, j]
    x <- x[2:length(x)]
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    data_organized[, j] <- x
  }
  if (categories[j] == 3) {
    x <- data_[, j]
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    data_organized[, j] <- x
  }
  if (categories[j] == 4) {
    x <- log(data_[, j])
    data_organized[, j] <- x[c(-1, -2)]
  }
  if (categories[j] == 5) {
    x <- log(data_[, j])
    x <- x[2:length(x)] - x[1:(length(x) - 1)]
    data_organized[, j] <- x[-1]
  }
  if (categories[j] == 1) {
    x <- data_[, j]
    data_organized[, j] <- x[3:length(x)]
  }
  #} #end for else
}


# record the dates of responses; from 3/1/1960 to 12/1/2014
date <- date[c(-1, -2)]

# record the response; from 3/1/1960 to 12/1/2014
#response <- data_organized[, 113]

##response id
res_name <- c("RPI", "INDPRO", "CUMFNS", "UNRATE", "PAYEMS", "CES0600000007", "NAPMNOI",
              "HOUST", "DPCERA3M086SBEA", "CMRMTSPLx", "FEDFUNDS", "T1YFFM", "T10YFFM", "BAAFFM",
              "EXUSUKx", "PPIFGS", "PPICMM", "CPIAUCSL", "PCEPI", "S.P.500") #20

id_response <- change_name_id(res_name,cname)
length(id_response) #check the length


response <- data_organized[, id_response]
# remove the response (percehtage chages of CPI) and series with NA from the design
# series 64, 66, 101, 130 contain NA at the beginning
na.series <- c("UMCSENTx", "TWEXMMTH", "ANDENOx", "ACOGNO")
id_na <- change_name_id(na.series, cname)
#data_organized <- data_organized[, -c(id_response, 64, 66, 101, 130)]
data_organized <- data_organized[, -c(id_response, id_na)]

# double check there are no NA in the data set
sum(is.na(data_organized))


##delete the high correlation X
X <- data_organized
#######test corrlation of X
n <- dim(X)[1]; p <- dim(X)[2]; 
corx = matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    corx[i,j] = as.numeric(cor.test(X[,i],X[,j],method="pearson")$estimate) 
  }
}
diag(corx) = 0
#high_cor <- matrix(0,p,p)
ini_index <- seq(1,p,1)
high_cor <- list()
for(i in 1:p){
  high_cor[[i]] <- which(corx[,i] > 0.9)
  #temp <-  which(corx[,i] > 0.9)
  # if(length(high_cor[[i]]) > 0){
  #   temp = high_cor[[i]]
  # }
}
high_cor

index_del <- c(3, 5, 15, 27, 31, 41, 42, 43,  45, 46, 72, 73, 75, 81, 82, 102  )
index_del <- index_del[order(index_del)]
X <- X[,-index_del]
dim(X) #658 94


data_organized <- X
##end of decorrlation  step
ar.fit <- AR.lag(4, data_organized, response)
X <- ar.fit$X; Y <- ar.fit$Y
dim(X) # 654 456
dim(Y) # 654 20
#delete the four factors
# X <- cbind(X, svd(X)$u[, 1:4])
# dim(X)


#X <- scale(X); Y <- scale(Y)

##standarized X and Y
X <- standarize_data(X)
Y <- standarize_data(Y)

names1 <- c(colnames(response), colnames(data_organized))
# Y, X are ready original Y 20 X 94
# X:  20 response t=1 lag + 94 predictors t = 1 lags  (from index 1 to 114)
#   + 20 response t=2 lag + 94 predictors t = 2 lags  (from index 115 to 228)
#   + 20 response t=3 lag + 94 predictors t = 3 lags  (from index 229 to 342)
#   + 20 response t=4 lag + 94 predictors t = 4 lags  (from index 343 to 456)
# return 20 Y and 456 X
# ---- the end of data process ---- 

write.csv(X, file = "final_XX.csv", row.names = F)
write.csv(Y, file = "final_YY.csv", row.names = F)

#############end of data process############


####### do prediction
X <- read.csv("final_XX.csv")
Y <- read.csv("final_YY.csv")
X <- as.matrix(X); Y <- as.matrix(Y)

#trian the data
X0 <- X; Y0 <- Y
sam <- 474
train_sam <- 1:sam
Y<- Y0[train_sam, ]; X <- X0[train_sam, ]
sam / dim(X0)[1]
# X <- X0; Y<-Y0
est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
#"STRS-MC" "STRS-DB" "SSTRS"
r = 2
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
fit5 <- secure.path(Y,X,nrank= r)



fit_all <- list()
fit_all[[1]] <- fit10
fit_all[[2]] <- fit30
fit_all[[3]] <- fit2
fit_all[[4]] <- fit3
fit_all[[5]] <- fit4
fit_all[[6]] <- fit5
loss <- rep(0,6)
for(i in 1:6){
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
  if(i ==6){
    C = fit1$C.est
    loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
    cat(i, loss, "\n")
  }
}
cat("adlasso", loss[1], "\n",
    "adglasso", loss[2], "\n",
    "rssvd", loss[3], "\n",
    "srrr", loss[4], "\n",
    "rrr", loss[5], "\n",
    "secure", loss[6], "\n" )





# 
# #trian the data
# X0 <- X; Y0 <- Y
# sam <- 474
# train_sam <- 1:sam
# Y<- Y0[train_sam, ]; X <- X0[train_sam, ]
# sam / dim(X0)[1]
# # X <- X0; Y<-Y0
# est_rank <- STRS(Y, X, type = "SSTRS", rep_MC = 200, rank_tol = 1e-8, C = 2.01); est_rank #4
# est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
# est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
# #"STRS-MC" "STRS-DB" "SSTRS"
# r = 2
# fit10 <- sofar(Y, X, ic.type = "GIC", nrank = r, 
#                control = list(methodA = "adlasso", methodB = "adlasso",
#                               nlam = 1000, lam.max.factor = 2000, lam.min.factor = 1e-10, 
#                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
#                ), screening = FALSE);  summary(fit10)
# fit30 <- sofar(Y, X, ic.type = "GIC", nrank = r, 
#                control = list(methodA = "adglasso", methodB = "adglasso",
#                               nlam = 1000, lam.max.factor = 3000, lam.min.factor = 1e-12, 
#                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
#                ), screening = FALSE);  summary(fit30)
# fit10$lam.id
# 
# fit30$lam.id
# 
# # > fit10$lam.id  nlam = 300, lam.max.factor = 2000, lam.min.factor = 1e-10, 
# # [1] 186 186 166 166
# # > fit30$lam.id   nlam = 300, lam.max.factor = 3000, lam.min.factor = 1e-12, 
# # [1] 200 193 181 185
# # fit1 = fit10
# # fit1 = fit30
# 
# fit2 <- rssvd(Y,X, nrank = r); summary(fit2)
# fit3 <- srrr(Y,X, method = "adglasso", ic.type = "GIC", nrank = r); summary(fit3)
# fit4 <- rrs.fit(Y, X,  nrank = r); summary(fit4)
# fit5 <- secure.path(Y,X,nrank= r)
# 
# 
# 
# fit_all <- list()
# fit_all[[1]] <- fit10
# fit_all[[2]] <- fit30
# fit_all[[3]] <- fit2
# fit_all[[4]] <- fit3
# fit_all[[5]] <- fit4
# fit_all[[6]] <- fit5
# loss <- rep(0,6)
# for(i in 1:6){
#   fit1 <- fit_all[[i]]
#   if(i==4){
#     Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
#     C = Ur%*%(Dr)%*%t(V)
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   if(i <=4 && i !=4){
#     Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
#     if(fit1$rank <2){C = Ur%*%t(V)*Dr
#     }else{
#       C = Ur%*%diag(Dr)%*%t(V)
#     }
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   
#   if(i ==5){
#     C = fit1$coef
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   if(i ==6){
#     C = fit1$C.est
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
# }
# cat("adlasso", loss[1], "\n",
#     "adglasso", loss[2], "\n",
#     "rssvd", loss[3], "\n",
#     "srrr", loss[4], "\n",
#     "rrr", loss[5], "\n",
#     "secure", loss[6], "\n" )
# loss1 = loss
# ##rank 2
# # adlasso 0.9211238 
# # adglasso 0.9353401 
# # rssvd 1.00214 
# # srrr 0.9679916 
# # rrr 1.775065 
# # secure 1.409139
# #################end of training
# 
# 
# 
# 
# 
# 
# #trian the data
# X0 <- X; Y0 <- Y
# sam <- 474
# train_sam <- 1:sam
# Y<- Y0[train_sam, ]; X <- X0[train_sam, ]
# 
# # X <- X0; Y<-Y0
# est_rank <- STRS(Y, X, type = "SSTRS", rep_MC = 200, rank_tol = 1e-8, C = 2.01); est_rank #4
# est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
# est_rank <- STRS(Y, X, type = "STRS-MC", rep_MC = 500, rank_tol = 1e-8, C = 2.01); est_rank 
# #"STRS-MC" "STRS-DB" "SSTRS"
# r = 3
# fit10 <- sofar(Y, X, ic.type = "BIC", nrank = r, 
#                control = list(methodA = "adlasso", methodB = "adlasso",
#                               nlam = 100, lam.max.factor = 1500, lam.min.factor = 1e-12, 
#                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
#                ), screening = FALSE);  summary(fit10)
# fit30 <- sofar(Y, X, ic.type = "GIC", nrank = r, 
#                control = list(methodA = "adglasso", methodB = "adglasso",
#                               nlam = 300, lam.max.factor = 3000, lam.min.factor = 1e-12, 
#                               penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-6
#                ), screening = FALSE);  summary(fit30)
# fit10$lam.id
# 
# fit30$lam.id
# 
# fit1 = fit10
# fit1 = fit30
# 
# fit2 <- rssvd(Y,X, nrank = r); summary(fit2)
# fit3 <- srrr(Y,X, method = "glasso", ic.type = "BIC", nrank = r); summary(fit3)
# fit4 <- rrs.fit(Y, X,  nrank = r); summary(fit4)
# fit5 <- secure.path(Y,X,nrank= r)
# 
# 
# 
# fit_all <- list()
# fit_all[[1]] <- fit10
# fit_all[[2]] <- fit30
# fit_all[[3]] <- fit2
# fit_all[[4]] <- fit3
# fit_all[[5]] <- fit4
# fit_all[[6]] <- fit5
# loss <- rep(0,6)
# for(i in 1:6){
#   fit1 <- fit_all[[i]]
#   if(i==4){
#     Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
#     C = Ur%*%(Dr)%*%t(V)
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   if(i <=4 && i !=4){
#     Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
#     if(fit1$rank <2){C = Ur%*%t(V)*Dr
#     }else{
#       C = Ur%*%diag(Dr)%*%t(V)
#     }
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   
#   if(i ==5){
#     C = fit1$coef
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
#   if(i ==6){
#     C = fit1$C.est
#     loss[i] <- norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
#     cat(i, loss, "\n")
#   }
# }
# cat("adlasso", loss[1], "\n",
#     "adglasso", loss[2], "\n",
#     "rssvd", loss[3], "\n",
#     "srrr", loss[4], "\n",
#     "rrr", loss[5], "\n",
#     "secure", loss[6], "\n" )
# loss1 = loss
# ##rank = 2
# # adlasso 1.971567 
# # adglasso 1.396153 
# # rssvd 1.38425 
# # srrr 1.33542 
# # rrr 2.14795 
# # secure 1.41327
# 
# C = fit5$C.est
# # adlasso 1.333522 
# # adglasso 1.337779 
# # rssvd 1.33222 
# # srrr 1.337489 
# # rrr 1.889156 
# # secure 1.41327 
# 
# # adlasso 1.336303 
# # adglasso 1.343324 
# # rssvd 1.33222 
# # srrr 1.337489 
# # rrr 1.889156 
# 
# fit1 = fit2
# 
# Ur = fit1$U; Dr = fit1$D; V = fit1$V; 
# C = Ur%*%t(V)*Dr
# C = Ur%*%diag(Dr)%*%t(V)
# 
# norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/((dim(X0)[1] - sam)*q)
# 
# 
# fit1 = fit3
# Ur = fit1$U; Dr = fit1$D; V = fit1$V; C = Ur%*%(Dr)%*%t(V)
# 
# C = fit4$coef
# 
# 
# #####estimate and inference from the all data#####
# X <- X0; Y<-Y0
# est_rank <- STRS(Y, X, type = "STRS-DB", rep_MC = 200, rank_tol = 1e-2, C = 2.01); est_rank 
# 
# n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2]
# r = est_rank #r = 4
# fit1 <- sofar(Y, X, ic.type = "BIC", nrank = r, 
#               control = list(methodA = "adlasso", methodB = "adlasso",
#                              nlam = 100, lam.max.factor = 1500, lam.min.factor = 1e-10, 
#                              penA = TRUE, penB = TRUE, penD = TRUE, epsilon = 1e-8
#               ), screening = FALSE);  summary(fit1)
# 
# fit1$lam.id
# 
# Ur = fit1$U; Dr = fit1$D; V = fit1$V;
# dU = Ur%*%diag(Dr)
# r = fit1$rank
# 
# C = Ur%*%diag(Dr)%*%t(V)
# E_est = Y - X%*%C
# 
# # norm(E_est,"F")/norm(Y,"F")
# # norm(Y0[-train_sam, ] - X0[-train_sam, ]%*%C, "F")^2/(n*q)
# 
# Sigmae_null = matrix(0, q, q) #no meaning
# obj_sigmae <-  sim(p = q, n = n, delta = 2, thresholding = 'al', data = E_est, cov_true = Sigmae_null)
# Sigmae <- obj_sigmae$cov_est
# 
# # X <- scale(X, center = TRUE, scale = TRUE)
# # Y <- scale(Y, center = TRUE, scale = TRUE)
# Xsigma <- t(X)%*%X/n
# #estimate the precision matrix
# #Theta <- solve(Xsigma)
# # ghat=cv.glasso1(X,0)
# # Theta=ghat$wi
# # initialization of ISEE parameters
# # regfactor = "sqrt"  # or "one", "sqrt"
# # npermu = 2         # or >= 2
# # sis.use = 0       # or 1, whether to use SIS for screening
# # bia.cor = 1        # or 0, whether to apply bias correction for ISEE
# # obj.c = isee(X, regfactor, npermu, sis.use, bia.cor = 1)
# # Theta <- obj.c$Omega.isee.c
# # #Theta=MBLasso(X,lambda=1,w.mb = rep(1,p))
# obj<- score.nodewiselasso(X, wantTheta = TRUE,
#                           verbose = FALSE,
#                           lambdaseq = "quantile",
#                           parallel = FALSE,
#                           ncores = 6,
#                           oldschool = FALSE,
#                           lambdatuningfactor = 1,
#                           cv.verbose = FALSE,
#                           do.ZnZ = TRUE)
# Theta = obj$out #0.01357721
# 
# Theta <- read.csv("theta.csv")
# Theta <- as.matrix(Theta)
# #write.csv(Theta, file = "Theta_node.csv", row.names = F)
# # > obj$bestlambda
# # [1] 0.0005060029
# norm(diag(p)- Xsigma%*%Theta, "M")
# # 
# # obj2<- score.nodewiselasso(X, wantTheta = TRUE,
# #                           verbose = FALSE,
# #                           lambdaseq = "quantile",
# #                           parallel = FALSE,
# #                           ncores = 6,
# #                           oldschool = FALSE,
# #                           lambdatuningfactor = 1,
# #                           cv.verbose = FALSE,
# #                           do.ZnZ = FALSE)
# # 
# # obj3<- score.nodewiselasso(X, wantTheta = TRUE,
# #                            verbose = FALSE,
# #                            lambdaseq = "quantile",
# #                            parallel = FALSE,
# #                            ncores = 6,
# #                            oldschool = FALSE,
# #                            lambdatuningfactor = "lambda.1se",
# #                            cv.verbose = FALSE,
# #                            do.ZnZ = FALSE)
# 
# # #Theta <- diag(p)*n
# #norm(diag(p)- Xsigma%*%Theta, "M")
# #Theta <- solve(Xsigma + 0.01*diag(p))
# Theta <- InverseLinfty(Xsigma, n,  resol=1.5, maxiter=50, threshold=1e-2)
# #fina = Theta
# norm(diag(p)- Xsigma%*%Theta, "M")
# # library(huge)
# # out.mb = huge(X, method = "glasso")
# # bbv <- huge.select(out.mb,criterion = "ebic")
# # Theta = bbv$opt.icov
# # norm(diag(p)- Xsigma%*%Theta, "M")
# # View(Theta)
# # dU[,c(1,2,3)] = dU[,c(1,3,2)]
# # V[,c(1,2,3)] = V[,c(1,3,2)]
# ##############start debias procedure
# uk0 = rep(0,p)
# dk0 = 0
# #de-biased 
# #the first layer
# k = 1
# u_wh = list()
# #ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u1_wh, u2_wh)
# ukdk1_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
# #ukdk1_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
# uk1 <- ukdk1_weak$u
# dk1 <- ukdk1_weak$d
# dk1
# 
# varobj1 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae,  S_u)
# varuk1 <- varobj1$varu
# mean(varuk1)
# vardk1 <- varobj1$vard
# 
# uk = uk1; varuk = varuk1
# index = NULL
# len2 = 0
# zz = 2.81
# #zz = 1.96
# low <- rep(0,p); up <- rep(0,p)
# for(j in 1:p){
#   low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n) #[which(dU[,1]==0)]
#   low[j] <- low_ci
#   up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
#   up[j] <- up_ci
#   len2 = len2 + abs(up_ci - low_ci)
#   if( up_ci < 0 || low_ci >0) {index = c(index, j)}
# }
# length(index)
# len2/p
# 
# # dev.new()
# # plot(uk1, ylim = c(-1,1), main='Confidence Intervals based on de-biased LASSO', ylab='', xlab = 'Coefficients');
# # #points(rep(0,p),col="blue");
# # lines(rep(0,p),col="blue");
# # #lines(c(up[i],low[i]),col="red");
# # lines(up,col="red");
# # lines(low,col="red");
# # legend('topright', legend=c('LASSO','de-biased LASSO','Ground-truth','Confidence Intervals'), col=c('black', 'blue','green','red'), pch=c(1,1,1,NA_integer_), lty = c(0,0,0,1))
# 
# #the second layer
# k = 2
# uk1_t <- uk1
# uk1_t[abs(uk1) < log(n)/sqrt(n)] = 0
# S_u1 = which(uk1_t != 0)
# S_u1
# S_u <- list(S_u1)
# u1_wh <- uk1_t
# u_wh <- list(u1_wh)
# ukdk2_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
# #ukdk2_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
# uk2 <- ukdk2_weak$u
# dk2 <- ukdk2_weak$d
# dk2
# varobj2 <- var_weakly_simple(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)
# #varobj2 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)
# #varobj2 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
# varuk2 <- varobj2$varu
# vardk2 <- varobj2$vard
# 
# uk = uk2; varuk = varuk2
# index = NULL
# len2 = 0
# for(j in 1:p){
#   low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n)
#   up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
#   len2 = len2 + abs(up_ci - low_ci)
#   if( up_ci < 0 || low_ci > 0) {index = c(index, j)}
# }
# length(index)
# len2/p
# 
# #the third layer
# k = 3
# uk1_t <- uk1
# uk1_t[abs(uk1) < log(n)/sqrt(n)] = 0 #quantile(abs(uk1),0.95)
# S_u1 = which(uk1_t != 0)
# S_u1
# #u1_wh <- -uk1_t
# 
# uk2_t <- uk2
# uk2_t[abs(uk2) < 1.2*log(n)/sqrt(n)] = 0
# #uk2_t[abs(uk2) < 1] = 0
# S_u2 = which(uk2_t != 0)
# S_u2
# S_u <- list(S_u1, S_u2)
# #u2_wh <- -uk2_t
# u_wh <- list(u1_wh, u2_wh)
# #u_wh <- list(rep(0,p), rep(0,p))
# ukdk3_weak <- uk_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, u_wh)
# #ukdk3_weak <- uk_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y)
# uk3 <- ukdk3_weak$u
# dk3 <- ukdk3_weak$d
# dk3
# 
# varobj3 <- var_weakly_simple(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)
# #varobj3 <- var_weakly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, S_u)
# #varobj3 <- var_nearly(k, r, p, Theta, dU, V, Xsigma, X, Y, Sigmae, b)
# varuk3 <- varobj3$varu
# vardk3 <- varobj3$vard
# 
# uk = uk3; varuk = varuk3
# index = NULL
# len2 = 0
# for(j in 1:p){
#   low_ci <-  uk[j]-sqrt(varuk[j])*zz/sqrt(n)
#   up_ci <-  uk[j] + sqrt(varuk[j])*zz/sqrt(n)
#   len2 = len2 + abs(up_ci - low_ci)
#   if( up_ci <= 0 || 0 >=up_ci) {index = c(index, j)}
# }
# length(index)
# len2/p
# 
