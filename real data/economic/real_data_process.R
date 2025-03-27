rm(list = ls())
setwd(".../economic")
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


file <- '.../2015-04.csv'
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
