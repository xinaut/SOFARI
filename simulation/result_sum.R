setwd(".../simulation")
source("func_nearly.R")
source("func_weakly.R")
source("sim.R")
source("STRS.R")
source("helpers.nodewise.R")
source("functions.R")
library(rrpack)
library(secure)
library(MASS)
library(glmnet)
library(foreach)
library(doParallel)
source("main of secure and sofar.R")
###parallel
start <- (proc.time())[3][[1]]
no_cores <- detectCores(no_cores)-4
cl <- makeCluster(no_cores)
registerDoParallel(cl)
#clusterEvalQ(cl, library(rrpack,MASS))
#set.seed(1235)
set.seed(12345)
iter = 1000
res <- foreach( ii = 1:iter, .combine = 'rbind', .errorhandling = "remove",
                .packages = c("rrpack", "MASS","glmnet","parallel"))%dopar% pare_sofar_nonsparse(ii)
stopCluster(cl)
end <- (proc.time())[3][[1]]
print(paste(' time = ', (end-start), 's', sep=''))

temp_res <- sum_res(res)

res_nonsofar <- temp_res

###single
start <- (proc.time())[3][[1]]
set.seed(1235)
res <- list()
iter <- 1000
for (ii in 1:iter) {
  #result <- pare_secure(ii)
  result <- pare_sofar_nonsparse(ii)
  res[[ii]] <- result
  print(ii)
}
end <- (proc.time())[3][[1]]
print(paste(' time = ', (end - start), 's', sep=''))



##function: sum_res (see the end of this file)
temp_res <- sum_res(res)


##save result
secure_low = res
sum_res(secure_low)

secure_high = res
sum_res(secure_high)

sofar_low = res
sum_res(sofar_low)

sofar_high = res
sum_res(sofar_high)



sum_res <- function(res){
  ##24 mean results for all
  result1 <- res[,1]
  ##the normal variables for u1, u2, u3 and d1, d2, d3
  re_u1 <- res[,2]; re_u2 <- res[,3]; re_u3 <- res[,4]
  re_d <- res[,5]
  ##individual CI and length for u1j, u2j, u3j
  ci_u1j <- res[,6]; len_u1j <- res[,7]
  ci_u2j <- res[,8]; len_u2j <- res[,9]
  ci_u3j <- res[,10]; len_u3j <- res[,11]
  #########first index
  iter = length(result1)
  print(iter)
  result1.matrix <- matrix(0, iter, 24)
  for (i in 1:iter) {
    result1.matrix[i,] = result1[[i]]
  }
  colnames(result1.matrix)<- c("U1covs0","U1covs1","U1covall","U1len_s0", "U1len_s1","U1lenall",
                               "U2covs0","U2covs1","U2covall","U2len_s0", "U2len_s1","U2lenall",
                               "U3covs0","U3covs1","U3covall","U3len_s0", "U3len_s1","U3lenall",
                               "d1covs0","d1len_s0",
                               "d2covs0","d2len_s0",
                               "d3covs0","d3len_s0")
  #cat(round(colMeans(result1.matrix),3), "\n")
  print(round(colMeans(result1.matrix),3))
  #######normal plot
  pp = length(re_u1[[1]])
  re_u1.m <- matrix(0, iter, pp); re_u2.m <- matrix(0, iter, pp); re_u3.m <- matrix(0, iter, pp)
  for (i in 1:iter) {
    re_u1.m[i,] = re_u1[[i]]
    re_u2.m[i,] = re_u2[[i]]
    re_u3.m[i,] = re_u3[[i]]
  }
  # hist(re_u2.m[,4], freq = F)
  # hist(re_u1.m[,1], xlim = c(-3,3),  freq = F)
  # hist(re_u2.m[,1], xlim = c(-3,3), breaks =  20, freq = F)
  # hist(re_u3.m[,1], xlim = c(-3,3), breaks =  10, freq = F)
  
  re_d.m <- matrix(0, iter, 3)
  for (i in 1:iter) {
    re_d.m[i,] = re_d[[i]]
  }
  # hist(re_d.m[,1],  breaks =  20, freq = F)
  # hist(re_d.m[,2], xlim = c(-3,3), breaks =  20, freq = F)
  # hist(re_d.m[,3], xlim = c(-3,3), breaks =  15, freq = F)
  
  
  ###CI and length
  ci_u1j.matrix <- matrix(0, iter, pp); ci_u2j.matrix <- matrix(0, iter, pp);  ci_u3j.matrix <- matrix(0, iter, pp)
  len_u1j.matrix <- matrix(0, iter, pp); len_u2j.matrix <- matrix(0, iter, pp);  len_u3j.matrix <- matrix(0, iter, pp)
  for (i in 1:iter) {
    ci_u1j.matrix[i,] = ci_u1j[[i]]
    ci_u2j.matrix[i,] = ci_u2j[[i]]
    ci_u3j.matrix[i,] = ci_u3j[[i]]
    len_u1j.matrix[i,] = len_u1j[[i]]
    len_u2j.matrix[i,] = len_u2j[[i]]
    len_u3j.matrix[i,] = len_u3j[[i]]
  }
  print(colMeans(ci_u1j.matrix))
  print(colMeans(ci_u2j.matrix))
  print(colMeans(ci_u3j.matrix))
  print(round(colMeans(len_u1j.matrix),3))
  print(round(colMeans(len_u2j.matrix),3))
  print(round(colMeans(len_u3j.matrix),3))
  # colMeans(ci_u1j.matrix)
  # colMeans(ci_u2j.matrix)
  # colMeans(ci_u3j.matrix)
  
  # cat(colMeans(ci_u1j.matrix), "\n", 
  #     colMeans(ci_u2j.matrix), "\n",
  #     colMeans(ci_u3j.matrix), "\n",
  #     round(colMeans(len_u1j.matrix),3), "\n",
  #     round(colMeans(len_u2j.matrix),3), "\n",
  #     round(colMeans(len_u3j.matrix),3), "\n"
  #     )
  
}




sum_res_single <- function(res){
  # ##24 mean results for all
  # result1 <- res[,1]
  # 
  # 
  # ##the normal variables for u1, u2, u3 and d1, d2, d3
  # re_u1 <- res[,2]; re_u2 <- res[,3]; re_u3 <- res[,4]
  # re_d <- res[,5]
  # ##individual CI and length for u1j, u2j, u3j
  # ci_u1j <- res[,6]; len_u1j <- res[,7]
  # ci_u2j <- res[,8]; len_u2j <- res[,9]
  # ci_u3j <- res[,10]; len_u3j <- res[,11]
  
  ####begin
  result_list <- vector("list", 11)
  for (j in 1:11) {
    result_list[[j]] <- do.call(rbind, lapply(res, function(x) x[[j]]))
  }
  #########first index
  result1 = result_list[[1]]
  iter = nrow(result1)
  print(iter)
  
  
  colnames(result1)<- c("U1covs0","U1covs1","U1covall","U1len_s0", "U1len_s1","U1lenall",
                               "U2covs0","U2covs1","U2covall","U2len_s0", "U2len_s1","U2lenall",
                               "U3covs0","U3covs1","U3covall","U3len_s0", "U3len_s1","U3lenall",
                               "d1covs0","d1len_s0",
                               "d2covs0","d2len_s0",
                               "d3covs0","d3len_s0")
  #cat(round(colMeans(result1.matrix),3), "\n")
  print(round(colMeans(result1),3))

  
}
