## This script is used to clean raw data and perform marginal screening
## The collected data is saved as "yeast_preprocess_data.RData"

### Set working directory to the `SOFARI_code` folder of the reproducibility_materials
### Please replace the path below with your own path containing the `SOFARI_code` folder, remove the `#`, and then run the following line
# setwd("~/SOFARI_code/real data/yeast")  ### <-- Replace this with your own path

load("yeast.rda")
marker <- t(yeast$marker)
n <- dim(marker)[1]
p <-dim(marker)[2]


#####Combine markers differ at most 1
group <- vector()
group[1] <- 1

#####selected markers in combined groups
markerid <- vector()
markerid[1] <- 1

j <- 1
i <- 1
while(i<p){
  while(i<p && sum(yeast$marker[i,]==yeast$marker[i+1,])>n-2){
    group[i+1] <- group[i]
    i <- i + 1
  }
  j <-  j + 1
  if(i<p) markerid[j] <- i+1
  if(i<p) group[i+1] <- j
  i <- i + 1
}

#####Selected genes on MAPK pathway
MAPK <- c("MFA1","MFA2","STE2","STE3","SLG1","WSC2","WSC3","MID2","SHO1",
          "SLN1","SHO1","RAS2","GPA1","STE18","YPD1","CDC42","CDC24","RHO1","FKS1","CDC42",
          "BEM1","STE20","BNI1","PKC1","STE20","SSK1","STE20","STE11","BCK1","STE11","SSK22","SSK2",
          "STE11","STE5","STE7","MKK1","MKK2","PBS2","STE7","FUS3","MSG5","SLT2","MLP1","MLP2","HOG1",
          "KSS1","DIG1","DIG2","STE12","FAR1","RLM1","SWI4","SWI6","MCM1","MSN2","MSN4","DIG1","DIG2",
          "STE12","MCM1","TEC1","FUS1","GSC2","GLO1","CTT1")

genenames <- rownames(yeast$exp.pos)
matchid <- match(MAPK,genenames)
geneid <- unique(na.omit(matchid))

####Complete data; including all the combined markers
Yall <- t(yeast$exp[geneid,])
colnames(Yall) <- genenames[geneid]
Xall <- t(yeast$marker[markerid,])
qall <- dim(Yall)[2]
pall <- dim(Xall)[2]


#####Select usuful pathways using marginal screening
pvalue.cutoff <- 0.05
markersig <- vector()
pvaluemat <- matrix(nrow=pall,ncol=qall,1)
for(ip in 1:pall){
  
  for(jq in 1:qall){
    pvaluemat[ip,jq] <- cor.test(Xall[,ip],Yall[,jq])$p.value
  }
  if(sum(pvaluemat[ip,]<pvalue.cutoff)>2) markersig <- c(markersig,ip)
}
markerinsig <- vector()
Y <- scale(Yall,center=FALSE,scale=FALSE)
X <- scale(Xall[,c(markersig,markerinsig)],center=FALSE,scale=FALSE)
colnames(Y) <- genenames[geneid]
n <- nrow(Y)
q <- ncol(Y)
p <- ncol(X)
save(X,Y,file = "yeast_preprocess_data.RData")
