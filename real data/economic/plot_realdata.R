setwd(".../economic")
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggfittext)

load(file = "realdata2.RData")
##datauu is the debiased estimate of matrix U
df = as.data.frame(cbind(as.matrix(1:456), datauu)) #add an index colomn 1:456
colnames(df) <- c("X", "F1", "F2", "F3")
names_lag1 =  names1; names_lag2 =  names1; names_lag3 =  names1; names_lag4 =  names1
# names1 is all the variable names, since we added four lag to these varibles, all variable names are names1 with lag
for(i in 1:length(names1)){
  names_lag1[i] <- paste(names1[i], sep = "", "_lag1")
  names_lag2[i] <- paste(names1[i], sep = "", "_lag2")
  names_lag3[i] <- paste(names1[i], sep = "", "_lag3")
  names_lag4[i] <- paste(names1[i], sep = "", "_lag4")
}
names_X <- c(names_lag1,names_lag2,names_lag3,names_lag4)
df[,1] = names_X
ff1 = which(df[,2]!=0); ff2 = which(df[,3]!=0); ff3 = which(df[,4]!=0)
f_sum1 = c(ff1, ff2, ff3)
f_sum = unique(f_sum1)
f_sum = f_sum[order(f_sum)]
df = df[f_sum,]
#ff1, ff2, ff3 are nonzero location of datauu U matrix
#the barplot shows all the nonzero factors

df = melt(df)                   
head(df)                       

p_1 = ggplot(df, aes(x=factor(X,levels =unique(X)), 
                   y=value, 
                   fill=factor(variable,levels = unique(variable)), 
))+
  labs(
    x="",   
    y="",   
    fill=""
  )

p_1 +  geom_bar(
  position="stack",
  stat="identity"
)  + theme_bw() + theme(panel.grid = element_blank(),
                       axis.text.x = element_text(size = 0.2,  color = "black", angle = 0, hjust = 0.5,
                                                  vjust = 0.5),
                     
                     legend.position = "top",
                  
                     text=element_text(size=12, family="serif"))   + geom_text_repel(data=subset(df, abs(df[,3]) > 3), aes(label=X),vjust = -0.01,  size = 3.5) + geom_line(
) + labs(x="",y="",) + scale_x_discrete(name = "",labels = ""
          )  +scale_fill_manual(values=c("#FF7256", "#66CD00", "#B23AEE")) 
##in the plot 
#the command `geom_text_repel(data=subset(df, abs(df[,3]) > 3), aes(label=X),vjust = -0.01,  size = 3.5)'
#add the name of variables whose absolute value are larger than 3
#These variables are treated as improtant 
ggsave("plot22.png", width = 11, height = 8, dpi=300)
ggsave("plot1006.pdf", width = 12, height = 8)
