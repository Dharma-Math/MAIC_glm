library(brglm2)

DataM <- read.csv('DataM_logit_R0.csv')
ys <- as.numeric(DataM[,1])
Xs <- cbind("Intercept" = 1, as.matrix(DataM[,-1]))

result <- glm(Y~ X2 + X3 + X4, family=binomial(link = "logit"), data=DataM)
source("MAIC_glm.R") 

MAIC <- MAIC_glm(result,Xs,ys)

