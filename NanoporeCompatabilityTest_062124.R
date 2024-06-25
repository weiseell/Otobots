### Testing Nanopore results with existing model

#load packages
library(tidyverse)
library(glmnet)
library(MASS)

#load data
load("Output/AgeLengthCpG_modelrun_021524.RData")
trim <- read.table("C:/Users/mvanw/Downloads/Barcode3_21-SS-HIP-5533_testtrim.bedMethyl",header = F,sep = "\t",stringsAsFactors = F)
notrim <- read.table("C:/Users/mvanw/Downloads/Barcode3_21-SS-HIP-5533_test.bedMethyl",header = F,sep = "\t",stringsAsFactors = F)

#selecting just percent mod and CpG position
pm1 <- trim[,c(1,3,11)]
colnames(pm1) <- c("CHROM","END","PERCMOD")
pm2 <- notrim[,c(1,3,11)]
colnames(pm2) <- c("CHROM","END","PERCMOD")
#reformatting to match Illumina individual data
pm1 <- pm1 %>% 
  mutate(CpG = paste(CHROM,END,sep = "."))
pm1 <- t(pm1[,c(4,3)])
colnames(pm1) <- pm1[1,]
pm1plot <- pm1[2,]

pm2 <- pm2 %>% 
  mutate(CpG = paste(CHROM,END,sep = "."))
pm2 <- t(pm2[,c(4,3)])
colnames(pm2) <- pm2[1,]
pm2plot <- pm1[2,]

best1 <- pm1[2,which(colnames(pm1) %in% CpGlist_best_comb$CpGlist_best)]
best2 <- pm2[2,which(colnames(pm2) %in% CpGlist_best_comb$CpGlist_best)]


plot(pm1plot,pm2plot)

#select data and run glmnet to select sites
y <- traindat$LogAge

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% names(best1))])
model1 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0.05,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))

bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])

#predict testdat
ytest <- c(testdat$LogAge,log(4),log(15))
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% names(best1))])
xtest <- rbind(xtest,as.numeric(best1),as.numeric(best2))

#predict and correlate train data
pred_logage <- predict(model1,newx = bestCpG_train,type = "response",s = model1$lambda.min)
train_cor <- cor.test(y,pred_logage)
#plot(y,pred_logage,col = traindat$ObsSex)

plot(model1, xvar = "dev", label = TRUE)

#predict and correlate test data
pred_logage_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.min)
test_cor <- cor.test(ytest,pred_logage_test)

plot(ytest,pred_logage_test)




