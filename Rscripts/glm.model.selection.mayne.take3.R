#### Load custom functions ####
## Log linear transformation
## From Horvath et al. 2022 PNAS 119:e2120887119, doi:10.1073/pnas.2120887119
logLin = Vectorize(function(x, maturity, ...) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

## Inverse log linear transformation
## From Horvath et al. 2022 PNAS 119:e2120887119, doi:10.1073/pnas.2120887119
logLin.inverse = Vectorize(function(y, maturity, ...) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

##set working directory
setwd("E:/MethylSeq_Otolith/Prelim.methyl.analysis/")

##load data
load("allCpGs_nosingletons.rda")
load("allCpGs_HPEISamples_nonas.rda")
sexdata <- read.csv("MethylWild_SexData.csv")

#load libraries
library(tidyverse)
library(glmnet)
library(MASS)


#remove CpGs with NAs for the model
nas <- colSums(is.na(allCpGs1))
allCpGs1 <- allCpGs1[,which(nas < 1)]
#write.csv(allCpGs1,"allCpGs_nonas.csv",append = F,quote = F)

#remove NW and mito data from data set
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NW",
                                    colnames(allCpGs1),perl = T)))]
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_009709.1",
                                    colnames(allCpGs1),perl = T)))]

#log age data
sexdata <- sexdata %>% 
  mutate(LogAge = log(Age))

#randomly select individuals for testing
testpick <- sexdata %>% 
  filter(DataType == "Wild") %>% 
  group_by(Age,ObsSex) %>% 
  slice_sample(n = 1)

#sort other individuals into training data
trainpick <- sexdata[which(!(sexdata$indivname %in% testpick$indivname)),]
trainpick <- trainpick %>% 
  filter(DataType == "Wild")

#select CpG sites for each group
testdat <- allCpGs1[which(allCpGs1$indivname %in% testpick$indivname),]
testdat <- full_join(testpick,testdat,by = "indivname")

traindat <- allCpGs1[which(!(allCpGs1$indivname %in% testpick$indivname)),]
traindat <- full_join(trainpick,traindat,by = "indivname")

#select data and run glmnet to select sites
y <- traindat$LogAge
x <- as.matrix(traindat[,-1:-5])

model1 <- cv.glmnet(x,y,family = "gaussian",alpha = 0.5,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.1se"))

bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% bestCpGnames)])
model2 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0,nfolds = 10)

#predict testdat
ytest <- testdat$LogAge
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% bestCpGnames)])

#predict and correlate train data
pred_logage <- predict(model2,newx = bestCpG_train,type = "response",s = model1$lambda.1se)
train_cor <- cor.test(y,pred_logage)
plot(y,pred_logage,col = traindat$ObsSex)

#predict and correlate test data
pred_logage_test <- predict(model2,newx = xtest,type = "response",s = model1$lambda.1se)
test_cor <- cor.test(ytest,pred_logage_test)
plot(ytest,pred_logage_test,col = testdat$ObsSex)

#uncertainty for train and test
train_MAE <- median(abs(exp(y) - exp(pred_logage[ ,1])))
test_MAE <- median(abs(exp(ytest) - exp(pred_logage_test[ ,1])))

##split clock by male and female (same train/test split)
#female clock
trainf <- traindat %>% filter(ObsSex == 2)
testf <- testdat %>% filter(ObsSex == 2)
xf <- as.matrix(trainf[,which(colnames(trainf) %in% bestCpGnames)])
yf <- trainf$LogAge

ytestf <- testf$LogAge
xtestf <- as.matrix(testf[,which(colnames(testf) %in% bestCpGnames)])

modelf <- cv.glmnet(xf,yf,family = "gaussian",alpha = 0,nfolds = 10)

#predict train and test correlation
pred_logage_f <- predict(modelf,newx = xf,type = "response",s = model1$lambda.1se)
train_cor_f <- cor.test(yf,pred_logage_f)
plot(yf,pred_logage_f)

#predict test and correlation
pred_logage_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.1se)
test_cor_f <- cor.test(ytestf,pred_logage_test_f)
plot(ytestf,pred_logage_test_f)

#uncertainty for train and test
train_MAE_f <- median(abs(exp(yf) - exp(pred_logage_f[ ,1])))
test_MAE_f <- median(abs(exp(ytestf) - exp(pred_logage_test_f[ ,1])))


#male clock
trainm <- traindat %>% filter(ObsSex == 1)
testm <- testdat %>% filter(ObsSex == 1)
xm <- as.matrix(trainm[,which(colnames(trainm) %in% bestCpGnames)])
ym <- trainm$LogAge

ytestm <- testm$LogAge
xtestm <- as.matrix(testm[,which(colnames(testm) %in% bestCpGnames)])

#run glmnet model with alpha=0 to include all sites
modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0,nfolds = 10)

#predict train and correlate
pred_logage_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.1se)
train_cor_m <- cor.test(ym,pred_logage_m)$estimate
plot(ym,pred_logage_m)

#predict test and correlate
pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.1se)
test_cor_m <- cor.test(ytestm,pred_logage_test_m)
plot(ytestm,pred_logage_test_m)

#uncertainty for clock
train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))

##save outputs
modelout[i,13] <- model1$lambda.1se
CpGlist[[i]] <- bestCpG_train

modelout[i,1] <- train_cor$estimate
modelout[i,2] <- test_cor$estimate

modelout[i,3] <- train_MAE
modelout[i,4] <- test_MAE

modelout[i,5] <- train_cor_f
modelout[i,6] <- test_cor_f

modelout[i,7] <- train_MAE_f
modelout[i,8] <- test_MAE_f

modelout[i,9] <- train_cor_m
modelout[i,10] <- test_cor_m

modelout[i,11] <- train_MAE_m
modelout[i,12] <- test_MAE_m

























