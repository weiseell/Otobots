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
library(caret)
library(glmnet)
library(MASS)

#subselect for female halibut
CPGs2 <- merge(sexdata, allCpGs1, by = "indivname")
CPGs2 <- subset(CPGs2,CPGs2$ObsSex == 2)
CPGs2 <- CPGs2[,-2:-3]

CpGsComb <- full_join(CPGs2,allCpGs1_HPEI)

#remove CpGs with NAs for the model
nas <- colSums(is.na(CpGsComb))
allCpGs_female <- CpGsComb[,which(nas < 1)]

#add age column
allCpGs_female <- merge(sexdata,allCpGs_female,by = "indivname")

allCpGs_female <- allCpGs_female %>% 
  mutate(LogAge=logLin(Age,5.8))

# Training data sets (70% of samples)
trainCpGs <- allCpGs_female %>% 
  group_by(Age) %>% 
  slice_sample(n = 4)

x.train <- trainCpGs[,-1:-3]

y.train <- trainCpGs$LogAge

#split data into CpG input and loglin predictor sets

# Test data sets (28 samples)
testCpGs <- allCpGs_female %>% 
  filter(!indivname %in% trainCpGs$indivname)

x.test <- as.matrix(testCpGs[,-1:-3])

y.test <- testCpGs$LogAge

### optimize with train function
#set train control variables
custom <- trainControl(method = "LOOCV",number = 5)


tuneGrid = expand.grid(alpha=seq(0,1,length=20),
                       lambda = seq(0.0001,1000,length=100))

#!# Run this line in Rgui:
#1. run this line: Rgui.exe --max-ppsize=500000
#in command prompt in the bin/x64 folder for R 4.3
#2. Run this line in the Rgui before running the top of the script:
#options(expressions = 5e5)

en <- train(LogAge~.,x.train,method = 'glmnet',tuneGrid = tuneGrid,
            trControl=custom)


x.train <- x.train[,-ncol(x.train)]

model1 <- glmnet(x=x.train,y=y.train,family = "gaussian",alpha = 0.05263158,lambda = 0.0001)

modelcoef <- as.matrix(coef(model1, s="lambda.1se"))

bestCpGs <- modelcoef[which(modelcoef > 0.0001),]
CpGnames <- names(bestCpGs)

#### running model with only best CpGs ####
#subset train and test for best CpGs
x.train.select <- x.train[,which(colnames(x.train) %in% CpGnames)]
x.test.select <- x.test[,which(colnames(x.test) %in% CpGnames)]

#optimize lambda and alpha again with new CpG set
x.train.select$LogAge <- trainCpGs$LogAge

en <- train(LogAge~.,x.train.select,method = 'glmnet',tuneGrid = tuneGrid,
            trControl=custom)

#fit optimized model
x.train.select <- x.train.select[,-ncol(x.train.select)]
model.fit <- glmnet(x = x.train.select,y = y.train, family = "gaussian", alpha = 0.5263158,lambda = 0.0001)

#run leave-one out to generate predicted values for train set
i <- 1
trainCpGs$MolClockLOO <- NA
for (i in 1:length(y.train)){
  epi.mat.loo <- x.train.select[-i,]
  ind.dat.loo <- y.train[-i]
  testSample <- as.matrix(x.train.select[i,])
  
  fit1 <-
    glmnet(
      x = epi.mat.loo,
      y = ind.dat.loo,
      lambda = 0.0001,
      alpha = 0.5263158,
      family = "gaussian")
  
  y.fitted <- predict(fit1,
                      newx = testSample,
                      type = "response")
  
  trainCpGs$MolClockLOO[i] <- y.fitted
  print(i)
}

#convert back to age from log age
trainCpGs <- trainCpGs %>% mutate(inverseAge = logLin.inverse(MolClockLOO, 5.8))

r.LOOCV <- cor(trainCpGs$Age, 
               trainCpGs$inverseAge, 
               method = "pearson")


# Plot results
ggplot(trainCpGs, aes(x = Age, y = inverseAge)) +
  ylim(0, 15) + xlim(0, 15) +
  geom_point(shape = 1, size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    fullrange = TRUE,
    colour = "gray30"
  ) +
  xlab("Known Age (years)") +
  ggtitle(paste("Training data; R = ", round(r.LOOCV, 4))) +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18)

#generate predicted values for test set
agepred <- data.frame(logAge=y.test,pred.logAge = predict(model.fit,newx = x.test.select,type = "response"))
agepred <- agepred %>% mutate(invAgeEst = logLin.inverse(s0, 5.8))
agepred <- agepred %>% mutate(invAgeEst = logLin.inverse(s0, 5.8))
agepred$Age <- testCpGs$Age

cor(agepred$Age,agepred$invAgeEst)

ggplot(agepred,aes(x = Age,
           y = invAgeEst)) +
  geom_point(
    shape = 1,
    size = 3) +
  ggtitle("Test data") +
  ylim(0, 25) + xlim(0, 25) +
  xlab("Known Age (years)") +
  ylab("Molecular clock Age (years)") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    fullrange = TRUE,
    color = "gray30"
  ) +
  theme_classic(base_size = 18)




