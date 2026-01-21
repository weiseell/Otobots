---
title: "ModelComparisonOptions"
author: "Ellie Weise"
date: "2026-01-20"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Elastic Net (current option)

```{r,message=F,warning=FALSE}
#load libraries
library(tidyverse)
library(glmnet)
library(MASS)
library(wesanderson)
library(ggforce)
###run clock with best CpG sites
###run clock with best CpG sites
CpGlist_best_old <- read.csv(file = "Output/CpGList_LengthAge.csv")
CpGlist_best_old$Set <- "OldSites"
CpGlist_best_old <- CpGlist_best_old %>% dplyr::select(CpGlist_best,n,Set)
CpGlist_best_4to15 <- read.csv(file = "Output/CpGList_LengthAge_Age4to15_121025.csv")
CpGlist_best_4to15$Set <- "Age4to15"
CpGlist_best_4to29 <- read.csv(file = "Output/CpGList_LengthAge_Age4to29_121025.csv")
CpGlist_best_4to29$Set <- "Age4to29"
CpGlist_all <- rbind(CpGlist_best_4to15,CpGlist_best_4to29)
#CpGlist_all <- CpGlist_best_old
##load data
# if doing clock from scratch
load("MethylKitInput/allCpGs_MinimalMissing.rda")
sexdata <- read.csv("Input/MethylWild_Sex_AgeData.csv")

```

### Preparing data for elastic net
```{r}
## Data prep ####
set.seed(1230)
# subset to age 4-15 wild data
sexdata <- sexdata %>% 
  filter(DataType == "Wild" & Age < 20)

#running a glmnet with all selected sites
#CpGlist_all <- CpGlist_all %>% filter(ModelType == "Age")
bestCpGs <- allCpGs1[,which(colnames(allCpGs1) %in% CpGlist_all$CpGlist_best)]

#remove CpGs with NAs for the model
nas <- colSums(is.na(bestCpGs))
bestCpGs <- bestCpGs[,which(nas < 1)]
bestCpGs <- cbind(indivname=allCpGs1$indivname,bestCpGs)

#write.csv(allCpGs1,"allCpGs_nonas.csv",append = F,quote = F)

#!# remove sex-det and inversion!!!
#allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_047158.1",
#                                    colnames(allCpGs1),perl = T)))]
#allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_047162.1",
#                                    colnames(allCpGs1),perl = T)))]
# remove positive control from HPEI set
#allCpGs1 <- allCpGs1[-125,]

# check nas across individuals
# should be zero!!
#nas_row <- rowSums(is.na(allCpGs1))

#log age data
sexdata <- sexdata %>% 
  mutate(LogAge = log(Age),LogLen = log(Length))

sexdata_f <- sexdata %>% 
  filter(ObsSex == 2)

sexdata_m <- sexdata %>% 
  filter(ObsSex == 1)
```


### Separating into Training and Testing Data Sets 
```{r}
#randomly select individuals for testing
testpick <- sexdata %>% 
  group_by(Age,ObsSex) %>% 
  slice_sample(n = 1)

#sort other individuals into training data
trainpick <- sexdata[which(!(sexdata$indivname %in% testpick$indivname)),]

#select CpG sites for each group
testdat <- bestCpGs[which(bestCpGs$indivname %in% testpick$indivname),]
testdat <- left_join(testpick,testdat,by = "indivname")

traindat <- bestCpGs[which(!(bestCpGs$indivname %in% testpick$indivname)),]
traindat <- left_join(trainpick,traindat,by = "indivname")
```


### Run Elastic Net Model on Training Set and Predict Test Set
```{r,warning=F}
## Full Data Model ####
#select data and run glmnet to select sites
y <- traindat$LogAge

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% CpGlist_all$CpGlist_best)])

table(colSums(is.na(bestCpG_train)))
nas <- colSums(is.na(bestCpG_train))
bestCpG_train <- bestCpG_train[,which(nas < 1)]

model1 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0.05,nfolds = 20)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))

bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])

#predict testdat
ytest <- testdat$LogAge
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% CpGlist_all$CpGlist_best)])
nas <- colSums(is.na(bestCpG_train))
xtest <- xtest[,which(nas < 1)]

#predict and correlate train data
pred_logage <- predict(model1,newx = bestCpG_train,type = "response",s = model1$lambda.1se)
colnames(pred_logage) <- "Epi_LogAge"
train_cor <- cor.test(y,pred_logage)
plot(y,pred_logage,col = traindat$ObsSex)

plot(model1, xvar = "dev", label = TRUE)

#predict and correlate test data
pred_logage_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.1se)
colnames(pred_logage_test) <- "Epi_LogAge"
test_cor <- cor.test(ytest,pred_logage_test)
plot(ytest,pred_logage_test,col = testdat$ObsSex)

#uncertainty for train and test
train_MAE <- median(abs(exp(y) - exp(pred_logage[ ,1])))
test_MAE <- median(abs(exp(ytest) - exp(pred_logage_test[ ,1])))

# plot prediction
#clock correlation plot
train_res <- data.frame(indivname=traindat$indivname,
                        ObsSex=traindat$ObsSex,
                        Length=traindat$Length,
                        LogAge=y,
                        Epi_LogAge=pred_logage,
                        type="Train",
                        sampletype=traindat$DataType)
test_res <- data.frame(indivname=testdat$indivname,
                       ObsSex=testdat$ObsSex,
                       Length=testdat$Length,
                       LogAge=ytest,
                       Epi_LogAge=pred_logage_test,
                       type="Test",
                       sampletype=testdat$DataType)
res_comb_all <- rbind(train_res,test_res)
res_comb_all <- res_comb_all %>% mutate(Age=exp(LogAge),EpiAge=exp(Epi_LogAge))
res_comb_all$ModelName <- "Combined"
ggplot(res_comb_all,aes(x = Age,y = EpiAge,color=as.character(ObsSex),shape = type))+
  geom_point(size = 3.5)+
  geom_jitter(height = 0,width = 0.1)+
  ylim(0, 30) + xlim(0, 30) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  ylab("Molecular clock Age (years)") +
  ggtitle("Correlation Between Otolith and Molecular Clock Age") +
  scale_color_manual(values = c(wes_palette("IsleofDogs1")[1],wes_palette("IsleofDogs1")[3]))+
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")
```

I've also tried this option with male and female separated model options, but since we're comparing overall methods I'm going to try these other options with the full data set first and then try it with male and female separated clocks once I have a better handle on the models. I also removed the hatchery data for now, it should help with the lack of 20+ age data, but they are all female and the age 23s look a little different from the wild data points around the same age, and I don't know if that's due to environmental difference or otolith error. Either way seems like a 'next problem' sort of thing.


