## Length CpG Selection

##load data
# if doing clock from scratch
load("Input/allCpGs_nosingletons.rda")
sexdata <- read.csv("Input/MethylWild_Sex_AgeData.csv")

#load libraries
library(tidyverse)
library(glmnet)
library(MASS)
library(wesanderson)

## Data prep ####
#remove CpGs with NAs for the model
nas <- colSums(is.na(allCpGs1))
allCpGs1 <- allCpGs1[,which(nas < 1)]
#write.csv(allCpGs1,"allCpGs_nonas.csv",append = F,quote = F)

#remove NW and mito data from data set
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NW",
                                    colnames(allCpGs1),perl = T)))]
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_009709.1",
                                    colnames(allCpGs1),perl = T)))]
#!# remove sex-det and inversion!!!
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_047158.1",
                                    colnames(allCpGs1),perl = T)))]
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_047162.1",
                                    colnames(allCpGs1),perl = T)))]

#log age data
sexdata <- sexdata %>% 
  mutate(LogAge = log(Age),LogLen = log(Length))

sexdata_f <- sexdata %>% 
  filter(ObsSex == 2)

sexdata_m <- sexdata %>% 
  filter(ObsSex == 1)

#create blank objects for storing outputs in loop
modelout <- data.frame(matrix(data = NA,nrow = 100, ncol = 13))
colnames(modelout) <- c("Train_cor","Test_cor","Train_MAE","Test_MAE",
                        "F_Train_cor","F_Test_cor","F_Train_MAE","F_Test_MAE",
                        "M_Train_cor","M_Test_cor","M_Train_MAE","M_Test_MAE","lambda.1se")
CpGlist <- list()
CpGlist_f <- list()
CpGlist_m <- list()

predout <- data.frame(indivname = sexdata$indivname,
                      LogLen = sexdata$LogLen)

predout_f <- data.frame(indivname = sexdata_f$indivname,
                        LogLen = sexdata_f$LogLen)

predout_m <- data.frame(indivname = sexdata_m$indivname,
                        LogLen = sexdata_m$LogLen)

##Loop for glmnet testing with all CpG sites ####
i <- 1
for (i in 1:100) {
  print(i)
  set.seed(Sys.time())
  
  #randomly select individuals for testing
  testpick <- sexdata %>% 
    group_by(Age,ObsSex) %>% 
    slice_sample(n = 1)
  
  #sort other individuals into training data
  trainpick <- sexdata[which(!(sexdata$indivname %in% testpick$indivname)),]
  
  #select CpG sites for each group
  testdat <- allCpGs1[which(allCpGs1$indivname %in% testpick$indivname),]
  testdat <- full_join(testpick,testdat,by = "indivname")
  
  traindat <- allCpGs1[which(!(allCpGs1$indivname %in% testpick$indivname)),]
  traindat <- full_join(trainpick,traindat,by = "indivname")
  
  #select data and run glmnet to select sites
  y <- traindat$LogLen
  x <- as.matrix(traindat[,-1:-6])
  
  model1 <- cv.glmnet(x,y,family = "gaussian",alpha = 0.5,nfolds = 10)
  coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))
  
  bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])
  
  #predict testdat
  ytest <- testdat$LogLen
  xtest <- as.matrix(testdat[,-1:-6])
  
  #predict and correlate train data
  pred_LogLen <- predict(model1,newx = x,type = "response",s = model1$lambda.min)
  train_cor <- cor.test(y,pred_LogLen)
  
  #predict and correlate test data
  pred_LogLen_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.min)
  test_cor <- cor.test(ytest,pred_LogLen_test)
  
  #uncertainty for train and test
  train_MAE <- median(abs(exp(y) - exp(pred_LogLen[ ,1])))
  test_MAE <- median(abs(exp(ytest) - exp(pred_LogLen_test[ ,1])))
  
  ##split clock by male and female (same train/test split)
  #female clock
  trainf <- traindat %>% filter(ObsSex == 2)
  testf <- testdat %>% filter(ObsSex == 2)
  xf <- as.matrix(trainf[,-1:-6])
  yf <- trainf$LogLen
  
  ytestf <- testf$LogLen
  xtestf <- as.matrix(testf[,-1:-6])
  
  modelf <- cv.glmnet(xf,yf,family = "gaussian",alpha = 0.5,nfolds = 10)
  
  coef.allCpG_f <- as.matrix(coef(modelf,s = "lambda.min"))
  
  bestCpGnames_f <- names(coef.allCpG_f[which(coef.allCpG_f > 0),])
  
  #predict train and test correlation
  pred_LogLen_f <- predict(modelf,newx = xf,type = "response",s = model1$lambda.min)
  train_cor_f <- cor.test(yf,pred_LogLen_f)
  
  #predict test and correlation
  pred_LogLen_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.min)
  test_cor_f <- cor.test(ytestf,pred_LogLen_test_f)
  
  #uncertainty for train and test
  train_MAE_f <- median(abs(exp(yf) - exp(pred_LogLen_f[ ,1])))
  test_MAE_f <- median(abs(exp(ytestf) - exp(pred_LogLen_test_f[ ,1])))
  
  
  #male clock
  trainm <- traindat %>% filter(ObsSex == 1)
  testm <- testdat %>% filter(ObsSex == 1)
  xm <- as.matrix(trainm[,-1:-6])
  ym <- trainm$LogLen
  
  ytestm <- testm$LogLen
  xtestm <- as.matrix(testm[,-1:-6])
  
  #run glmnet model with alpha=0 to include all sites
  modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0.5,nfolds = 10)
  
  coef.allCpG_m <- as.matrix(coef(modelm,s = "lambda.min"))
  
  bestCpGnames_m <- names(coef.allCpG_m[which(coef.allCpG_m > 0),])
  
  #predict train and correlate
  pred_LogLen_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.min)
  train_cor_m <- cor.test(ym,pred_LogLen_m)
  
  #predict test and correlate
  pred_LogLen_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.min)
  test_cor_m <- cor.test(ytestm,pred_LogLen_test_m)
  
  #uncertainty for clock
  train_MAE_m <- median(abs(exp(ym) - exp(pred_LogLen_m[ ,1])))
  test_MAE_m <- median(abs(exp(ytestm) - exp(pred_LogLen_test_m[ ,1])))
  
  ##test male/female clocks on opposite glmnet
  #female clock with male test set
  xm_all <- as.matrix(rbind(as.data.frame(xm),as.data.frame(xtestm)))
  ym_all <- c(ym,ytestm)
  #predict test and correlate
  pred_LogLen_m_all <- predict(modelf,newx = xm_all,type = "response",s = model1$lambda.min)
  test_cor_m_all <- cor.test(ym_all,pred_LogLen_m_all)
  #male clock with female test set
  xf_all <- as.matrix(rbind(as.data.frame(xf),as.data.frame(xtestf)))
  yf_all <- c(yf,ytestf)
  #predict test and correlate
  pred_LogLen_f_all <- predict(modelm,newx = xf_all,type = "response",s = model1$lambda.min)
  test_cor_f_all <- cor.test(yf_all,pred_LogLen_f_all)
  
  ##save outputs
  modelout[i,13] <- model1$lambda.min
  CpGlist[[i]] <- bestCpGnames
  CpGlist_f[[i]] <- bestCpGnames_f
  CpGlist_m[[i]] <- bestCpGnames_m
  
  modelout[i,1] <- train_cor$estimate
  modelout[i,2] <- test_cor$estimate
  
  modelout[i,3] <- train_MAE
  modelout[i,4] <- test_MAE
  
  modelout[i,5] <- train_cor_f$estimate
  modelout[i,6] <- test_cor_f$estimate
  
  modelout[i,7] <- train_MAE_f
  modelout[i,8] <- test_MAE_f
  
  modelout[i,9] <- train_cor_m$estimate
  modelout[i,10] <- test_cor_m$estimate
  
  modelout[i,11] <- train_MAE_m
  modelout[i,12] <- test_MAE_m
  
  #save predictions
  #full model
  trainpred <- data.frame(indivname = trainpick$indivname,
                          LogLen = trainpick$LogLen,
                          predlen = pred_LogLen)
  testpred <- data.frame(indivname = testpick$indivname,
                         LogLen = testpick$LogLen,
                         predlen = pred_LogLen_test)
  predout_tmp <- rbind(trainpred,testpred)
  predout <- full_join(predout,predout_tmp, by = c("indivname","LogLen"))
  
  #female model
  trainpred_f <- data.frame(indivname = trainf$indivname,
                            LogLen = trainf$LogLen,
                            predlen = pred_LogLen_f)
  testpred_f <- data.frame(indivname = testf$indivname,
                           LogLen = testf$LogLen,
                           predlen = pred_LogLen_test_f)
  predout_tmp <- rbind(trainpred_f,testpred_f)
  predout_f <- full_join(predout_f,predout_tmp, by = c("indivname","LogLen"))
  
  #male model
  trainpred_m <- data.frame(indivname = trainm$indivname,
                            LogLen = trainm$LogLen,
                            predlen = pred_LogLen_m)
  testpred_m <- data.frame(indivname = testm$indivname,
                           LogLen = testm$LogLen,
                           predlen = pred_LogLen_test_f)
  predout_tmp <- rbind(trainpred_m,testpred_m)
  predout_m <- full_join(predout_m,predout_tmp, by = c("indivname","LogLen"))
  
}

hist(modelout$Train_cor)
hist(modelout$Test_cor)
hist(modelout$F_Train_cor)
hist(modelout$M_Train_cor)
hist(modelout$F_Test_cor)
hist(modelout$M_Test_cor)

#concatenate all the best CpGs
CpGlist_best_len <- CpGlist[[1]]
for (i in 2:100) {
  CpGlist_best_len <- c(CpGlist_best_len,CpGlist[[i]])
}

counting <- data.frame(CpGlist_best_len) %>% 
  group_by_all() %>% 
  count() %>% 
  arrange(desc(n))
table(counting$n)

top_CpG_group_len <- counting %>% 
  arrange(desc(n)) %>% 
  filter(n > 20) %>% 
  rename(CpGlist_best = CpGlist_best_len)

CpGlist_best_comb <- unique(rbind(top_CpG_group,top_CpG_group_len))

write.csv(CpGlist_best_comb,file = "Output/CpGList_LengthAge.csv",row.names = F)
