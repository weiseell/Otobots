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
library(wesanderson)

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

#generate adjusted age data based on catch month

#log age data
sexdata <- sexdata %>% 
  filter(DataType == "Wild") %>% 
  mutate(LogAge = log(Age))

sexdata_f <- sexdata %>% 
  filter(ObsSex == 2)

sexdata_m <- sexdata %>% 
  filter(ObsSex == 1)

#asymptotic age transformation
#sexdata_tmp1 <- sexdata %>% 
#  filter(DataType == "Wild") %>% 
#  filter(ObsSex == 1) %>% 
#  mutate(LogAge = logLin(x = Age, maturity = 8))

#sexdata_tmp2 <- sexdata %>% 
#  filter(DataType == "Wild") %>% 
#  filter(ObsSex == 2) %>% 
#  mutate(LogAge = logLin(x = Age, maturity = 10))

#sexdata <- rbind(sexdata_tmp1,sexdata_tmp2)

#create blank objects for storing outputs in loop
modelout <- data.frame(matrix(data = NA,nrow = 100, ncol = 13))
colnames(modelout) <- c("Train_cor","Test_cor","Train_MAE","Test_MAE",
                        "F_Train_cor","F_Test_cor","F_Train_MAE","F_Test_MAE",
                        "M_Train_cor","M_Test_cor","M_Train_MAE","M_Test_MAE","lambda.1se")
CpGlist <- list()
CpGlist_f <- list()
CpGlist_m <- list()

predout <- data.frame(indivname = sexdata$indivname,
                      LogAge = sexdata$LogAge)

predout_f <- data.frame(indivname = sexdata_f$indivname,
                      LogAge = sexdata_f$LogAge)

predout_m <- data.frame(indivname = sexdata_m$indivname,
                      LogAge = sexdata_m$LogAge)

##Loop for glmnet testing with all CpG sites ####
i <- 1
for (i in 1:100) {
  print(i)
  set.seed(Sys.time())
  
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
  coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))
  
  bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])
  
  #predict testdat
  ytest <- testdat$LogAge
  xtest <- as.matrix(testdat[,-1:-5])
  
  #predict and correlate train data
  pred_logage <- predict(model1,newx = x,type = "response",s = model1$lambda.min)
  train_cor <- cor.test(y,pred_logage)
  
  #predict and correlate test data
  pred_logage_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.min)
  test_cor <- cor.test(ytest,pred_logage_test)
  
  #uncertainty for train and test
  train_MAE <- median(abs(exp(y) - exp(pred_logage[ ,1])))
  test_MAE <- median(abs(exp(ytest) - exp(pred_logage_test[ ,1])))
  
  ##split clock by male and female (same train/test split)
  #female clock
  trainf <- traindat %>% filter(ObsSex == 2)
  testf <- testdat %>% filter(ObsSex == 2)
  xf <- as.matrix(trainf[,-1:-5])
  yf <- trainf$LogAge
  
  ytestf <- testf$LogAge
  xtestf <- as.matrix(testf[,-1:-5])
  
  modelf <- cv.glmnet(xf,yf,family = "gaussian",alpha = 0.5,nfolds = 10)
  
  coef.allCpG_f <- as.matrix(coef(modelf,s = "lambda.min"))
  
  bestCpGnames_f <- names(coef.allCpG_f[which(coef.allCpG_f > 0),])
  
  #predict train and test correlation
  pred_logage_f <- predict(modelf,newx = xf,type = "response",s = model1$lambda.min)
  train_cor_f <- cor.test(yf,pred_logage_f)
  
  #predict test and correlation
  pred_logage_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.min)
  test_cor_f <- cor.test(ytestf,pred_logage_test_f)
  
  #uncertainty for train and test
  train_MAE_f <- median(abs(exp(yf) - exp(pred_logage_f[ ,1])))
  test_MAE_f <- median(abs(exp(ytestf) - exp(pred_logage_test_f[ ,1])))
  
  
  #male clock
  trainm <- traindat %>% filter(ObsSex == 1)
  testm <- testdat %>% filter(ObsSex == 1)
  xm <- as.matrix(trainm[,-1:-5])
  ym <- trainm$LogAge
  
  ytestm <- testm$LogAge
  xtestm <- as.matrix(testm[,-1:-5])
  
  #run glmnet model with alpha=0 to include all sites
  modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0.5,nfolds = 10)
  
  coef.allCpG_m <- as.matrix(coef(modelm,s = "lambda.min"))
  
  bestCpGnames_m <- names(coef.allCpG_m[which(coef.allCpG_m > 0),])
  
  #predict train and correlate
  pred_logage_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.min)
  train_cor_m <- cor.test(ym,pred_logage_m)
  
  #predict test and correlate
  pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.min)
  test_cor_m <- cor.test(ytestm,pred_logage_test_m)
  
  #uncertainty for clock
  train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
  test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))
  
  ##test male/female clocks on opposite glmnet
  #female clock with male test set
  xm_all <- as.matrix(rbind(as.data.frame(xm),as.data.frame(xtestm)))
  ym_all <- c(ym,ytestm)
  #predict test and correlate
  pred_logage_m_all <- predict(modelf,newx = xm_all,type = "response",s = model1$lambda.min)
  test_cor_m_all <- cor.test(ym_all,pred_logage_m_all)
  #male clock with female test set
  xf_all <- as.matrix(rbind(as.data.frame(xf),as.data.frame(xtestf)))
  yf_all <- c(yf,ytestf)
  #predict test and correlate
  pred_logage_f_all <- predict(modelm,newx = xf_all,type = "response",s = model1$lambda.min)
  test_cor_f_all <- cor.test(yf_all,pred_logage_f_all)
  
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
                          LogAge = trainpick$LogAge,
                          predage = pred_logage)
  testpred <- data.frame(indivname = testpick$indivname,
                          LogAge = testpick$LogAge,
                          predage = pred_logage_test)
  predout_tmp <- rbind(trainpred,testpred)
  predout <- full_join(predout,predout_tmp, by = c("indivname","LogAge"))
  
  #female model
  trainpred_f <- data.frame(indivname = trainf$indivname,
                          LogAge = trainf$LogAge,
                          predage = pred_logage_f)
  testpred_f <- data.frame(indivname = testf$indivname,
                         LogAge = testf$LogAge,
                         predage = pred_logage_test_f)
  predout_tmp <- rbind(trainpred_f,testpred_f)
  predout_f <- full_join(predout_f,predout_tmp, by = c("indivname","LogAge"))
  
  #male model
  trainpred_m <- data.frame(indivname = trainm$indivname,
                          LogAge = trainm$LogAge,
                          predage = pred_logage_m)
  testpred_m <- data.frame(indivname = testm$indivname,
                         LogAge = testm$LogAge,
                         predage = pred_logage_test_f)
  predout_tmp <- rbind(trainpred_m,testpred_m)
  predout_m <- full_join(predout_m,predout_tmp, by = c("indivname","LogAge"))
  
}

hist(modelout$Train_cor)
hist(modelout$Test_cor)
hist(modelout$F_Train_cor)
hist(modelout$M_Train_cor)
hist(modelout$F_Test_cor)
hist(modelout$M_Test_cor)
#concatenate all the best CpGs
CpGlist_best <- CpGlist[[1]]
for (i in 2:100) {
  CpGlist_best <- c(CpGlist_best,CpGlist[[i]])
}

counting <- data.frame(CpGlist_best) %>% 
  group_by_all() %>% 
  count() %>% 
  arrange(desc(n))
table(counting$n)

top_CpG_group <- counting %>% 
  arrange(desc(n)) %>% 
  filter(n > 20)

#best Cpgs for female clock
CpGlist_best_f <- CpGlist_f[[1]]
for (i in 2:100) {
  CpGlist_best_f <- c(CpGlist_best_f,CpGlist[[i]])
}

counting <- data.frame(CpGlist_best_f) %>% 
  group_by_all() %>% 
  count() %>% 
  arrange(desc(n))
table(counting$n)

top_CpG_group_f <- counting %>% 
  arrange(desc(n)) %>% 
  filter(n > 20)

#best Cpgs for male clock
CpGlist_best_m <- CpGlist_m[[1]]
for (i in 2:100) {
  CpGlist_best_m <- c(CpGlist_best_m,CpGlist[[i]])
}

counting <- data.frame(CpGlist_best_m) %>% 
  group_by_all() %>% 
  count() %>% 
  arrange(desc(n))
table(counting$n)

top_CpG_group_m <- counting %>% 
  arrange(desc(n)) %>% 
  filter(n > 20)

top_CpG_group_f$CpGlist_best_f %in% top_CpG_group_m$CpGlist_best_m
top_CpG_group$CpGlist_best %in% top_CpG_group_m$CpGlist_best_m

#average predictions from looped model ####
reps <- predout[,-1:-2]
predout$AvePredLogAge <- rowMeans(reps)
cor.test(predout$LogAge,predout$AvePredLogAge)

plot(predout$LogAge,predout$AvePredLogAge)

predout$AvePredAge <- exp(predout$AvePredLogAge)

#repeat for M/F separated groups
reps <- predout_f[,-1:-2]
predout_f$AvePredLogAge <- rowMeans(reps)
cor.test(predout_f$LogAge,predout_f$AvePredLogAge)

plot(predout_f$LogAge,predout_f$AvePredLogAge)

predout_f$AvePredAge <- exp(predout_f$AvePredLogAge)

reps <- predout_m[,-1:-2]
predout_m$AvePredLogAge <- rowMeans(reps)
cor.test(predout_m$LogAge,predout_m$AvePredLogAge)

plot(predout_m$LogAge,predout_m$AvePredLogAge)

predout_m$AvePredAge <- exp(predout_m$AvePredLogAge)

#calculate MAE
median(abs(exp(predout$LogAge) - exp(predout$AvePredLogAge)),na.rm = T)
median(abs(exp(predout_f$LogAge) - exp(predout_f$AvePredLogAge)),na.rm = T)
median(abs(exp(predout_m$LogAge) - exp(predout_m$AvePredLogAge)),na.rm = T)
##make plots of correlations for ICES
predout1 <- predout %>% 
  dplyr::select(indivname,AvePredAge)
predout_all <- merge(predout1,sexdata)

predout_f1 <- predout_f %>% 
  dplyr::select(indivname,AvePredAge)
predout_f1 <- merge(predout_f1,sexdata_f)

predout_m1 <- predout_m %>% 
  dplyr::select(indivname,AvePredAge)
predout_m1 <- merge(predout_m1,sexdata_m)

ggplot(predout_f1,aes(x = Age,y = AvePredAge,color=ObsSex))+
  geom_point()+
  ylim(0, 15) + xlim(0, 15) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  #ggtitle(paste("Female Halibut - Training data; R = ", round(corr3$estimate, 4))) +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18)

##Clock with only best CpG sites ####
###run clock with best CpG sites
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

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% top_CpG_group$CpGlist_best)])
model1 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0.5,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))

#predict testdat
ytest <- testdat$LogAge
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% top_CpG_group$CpGlist_best)])

#predict and correlate train data
pred_logage <- predict(model1,newx = bestCpG_train,type = "response",s = model1$lambda.min)
train_cor <- cor.test(y,pred_logage)
plot(y,pred_logage,col = traindat$ObsSex)

#predict and correlate test data
pred_logage_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.min)
test_cor <- cor.test(ytest,pred_logage_test)
plot(ytest,pred_logage_test,col = testdat$ObsSex)

#uncertainty for train and test
train_MAE <- median(abs(exp(y) - exp(pred_logage[ ,1])))
test_MAE <- median(abs(exp(ytest) - exp(pred_logage_test[ ,1])))

##Sex-separated clocks
#female clock
trainf <- traindat %>% filter(ObsSex == 2)
testf <- testdat %>% filter(ObsSex == 2)
xf <- as.matrix(trainf[,which(colnames(trainf) %in% top_CpG_group$CpGlist_best)])
yf <- trainf$LogAge

ytestf <- testf$LogAge
xtestf <- as.matrix(testf[,which(colnames(testf) %in% top_CpG_group$CpGlist_best)])

modelf <- cv.glmnet(xf,yf,family = "gaussian",alpha = 0,nfolds = 10)

#predict train and test correlation
pred_logage_f <- predict(modelf,newx = xf,type = "response",s = model1$lambda.min)
train_cor_f <- cor.test(yf,pred_logage_f)
plot(yf,pred_logage_f)

#predict test and correlation
pred_logage_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.min)
test_cor_f <- cor.test(ytestf,pred_logage_test_f)
plot(ytestf,pred_logage_test_f)

#uncertainty for train and test
train_MAE_f <- median(abs(exp(yf) - exp(pred_logage_f[ ,1])))
test_MAE_f <- median(abs(exp(ytestf) - exp(pred_logage_test_f[ ,1])))


#male clock
trainm <- traindat %>% filter(ObsSex == 1)
testm <- testdat %>% filter(ObsSex == 1)
xm <- as.matrix(trainm[,which(colnames(trainm) %in% top_CpG_group$CpGlist_best)])
ym <- trainm$LogAge

ytestm <- testm$LogAge
xtestm <- as.matrix(testm[,which(colnames(testm) %in% top_CpG_group$CpGlist_best)])

#run glmnet model with alpha=0 to include all sites
modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0,nfolds = 10)

#predict train and correlate
pred_logage_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.min)
train_cor_m <- cor.test(ym,pred_logage_m)$estimate
plot(ym,pred_logage_m)

#predict test and correlate
pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.min)
test_cor_m <- cor.test(ytestm,pred_logage_test_m)
plot(ytestm,pred_logage_test_m)

#uncertainty for clock
train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))

##Plots ####
#clock correlation plot
train_res <- data.frame(indivname=traindat$indivname,ObsSex=traindat$ObsSex,LogAge=y,LogEpiAge=pred_logage,type="Train")
test_res <- data.frame(indivname=testdat$indivname,ObsSex=testdat$ObsSex,LogAge=ytest,LogEpiAge=pred_logage_test,type="Test")
res_comb_all <- rbind(train_res,test_res)
res_comb_all <- res_comb_all %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
ggplot(res_comb_all,aes(x = Age,y = EpiAge,color=ObsSex,shape = type))+
  geom_point()+
  geom_jitter(height = 0,width = 0.1)+
  ylim(0, 17) + xlim(0, 17) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  #ggtitle(paste("Female Halibut - Training data; R = ", round(corr3$estimate, 4))) +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18)

#comparison of error between otolith and molecular estimates
oto_df <- data.frame(OtoAge = c(4,6,8,10,12,15),
                     mean = c(4.18,6.25,8.04,9.90,11.67,14.67),
                     sd = c(0.40,0.97,0.96,1.22,1.56,2.16))

ggplot(oto_df,aes(x = OtoAge,y=mean))+
  geom_point()+
  geom_errorbar(aes(xmin = mean-sd, xmax = mean+sd))


res_summ <- res_comb_all %>% 
  group_by(Age) %>% 
  summarise(meanEpi = mean(EpiAge),
            sdEpi = sd(EpiAge)) %>% 
  rename(OtoAge = Age) %>% 
  mutate(meanOto = oto_df$mean, sdOto = oto_df$sd)

ggplot(res_summ, aes(x = meanOto, y = meanEpi, color = as.character(OtoAge))) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEpi - sdEpi, ymax = meanEpi + sdEpi)) +
  geom_errorbarh(aes(xmin = meanOto-sdOto, xmax = meanOto+sdOto))+
  ylim(0,17) + xlim(0,17) +
  geom_abline(slope = 1,intercept = 0) +
  theme_classic(base_size = 18) +
  xlab("Otolith Age (years)") +
  ylab("Molecular clock Age (years)") +
  ggtitle("Comparison of Uncertainty Between \n Otolith and Epigenetic Aging Methods")+
  labs(color = "Otolith Age")


#combining all three models in a plot
res_comb_all$ModelName <- "Combined"

train_res_f <- data.frame(indivname=trainf$indivname,ObsSex=trainf$ObsSex,LogAge=yf,LogEpiAge=pred_logage_f,type="Train")
test_res_f <- data.frame(indivname=testf$indivname,ObsSex=testf$ObsSex,LogAge=ytestf,LogEpiAge=pred_logage_test_f,type="Test")
res_comb_all_f <- rbind(train_res_f,test_res_f)
res_comb_all_f <- res_comb_all_f %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
res_comb_all_f$ModelName <- "Female Only"

train_res_m <- data.frame(indivname=trainm$indivname,ObsSex=trainm$ObsSex,LogAge=ym,LogEpiAge=pred_logage_m,type="Train")
test_res_m <- data.frame(indivname=testm$indivname,ObsSex=testm$ObsSex,LogAge=ytestm,LogEpiAge=pred_logage_test_m,type="Test")
res_comb_all_m <- rbind(train_res_m,test_res_m)
res_comb_all_m <- res_comb_all_m %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
res_comb_all_m$ModelName <- "Male Only"

res_merged <- rbind(res_comb_all,res_comb_all_f,res_comb_all_m)

ggplot(res_merged,aes(x = Age,y = EpiAge,color=as.character(ObsSex),shape = type))+
  facet_grid(~ModelName) +
  geom_point(size = 3.5)+
  geom_jitter(height = 0,width = 0.1)+
  ylim(0, 17) + xlim(0, 17) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  ylab("Molecular clock Age (years)") +
  ggtitle("Correlation Between Otolith and Molecular Clock Age") +
  scale_color_manual(values = c(wes_palette("IsleofDogs1")[1],wes_palette("IsleofDogs1")[3]))+
  theme_classic(base_size = 18)

#MAE values
res_merged$AbsErr <- abs(res_merged$Age - res_merged$EpiAge)

ggplot(res_merged,aes(x = type, y = AbsErr)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.character(ObsSex)), size = 3) +
  facet_wrap(~ModelName) +
  scale_color_manual(values = c(wes_palette("IsleofDogs1")[1],wes_palette("IsleofDogs1")[3])) +
  theme_classic(base_size = 18) + 
  labs(color = "Sex", x = element_blank(), y = "Absolute Error")









