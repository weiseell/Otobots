
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

modelout <- data.frame(matrix(data = NA,nrow = 100, ncol = 13))
colnames(modelout) <- c("Train_cor","Test_cor","Train_MAE","Test_MAE",
                        "F_Train_cor","F_Test_cor","F_Train_MAE","F_Test_MAE",
                        "M_Train_cor","M_Test_cor","M_Train_MAE","M_Test_MAE","lambda.1se")
CpGlist <- list()

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
  
  #predict and correlate test data
  pred_logage_test <- predict(model2,newx = xtest,type = "response",s = model1$lambda.1se)
  test_cor <- cor.test(ytest,pred_logage_test)
  
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
  
  #predict test and correlation
  pred_logage_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.1se)
  test_cor_f <- cor.test(ytestf,pred_logage_test_f)
  
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
  train_cor_m <- cor.test(ym,pred_logage_m)
  
  #predict test and correlate
  pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.1se)
  test_cor_m <- cor.test(ytestm,pred_logage_test_m)
  
  #uncertainty for clock
  train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
  test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))
  
  ##save outputs
  modelout[i,13] <- model1$lambda.1se
  CpGlist[[i]] <- bestCpGnames
  
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
}

hist(modelout$Train_cor)
hist(modelout$F_Train_cor)
hist(modelout$M_Train_cor)

#concatenate all the best CpGs
CpGlist_best <- CpGlist[[1]]
for (i in 2:100) {
  CpGlist_best <- c(CpGlist_best,CpGlist[[i]])
}
unique(CpGlist_best)

counting <- data.frame(CpGlist_best) %>% 
  group_by_all() %>% 
  count()
table(counting$n)

top_CpG_group <- counting %>% 
  arrange(desc(n)) %>% 
  filter(n > 20)

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
model1 <- cv.glmnet(x,y,family = "gaussian",alpha = 0.5,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.1se"))

model2 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0,nfolds = 10)

#predict testdat
ytest <- testdat$LogAge
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% top_CpG_group$CpGlist_best)])

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
xf <- as.matrix(trainf[,which(colnames(trainf) %in% top_CpG_group$CpGlist_best)])
yf <- trainf$LogAge

ytestf <- testf$LogAge
xtestf <- as.matrix(testf[,which(colnames(testf) %in% top_CpG_group$CpGlist_best)])

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
xm <- as.matrix(trainm[,which(colnames(trainm) %in% top_CpG_group$CpGlist_best)])
ym <- trainm$LogAge

ytestm <- testm$LogAge
xtestm <- as.matrix(testm[,which(colnames(testm) %in% top_CpG_group$CpGlist_best)])

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

##make plots of correlations for ICES
train_res <- data.frame(indivname=traindat$indivname,ObsSex=traindat$ObsSex,LogAge=y,LogEpiAge=pred_logage,type="Train")
test_res <- data.frame(indivname=testdat$indivname,ObsSex=testdat$ObsSex,LogAge=ytest,LogEpiAge=pred_logage_test,type="Test")
res_comb_all <- rbind(train_res,test_res)
ggplot(res_comb_all,aes(x = LogAge,y = s1,color=ObsSex,shape = type))+
  geom_point()+
  #ylim(0, 15) + xlim(0, 15) +
  geom_smooth(method = "lm",
              se = FALSE,
              fullrange = TRUE,
              colour = "gray30")+
  xlab("Otolith Age (years)") +
  #ggtitle(paste("Female Halibut - Training data; R = ", round(corr3$estimate, 4))) +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18)


train_res_f <- data.frame(indivname=trainf$indivname,ObsSex=trainf$ObsSex,LogAge=yf,LogEpiAge=pred_logage_f,type="Train")
test_res_f <- data.frame(indivname=testf$indivname,ObsSex=testf$ObsSex,LogAge=ytestf,LogEpiAge=pred_logage_test_f,type="Test")
res_comb_all_f <- rbind(train_res_f,test_res_f)
res_comb_all_f <- res_comb_all_f %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))

train_res_m <- data.frame(indivname=trainm$indivname,ObsSex=trainm$ObsSex,LogAge=ym,LogEpiAge=pred_logage_m,type="Train")
test_res_m <- data.frame(indivname=testm$indivname,ObsSex=testm$ObsSex,LogAge=ytestm,LogEpiAge=pred_logage_test_m,type="Test")
res_comb_all_m <- rbind(train_res_m,test_res_m)
res_comb_all_m <- res_comb_all_m %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))

res_all <- rbind(res_comb_all_f,res_comb_all_m)
ggplot(res_all,aes(x = Age,y = EpiAge,color=ObsSex,shape = type))+
  geom_point()+
  #geom_abline(slope = 1,intercept = 0)+
  ylim(0, 16) + xlim(0, 16) +
  xlab("Otolith Age (years)") +
  #ggtitle(paste("Female Halibut - Training data; R = ", round(corr3$estimate, 4))) +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18)


##plot %methylation vs age at each site
plot(y,bestCpG_train[,1],col = traindat$ObsSex)  

##try asymptotic transformation/non-linear model



