##Clock with only best CpG sites ####

## load model outputs from Age clock
load("ModelResult_Workspace_111623.RData")
sexdata <- read.csv("MethylWild_Sex_AgeData.csv")

sexdata$LogLen <- log(sexdata$Length)
###run clock with best CpG sites
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
x <- as.matrix(traindat[,-1:-5])

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% top_CpG_group$CpGlist_best)])
model1 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))

#predict testdat
ytest <- testdat$LogLen
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% top_CpG_group$CpGlist_best)])

#predict and correlate train data
pred_logage <- predict(model1,newx = bestCpG_train,type = "response",s = model1$lambda.min)
train_cor <- cor.test(y,pred_logage)
plot(y,pred_logage,col = traindat$ObsSex)

plot(model1, xvar = "dev", label = TRUE)

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
yf <- trainf$LogLen

ytestf <- testf$LogLen
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
ym <- trainm$LogLen

ytestm <- testm$LogLen
xtestm <- as.matrix(testm[,which(colnames(testm) %in% top_CpG_group$CpGlist_best)])

#run glmnet model with alpha=0 to include all sites
modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0,nfolds = 10)

#predict train and correlate
pred_logage_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.min)
train_cor_m <- cor.test(ym,pred_logage_m)
plot(ym,pred_logage_m)

#predict test and correlate
pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.min)
test_cor_m <- cor.test(ytestm,pred_logage_test_m)
plot(ytestm,pred_logage_test_m)

#uncertainty for clock
train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))

##Plots ####

#combining all three models in a plot
train_res <- data.frame(indivname=traindat$indivname,ObsSex=traindat$ObsSex,Age=traindat$Age,LogLen=y,LogEpiLen=pred_logage,type="Train")
test_res <- data.frame(indivname=testdat$indivname,ObsSex=testdat$ObsSex,Age=testdat$Age,LogLen=ytest,LogEpiLen=pred_logage_test,type="Test")
res_comb_all <- rbind(train_res,test_res)
res_comb_all <- res_comb_all %>% mutate(Length=exp(LogLen),EpiLen=exp(s1))
res_comb_all$ModelName <- "Combined"

train_res_f <- data.frame(indivname=trainf$indivname,ObsSex=trainf$ObsSex,Age=trainf$Age,LogLen=yf,LogEpiLen=pred_logage_f,type="Train")
test_res_f <- data.frame(indivname=testf$indivname,ObsSex=testf$ObsSex,Age=testf$Age,LogLen=ytestf,LogEpiLen=pred_logage_test_f,type="Test")
res_comb_all_f <- rbind(train_res_f,test_res_f)
res_comb_all_f <- res_comb_all_f %>% mutate(Length=exp(LogLen),EpiLen=exp(s1))
res_comb_all_f$ModelName <- "Female Only"

train_res_m <- data.frame(indivname=trainm$indivname,ObsSex=trainm$ObsSex,Age=trainm$Age,LogLen=ym,LogEpiLen=pred_logage_m,type="Train")
test_res_m <- data.frame(indivname=testm$indivname,ObsSex=testm$ObsSex,Age=testm$Age,LogLen=ytestm,LogEpiLen=pred_logage_test_m,type="Test")
res_comb_all_m <- rbind(train_res_m,test_res_m)
res_comb_all_m <- res_comb_all_m %>% mutate(Length=exp(LogLen),EpiLen=exp(s1))
res_comb_all_m$ModelName <- "Male Only"

res_merged <- rbind(res_comb_all,res_comb_all_f,res_comb_all_m)

ggplot(res_merged,aes(x = Length,y = EpiLen,color=Age,shape = type))+
  facet_grid(~ModelName) +
  geom_point(size = 3.5)+
  geom_jitter(height = 0,width = 0.1)+
  ylim(0, 175) + xlim(0, 175) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Measured Length (cm)") +
  ylab("Predicted Length (cm)") +
  ggtitle("Correlation Between Measure and Epigenetic Predicted Length") +
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









