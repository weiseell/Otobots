##Clock with only best CpG sites ####
###run clock with best CpG sites
CpGlist_best_comb <- read.csv(file = "Output/CpGList_LengthAge.csv")
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
y <- traindat$LogAge
x <- as.matrix(traindat[,-1:-5])

#running a glmnet with all selected sites
bestCpG_train <- as.matrix(traindat[,which(colnames(traindat) %in% CpGlist_best_comb$CpGlist_best)])
model1 <- cv.glmnet(bestCpG_train,y,family = "gaussian",alpha = 0.05,nfolds = 10)
coef.allCpG <- as.matrix(coef(model1,s = "lambda.min"))

bestCpGnames <- names(coef.allCpG[which(coef.allCpG > 0),])

#predict testdat
ytest <- testdat$LogAge
xtest <- as.matrix(testdat[,which(colnames(testdat) %in% CpGlist_best_comb$CpGlist_best)])

#predict and correlate train data
pred_logage <- predict(model1,newx = bestCpG_train,type = "response",s = model1$lambda.min)
train_cor <- cor.test(y,pred_logage)
#plot(y,pred_logage,col = traindat$ObsSex)

plot(model1, xvar = "dev", label = TRUE)



#predict and correlate test data
pred_logage_test <- predict(model1,newx = xtest,type = "response",s = model1$lambda.min)
test_cor <- cor.test(ytest,pred_logage_test)
#plot(ytest,pred_logage_test,col = testdat$ObsSex)

#uncertainty for train and test
train_MAE <- median(abs(exp(y) - exp(pred_logage[ ,1])))
test_MAE <- median(abs(exp(ytest) - exp(pred_logage_test[ ,1])))

##Sex-separated clocks
#female clock
trainf <- traindat %>% filter(ObsSex == 2)
testf <- testdat %>% filter(ObsSex == 2)
xf <- as.matrix(trainf[,which(colnames(trainf) %in% CpGlist_best_comb$CpGlist_best)])
yf <- trainf$LogAge

ytestf <- testf$LogAge
xtestf <- as.matrix(testf[,which(colnames(testf) %in% CpGlist_best_comb$CpGlist_best)])

modelf <- cv.glmnet(xf,yf,family = "gaussian",alpha = 0.05,nfolds = 10)

#predict train and test correlation
pred_logage_f <- predict(modelf,newx = xf,type = "response",s = model1$lambda.min)
train_cor_f <- cor.test(yf,pred_logage_f)
#plot(yf,pred_logage_f)

#predict test and correlation
pred_logage_test_f <- predict(modelf,newx = xtestf,type = "response",s = model1$lambda.min)
test_cor_f <- cor.test(ytestf,pred_logage_test_f)
#plot(ytestf,pred_logage_test_f)

#uncertainty for train and test
train_MAE_f <- median(abs(exp(yf) - exp(pred_logage_f[ ,1])))
test_MAE_f <- median(abs(exp(ytestf) - exp(pred_logage_test_f[ ,1])))


#male clock
trainm <- traindat %>% filter(ObsSex == 1)
testm <- testdat %>% filter(ObsSex == 1)
xm <- as.matrix(trainm[,which(colnames(trainm) %in% CpGlist_best_comb$CpGlist_best)])
ym <- trainm$LogAge

ytestm <- testm$LogAge
xtestm <- as.matrix(testm[,which(colnames(testm) %in% CpGlist_best_comb$CpGlist_best)])

#run glmnet model with alpha=0 to include all sites
modelm <- cv.glmnet(xm,ym,family = "gaussian",alpha = 0.05,nfolds = 10)

#predict train and correlate
pred_logage_m <- predict(modelm,newx = xm,type = "response",s = model1$lambda.min)
train_cor_m <- cor.test(ym,pred_logage_m)$estimate
#plot(ym,pred_logage_m)

#predict test and correlate
pred_logage_test_m <- predict(modelm,newx = xtestm,type = "response",s = model1$lambda.min)
test_cor_m <- cor.test(ytestm,pred_logage_test_m)
#plot(ytestm,pred_logage_test_m)

#uncertainty for clock
train_MAE_m <- median(abs(exp(ym) - exp(pred_logage_m[ ,1])))
test_MAE_m <- median(abs(exp(ytestm) - exp(pred_logage_test_m[ ,1])))

# save Rdata input for consistent plots/downstream analysis

##Plots ####
## Load manuscript clock results for consistency 

## make results into dataframes
#clock correlation plot
train_res <- data.frame(indivname=traindat$indivname,ObsSex=traindat$ObsSex,Length=traindat$Length,LogAge=y,LogEpiAge=pred_logage,type="Train")
test_res <- data.frame(indivname=testdat$indivname,ObsSex=testdat$ObsSex,Length=testdat$Length,LogAge=ytest,LogEpiAge=pred_logage_test,type="Test")
res_comb_all <- rbind(train_res,test_res)
res_comb_all <- res_comb_all %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
res_comb_all$ModelName <- "Combined"

train_res_f <- data.frame(indivname=trainf$indivname,ObsSex=trainf$ObsSex,Length=trainf$Length,LogAge=yf,LogEpiAge=pred_logage_f,type="Train")
test_res_f <- data.frame(indivname=testf$indivname,ObsSex=testf$ObsSex,Length=testf$Length,LogAge=ytestf,LogEpiAge=pred_logage_test_f,type="Test")
res_comb_all_f <- rbind(train_res_f,test_res_f)
res_comb_all_f <- res_comb_all_f %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
res_comb_all_f$ModelName <- "Female Only"

train_res_m <- data.frame(indivname=trainm$indivname,ObsSex=trainm$ObsSex,Length=trainm$Length,LogAge=ym,LogEpiAge=pred_logage_m,type="Train")
test_res_m <- data.frame(indivname=testm$indivname,ObsSex=testm$ObsSex,Length=testm$Length,LogAge=ytestm,LogEpiAge=pred_logage_test_m,type="Test")
res_comb_all_m <- rbind(train_res_m,test_res_m)
res_comb_all_m <- res_comb_all_m %>% mutate(Age=exp(LogAge),EpiAge=exp(s1))
res_comb_all_m$ModelName <- "Male Only"

res_merged <- rbind(res_comb_all,res_comb_all_f,res_comb_all_m)

tiff(filename = "Figures/AgingClock_Combined_AgeLengthCpGs.tiff",width = 300,height = 150,units = "mm",res = 400)
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
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")
dev.off()
#MAE values
res_merged$AbsErr <- abs(res_merged$Age - res_merged$EpiAge)

tiff(filename = "Figures/AgingClock_MAE.tiff",width = 300,height = 150,units = "mm",res = 400)
ggplot(res_merged,aes(x = type, y = AbsErr)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.character(ObsSex)), size = 3) +
  facet_wrap(~ModelName) +
  scale_color_manual(values = c(wes_palette("IsleofDogs1")[1],wes_palette("IsleofDogs1")[3])) +
  theme_classic(base_size = 18) + 
  labs(color = "Sex", x = element_blank(), y = "Absolute Error")
dev.off()

# Sorting clock by length percentile
plotdata <- res_merged[,c(-1,-4:-5)]

tiff(filename = "Figures/AgingClock_LengthPercentile.tiff",width = 300,height = 150,units = "mm",res = 400)
plotdata %>% 
  group_by(ModelName,Age) %>% 
  arrange(Length) %>% 
  mutate(rank=rank(Length)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Age,y = EpiAge,color=rank,shape = as.character(ObsSex)))+
  facet_grid(~ModelName) +
  geom_point(size = 3.5) +
  scale_color_gradient(low = "grey", high = "darkblue") +
  geom_jitter(height = 0,width = 0.1) +
  ylim(0, 17) + xlim(0, 17) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  ylab("Molecular clock Age (years)") +
  ggtitle("Correlation Between Otolith and Molecular Clock Age") +
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")
dev.off()

## Plot by length, age color ellipses
plotdata1 <- res_merged[,c(-1,-4,-8)]
tiff(filename = "Figures/AgingClock_LengthPlot_AgeEllipses.tiff",width = 200,height = 300,units = "mm",res = 400)
plotdata1 %>% 
  ggplot(aes(x = Length,y = s1,color=as.factor(Age),group=as.factor(Age)))+
  facet_wrap(~ModelName,nrow=3,scales = "free_x") +
  geom_point(size = 3.5) +
  geom_mark_ellipse()+
  geom_jitter(height = 0,width = 0.1) +
  ylim(1, 3)+
  xlab("Length (cm)") +
  ylab("Molecular clock Age (years)") +
  ggtitle("Correlation Between Otolith and Molecular Clock Age") +
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")
dev.off()


