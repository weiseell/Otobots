### Best CpG correlation and colinearity testing
#load libraries
library(tidyverse)
library(olsrr)

##load data
load("MethylKitInput/allCpGs_MinimalMissing.rda")
sexdata <- read.csv("Input/MethylWild_Sex_AgeData.csv")
CpGRegion <- read.csv(file = "Output/BestCpG_RegionGroups.csv")

###run clock with best CpG sites
CpGlist_best_old <- read.csv(file = "Output/CpGList_LengthAge.csv")
CpGlist_best_old$Set <- "OldSites"
CpGlist_best_old <- CpGlist_best_old %>% dplyr::select(CpGlist_best,n,Set)
CpGlist_best_4to15 <- read.csv(file = "Output/CpGList_LengthAge_Age4to15_121025.csv")
CpGlist_best_4to15$Set <- "Age4to15"
CpGlist_best_4to19 <- read.csv(file = "Output/CpGList_LengthAge_Age4to19_121025.csv")
CpGlist_best_4to19$Set <- "Age4to19"
CpGlist_best_4to29 <- read.csv(file = "Output/CpGList_LengthAge_Age4to29_121025.csv")
CpGlist_best_4to29$Set <- "Age4to29"
CpGlist_all <- rbind(CpGlist_best_4to15,CpGlist_best_4to29,CpGlist_best_old)

# plot missing data levels and subset to minimize missing data
nas <- colSums(is.na(allCpGs1))
hist(nas,breaks = 80)
#!# for this data set, looks like low levels of increasing missing data
#!# no jumps at numbers where it would be data-set based, as we're combining
#!# three different sequencing runs

## grab all potentially important CpGs
CpGlist_all <- CpGlist_all %>% 
  mutate(tmp=CpGlist_best) %>% 
  separate_wider_delim(tmp,delim = ".",names = c("CHROM","ones","MapInfo"),too_few = "align_start") %>% 
  unite(col=CHROM,CHROM,ones,sep = ".")

CpGlist <- data.frame(CpG=unique(CpGlist_all$CpGlist_best))

# subset to just wild data  for testing
# goal is to prevent low sample size and environmental bias
# for initial model, just to see if it works before adding more problems
RF_CpGs <- sexdata %>% 
  filter(DataType == "Wild") %>% 
  # subselect a training data set
  group_by(Age,ObsSex)

tmpCpGs <- allCpGs1[,which(colnames(allCpGs1) %in% CpGlist$CpG)]
CpGsbest_allsets <- data.frame(indivname=allCpGs1$indivname,tmpCpGs)

testdf <- left_join(RF_CpGs,CpGsbest_allsets)

# run cor.test for each CpG with age and with length
cor_age4to29 <- testdf %>% 
  pivot_longer(cols = starts_with("NC"),
               names_to = "CpG",
               values_to = "PercMethyl") %>% 
  group_by(CpG) %>% 
  summarise(CorAge4to29=cor.test(x=Age,y=PercMethyl)[["estimate"]],
            CorAge4to29_LCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][1],
            CorAge4to29_HCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][2],
            CorLength4to29=cor.test(x=Length,y=PercMethyl)[["estimate"]],
            CorLength4to29_LCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][1],
            CorLength4to29_HCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][2])

cor_age4to19 <- testdf %>% 
  pivot_longer(cols = starts_with("NC"),
               names_to = "CpG",
               values_to = "PercMethyl") %>% 
  filter(Age < 20) %>% 
  group_by(CpG) %>% 
  summarise(CorAge4to19=cor.test(x=Age,y=PercMethyl)[["estimate"]],
            CorAge4to19_LCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][1],
            CorAge4to19_HCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][2],
            CorLength4to19=cor.test(x=Length,y=PercMethyl)[["estimate"]],
            CorLength4to19_LCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][1],
            CorLength4to19_HCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][2])

cor_age4to15 <- testdf %>% 
  pivot_longer(cols = starts_with("NC"),
               names_to = "CpG",
               values_to = "PercMethyl") %>% 
  filter(Age < 16) %>% 
  group_by(CpG) %>% 
  summarise(CorAge4to15=cor.test(x=Age,y=PercMethyl)[["estimate"]],
            CorAge4to15_LCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][1],
            CorAge4to15_HCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][2],
            CorLength4to15=cor.test(x=Length,y=PercMethyl)[["estimate"]],
            CorLength4to15_LCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][1],
            CorLength4to15_HCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][2])

cor_age4to12 <- testdf %>% 
  pivot_longer(cols = starts_with("NC"),
               names_to = "CpG",
               values_to = "PercMethyl") %>% 
  filter(Age < 13) %>% 
  group_by(CpG) %>% 
  summarise(CorAge4to12=cor.test(x=Age,y=PercMethyl)[["estimate"]],
            CorAge4to12_LCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][1],
            CorAge4to12_HCI=cor.test(x=Age,y=PercMethyl)[["conf.int"]][2],
            CorLength4to12=cor.test(x=Length,y=PercMethyl)[["estimate"]],
            CorLength4to12_LCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][1],
            CorLength4to12_HCI=cor.test(x=Length,y=PercMethyl)[["conf.int"]][2])

#!# correlations go up as the older individuals are excluded
#!# improvements happen as 20+, 16+ and 15+ individuals are removed
#!# is this otolith error increasing or methyl patterns getting messier? Probably both

## use age 4-12 as a test set and grab CpGs with decent linear correlation
cor_age4to19 <- testdf

cor_age4to12$CpG[which(cor_age4to12$CorAge4to12 > 0.5)]

testdf_cor0.5 <- testdf[,colnames(testdf) %in% cor_age4to29$CpG[which(cor_age4to29$CorAge4to29 > 0.6)]]
#!# upping to 0.6 cor minimum eliminated the perfectly linear covariates
testdf_cor0.5 <- data.frame(Age=testdf$Age,testdf_cor0.5)
#testdf_cor0.5 <- testdf_cor0.5 %>% filter(Age < 20)

## run a basic linear model
model <- lm(Age~.,data = testdf_cor0.5)

# looking for perfectly colinear variables
ld.vars <- attributes(alias(model)$Complete)$dimnames[[1]]

reduceddf <- testdf_cor0.5[,which(!(colnames(testdf_cor0.5) %in% ld.vars))]

# rerun model without perfect colinear vars if needed
model <- lm(Age~.,data = reduceddf)

## get VIF colinearity tests
ols_vif_tol(model) %>% 
  arrange(desc(VIF))

# with cor < 0.6 eliminated, colinearity is below 10?
# not what I was expecting but I will take it

preds <- predict.lm(model)
plot(testdf_cor0.5$Age,preds)
cor.test(testdf_cor0.5$Age,preds)
## Age 4-12 has cor > 0.99 and good predictions
## Age 4-29 has cor~0.85 and underestimates the 20+ group but that 
## could very well be a sample size artifact

## plot CpG relationships and predictions
reduceddf %>% 
  pivot_longer(cols = starts_with("NC"),
               names_to = "CpG",
               values_to = "PercMethyl") %>% 
  ggplot(aes(x=Age,y=PercMethyl))+
    geom_point()+
    facet_wrap(~CpG,scales = "free_y")

data.frame(OtoAge=testdf_cor0.5$Age,EpiAge=preds) %>% 
  ggplot(aes(x = OtoAge,y = EpiAge))+
  geom_point(size = 3.5)+
  geom_jitter(height = 0,width = 0.1)+
  ylim(0, 30) + xlim(0, 30) +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Otolith Age (years)") +
  ylab("Molecular clock Age (years)") +
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

MAE = median(abs(preds - testdf_cor0.5$Age))

