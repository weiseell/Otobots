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
CpGlist_all <- rbind(CpGlist_best_4to15,CpGlist_best_4to29)

# plot missing data levels and subset to minimize missing data
nas <- colSums(is.na(allCpGs1))
hist(nas,breaks = 80)
#!# for this data set, looks like low levels of increasing missing data
#!# no jumps at numbers where it would be data-set based, as we're combining
#!# three different sequencing runs

## grab all potentially important CpGs
CpGlist_all <- CpGlist_all %>% 
  mutate(tmp=CpGlist_best) %>% 
  filter(CpGlist_best != "(Intercept)") %>% 
  separate_wider_delim(tmp,delim = ".",names = c("CHROM","ones","MapInfo"),too_few = "align_start") %>% 
  unite(col=CHROM,CHROM,ones,sep = ".")

CpGlist <- data.frame(CpG=unique(CpGlist_all$CpGlist_best))

bestCpGs <- allCpGs1[,which(colnames(allCpGs1) %in% CpGlist_all$CpGlist_best)]
CpGs_covars <- left_join(sexdata,bestCpGs)
CpGs_covars <- CpGs_covars %>% filter(DataType == "Wild")
## plot CpGs against age and length for each chromosome, print to pdf, to see all the CpGs big enough
## and to make it easier to compare regions
pdf(file = "CpGPlotsByChrom_AgeLength.pdf",width = 14,height = 10)
for (i in 1:length(unique(CpGlist_all$CHROM))) {
  CHROMtmp <- unique(CpGlist_all$CHROM)[i]
  # select covars and CpGs from chromosome
  plot1 <- CpGs_covars %>% 
    dplyr::select(indivname,Age,ObsSex,Length,starts_with(CHROMtmp)) %>% 
    pivot_longer(cols = !c(indivname,Age,ObsSex,Length),names_to = "CpG",values_to = "PercMethyl") %>% 
    ggplot(aes(x = Age,y = PercMethyl,colour = ObsSex))+
      geom_point()+
    facet_wrap(~CpG,scales = "free_y")+
    ggtitle(paste0("Otolith Age vs. Percent Methylation - ",CHROMtmp))
  
  plot2 <- CpGs_covars %>% 
    dplyr::select(indivname,Age,ObsSex,Length,starts_with(CHROMtmp)) %>% 
    pivot_longer(cols = !c(indivname,Age,ObsSex,Length),names_to = "CpG",values_to = "PercMethyl") %>% 
    ggplot(aes(x = Length,y = PercMethyl,colour = ObsSex))+
    geom_point()+
    facet_wrap(~CpG,scales = "free")+
    ggtitle(paste0("Length vs. Percent Methylation - ",CHROMtmp))
  print(plot1)
  print(plot2)
}

dev.off()


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

cor_age4to12$CpG[which(cor_age4to12$CorAge4to12 > 0.6)]

testdf_cor0.5 <- testdf[,colnames(testdf) %in% cor_age4to29$CpG[which(cor_age4to29$CorAge4to29 > 0.5)]]
#!# upping to 0.6 cor minimum eliminated the perfectly linear covariates
testdf_cor0.5 <- data.frame(Age=testdf$Age,testdf_cor0.5)
testdf_cor0.5 <- testdf_cor0.5 %>% filter(Age < 20)

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
# eliminate colinar factors one by one until they are all less than 10
reduceddf <- reduceddf[,which(colnames(reduceddf)!="NC_047172.1.26308364")]

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

# try an additive model 
# try a single CpG as a test
library(mgcv)
gam_model <- gam(Age ~ s(NC_047159.1.29170063),data = reduceddf)
summary(gam_model)
plot(reduceddf$Age,reduceddf$NC_047159.1.29170063)

# try the reduceddf set to compare directly to linear model
gam_model1 <- gam(Age ~ s(NC_047151.1.10785,bs = "cr") +
                    s(NC_047151.1.11792,bs = "cr") +
                    s(NC_047151.1.13802,bs = "cr") +
                    s(NC_047151.1.9880,bs = "cr") +
                    s(NC_047152.1.5986,bs = "cr") +
                    s(NC_047153.1.32669800,bs = "cr") +
                    s(NC_047155.1.30745806,bs = "cr") +
                    s(NC_047155.1.30745807,bs = "cr") +
                    s(NC_047155.1.30747359,bs = "cr") +
                    s(NC_047155.1.30748169,bs = "cr") +
                    s(NC_047155.1.30750614,bs = "cr") +
                    s(NC_047155.1.30751619,bs = "cr") +
                    s(NC_047159.1.29170063,bs = "cr") +
                    s(NC_047169.1.16885311,bs = "cr") +
                    s(NC_047172.1.26293714,bs = "cr") +
                    s(NC_047172.1.26293982,bs = "cr") +
                    s(NC_047172.1.26294006,bs = "cr") +
                    s(NC_047172.1.26303009,bs = "cr") +
                    s(NC_047172.1.26305360,bs = "cr") +
                    s(NC_047172.1.26307361,bs = "cr") +
                    s(NC_047174.1.21431948,bs = "cr"),data = reduceddf)

summary(gam_model1)
gam.check(gam_model1, k.rep = 1000)
plot.gam(gam_model1)

# these predictions look WAY too nice, bias from overfit probably lol
predictions = predict(gam_model1,type = 'response',se = TRUE)

df_preds = data.frame(reduceddf$Age, predictions) %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit)

ggplot(aes(x = reduceddf.Age, y = fit), data = df_preds) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray92') +
  geom_line(color = '#56B4E9') +
  geom_point(size = 3.5)

# try a gam with one per genome section represented
#!# this is based on the age vs. otolith age graphs not any automatic selection
# set options as generic 's' to see what the fit types should be
gam_model2 <- gam(Age ~ s(NC_047151.1.9880) +
                    s(NC_047163.1.29069872) +
                    s(NC_047153.1.32666785) +
                    s(NC_047155.1.30745806) +
                    s(NC_047159.1.29166056) +
                    s(NC_047164.1.2171813) +
                    s(NC_047174.1.21430655) +
                    s(NC_047172.1.26303009) +
                    s(NC_047165.1.15408994) +
                    s(NC_047165.1.910853) +
                    s(NC_047167.1.11015458) +
                    s(NC_047170.1.1470) +
                    s(NC_047171.1.15408463) +
                    s(NC_047171.1.25040373) +
                    s(NC_047152.1.5986) +
                    s(NC_047168.1.480431) +
                    s(NC_047157.1.604941) +
                    s(NC_047161.1.22201614) +
                    s(NC_047169.1.16885311),data = testdf, method = "REML")

# look at linearity and smooth factors
basic_summary <- summary(gam_model2)
# k values are below edf so I think that's okay
gam.check(gam_model2)

# look at edf values to see linear vs. nonlinear
basic_summary$s.table

par(mfrow = c(1, 1))
plot(gam_model2, all.terms = TRUE)


### try a gam with one per 'smooth type' within the best CpG set

# look at coefficients within 'best set' of elastic net, which ones didn't shrink?
coef_impor <- names(coef.allCpG[which(coef.allCpG > 0),])
CpG_impor <- testdf[,which(colnames(testdf) %in% coef_impor)]

## try a simple AIC add-in model
# try each factor to see if it improves the model
# if not skip
for (i in vector) {
  
}

## other ways to shrink these variables to something reasonable for a gam?
## AIC table?


## Trying to mess with age variation
#!# this reference collection defines agreement as within two years
#!# for some context to make myself feel better
refages <- read.csv(file = "Input/age_comp_ref_coll_June21_2022_Rinput.csv")
means_ref <- refages %>% 
  group_by(RefAge) %>% 
  summarise(totN=sum(N),
            meanAge=sum(TestAge*N)/totN,
            stddevAge=sqrt(sum((TestAge-meanAge)^2)/totN)) %>% 
  filter(RefAge >= 4 & RefAge < 20) %>% 
  mutate(lower = meanAge - 1.96 * stddevAge,
         upper = meanAge + 1.96 * stddevAge)

ggplot(aes(x = RefAge, y = meanAge), data = means_ref) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray92') +
  geom_line(color = '#56B4E9') +
  geom_point(size = 2) +
  theme_minimal()

#!# If we do this officially, smooth mean/sd to account for sample size better
#!# should be similar to Jonathan's spline functions for the length/age key should work
#!# or actually just exactly those

## test age noise impact
model_reps <- list()
reppreds <- c()
corrs <- c()
nrep = 1
#loop to select ages with rnorm
while (nrep < 101) {
  i <- 1
  testage <- c()
  for (i in 1:length(unique(testdf$Age))) {
    row <- means_ref[which(means_ref$RefAge == unique(testdf$Age)[i]),]
    ages <- testdf %>% 
      dplyr::select(Age) %>% 
      filter(Age == unique(testdf$Age)[i])
    samp <- rnorm(length(ages$Age), mean = row$meanAge, sd = row$stddevAge)
    names(samp) <- ages$Age
    testage <- c(testage,samp)
  }
  
  # test linear model (simple) with the different ages
  # switch ages in reduceddf
  reduceddf_noise <- reduceddf %>% dplyr::select(-Age) %>% mutate(NoiseAge=testage)
  # run a linear model
  model <- lm(NoiseAge~.,data = reduceddf_noise)
  model_reps[[nrep]] <- model
  
  #run predictions
  preds <- predict.lm(model)
  tempreds <- data.frame(NoiseAge=reduceddf_noise$NoiseAge,
                                 PredAge=preds,
                                 rep=nrep)
  
  reppreds <- rbind(reppreds,tempreds)
  
  # test a correlation
  corrs <- c(corrs,cor.test(reduceddf_noise$NoiseAge,preds)[["estimate"]])
  
  nrep = nrep + 1
}

hist(corrs)

reppreds %>% 
  group_by(rep) %>% 
  summarise(MAE=median(abs(PredAge - NoiseAge))) %>% 
  ggplot(aes(x=MAE))+
  geom_boxplot()+
  theme_minimal()

reppreds %>% 
  ggplot(aes(x=NoiseAge,y=PredAge))+
  geom_line(aes(group=rep),alpha = 0.1,color="forestgreen")+
  geom_abline(slope = 1,intercept = 0) +
  theme_classic()+
  ylim(0, 30) + xlim(0, 30)







