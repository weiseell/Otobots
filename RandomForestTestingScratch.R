### Random Forest Testing hooray!!
##load data
load("MethylKitInput/allCpGs_MinimalMissing.rda")
sexdata <- read.csv("Input/MethylWild_Sex_AgeData.csv")
CpGRegion <- read.csv(file = "Output/BestCpG_RegionGroups.csv")

CpGlist_best_Age4to15 <- read.csv(file = "Output/CpGList_LengthAge_Age4to15_121025.csv")
CpGlist_best_Age16to29 <- read.csv(file = "Output/CpGList_LengthAge_Age4to29_121025.csv")

#load libraries
library(tidyverse)
library(ranger)

# plot missing data levels and subset to minimize missing data
nas <- colSums(is.na(allCpGs1))
hist(nas,breaks = 80)
#!# for this data set, looks like low levels of increasing missing data
#!# no jumps at numbers where it would be data-set based, as we're combining
#!# three different sequencing runs

# subset to just wild data from ages 4-19 for testing
# goal is to prevent low sample size and environmental bias
# for initial model, just to see if it works before adding more problems
RF_CpGs <- sexdata %>% 
  filter(Age < 20 & DataType == "Wild") %>% 
  # subselect a training data set
  group_by(Age,ObsSex) %>% 
  slice_sample(prop = 0.8) %>% 
  left_join(allCpGs1)

nas <- colSums(is.na(RF_CpGs))
hist(nas,breaks = 80)
#!# this one has a right-skewed bump at around 80 missing points
#!# likely due to previous pre-filtering
#!# removing CpGs with missing data to test
RF_CpGs <- RF_CpGs[,which(nas < 1)]

## subsetting RF_CpGs based on elastic net identified regions
CpGnames <- data.frame(CpG=colnames(RF_CpGs)) %>% 
  mutate(tmp=CpG) %>% 
  separate_wider_delim(tmp,delim = ".",names = c("CHROM","ones","MapInfo"),too_few = "align_start") %>% 
  unite(col=CHROM,CHROM,ones,sep = ".") %>% 
  dplyr::select(CHROM,MapInfo) %>% 
  mutate(MapInfo=as.numeric(MapInfo)) %>% 
  arrange(CHROM,MapInfo)

RegCpGs <- NA
for (i in 1:nrow(CpGRegion)) {
  regtmp <- CpGRegion[i,]
  nametmp <- CpGnames[which(CpGnames$CHROM == regtmp$CHROM &
          CpGnames$MapInfo >= regtmp$StartPOS &
          CpGnames$MapInfo <= regtmp$EndPOS),]
  RegCpGs <- rbind(RegCpGs,nametmp)
}

RegCpGs <- RegCpGs %>% 
  unite(col=names,CHROM,MapInfo,sep = ".")

RF_CpGs_reg <- RF_CpGs[,which(colnames(RF_CpGs) %in% RegCpGs$names)]
# removing columns except for age and CpGs and running through a test RF model
# this one uses all possible options
RF_Ages <- as.matrix(RF_CpGs$Age) 

# test variable option
RF_test <- ranger(x=RF_CpGs_reg,y = RF_Ages,
                  mtry = floor(ncol(RF_CpGs_reg/3)),num.trees = (nrow(RF_CpGs_reg)*100))
#!# performance is still pretty poor with regions

# trying hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(ncol(RF_CpGs_reg) * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)

for(i in 1:nrow(hyper_grid)) {
  print(i)
  # fit model for ith hyperparameter combination
  fit <- ranger(x=RF_CpGs_reg,y = RF_Ages,
                num.trees       = ncol(RF_CpGs_reg) * 10,
                mtry            = hyper_grid$mtry[i],
                min.node.size   = hyper_grid$min.node.size[i],
                replace         = hyper_grid$replace[i],
                sample.fraction = hyper_grid$sample.fraction[i],
                verbose         = FALSE,
                respect.unordered.factors = 'order')
  # report prediction error for each parameter set
  hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

default_rmse <- min(hyper_grid$rmse)
hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100)
#!# these are all quite tight, 30 models with <2% difference in performance

# trying one of the best model options to plot
#!# trying with way more trees too
#!# tried log ages since that helps with EN
fit_best <- ranger(x=RF_CpGs_reg,y = RF_Ages,
              num.trees       = ncol(RF_CpGs_reg) * 100,
              mtry            = hyper_grid$mtry[117],
              min.node.size   = hyper_grid$min.node.size[117],
              replace         = hyper_grid$replace[117],
              sample.fraction = hyper_grid$sample.fraction[117],
              verbose         = FALSE,
              respect.unordered.factors = 'order')

plot(RF_Ages,fit_best[["predictions"]])
#!# still not great haha

## trying with only BestCpG SITES id'd from EN models, rather than the regions
RF_CpGs_best <- RF_CpGs[,which(colnames(RF_CpGs) %in% CpGlist_best_Age4to15$CpGlist_best)]

RF_test2 <- ranger(x=RF_CpGs_best,y = RF_Ages,
                  mtry = floor(ncol(RF_CpGs_best/3)),num.trees = (nrow(RF_CpGs_best)*100))
#!# THIS MADE IT WAY WORSE UGH


#!# trying a totally different vibe, running a quick colinearity test on the 'Best CpG' set

