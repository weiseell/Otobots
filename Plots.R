## load in model predictions from all three options
preds_comb <- read.csv("PredsComb_BestModels.csv")
preds_fem <- read.csv("PredsFemaleOnly_BestModels.csv")
preds_mal <- read.csv("PredsMaleOnly_BestModels.csv")
allouts <- rbind(preds_comb,preds_fem,preds_mal)

# 

## Model Prediction Plots
allouts %>% 
  ggplot(aes(x = Age, y = EpiAge, colour = as.factor(ObsSex),shape = type))+
    facet_wrap(SampleSet~ModelType)+
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

## Absolute Error Plots
allouts %>% 
  mutate(AbsErr = abs(Age - EpiAge)) %>% 
  ggplot(aes(y = AbsErr, x = ModelType, fill = SampleSet))+
    geom_boxplot(alpha = 0.7)+
    theme_bw()+
    scale_fill_manual(values = c(wes_palette("IsleofDogs1")[6],
                                  wes_palette("IsleofDogs1")[1],
                                  wes_palette("IsleofDogs1")[3]))+
  xlab("Model Option") +
  ylab("Absolute Error") +
  ggtitle("Absolute Error")+
  theme_classic(base_size = 15)
  
# plot intron regions and sizes
library(tidyverse)
library(GenomicRanges)

CpGlist_all <- read.csv("Output/CpGList_LengthAge_Regions.csv")
ChromSwap <- read.csv("Input/ChromosomeSwap.csv")







#make SNP panel plot####
chrom <- CpGlist_all %>% 
  dplyr::select(Intron,POS,AgeRange) %>% 
  rename(CHR=Intron,MapInfo=POS)
chrom$CHR <- as.integer(chrom$CHR)
chrom$MapInfo <- as.integer(chrom$MapInfo)

chrom <- chrom %>% arrange(CHR,MapInfo)





#tiff(filename = "Figures/CpGIntronVisualization.tiff",width = 84,height = 150,units = "mm",res = 400)
#chrompos <- prepareGenomePlot(chrom, cols = "grey50",bleach = 0, topspace = 1, sexChromosomes = FALSE)
#points(chrompos[,2],chrompos[,1]+0.05, pch="|", cex = 0.75, col="deepskyblue4")
#dev.off()





