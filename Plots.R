library(tidyverse)
library(wesanderson)
## load in model predictions from all three options
preds_comb <- read.csv("PredsComb_BestModels.csv")
preds_fem <- read.csv("PredsFemaleOnly_BestModels.csv")
preds_mal <- read.csv("PredsMaleOnly_BestModels.csv")
allouts <- rbind(preds_comb,preds_fem,preds_mal)

# 

## Model Prediction Plots
jpeg(file = "Figures/CpGModelComparisonLinearvsNonLinear.jpeg",width = 800, height = 800)
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
dev.off()
## Absolute Error Plots
jpeg(file = "Figures/AbsError_ComparisonPlots.jpeg",width = 500, height = 500)
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
dev.off()
# plot intron regions and sizes
library(tidyverse)
library(GenomicRanges)

CpGlist_all <- read.csv("Output/CpGList_LengthAge_Regions.csv")
ChromSwap <- read.csv("Input/ChromosomeSwap.csv")

CpGlist_all$CHROM1 <- NA
for (i in 1:length(unique(CpGlist_all$CHROM))) {
  CHRtmp <- ChromSwap[which(ChromSwap$NCseq == unique(CpGlist_all$CHROM)[i]),]
  CpGlist_all$CHROM1[which(CpGlist_all$CHROM == CHRtmp$NCseq)] <- CHRtmp$Einfeldt
}
## remove the regions with one or two CpGs
# focusing on regions for now
CpGlist_plot <- CpGlist_all %>% 
  add_count(Intron) %>% 
  group_by(Intron) %>% 
  mutate(regmin=min(POS),
             regmax=max(POS)) %>% 
  filter(n > 4)

labels <- CpGlist_plot %>% 
  group_by(Intron) %>% 
  summarise(Length = max(POS)-min(POS),
            LenName = paste("Length = \n",Length))

jpeg(file = "Figures/CpGRegion_VisualbyGenomeRegion.jpeg",
     width = 500, height = 900)
CpGlist_plot %>% 
  mutate(x1=POS-0.2,x2=POS+0.2,y1=0,y2=1,
         IntronName=paste0(CHROM1, ":",regmin,"-\n",regmax)) %>% 
  arrange(IntronName) %>% 
  ggplot()+
    facet_wrap(~IntronName,ncol = 1,scales = "free_x",strip.position = "left")+
    geom_rect(aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,colour = ModelType, fill = ModelType))+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(size = 10),
          axis.text.x=element_text(angle=-90),
          title = element_text("Location of Age and Length Important CpGs within Genomic Regions"))
dev.off()



