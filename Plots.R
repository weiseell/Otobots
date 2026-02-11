library(tidyverse)
library(wesanderson)
library(karyoploteR)
#!# install with: BiocManager::install("karyoploteR")
## load in model predictions from all three options
preds_comb <- read.csv("PredsComb_BestModels.csv")
preds_fem <- read.csv("PredsFemaleOnly_BestModels.csv")
preds_mal <- read.csv("PredsMaleOnly_BestModels.csv")
allouts <- rbind(preds_comb,preds_fem,preds_mal)

# 
## load in model predictions from all three options
preds_comb <- read.csv("PredsComb_BestModels.csv")
preds_fem <- read.csv("PredsFemaleOnly_BestModels.csv")
preds_mal <- read.csv("PredsMaleOnly_BestModels.csv")
allouts <- rbind(preds_comb,preds_fem,preds_mal)

# load libraries
library(tidyverse)
library(wesanderson)
## Model Prediction Plots
pdf(file = "Figures/ModelOptions_012926.pdf",
    height = 10, width = 10)
allouts %>% 
  ggplot(aes(x = Age, y = EpiAge, colour = as.factor(ObsSex),shape = type))+
  facet_grid(SampleSet~ModelType,switch = "y")+
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
pdf(file = "Figures/ModelMAE_012926.pdf",
    height = 10, width = 8)
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
  theme_classic(base_size = 15)+
  theme(legend.position = "none")
dev.off()

# calculate per age MAE and check for a trend
pdf(file = "Figures/ModelMAETrends_012926.pdf",
    height = 10, width = 10)
allouts %>% 
  group_by(SampleSet,ModelType,Age) %>% 
  summarise(MAE_Age = mean(abs(Age - EpiAge))) %>% 
  ungroup() %>% 
  ggplot(aes(x = Age, y = MAE_Age))+
  facet_grid(SampleSet~ModelType,switch = "y",scales = "free_y")+
  geom_point()+
  theme_bw()+
  xlab("Otolith Age") +
  ylab("Mean Absolute Error - Age") +
  ggtitle("Absolute Error Per Age")
dev.off()
# plot intron regions and sizes
CpGlist_all <- read.csv("Output/CpGList_LengthAge_Regions.csv")
ChromSwap <- read.csv("Input/ChromosomeSwap.csv")
GenInfo <- read.delim("~/Downloads/sequence_report.tsv")

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
             regmax=max(POS),
         CpGName=paste0(CHROM1, ":",POS))

# add colors to make both plots consistent
CpGlist_plot$col <- NA
CpGlist_plot$col[which(CpGlist_plot$ModelType == "Length")] <- wes_palette("AsteroidCity1")[1]
CpGlist_plot$col[which(CpGlist_plot$ModelType == "Age")] <- wes_palette("AsteroidCity1")[3]
CpGlist_plot$col[which(CpGlist_plot$ModelType == "Both")] <- wes_palette("AsteroidCity1")[2]


labels <- CpGlist_plot %>% 
  filter(n > 4) %>% 
  group_by(Intron) %>% 
  summarise(Length = max(POS)-min(POS),
            LenName = paste("Length = \n",Length))
labels$RegName = paste0("Region ",LETTERS[1:length(labels$Intron)])
CpGlist_plot$RegName <- NA
for(i in 1:length(labels$Intron)){
  CpGlist_plot$RegName[which(CpGlist_plot$Intron == labels$Intron[i])] <- labels$RegName[i]
}

jpeg(file = "Figures/CpGRegion_VisualbyGenomeRegion.jpeg",
     width = 500, height = 700,quality = 1.0)
CpGlist_plot %>% 
  filter(n > 4) %>% 
  mutate(x1=POS-0.2,x2=POS+0.2,y1=0,y2=1,
         #IntronName=paste0(CHROM1, ":",regmin,"-\n",regmax),
         Length = max(POS)-min(POS),
         LenName = paste("Length = \n",Length)) %>% 
  arrange(RegName) %>% 
  ggplot()+
    facet_wrap(~RegName,ncol = 2,
               scales = "free_x",
               strip.position = "left",axis.labels = "all")+
    geom_rect(aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,colour = ModelType, fill = ModelType))+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    scale_color_manual(values = c(wes_palette("AsteroidCity1")[1],
                                wes_palette("AsteroidCity1")[2],
                                wes_palette("AsteroidCity1")[3]))+
    scale_fill_manual(values = c(wes_palette("AsteroidCity1")[1],
                                wes_palette("AsteroidCity1")[2],
                                wes_palette("AsteroidCity1")[3]))+
    labs(x = "Chromosome Position (base)")+
    ggtitle("Location of Age and Length Important CpGs within Genomic Regions") %>% 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(size = 14),
          axis.text.x=element_text(angle=90,size = 10),
          axis.title.x = element_text(size = 12),
          panel.spacing.x = unit(2,"lines"),
          panel.spacing.y = unit(1,"lines"),
          legend.position = "none")
dev.off()

# make a compact genome plot with intron and CpG location labels
# lengths of chromosomes
RegMinMax <- CpGlist_plot %>% 
  group_by(Intron) %>% 
  dplyr::slice(1)

GenInfo <- GenInfo %>% 
  filter(Chromosome.name != "Un" & Chromosome.name != "MT") %>% 
  filter(Chromosome.name %in% CpGlist_plot$CHROM1)

chrs <- toGRanges(data.frame(chr = GenInfo$Chromosome.name,
                             start = 1,end = GenInfo$Seq.length))

CpGs <- toGRanges(data.frame(chr = CpGlist_plot$CHROM1,
                             start = CpGlist_plot$POS,end = CpGlist_plot$POS+1,
                             labels = CpGlist_plot$CpGName),col = CpGlist_plot$col)

regs <- toGRanges(data.frame(chr = RegMinMax$CHROM1,
                             start = RegMinMax$regmin,end = RegMinMax$regmax,
                             labels = paste0(LETTERS[1:length(RegMinMax$CHROM1)]," (n=",RegMinMax$n,")")))
plotpars <- getDefaultPlotParams(plot.type = 2)
plotpars$topmargin <- 400
plotpars$data1height <- 250
plotpars$data2height <- 250
jpeg(file = "Figures/CpGRegion_GenomeWithCpGAnnotation.jpeg",
     width = 300, height = 700)
kp <- plotKaryotype(genome = chrs, 
                    plot.type = 2,
                    plot.params = plotpars,
                    main = "Selected CpG (n = 151) Chromosome Positions")
kpPlotMarkers(kp,regs,text.orientation = "horizontal",adjust.label.position = T)
kpArrows(kp,CpGs,length = 0.1,y0 = -0.4,y1 = -0.28,col = CpGs$col)
dev.off()




