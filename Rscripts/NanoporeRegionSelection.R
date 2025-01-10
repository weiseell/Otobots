#load libraries
library(quantsmooth)
library(tidyverse)
library(stringr)

#load data
df <- read.csv("Output/CpGList_LengthAge.csv")
ChromSwap <- read.csv("Input/ChromosomeSwap.csv")
load("Input/allCpGs_nosingletons.rda")
CpGlist_best_comb <- read.csv(file = "Output/CpGList_LengthAge.csv")

##Dataset cleaning ####
#remove CpGs with NAs for the model
nas <- colSums(is.na(allCpGs1))
allCpGs1 <- allCpGs1[,which(nas < 1)]
#write.csv(allCpGs1,"allCpGs_nonas.csv",append = F,quote = F)

#remove NW and mito data from data set
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NW",
                                    colnames(allCpGs1),perl = T)))]
allCpGs1 <- allCpGs1[,which(!(grepl(pattern = "NC_009709.1",
                                    colnames(allCpGs1),perl = T)))]

# making CpG file
CpGnames <- colnames(allCpGs1)
tmp <- data.frame(str_split_fixed(string = CpGnames,pattern = "\\.",n = 3))
CpGlocs <- tmp %>% 
  mutate(CHROM = paste(X1,X2,sep = "."),Type = "CpG") %>% 
  rename(MapInfo=X3) %>% 
  dplyr::select(CHROM,MapInfo,Type)

#alternative file with only best CpGs from model
CpGlocs <- df %>% 
  dplyr::select(CHROM,POS,ModelType) %>% 
  rename(MapInfo = POS,Type = ModelType)

#make SNP panel plot####
CpGlocs$CHR<- NA
for (i in 1:length(ChromSwap$NCseq)) {
  CpGlocs$CHR[which(CpGlocs$CHROM == ChromSwap$NCseq[i])] <- ChromSwap$Einfeldt[i]
}

CpGlocs1 <- CpGlocs %>% dplyr::select(CHR,MapInfo,Type)
chrom <- ChromSwap[,3:4]
chrom <- chrom %>% rename(CHR = Einfeldt,MapInfo = Length) %>% mutate(Type = "ChrLength")


chrom <- rbind(chrom,CpGlocs1)
chrom$MapInfo <- as.integer(chrom$MapInfo)

## Making plot to show all CpGs and their location on the genome
tiff(filename = "Figures/CpGVisualization_032224.tiff",width = 84,height = 150,units = "mm",res = 400)
chrompos <- prepareGenomePlot(chrompos = chrom[,1:2], cols = "grey50",bleach = 0, topspace = 1, sexChromosomes = FALSE)
points(chrompos[25:length(chrom$CHR),2],chrompos[25:length(chrom$CHR),1]+0.05, , pch="|", cex = 0.75, col="deepskyblue4")
title("Atlantic Halibut CpG \nPositions by Chromosome")
dev.off()

## ID CpG-rich regions (less than 1000b apart) ####
CpGRegion <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(CpGRegion) <- c("CHROM","StartPOS","EndPOS")

for (k in 1:length(unique(CpGlocs$CHROM))) {
  currentCHR <- unique(CpGlocs$CHROM)[k]
  CpG_Chr <- subset(CpGlocs,CpGlocs$CHROM == currentCHR)
  
  #need at least 2 CpGs in the chromosome for this to work (different loop for chromosomes with one site)
  if(length(CpG_Chr$CHROM) >= 2){
    CpG_Chr$MapInfo <- as.integer(CpG_Chr$MapInfo)
    CpG_Chr <- CpG_Chr %>% arrange(CHROM,MapInfo)
    
    #setting the first two CpGs in Chr to be the beginning of the first region
    #start 1
    mins <- c()
    mins <- CpG_Chr[1,2]
    #start 2
    maxes <- c()
    if(CpG_Chr[2,2] - CpG_Chr[1,2] >= 10000){
      maxes <- CpG_Chr[2,2]
    }
    tmpmax <- CpG_Chr[2,2]
    
    ## loop to look for islands
    for (i in 3:length(CpG_Chr$CHROM)) {
      #set current max
      #if previous loop set a new maximum, set a new miniumum
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] >= 10000){
        mins <- c(mins,CpG_Chr[i,2])
        #tmpmax <- CpG_Chr[i,2]
      }
      
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] >= 10000 & CpG_Chr[i,2] - CpG_Chr[i-1,2] >= 10000){
        maxes <- c(maxes,CpG_Chr[i,2])
      }
      
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] < 10000){
        if(CpG_Chr[i,2] - CpG_Chr[i-1,2] < 10000){
          tmpmax <- CpG_Chr[i,2]
        }
        
        if(CpG_Chr[i,2] - CpG_Chr[i-1,2] >= 10000){
          maxes <- c(maxes,tmpmax)
          
        }
      }
    }
    
    if(CpG_Chr[length(CpG_Chr$CHROM),2] - CpG_Chr[length(CpG_Chr$CHROM)-1,2] >= 10000){
      mins <- c(mins,CpG_Chr[length(CpG_Chr$CHROM),2])
    }
    
    maxes <- c(maxes,CpG_Chr[length(CpG_Chr$CHROM),2])
    ## making data frame out of mins and maxes
    df <- data.frame(matrix(nrow = length(mins),ncol = 3))
    colnames(df) <- c("CHROM","StartPOS","EndPOS")
    
    df[,1] <- currentCHR
    df[,2] <- mins
    df[,3] <- maxes
    
    CpGRegion <- rbind(CpGRegion,df)
  }
  
  if(length(CpG_Chr$CHROM) == 1){
    df <- data.frame(matrix(nrow = 1,ncol = 3))
    colnames(df) <- c("CHROM","StartPOS","EndPOS")
    df[1,1] <- CpG_Chr$CHROM
    df[1,2] <- CpG_Chr$MapInfo
    df[1,3] <- CpG_Chr$MapInfo
    
    CpGRegion <- rbind(CpGRegion,df)
  }
}
## Adjusting the regions so they cover at least 10000 bp
for (i in 1:length(CpGRegion$CHROM)) {
  dist <- CpGRegion$EndPOS[i] - CpGRegion$StartPOS[i]
  if(dist < 10000){
    addto <- round((10000 - dist)/2,digits = 0)
    newmin <- CpGRegion$StartPOS[i] - addto
    newmax <- CpGRegion$EndPOS[i] + addto
    chrmax <- ChromSwap$Length[which(ChromSwap$NCseq == CpGRegion$CHROM[i])]
    if(newmin >= 0 & newmax <= chrmax){
      CpGRegion$StartPOS[i] <- newmin
      CpGRegion$EndPOS[i] <- newmax
    }
    
    if(newmin < 0){
      newdist <- 10000 - CpGRegion$StartPOS[i] + 1
      newmax <- CpGRegion$EndPOS[i] + newdist
      CpGRegion$StartPOS[i] <- 1
      CpGRegion$EndPOS[i] <- newmax
    }
    
    if(newmax > chrmax){
      newdist <- 10000 - (chrmax - CpGRegion$EndPOS)
      CpGRegion$StartPOS[i] <- CpGRegion$StartPOS - newdist
      CpGRegion$EndPOS[i] <- chrmax
    }
  }
}

#size of all CpG regions to be amplified
CpGRegion$diff <- CpGRegion$EndPOS - CpGRegion$StartPOS

#comparing total size compared to full genome sequence (750MB)
sum(CpGRegion$diff)/750000000

## Removing sections related to inversion and then adding the full inversion region with a buffer
CpGregion_noinv <- CpGRegion[which(!(CpGRegion$CHROM == "NC_047158.1" & CpGRegion$StartPOS >= 4000000 & CpGRegion$EndPOS <= 13000000)),]

# inversion
colnames(CpGRegion)
inv <- c("NC_047158.1",as.integer(4000000),as.integer(13000000),as.integer(13000000-4000000))
length(CpGregion_noinv$CHROM)
#add a blank row
CpGregion_noinv[264,] <- c("NC_047158.1",as.integer(4000000),as.integer(13000000),as.integer(13000000-4000000))
CpGregion_addinv <- CpGregion_noinv %>% arrange(CHROM,as.integer(StartPOS))

sum(as.integer(CpGregion_addinv$diff))/750000000
write.table(CpGregion_addinv[,1:3],file = "Output/NanoporeRegionSelection_040424.txt",quote = F,sep = "\t",row.names = F,col.names = F)






