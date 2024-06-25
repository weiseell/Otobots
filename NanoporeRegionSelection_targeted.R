#Scripts for additional figures in SNP panel manuscript

#1. Visualization of SNP panel

#load libraries
library(quantsmooth)
library(tidyverse)
library(stringr)

#load data
ChromSwap <- read.csv("Input/ChromosomeSwap.csv")
CpGlist_best_comb <- read.csv(file = "Output/CpGList_LengthAge.csv")

##Dataset cleaning ####
#alternative file with only best CpGs from model
CpGlocs <- CpGlist_best_comb %>% 
  select(CHROM,POS,ModelType) %>% 
  rename(MapInfo = POS,Type = ModelType)

## ID CpG-rich regions (less than 1000b apart) ####
CpGRegion <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(CpGRegion) <- c("CHROM","StartPOS","EndPOS")

for (k in 1:length(unique(CpGlocs$CHROM))) {
  currentCHR <- unique(CpGlocs$CHROM)[k]
  CpG_Chr <- subset(CpGlocs,CpGlocs$CHROM == currentCHR)
  
  if(nrow(CpG_Chr) == 1){
    #make a data frame
    df <- data.frame(matrix(nrow = 1,ncol = 3))
    colnames(df) <- c("CHROM","StartPOS","EndPOS")
    
    #add the info from the one CpG
    df[1,1] <- CpG_Chr$CHROM
    df[1,2] <- CpG_Chr$MapInfo
    df[1,3] <- CpG_Chr$MapInfo
  }
  
  if(nrow(CpG_Chr) == 2){
    if(abs(CpG_Chr$MapInfo[2] - CpG_Chr$MapInfo[1]) <= 5000){
      #make a data frame
      df <- data.frame(matrix(nrow = 1,ncol = 3))
      colnames(df) <- c("CHROM","StartPOS","EndPOS")
      
      #add the info from the combined CpG group
      df[1,1] <- CpG_Chr$CHROM[1]
      df[1,2] <- CpG_Chr$MapInfo[1]
      df[1,3] <- CpG_Chr$MapInfo[2]
    }
    if(abs(CpG_Chr$MapInfo[2] - CpG_Chr$MapInfo[1]) > 5000){
      #make two separate regions since the CpGs are too far apart to be combines
      #make a data frame
      df <- data.frame(matrix(nrow = 2,ncol = 3))
      colnames(df) <- c("CHROM","StartPOS","EndPOS")
      
      #add the info from the two CpGs
      df[1,1] <- CpG_Chr$CHROM[1]
      df[1,2] <- CpG_Chr$MapInfo[1]
      df[1,3] <- CpG_Chr$MapInfo[1]
      
      df[2,1] <- CpG_Chr$CHROM[2]
      df[2,2] <- CpG_Chr$MapInfo[2]
      df[2,3] <- CpG_Chr$MapInfo[2]
    }
  }
  if(nrow(CpG_Chr) > 2){
    #if there's more than 3 CpGs we have to check all of them like we did for the big group
    #setting the first two CpGs in Chr to be the beginning of the first region
    #start 1
    mins <- c()
    mins <- CpG_Chr[1,2]
    #start 2
    maxes <- c()
    if(CpG_Chr[2,2] - CpG_Chr[1,2] >= 5000){
      maxes <- CpG_Chr[1,2]
    }
    tmpmax <- CpG_Chr[2,2]
    
    ## loop to look for islands
    for (i in 3:length(CpG_Chr$CHROM)) {
      #set current max
      #if previous loop set a new maximum, set a new miniumum
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] >= 5000){
        mins <- c(mins,CpG_Chr[i,2])
      }
      
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] >= 5000 & CpG_Chr[i,2] - CpG_Chr[i-1,2] >= 5000){
        maxes <- c(maxes,CpG_Chr[i,2])
      }
      
      if(CpG_Chr[i-1,2] - CpG_Chr[i-2,2] < 5000){
        if(CpG_Chr[i,2] - CpG_Chr[i-1,2] < 5000){
          tmpmax <- CpG_Chr[i,2]
        }
        
        if(CpG_Chr[i,2] - CpG_Chr[i-1,2] >= 5000){
          maxes <- c(maxes,tmpmax)
          
        }
      }
    }
    
    if(CpG_Chr[length(CpG_Chr$CHROM),2] - CpG_Chr[length(CpG_Chr$CHROM)-1,2] >= 5000){
      mins <- c(mins,CpG_Chr[length(CpG_Chr$CHROM),2])
    }
    
    maxes <- c(maxes,CpG_Chr[length(CpG_Chr$CHROM),2])
    ## making data frame out of mins and maxes
    df <- data.frame(matrix(nrow = length(mins),ncol = 3))
    colnames(df) <- c("CHROM","StartPOS","EndPOS")
    
    df[,1] <- currentCHR
    df[,2] <- mins
    df[,3] <- maxes
    
  }
  CpGRegion <- rbind(CpGRegion,df)
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

sum(as.integer(CpGregion_addinv$diff))/750000000
write.table(CpGRegion[,1:3],file = "Output/NanoporeRegionSelection_040424.txt",quote = F,sep = "\t",row.names = F,col.names = F)


