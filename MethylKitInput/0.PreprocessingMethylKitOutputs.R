
##pre-models: concatenate all the methyl files together and remove missing data
setwd("E:/ClonedRepositories/Otobots/MethylKitInput/")

library(tidyverse)

#read file names in the folder
filenames <- list.files(pattern = "MethylWild")
filenames <- filenames[grep(x = filenames,pattern = "CpG.txt")]
filenames1 <- list.files(pattern = "HPEI")
filenames1 <- filenames1[grep(x = filenames1,pattern = "CpG.txt")]
filenames <- c(filenames,filenames1)

tmpname <- unlist(strsplit(filenames,split = "2023_"))
tmpname <- unlist(strsplit(tmpname,split = "2025_"))
tmpname <- matrix(unlist(strsplit(tmpname,split = "HPEI_")),nrow = 2)[2,]
tmpname <- matrix(unlist(strsplit(tmpname,split = "_CpG")),nrow = 2)[1,]

#process the first file in the list to start the data frame
indiv1 <- read.table(filenames[1],header = T,sep = "\t")
indiv1 <- indiv1 %>% 
  filter(coverage >= 30) %>% 
  slice(grep(pattern = "NW",indiv1$chrBase,perl = T,invert=T)) %>% 
  slice(grep(pattern = "NC_009709.1",indiv1$chrBase,perl = T,invert=T)) %>% 
  dplyr::select(chrBase,freqT,coverage) %>% 
  mutate(indivname=tmpname[1],filename=filenames[1])

allCpGs <- indiv1

#merge all the %T sites together
# makes a very long skinny data frame that can be spread after it's made
for (i in 2:125) {
  print(i)
  indiv1 <- read.table(filenames[i],header = T,sep = "\t")
  indiv1 <- indiv1 %>% 
    filter(coverage >= 30) %>% 
    slice(grep(pattern = "NW",indiv1$chrBase,perl = T,invert=T)) %>% 
    slice(grep(pattern = "NC_009709.1",indiv1$chrBase,perl = T,invert=T)) %>% 
    dplyr::select(chrBase,freqT,coverage) %>% 
    mutate(indivname=tmpname[i],filename=filenames[i])
  
  allCpGs <- rbind(allCpGs,indiv1)
}

#id number of indivs per CpG
CpGcount <- allCpGs %>% 
  group_by(chrBase) %>% 
  summarise(counts=n()) %>% 
  # remove singletons
  filter(counts > 1)

hist(CpGcount$counts)

CpGcount <- CpGcount %>% 
  # filter out mostly missing CpGs
  # include only CpGs across all data sets to start
  filter(counts > 60)
# remove singletons in long data frame and then spread with 
# indivs in rows and CpGs in columns
allCpGs1 <- allCpGs[which(allCpGs$chrBase %in% CpGcount$chrBase),]

allCpGs1 <- allCpGs1 %>%
  dplyr::select(chrBase,freqT,indivname) %>% 
  spread(key = chrBase,value = freqT)

allCpGs1_30X <- allCpGs1
# save file with no missing data
save(allCpGs_30X, file = "allCpGs_AllReads_30Xcov.rda")
save(allCpGs1_30X,file = "allCpGs_MinimalMissing_30Xcov.rda")

# look at average depth per individual for called CpGs
allCpGs <- allCpGs %>% 
  mutate(tmp = chrBase) %>% 
  separate_wider_delim(tmp,delim = ".",names = c("CHROM","ones","POS")) %>% 
  unite(col=CHROM,CHROM,ones,sep = ".")
  
indivdepth <- allCpGs %>% 
  group_by(indivname,CHROM) %>% 
  summarise(avg = mean(coverage),
            meds = median(coverage))

# are there depth differences across the batches?
#!# doesn't look like it at first glance

# run PCA to see if there's batch issues/environmental impacts
















 