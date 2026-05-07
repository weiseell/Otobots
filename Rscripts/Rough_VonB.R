## Rough VonB curves based on current otolith data

## load libraries
library(tidyverse)
library(FSA)
library(rsample)
library(car)
library(ggpubr)

## load data
dat <- read.csv("Input/SampleMetadataExport_WithAge_2025_Redit.csv")
episamps <- read.csv("Input/MethylWild_Sex_AgeData.csv")
load("Input/SCALmodeldata2022.rdata")
table(dat$OtoAge)

# eliminate data with small sample size
vb <- makeGrowthFun(type = "von Bertalanffy")

#age range is 3 and 15
dat <- data.frame(indivname=dat$SampleID,ObsSex = dat$ObsSex, Length = dat$Length, OtoAge = dat$OtoAge)

# separating data between males and females
fem <- dat %>% filter(ObsSex == 2)
male <- dat %>% filter(ObsSex == 1)

# identify initial parameters for VonB curves
parainit_f <- findGrowthStarts(Length~OtoAge,fem,type = "von Bertalanffy")
parainit_m <- findGrowthStarts(Length~OtoAge,male,type = "von Bertalanffy")

# estimate vonb parameters
fit_f <- nls(Length~vb(OtoAge,Linf,K,t0),data=fem,start=parainit_f)
fit_m <- nls(Length~vb(OtoAge,Linf,K,t0),data=male,start=parainit_m)
# coefficients
coef(fit_f)
coef(fit_m)

# bootstrapping for uncertainty
boot_f <- car::Boot(fit_f) 
confint(boot_f)

boot_m <- car::Boot(fit_m)
confint(boot_m)

## plotting data for female clock
# predict mean length at age for modeled ages
len_pred <- predict(fit_f,newdata = data.frame(OtoAge=3:29))
predict2 <- function(x) predict(x,data.frame(OtoAge=ages))

ages <- seq(3,29,by=1)

f.boot2 <- car::Boot(fit_f,f=predict2)

preds1 <- data.frame(ages,
                     predict(fit_f,data.frame(OtoAge=ages)),
                     confint(f.boot2,type = "perc"))
names(preds1) <- c("age","fit","LCI","UCI")
headtail(preds1)

## using parameters from the Halibut SCAL model
para_fSCAL <- list(Linf = BiologicalParameters$lInf_f,
                   K = BiologicalParameters$vonK_f,
                   t0 = BiologicalParameters$wt_a)
preds_fSCAL <- data.frame(1:30,
                        vb(1:30,para_fSCAL$Linf,para_fSCAL$K,para_fSCAL$t0))
names(preds_fSCAL) <- c("age","fit")
headtail(preds_fSCAL)

## get data for episamps
genf <- episamps %>% 
  filter(ObsSex == 2)

genm <- episamps %>% 
  filter(ObsSex == 1)
# plot curve
fplot <- ggplot() + 
  geom_point(data=fem,aes(y=Length,x=OtoAge),size=2,alpha=0.1) +
  geom_point(data=genf,aes(y=Length,x=Age),size=2,alpha=0.8,col="darkblue") +
  #geom_ribbon(data=preds1,aes(x=age,ymin=LCI,ymax=UCI),fill="gray90") +
  #geom_line(data=preds1,aes(y=fit,x=age),size=1,linetype=2) +
  #geom_line(data=preds_fSCAL,aes(y=fit,x=age),size=1,linetype=3) +
  scale_y_continuous(name="Total Length (cm)",limits=c(0,240),expand=c(0,0)) +
  scale_x_continuous(name="Age (years)",expand=c(0,0),
                     limits=c(0,30),breaks=seq(0,30,2)) +
  theme_bw() +
  theme(panel.grid=element_blank())+
  ggtitle("Length-Age Genetic Females")

## plotting data for male clock
# predict mean length at age for modeled ages
len_pred <- predict(fit_m,newdata = data.frame(OtoAge=3:29))

ages <- seq(3,29,by=1)

boot_m2 <- car::Boot(fit_m,f=predict2)

preds2 <- data.frame(ages,
                     predict(fit_m,data.frame(OtoAge=ages)),
                     confint(boot_m2,type = "perc"))
names(preds2) <- c("age","fit","LCI","UCI")
headtail(preds2)

## using parameters from the Halibut SCAL model
para_mSCAL <- list(Linf = BiologicalParameters$lInf_m,
                   K = BiologicalParameters$vonK_m,
                   t0 = BiologicalParameters$wt_a)
preds_mSCAL <- data.frame(1:30,
                          vb(1:30,para_mSCAL$Linf,para_mSCAL$K,para_mSCAL$t0))
names(preds_mSCAL) <- c("age","fit")
headtail(preds_mSCAL)

# plot curve
mplot <- ggplot() + 
  #geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI),fill="gray90") +
  geom_point(data=male,aes(y=Length,x=OtoAge),size=2,alpha=0.1) +
  geom_point(data=genm,aes(y=Length,x=Age),size=2,alpha=0.8,col="darkblue") +
  #geom_line(data=preds2,aes(y=fit,x=age),size=1,linetype=2) +
  #geom_line(data=preds_mSCAL,aes(y=fit,x=age),size=1,linetype=3) +
  scale_y_continuous(name="Total Length (cm)",limits=c(0,200),expand=c(0,0)) +
  scale_x_continuous(name="Age (years)",expand=c(0,0),
                     limits=c(0,30),breaks=seq(0,30,2)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  ggtitle("Length-Age for Genetic Males")

## combine plots into two subplots
tiff(filename = "Figures/LengthAge_EpiAge_042226.tiff",width = 300,height = 150,units = "mm",res = 400)
ggarrange(fplot,mplot,labels = c("A","B"),ncol = 2,nrow = 1)
dev.off()
