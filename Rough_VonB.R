## Rough VonB curves based on current otolith data

## load libraries
library(tidyverse)
library(FSA)
library(FSAtools)
library(car)
library(ggpubr)

## load data
dat <- read.csv("Input/SampleMetadataExport_WithAge_9Feb2023_Redit.csv")

table(dat$OtoAge)

load("Input/")
# eliminate data with small sample size
vb <- vbFuns(param="Typical")

#age range is 3 and 15
dat <- dat %>% 
  filter(OtoAge > 2 & OtoAge < 16)
dat <- data.frame(ObsSex = dat$ObsSex, Length = dat$Length, OtoAge = dat$OtoAge)

# separating data between males and females
fem <- dat %>% filter(ObsSex == 2)
male <- dat %>% filter(ObsSex == 1)

# identify initial parameters for VonB curves
parainit_f <- vbStarts(Length~OtoAge,fem)
parainit_m <- vbStarts(Length~OtoAge,male)

# estimate vonb parameters
fit_f <- nls(Length~vb(OtoAge,Linf,K,t0),data=fem,start=parainit_f)
fit_m <- nls(Length~vb(OtoAge,Linf,K,t0),data=male,start=parainit_m)
# coefficients
coef(fit_f)
coef(fit_m)

# bootstrapping for uncertainty
boot_f <- Boot(fit_f)
confint(boot_f)

boot_m <- Boot(fit_m)
confint(boot_m)

## plotting data for female clock
# predict mean length at age for modeled ages
len_pred <- predict(fit_f,newdata = data.frame(OtoAge=3:15))
predict2 <- function(x) predict(x,data.frame(OtoAge=ages))

ages <- seq(3,15,by=1)

f.boot2 <- Boot(fit_f,f=predict2)

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

# plot curve
fplot <- ggplot() + 
  geom_point(data=fem,aes(y=Length,x=OtoAge),size=2,alpha=0.1) +
  geom_ribbon(data=preds1,aes(x=age,ymin=LCI,ymax=UCI),fill="gray90") +
  geom_line(data=preds1,aes(y=fit,x=age),size=1,linetype=2) +
  geom_line(data=preds_fSCAL,aes(y=fit,x=age),size=1,linetype=3) +
  scale_y_continuous(name="Total Length (cm)",limits=c(0,220),expand=c(0,0)) +
  scale_x_continuous(name="Age (years)",expand=c(0,0),
                     limits=c(0,30),breaks=seq(0,30,2)) +
  theme_bw() +
  theme(panel.grid=element_blank())+
  ggtitle("VonB Length-Age Curve for Genetic Females")

## plotting data for male clock
# predict mean length at age for modeled ages
len_pred <- predict(fit_m,newdata = data.frame(OtoAge=3:15))

ages <- seq(3,15,by=1)

boot_m2 <- Boot(fit_m,f=predict2)

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
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI),fill="gray90") +
  geom_point(data=male,aes(y=Length,x=OtoAge),size=2,alpha=0.1) +
  geom_line(data=preds2,aes(y=fit,x=age),size=1,linetype=2) +
  geom_line(data=preds_mSCAL,aes(y=fit,x=age),size=1,linetype=3) +
  scale_y_continuous(name="Total Length (cm)",limits=c(0,200),expand=c(0,0)) +
  scale_x_continuous(name="Age (years)",expand=c(0,0),
                     limits=c(0,30),breaks=seq(0,30,2)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  ggtitle("VonB Length-Age Curve for Genetic Males")

## combine plots into two subplots
tiff(filename = "Figures/VonB_GenAndSCAL_042424.tiff",width = 300,height = 150,units = "mm",res = 400)
ggarrange(fplot,mplot,labels = c("A","B"),ncol = 2,nrow = 1)
dev.off()
