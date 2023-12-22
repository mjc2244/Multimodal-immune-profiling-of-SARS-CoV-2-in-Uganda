#Figure 4e and Table S10b - Associations between immune cell abundances and Delta COVID via linear regression

#Clear R environment 
rm(list = ls())

#Import whole dataset 
combined <- read.csv(file.choose(), header=TRUE)
covid <- subset(combined, covid==1)

#Convert all PIDs to fixed row identifiers 
rownames(covid) = covid$pid
covid$pid = NULL

#Include 4 patients with unknown HIV status with HIV neg
covid$hivrdtresult[is.na(covid$hivrdtresult)] = 0

#Linear regressions looking at relationship between COVID vs. non COVID adjusted for age, sex, HIV status, steroid exposure, and WHO clinical severity
library(performance)
library(see)
library(patchwork)
library(jtools)

#Recode 0 to 0.0001 for immune cell abundances to allow log transformation
library(dplyr)   
covid <- covid %>% mutate(bcellsnaive = replace(bcellsnaive, bcellsnaive == 0, 0.0001))
covid <- covid %>% mutate(bcellsmemory = replace(bcellsmemory, bcellsmemory == 0, 0.0001))
covid <- covid %>% mutate(plasmacells = replace(plasmacells, plasmacells == 0, 0.0001))
covid <- covid %>% mutate(tcellscd8 = replace(tcellscd8, tcellscd8 == 0, 0.0001))
covid <- covid %>% mutate(tcellscd4naive = replace(tcellscd4naive, tcellscd4naive == 0, 0.0001))
covid <- covid %>% mutate(tcellscd4memoryresting = replace(tcellscd4memoryresting, tcellscd4memoryresting == 0, 0.0001))
covid <- covid %>% mutate(tcellscd4memoryactivated = replace(tcellscd4memoryactivated, tcellscd4memoryactivated == 0, 0.0001))
covid <- covid %>% mutate(nkcellsresting = replace(nkcellsresting, nkcellsresting == 0, 0.0001))
covid <- covid %>% mutate(monocytes = replace(monocytes, monocytes == 0, 0.0001))
covid <- covid %>% mutate(dendriticcellsactivated = replace(dendriticcellsactivated , dendriticcellsactivated  == 0, 0.0001))
covid <- covid %>% mutate(mastcellsresting = replace(mastcellsresting , mastcellsresting  == 0, 0.0001))
covid <- covid %>% mutate(eosinophils = replace(eosinophils , eosinophils  == 0, 0.0001))
covid <- covid %>% mutate(neutrophils = replace(neutrophils , neutrophils  == 0, 0.0001))

#Now log10 transform immune cell abundances
covid$logbcellsnaive = log10(covid$bcellsnaive)
covid$logbcellsmemory = log10(covid$bcellsmemory)
covid$logplasmacells = log10(covid$plasmacells)
covid$logtcellscd8 = log10(covid$tcellscd8)
covid$logtcellscd4naive = log10(covid$tcellscd4naive)
covid$logtcellscd4memoryresting = log10(covid$tcellscd4memoryresting)
covid$logtcellscd4memoryactivated = log10(covid$tcellscd4memoryactivated)
covid$lognkcellsresting = log10(covid$nkcellsresting)
covid$logmonocytes = log10(covid$monocytes)
covid$logdendriticcellsactivated = log10(covid$dendriticcellsactivated)
covid$logmastcellsresting = log10(covid$mastcellsresting)
covid$logeosinophils = log10(covid$eosinophils)
covid$logneutrophils = log10(covid$neutrophils)

#Naive B-cells
lmbcellsnaivecovid = lm(logbcellsnaive ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmbcellsnaivecovid)
check_model(lmbcellsnaivecovid)
summ(lmbcellsnaivecovid, confint = TRUE, digits=3)

#Memory B-cells
lmbcellsmemorycovid = lm(logbcellsmemory ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmbcellsmemorycovid)
check_model(lmbcellsmemorycovid)
summ(lmbcellsmemorycovid, confint = TRUE, digits=3)

#Plasma cells
lmplasmacellscovid = lm(logplasmacells ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmplasmacellscovid)
check_model(lmplasmacellscovid)
summ(lmplasmacellscovid, confint = TRUE, digits=3)

#CD8
lmtcellscd8covid = lm(logtcellscd8 ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmtcellscd8covid)
check_model(lmtcellscd8covid)
summ(lmtcellscd8covid, confint = TRUE, digits=3)

#Naive CD4
lmtcellscd4naivecovid = lm(logtcellscd4naive ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmtcellscd4naivecovid)
check_model(lmtcellscd4naivecovid)
summ(lmtcellscd4naivecovid, confint = TRUE, digits=3)

#Resting CD4
lmtcellscd4memoryrestingcovid = lm(logtcellscd4memoryresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmtcellscd4memoryrestingcovid)
check_model(lmtcellscd4memoryrestingcovid)
summ(lmtcellscd4memoryrestingcovid, confint = TRUE, digits=3)

#Activated CD4
lmtcellscd4memoryactivatedcovid = lm(logtcellscd4memoryactivated ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmtcellscd4memoryactivatedcovid)
check_model(lmtcellscd4memoryactivatedcovid)
summ(lmtcellscd4memoryactivatedcovid, confint = TRUE, digits=3)

#Resting NK
lmnkcellsrestingcovid = lm(lognkcellsresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmnkcellsrestingcovid)
check_model(lmnkcellsrestingcovid)
summ(lmnkcellsrestingcovid, confint = TRUE, digits=3)

#Monocytes
lmmonocytescovid = lm(logmonocytes ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmmonocytescovid)
check_model(lmmonocytescovid)
summ(lmmonocytescovid, confint = TRUE, digits=3)

#Activated DCs
lmdendriticcellsactivatedcovid = lm(logdendriticcellsactivated ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmdendriticcellsactivatedcovid)
check_model(lmdendriticcellsactivatedcovid)
summ(lmdendriticcellsactivatedcovid, confint = TRUE, digits=3)

#Resting mast cells 
lmmastcellsrestingcovid = lm(logmastcellsresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmmastcellsrestingcovid)
check_model(lmmastcellsrestingcovid)
summ(lmmastcellsrestingcovid, confint = TRUE, digits=3)

#Eosinophils
lmeosinophilscovid = lm(logeosinophils ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmeosinophilscovid)
check_model(lmeosinophilscovid)
summ(lmeosinophilscovid, confint = TRUE, digits=3)

#Neutrophils
lmneutrophilscovid = lm(logneutrophils ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) 
summary(lmneutrophilscovid)
check_model(lmneutrophilscovid)
summ(lmneutrophilscovid, confint = TRUE, digits=3)

#Plot Delta phase coefficients using forest plots
library(jtools)
library(ggplot2)
library(ggstance)
library(broom.mixed)
library(RColorBrewer)
library(jcolors)

library(modelsummary)
models <- list(
  "Naive B cells" = lm(logbcellsnaive ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) ,
  "Memory B cells" = lm(logbcellsmemory ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid) ,
  "Plasma cells" = lm(logplasmacells ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "CD8+" = lm(logtcellscd8 ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid), 
    "Naive CD4+" = lm(logtcellscd4naive ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Resting CD4+" = lm(logtcellscd4memoryresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Activated CD4+" = lm(logtcellscd4memoryactivated ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Resting NK cells" = lm(lognkcellsresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Monocytes" = lm(logmonocytes ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Activated DCs" = lm(logdendriticcellsactivated ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Resting mast cells" = lm(logmastcellsresting ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Eosinophils" = lm(logeosinophils ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid),
    "Neutrophils" = lm(logneutrophils ~ phasedelta + age + sex + hivrdtresult + steroids + whoseverity_imp, data = covid))

modelsummary(models, statistic = 'conf.int', "p.value")
modelplot(models, coef_omit = c(1,3,4,5,6,7))

library(ggplot2)
library(pals)
library(pals)

manualcolors2<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

b <- list(geom_vline(xintercept = 0, color = 'black', linetype=2))

regressionplotcells <- modelplot(models, background = b, coef_omit = c(1,3,4,5,6,7)) +
  labs(x = "Adjusted \u03B2 with 95% Confidence Interval", 
       y = "",
       title = "") + theme_bw() +
  scale_color_manual(values = manualcolors2) + theme(text=element_text(size = 12)) 
regressionplotcells <- regressionplotcells + theme(axis.title=element_text(size = 12), axis.text.y=element_blank(), legend.text = element_text(size = 12), axis.text.x = element_text(size = 12))
regressionplotcells <- regressionplotcells + guides(color = guide_legend(reverse=TRUE)) + theme(legend.title=element_blank())
regressionplotcells