#Figure 6a - Mediators in PLWH vs. those without via linear regression

#Clear R environment 
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid patients
covid <- subset(combined, pathogencode == 1)

#Linear regressions for relationship between biomarkers and HIV status, adjusted for age, sex, WHO clinical severity
library(performance)
library(see)
library(patchwork)
library(sjPlot)
library(jtools)

#sCD40L
lmscd40lcovid = lm(logscd40l ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmscd40lcovid)
check_model(lmscd40lcovid)
summ(lmscd40lcovid, confint = TRUE, digits=3)

#EGF
lmegfcovid = lm(logegf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmegfcovid)
check_model(lmegfcovid)
summ(lmegfcovid, confint = TRUE, digits=3)

#Eotaxin
lmeotaxincovid = lm(logeotaxin ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmeotaxincovid)
check_model(lmeotaxincovid)
summ(lmeotaxincovid, confint = TRUE, digits=3)

#FGF-2
lmfgf2covid = lm(logfgf2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmfgf2covid)
check_model(lmfgf2covid)
summ(lmfgf2covid, confint = TRUE, digits=3)

#FLT-3L
lmflt3lcovid = lm(logflt3l ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmflt3lcovid)
check_model(lmflt3lcovid)
summ(lmflt3lcovid, confint = TRUE, digits=3)

#Fractalkine
lmfractalkinecovid = lm(logfractalkine ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmfractalkinecovid)
check_model(lmfractalkinecovid)
summ(lmfractalkinecovid, confint = TRUE, digits=3)

#G-CSF
lmgcsfcovid = lm(loggcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmgcsfcovid)
check_model(lmgcsfcovid)
summ(lmgcsfcovid, confint = TRUE, digits=3)

#GM-CSF
lmgmcsfcovid = lm(loggmcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmgmcsfcovid)
check_model(lmgmcsfcovid)
summ(lmgmcsfcovid, confint = TRUE, digits=3)

#GRO-alpha/CXCL1
lmgroalphacovid = lm(loggroalpha ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmgroalphacovid)
check_model(lmgroalphacovid)
summ(lmgroalphacovid, confint = TRUE, digits=3)

#IFN-a2
lmifna2covid = lm(logifna2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmifna2covid)
check_model(lmifna2covid)
summ(lmifna2covid, confint = TRUE, digits=3)

#IFN-gamma
lmifngammacovid = lm(logifngamma ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmifngammacovid)
check_model(lmifngammacovid)
summ(lmifngammacovid, confint = TRUE, digits=3)

#IL-1a
lmil1acovid = lm(logil1a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil1acovid)
check_model(lmil1acovid)
summ(lmil1acovid, confint = TRUE, digits=3)

#IL-1b
lmil1bcovid = lm(logil1b ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil1bcovid)
check_model(lmil1bcovid)
summ(lmil1bcovid, confint = TRUE, digits=3)

#IL-1Ra
lmil1racovid = lm(logil1ra ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil1racovid)
check_model(lmil1racovid)
summ(lmil1racovid, confint = TRUE, digits=3)

#IL-2
lmil2covid = lm(logil2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil2covid)
check_model(lmil2covid)
summ(lmil2covid, confint = TRUE, digits=3)

#IL-3
lmil3covid = lm(logil3 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil3covid)
check_model(lmil3covid)
summ(lmil3covid, confint = TRUE, digits=3)

#IL-4
lmil4covid = lm(logil4 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil4covid)
check_model(lmil4covid)
summ(lmil4covid, confint = TRUE, digits=3)

#IL-5
lmil5covid = lm(logil5 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil5covid)
check_model(lmil5covid)
summ(lmil5covid, confint = TRUE, digits=3)

#IL-6
lmil6covid = lm(logil6 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil6covid)
check_model(lmil6covid)
summ(lmil6covid, confint = TRUE, digits=3)

#IL-7
lmil7covid = lm(logil7 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil7covid)
check_model(lmil7covid)
summ(lmil7covid, confint = TRUE, digits=3)

#IL-8
lmil8covid = lm(logil8 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil8covid)
check_model(lmil8covid)
summ(lmil8covid, confint = TRUE, digits=3)

#IL-9
lmil9covid = lm(logil9 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil9covid)
check_model(lmil9covid)
summ(lmil9covid, confint = TRUE, digits=3)

#IL-10
lmil10covid = lm(logil10 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil10covid)
check_model(lmil10covid)
summ(lmil10covid, confint = TRUE, digits=3)

#IL-12p40
lmil12p40covid = lm(logil12p40 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil12p40covid)
check_model(lmil12p40covid)
summ(lmil12p40covid, confint = TRUE, digits=3)

#IL-12p70
lmil12p70covid = lm(logil12p70 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil12p70covid)
check_model(lmil12p70covid)
summ(lmil12p70covid, confint = TRUE, digits=3)

#IL-13
lmil13covid = lm(logil13 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil13covid)
check_model(lmil13covid)
summ(lmil13covid, confint = TRUE, digits=3)

#IL-15
lmil15covid = lm(logil15 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil15covid)
check_model(lmil15covid)
summ(lmil15covid, confint = TRUE, digits=3)

#IL-17A
lmil17acovid = lm(logil17a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil17acovid)
check_model(lmil17acovid)
summ(lmil17acovid, confint = TRUE, digits=3)

#IL-17E/IL-25
lmil17eil25covid = lm(logil17eil25 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil17eil25covid)
check_model(lmil17eil25covid)
summ(lmil17eil25covid, confint = TRUE, digits=3)

#IL-17F
lmil17fcovid = lm(logil17f ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil17fcovid)
check_model(lmil17fcovid)
summ(lmil17fcovid, confint = TRUE, digits=3)

#IL-18
lmil18covid = lm(logil18 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil18covid)
check_model(lmil18covid)
summ(lmil18covid, confint = TRUE, digits=3)

#IL-22
lmil22covid = lm(logil22 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil22covid)
check_model(lmil22covid)
summ(lmil22covid, confint = TRUE, digits=3)

#IL-27
lmil27covid = lm(logil27 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmil27covid)
check_model(lmil27covid)
summ(lmil27covid, confint = TRUE, digits=3)

#IP-10/CXCL10
lmip10covid = lm(logip10 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmip10covid)
check_model(lmip10covid)
summ(lmip10covid, confint = TRUE, digits=3)

#MCP-1
lmmcp1covid = lm(logmcp1 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmcp1covid)
check_model(lmmcp1covid)
summ(lmmcp1covid, confint = TRUE, digits=3)

#MCP-3
lmmcp3covid = lm(logmcp3 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmcp3covid)
check_model(lmmcp3covid)
summ(lmmcp3covid, confint = TRUE, digits=3)

#M-CSF
lmmcsfcovid = lm(logmcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmcsfcovid)
check_model(lmmcsfcovid)
summ(lmmcsfcovid, confint = TRUE, digits=3)

#MDC
lmmdccovid = lm(logmdc ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmdccovid)
check_model(lmmdccovid)
summ(lmmdccovid, confint = TRUE, digits=3)

#MIG/CXCL9
lmmigcxcl9covid = lm(logmigcxcl9 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmigcxcl9covid)
check_model(lmmigcxcl9covid)
summ(lmmigcxcl9covid, confint = TRUE, digits=3)

#MIP-1a/CCL3
lmmip1acovid = lm(logmip1a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmip1acovid)
check_model(lmmip1acovid)
summ(lmmip1acovid, confint = TRUE, digits=3)

#MIP-1b/CCL4
lmmip1bcovid = lm(logmip1b ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmmip1bcovid)
check_model(lmmip1bcovid)
summ(lmmip1bcovid, confint = TRUE, digits=3)

#PDGF-AA
lmpdgfaacovid = lm(logpdgfaa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmpdgfaacovid)
check_model(lmpdgfaacovid)
summ(lmpdgfaacovid, confint = TRUE, digits=3)

#PDGF-AB/BB
lmpdgfabbbcovid = lm(logpdgfabbb ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmpdgfabbbcovid)
check_model(lmpdgfabbbcovid)
summ(lmpdgfabbbcovid, confint = TRUE, digits=3)

#RANTES
lmrantescovid = lm(lograntes ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmrantescovid)
check_model(lmrantescovid)
summ(lmrantescovid, confint = TRUE, digits=3)

#TGF-alpha
lmtgfacovid = lm(logtgfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmtgfacovid)
check_model(lmtgfacovid)
summ(lmtgfacovid, confint = TRUE, digits=3)

#TNF-alpha
lmtnfacovid = lm(logtnfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmtnfacovid)
check_model(lmtnfacovid)
summ(lmtnfacovid, confint = TRUE, digits=3)

#TNF-beta/LT-alpha
lmtnfbcovid = lm(logtnfb ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmtnfbcovid)
check_model(lmtnfbcovid)
summ(lmtnfbcovid, confint = TRUE, digits=3)

#VEGF-A
lmvegfacovid = lm(logvegfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) 
summary(lmvegfacovid)
check_model(lmvegfacovid)
summ(lmvegfacovid, confint = TRUE, digits=3)

#Plot these coefficients using forest plots
library(jtools)
library(ggplot2)
library(ggstance)
library(broom.mixed)
library(RColorBrewer)
library(jcolors)

#Generate table of all regression models to extract HIV coefficients, 95% CI, p-values
library(modelsummary)

allmodels <- list(
  "sCD40L" = lm(logscd40l ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "EGF" = lm(logegf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "Eotaxin" = lm(logeotaxin ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid), 
  "FGF2" = lm(logfgf2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "FLT-3L" = lm(logflt3l ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),  
  "Fractalkine" = lm(logfractalkine ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "G-CSF" = lm(loggcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "GM-CSF" = lm(loggmcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "GROalpha" = lm(loggroalpha ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IFN-a2" = lm(logifna2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IFN-gamma" = lm(logifngamma ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-1a" = lm(logil1a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-1b" = lm(logil1b ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-1ra" = lm(logil1ra ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-2" = lm(logil2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-3" = lm(logil3 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-4" = lm(logil4 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-5" = lm(logil5 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-6" = lm(logil6 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-7" = lm(logil7 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-8" = lm(logil8 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-9" = lm(logil9 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-10" = lm(logil10 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-12p40" = lm(logil12p40 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-12p70" = lm(logil12p70 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-13" = lm(logil13 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-15" = lm(logil15 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-17A" = lm(logil17a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-17E/IL-25" = lm(logil17eil25 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid), 
  "IL-17F" = lm(logil17f ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-18" = lm(logil18 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-22" = lm(logil22 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "IL-27" = lm(logil27 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "CXCL10" = lm(logip10 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "MCP-1" = lm(logmcp1 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "MCP-3" = lm(logmcp3 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "M-CSF" = lm(logmcsf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "MDC" = lm(logmdc ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "CXCL9" = lm(logmigcxcl9 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "MIP1a" = lm(logmip1a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "MIP1b" = lm(logmip1b ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "PDGF-AA" = lm(logpdgfaa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "PDGF-AB/BB" = lm(logpdgfabbb ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "RANTES/CCL5" = lm(lograntes ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid) ,
  "TGF-alpha" = lm(logtgfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "TNF-alpha" = lm(logtnfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "LT-alpha" = lm(logtnfb ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "VEGF-A" = lm(logvegfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid))

modelsummary(allmodels, statistic = c('conf.int', "p.value"))

#Plot significant coefficients
library(modelsummary)
models <- list(
  "EGF" = lm(logegf ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "FGF2" = lm(logfgf2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "CXCL1" = lm(loggroalpha ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IFN-\u03b12" = lm(logifna2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-1\u03b2" = lm(logil1b ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-1Ra" = lm(logil1ra ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-2" = lm(logil2 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-3" = lm(logil3 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-5" = lm(logil5 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-12p40" = lm(logil12p40 ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "IL-17A" = lm(logil17a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "CCL3" = lm(logmip1a ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "PDGF-AB/BB" = lm(logpdgfabbb ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid),
  "TGF-\u03b1" = lm(logtgfa ~ hivrdtresult + age + sex  + whoseverity_imp, data = covid))

modelsummary(models, statistic = 'conf.int')
modelplot(models, coef_omit = c(1,3,4,5))

library(ggplot2)
library(pals)
library(pals)

manualcolors <- c("#F3DF6C", "#CEAB07","#D5D5D3","#24281A","#798E87","#C27D38","#CCC591","#29211F",
                  "#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B", "#D8B70A", "#02401B", "#A2A475", 
                  "#81A88D", "#972D15", "#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", 
                  "#D8A499", "#7294D4", "#9986A5", "#79402E", "#CCBA72") 

manualcolors2<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                 'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                 'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                 "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                 'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                 "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

b <- list(geom_vline(xintercept = 0, color = 'black', linetype=2))

regressionplotbiomarkers <- modelplot(models, background = b, coef_omit = c(1,3,4,5)) +
  labs(x = "Adjusted \u03B2 with 95% Confidence Interval", 
       y = "",
       title = "") + theme_bw() +
  scale_color_manual(values = manualcolors2) + theme(text=element_text(size=12)) 
regressionplotbiomarkers <- regressionplotbiomarkers + theme(axis.title=element_text(size=12), axis.text.y=element_blank(), legend.text = element_text(size = 12), axis.text.x = element_text(size = 12))
regressionplotbiomarkers <- regressionplotbiomarkers + guides(color = guide_legend(reverse=TRUE)) + theme(legend.title=element_blank())
regressionplotbiomarkers

#Combine with GSEA plot
library(ggpubr)
hivbiomarkergseaplot <- ggarrange(regressionplotbiomarkers, immunometabolichivbarplot, ncol = 1, nrow=2, labels = c("a", "b"))
hivbiomarkergseaplot











