#Figure 7b - Associations between mediators and COVID status via linear regression

#Clear R environment 
rm(list = ls())

#Import master dataset
combined <- read.csv(file.choose(), header=TRUE)
covidflusari <- subset(combined, min_no_sym==0 | is.na(min_no_sym))

#PIDs to fixed row identifiers 
rownames(covidflusari) = covidflusari$pid
covidflusari$pid = NULL

#Relabel pathogen code to binary 1/0
covidflusari[covidflusari$pathogencode == 1, "pathogencode2"] <- 1
covidflusari[covidflusari$pathogencode == 2 | covidflusari$pathogencode == 3 | covidflusari$pathogencode == 4, "pathogencode2"] <- 0

#Include 4 patients not known to be living with HIV (but with missing RDT) as with HIV neg
covidflusari$hivrdtresult[is.na(covidflusari$hivrdtresult)] = 0

#Linear regressions looking at relationship between COVID status and soluble mediators, adjusted for age, sex, HIV status, and WHO clinical severity
library(performance)
library(see)
library(patchwork)
library(jtools)

#sCD40L
lmscd40lcovidflusari = lm(logscd40l ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmscd40lcovidflusari)
check_model(lmscd40lcovidflusari)
summ(lmscd40lcovidflusari,  confint = TRUE, digits=3)

#EGF
lmegfcovidflusari = lm(logegf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmegfcovidflusari)
check_model(lmegfcovidflusari)
summ(lmegfcovidflusari,  confint = TRUE, digits=3)

#Eotaxin
lmeotaxincovidflusari = lm(logeotaxin ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmeotaxincovidflusari)
check_model(lmeotaxincovidflusari)
summ(lmeotaxincovidflusari,  confint = TRUE, digits=3)

#FGF-2
lmfgf2covidflusari = lm(logfgf2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmfgf2covidflusari)
check_model(lmfgf2covidflusari)
summ(lmfgf2covidflusari,  confint = TRUE, digits=3)

#FLT-3L
lmflt3lcovidflusari = lm(logflt3l ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmflt3lcovidflusari)
check_model(lmflt3lcovidflusari)
summ(lmflt3lcovidflusari,  confint = TRUE, digits=3)

#Fractalkine
lmfractalkinecovidflusari = lm(logfractalkine ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmfractalkinecovidflusari)
check_model(lmfractalkinecovidflusari)
summ(lmfractalkinecovidflusari,  confint = TRUE, digits=3)

#G-CSF
lmgcsfcovidflusari = lm(loggcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmgcsfcovidflusari)
check_model(lmgcsfcovidflusari)
summ(lmgcsfcovidflusari,  confint = TRUE, digits=3)

#GM-CSF
lmgmcsfcovidflusari = lm(loggmcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmgmcsfcovidflusari)
check_model(lmgmcsfcovidflusari)
summ(lmgmcsfcovidflusari,  confint = TRUE, digits=3)

#GRO-alpha/CXCL1
lmgroalphacovidflusari = lm(loggroalpha ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmgroalphacovidflusari)
check_model(lmgroalphacovidflusari)
summ(lmgroalphacovidflusari,  confint = TRUE, digits=3)

#IFN-a2
lmifna2covidflusari = lm(logifna2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmifna2covidflusari)
check_model(lmifna2covidflusari)
summ(lmifna2covidflusari,  confint = TRUE, digits=3)

#IFN-gamma
lmifngammacovidflusari = lm(logifngamma ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmifngammacovidflusari)
check_model(lmifngammacovidflusari)
summ(lmifngammacovidflusari,  confint = TRUE, digits=3)

#IL-1a
lmil1acovidflusari = lm(logil1a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil1acovidflusari)
check_model(lmil1acovidflusari)
summ(lmil1acovidflusari,  confint = TRUE, digits=3)

#IL-1b
lmil1bcovidflusari = lm(logil1b ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil1bcovidflusari)
check_model(lmil1bcovidflusari)
summ(lmil1bcovidflusari,  confint = TRUE, digits=3)

#IL-1Ra
lmil1racovidflusari = lm(logil1ra ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil1racovidflusari)
check_model(lmil1racovidflusari)
summ(lmil1racovidflusari,  confint = TRUE, digits=3)

#IL-2
lmil2covidflusari = lm(logil2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil2covidflusari)
check_model(lmil2covidflusari)
summ(lmil2covidflusari,  confint = TRUE, digits=3)

#IL-3
lmil3covidflusari = lm(logil3 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil3covidflusari)
check_model(lmil3covidflusari)
summ(lmil3covidflusari,  confint = TRUE, digits=3)

#IL-4
lmil4covidflusari = lm(logil4 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil4covidflusari)
check_model(lmil4covidflusari)
summ(lmil4covidflusari,  confint = TRUE, digits=3)

#IL-5
lmil5covidflusari = lm(logil5 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil5covidflusari)
check_model(lmil5covidflusari)
summ(lmil5covidflusari,  confint = TRUE, digits=3)

#IL-6
lmil6covidflusari = lm(logil6 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil6covidflusari)
check_model(lmil6covidflusari)
summ(lmil6covidflusari,  confint = TRUE, digits=3)

#IL-7
lmil7covidflusari = lm(logil7 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil7covidflusari)
check_model(lmil7covidflusari)
summ(lmil7covidflusari,  confint = TRUE, digits=3)

#IL-8
lmil8covidflusari = lm(logil8 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil8covidflusari)
check_model(lmil8covidflusari)
summ(lmil8covidflusari,  confint = TRUE, digits=3)

#IL-9
lmil9covidflusari = lm(logil9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil9covidflusari)
check_model(lmil9covidflusari)
summ(lmil9covidflusari,  confint = TRUE, digits=3)

#IL-10
lmil10covidflusari = lm(logil10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil10covidflusari)
check_model(lmil10covidflusari)
summ(lmil10covidflusari,  confint = TRUE, digits=3)

#IL-12p40
lmil12p40covidflusari = lm(logil12p40 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil12p40covidflusari)
check_model(lmil12p40covidflusari)
summ(lmil12p40covidflusari,  confint = TRUE, digits=3)

#IL-12p70
lmil12p70covidflusari = lm(logil12p70 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil12p70covidflusari)
check_model(lmil12p70covidflusari)
summ(lmil12p70covidflusari,  confint = TRUE, digits=3)

#IL-13
lmil13covidflusari = lm(logil13 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil13covidflusari)
check_model(lmil13covidflusari)
summ(lmil13covidflusari,  confint = TRUE, digits=3)

#IL-15
lmil15covidflusari = lm(logil15 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil15covidflusari)
check_model(lmil15covidflusari)
summ(lmil15covidflusari,  confint = TRUE, digits=3)

#IL-17A
lmil17acovidflusari = lm(logil17a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil17acovidflusari)
check_model(lmil17acovidflusari)
summ(lmil17acovidflusari,  confint = TRUE, digits=3)

#IL-17E/IL-25
lmil17eil25covidflusari = lm(logil17eil25 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil17eil25covidflusari)
check_model(lmil17eil25covidflusari)
summ(lmil17eil25covidflusari,  confint = TRUE, digits=3)

#IL-17F
lmil17fcovidflusari = lm(logil17f ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil17fcovidflusari)
check_model(lmil17fcovidflusari)
summ(lmil17fcovidflusari,  confint = TRUE, digits=3)

#IL-18
lmil18covidflusari = lm(logil18 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil18covidflusari)
check_model(lmil18covidflusari)
summ(lmil18covidflusari,  confint = TRUE, digits=3)

#IL-22
lmil22covidflusari = lm(logil22 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil22covidflusari)
check_model(lmil22covidflusari)
summ(lmil22covidflusari,  confint = TRUE, digits=3)

#IL-27
lmil27covidflusari = lm(logil27 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmil27covidflusari)
check_model(lmil27covidflusari)
summ(lmil27covidflusari,  confint = TRUE, digits=3)

#IP-10/CXCL10
lmip10covidflusari = lm(logip10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmip10covidflusari)
check_model(lmip10covidflusari)
summ(lmip10covidflusari,  confint = TRUE, digits=3)

#MCP-1
lmmcp1covidflusari = lm(logmcp1 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmcp1covidflusari)
check_model(lmmcp1covidflusari)
summ(lmmcp1covidflusari,  confint = TRUE, digits=3)

#MCP-3
lmmcp3covidflusari = lm(logmcp3 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmcp3covidflusari)
check_model(lmmcp3covidflusari)
summ(lmmcp3covidflusari,  confint = TRUE, digits=3)

#M-CSF
lmmcsfcovidflusari = lm(logmcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmcsfcovidflusari)
check_model(lmmcsfcovidflusari)
summ(lmmcsfcovidflusari,  confint = TRUE, digits=3)

#MDC
lmmdccovidflusari = lm(logmdc ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmdccovidflusari)
check_model(lmmdccovidflusari)
summ(lmmdccovidflusari,  confint = TRUE, digits=3)

#MIG/CXCL9
lmmigcxcl9covidflusari = lm(logmigcxcl9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmigcxcl9covidflusari)
check_model(lmmigcxcl9covidflusari)
summ(lmmigcxcl9covidflusari,  confint = TRUE, digits=3)

#MIP-1a/CCL3
lmmip1acovidflusari = lm(logmip1a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmip1acovidflusari)
check_model(lmmip1acovidflusari)
summ(lmmip1acovidflusari,  confint = TRUE, digits=3)

#MIP-1b/CCL4
lmmip1bcovidflusari = lm(logmip1b ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmmip1bcovidflusari)
check_model(lmmip1bcovidflusari)
summ(lmmip1bcovidflusari,  confint = TRUE, digits=3)

#PDGF-AA
lmpdgfaacovidflusari = lm(logpdgfaa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmpdgfaacovidflusari)
check_model(lmpdgfaacovidflusari)
summ(lmpdgfaacovidflusari,  confint = TRUE, digits=3)

#PDGF-AB/BB
lmpdgfabbbcovidflusari = lm(logpdgfabbb ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmpdgfabbbcovidflusari)
check_model(lmpdgfabbbcovidflusari)
summ(lmpdgfabbbcovidflusari,  confint = TRUE, digits=3)

#RANTES
lmrantescovidflusari = lm(lograntes ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmrantescovidflusari)
check_model(lmrantescovidflusari)
summ(lmrantescovidflusari,  confint = TRUE, digits=3)

#TGF-alpha
lmtgfacovidflusari = lm(logtgfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmtgfacovidflusari)
check_model(lmtgfacovidflusari)
summ(lmtgfacovidflusari,  confint = TRUE, digits=3)

#TNF-alpha
lmtnfacovidflusari = lm(logtnfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmtnfacovidflusari)
check_model(lmtnfacovidflusari)
summ(lmtnfacovidflusari,  confint = TRUE, digits=3)

#TNF-beta/LT-alpha
lmtnfbcovidflusari = lm(logtnfb ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmtnfbcovidflusari)
check_model(lmtnfbcovidflusari)
summ(lmtnfbcovidflusari,  confint = TRUE, digits=3)

#VEGF-A
lmvegfacovidflusari = lm(logvegfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) 
summary(lmvegfacovidflusari)
check_model(lmvegfacovidflusari)
summ(lmvegfacovidflusari,  confint = TRUE, digits=3)

#Generate table of all regression models to extract COVID coefficients, 95% CI, p-values
library(modelsummary)

allmodels <- list(
  "sCD40L" = lm(logscd40l ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "EGF" = lm(logegf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "Eotaxin" = lm(logeotaxin ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari), 
  "FGF2" = lm(logfgf2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "FLT-3L" = lm(logflt3l ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),  
  "Fractalkine" = lm(logfractalkine ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "G-CSF" = lm(loggcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "GM-CSF" = lm(loggmcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "GROalpha" = lm(loggroalpha ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IFN-a2" = lm(logifna2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IFN-gamma" = lm(logifngamma ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-1a" = lm(logil1a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-1b" = lm(logil1b ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-1ra" = lm(logil1ra ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-2" = lm(logil2 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-3" = lm(logil3 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-4" = lm(logil4 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-5" = lm(logil5 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-6" = lm(logil6 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-7" = lm(logil7 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-8" = lm(logil8 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-9" = lm(logil9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-10" = lm(logil10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-12p40" = lm(logil12p40 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-12p70" = lm(logil12p70 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-13" = lm(logil13 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-15" = lm(logil15 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-17A" = lm(logil17a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-17E/IL-25" = lm(logil17eil25 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari), 
  "IL-17F" = lm(logil17f ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-18" = lm(logil18 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-22" = lm(logil22 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-27" = lm(logil27 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "CXCL10" = lm(logip10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "MCP-1" = lm(logmcp1 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "MCP-3" = lm(logmcp3 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "M-CSF" = lm(logmcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "MDC" = lm(logmdc ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "CXCL9" = lm(logmigcxcl9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "MIP1a" = lm(logmip1a ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "MIP1b" = lm(logmip1b ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "PDGF-AA" = lm(logpdgfaa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "PDGF-AB/BB" = lm(logpdgfabbb ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),
  "RANTES/CCL5" = lm(lograntes ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "TGF-alpha" = lm(logtgfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),
  "TNF-alpha" = lm(logtnfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),
  "LT-alpha" = lm(logtnfb ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),
  "VEGF-A" = lm(logvegfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari))

modelsummary(allmodels, statistic = c('conf.int', "p.value"))

#Plot significant regression coefficients using forest plots
library(modelsummary)
library(ggplot2)
library(pals)

library(modelsummary)
models <- list(
  "EGF" = lm(logegf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "Eotaxin" = lm(logeotaxin ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari), 
  "FLT-3L" = lm(logflt3l ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),  
  "IL-6" = lm(logil6 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-8" = lm(logil8 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-9" = lm(logil9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-10" = lm(logil10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-17E/IL-25" = lm(logil17eil25 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari), 
  "IL-18" = lm(logil18 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "IL-27" = lm(logil27 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "CXCL10" = lm(logip10 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "M-CSF" = lm(logmcsf ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "CXCL9" = lm(logmigcxcl9 ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "PDGF-AA" = lm(logpdgfaa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
  "PDGF-AB/BB" = lm(logpdgfabbb ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari),
"RANTES/CCL5" = lm(lograntes ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari) ,
"TGF-\u03b1" = lm(logtgfa ~ pathogencode2 + age + sex + hivrdtresult + whoseverity_imp, data = covidflusari))

modelsummary(models, statistic = 'conf.int')
modelplot(models, coef_omit = c(1,3,4,5,6))

manualcolors2<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                 'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                 'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                 "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                 'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                 "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

b <- list(geom_vline(xintercept = 0, color = 'black', linetype=2))

regressionplot <- modelplot(models, background = b, coef_omit = c(1,3,4,5,6)) +
  labs(x = "Adjusted \u03B2 with 95% Confidence Interval", 
       y = "",
       title = "") + theme_bw() +
  scale_color_manual(values = manualcolors2) + theme(text=element_text(size=12)) 
regressionplot <- regressionplot + theme(axis.title=element_text(size=12), axis.text.y=element_blank(), legend.text = element_text(size = 12), axis.text.x = element_text(size = 12))
regressionplot <- regressionplot + guides(color = guide_legend(reverse=TRUE)) + theme(legend.title=element_blank())
regressionplot 

