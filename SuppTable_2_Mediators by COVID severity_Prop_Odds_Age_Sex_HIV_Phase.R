#Association between mediators and COVID severity via proportional odds regression

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, make row ID fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients 
covid <- subset(combined, pathogencode == 1)

#Include 4 patients not known to be living with HIV but with missing RDT as HIV neg
covid$hivrdtresult[is.na(covid$hivrdtresult)] = 0

#SARS-CoV-2 variant phase to factor
covid$phase3 <- as.factor(covid$phase3)

#Biomarker concentrations across illness severity groups in multivariable proportional odds model adjusted for age, sex, HIv co-infection, and SARS-CoV-2 variant phase
library(ordinal)

#sCD40L
scd40l = clm(as.factor(whoseverity_imp) ~ logscd40l + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(scd40l)
confint(scd40l, type = "Wald")
exp(coef(scd40l))

#EGF
egf = clm(as.factor(whoseverity_imp) ~ logegf + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(egf)
confint(egf, type = "Wald")
exp(coef(egf))

#Eotaxin
eotaxin = clm(as.factor(whoseverity_imp) ~ logeotaxin + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(eotaxin)
confint(eotaxin, type = "Wald")
exp(coef(eotaxin))

#FGF-2
fgf2 = clm(as.factor(whoseverity_imp) ~ logfgf2 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(fgf2)
confint(fgf2, type = "Wald")
exp(coef(fgf2))

#FLT3L
flt3l = clm(as.factor(whoseverity_imp) ~ logflt3l + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(flt3l)
confint(flt3l, type = "Wald")
exp(coef(flt3l))

#Fractalkine
fractalkine = clm(as.factor(whoseverity_imp) ~ logfractalkine + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(fractalkine)
confint(fractalkine, type = "Wald")
exp(coef(fractalkine))

#G-CSF
gcsf = clm(as.factor(whoseverity_imp) ~ loggcsf + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(gcsf)
confint(gcsf, type = "Wald")
exp(coef(gcsf))

#GM-CSF
gmcsf = clm(as.factor(whoseverity_imp) ~ loggmcsf + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(gmcsf)
confint(gmcsf, type = "Wald")
exp(coef(gmcsf))

#GRO-alpha/CXCL1
groalpha = clm(as.factor(whoseverity_imp) ~ loggroalpha + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(groalpha)
confint(groalpha, type = "Wald")
exp(coef(groalpha))

#IFN-a2
ifna2 = clm(as.factor(whoseverity_imp) ~ logifna2 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(ifna2)
confint(ifna2, type = "Wald")
exp(coef(ifna2))

#IFN-gamma
ifngamma = clm(as.factor(whoseverity_imp) ~ logifngamma + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(ifngamma)
confint(ifngamma, type = "Wald")
exp(coef(ifngamma))

#IL-1a
il1a = clm(as.factor(whoseverity_imp) ~ logil1a + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il1a)
confint(il1a, type = "Wald")
exp(coef(il1a))

#IL-1b
il1b = clm(as.factor(whoseverity_imp) ~ logil1b + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il1b)
confint(il1b, type = "Wald")
exp(coef(il1b))

#IL-1Ra
il1ra = clm(as.factor(whoseverity_imp) ~ logil1ra + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il1ra)
confint(il1ra, type = "Wald")
exp(coef(il1ra))

#IL-2
il2 = clm(as.factor(whoseverity_imp) ~ logil2 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il2)
confint(il2, type = "Wald")
exp(coef(il2))

#IL-3
il3 = clm(as.factor(whoseverity_imp) ~ logil3 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il3)
confint(il3, type = "Wald")
exp(coef(il3))

#IL-4
il4 = clm(as.factor(whoseverity_imp) ~ logil4 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il4)
confint(il4, type = "Wald")
exp(coef(il4))

#IL-5
il5 = clm(as.factor(whoseverity_imp) ~ logil5 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il5)
confint(il5, type = "Wald")
exp(coef(il5))

#IL-6
il6 = clm(as.factor(whoseverity_imp) ~ logil6 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il6)
confint(il6, type = "Wald")
exp(coef(il6))

#IL-7
il7 = clm(as.factor(whoseverity_imp) ~ logil7 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il7)
confint(il7, type = "Wald")
exp(coef(il7))

#IL-8
il8 = clm(as.factor(whoseverity_imp) ~ logil8 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il8)
confint(il8, type = "Wald")
exp(coef(il8))

#IL-9
il9 = clm(as.factor(whoseverity_imp) ~ logil9 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il9)
confint(il9, type = "Wald")
exp(coef(il9))

#IL-10
il10 = clm(as.factor(whoseverity_imp) ~ logil10 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il10)
confint(il10, type = "Wald")
exp(coef(il10))

#IL-12p40
il12p40 = clm(as.factor(whoseverity_imp) ~ logil12p40 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il12p40)
confint(il12p40, type = "Wald")
exp(coef(il12p40))

#IL-12p70
il12p70 = clm(as.factor(whoseverity_imp) ~ logil12p70 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il12p70)
confint(il12p70, type = "Wald")
exp(coef(il12p70))

#IL-13
il13 = clm(as.factor(whoseverity_imp) ~ logil13 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il13)
confint(il13, type = "Wald")
exp(coef(il13))

#IL-15
il15 = clm(as.factor(whoseverity_imp) ~ logil15 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il15)
confint(il15, type = "Wald")
exp(coef(il15))

#IL-17a
il17a = clm(as.factor(whoseverity_imp) ~ logil17a + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il17a)
confint(il17a, type = "Wald")
exp(coef(il17a))

#IL-17E/IL-25
il17eil25 = clm(as.factor(whoseverity_imp) ~ logil17eil25 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il17eil25)
confint(il17eil25, type = "Wald")
exp(coef(il17eil25))

#IL-17F
il17f = clm(as.factor(whoseverity_imp) ~ logil17f + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il17f)
confint(il17f, type = "Wald")
exp(coef(il17f))

#IL-18
il18 = clm(as.factor(whoseverity_imp) ~ logil18 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il18)
confint(il18, type = "Wald")
exp(coef(il18))

#IL-22
il22 = clm(as.factor(whoseverity_imp) ~ logil22 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il22)
confint(il22, type = "Wald")
exp(coef(il22))

#IL-27
il27 = clm(as.factor(whoseverity_imp) ~ logil27 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(il27)
confint(il27, type = "Wald")
exp(coef(il27))

#IP-10/CXCL10
ip10 = clm(as.factor(whoseverity_imp) ~ logip10 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(ip10)
confint(ip10, type = "Wald")
exp(coef(ip10))

#MCP-1
mcp1 = clm(as.factor(whoseverity_imp) ~ logmcp1 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mcp1)
confint(mcp1, type = "Wald")
exp(coef(mcp1))

#MCP-3
mcp3 = clm(as.factor(whoseverity_imp) ~ logmcp3 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mcp3)
confint(mcp3, type = "Wald")
exp(coef(mcp3))

#M-CSF
mcsf = clm(as.factor(whoseverity_imp) ~ logmcsf + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mcsf)
confint(mcsf, type = "Wald")
exp(coef(mcsf))

#MDC
mdc = clm(as.factor(whoseverity_imp) ~ logmdc + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mdc)
confint(mdc, type = "Wald")
exp(coef(mdc))

#MIG/CXCL9
migcxcl9 = clm(as.factor(whoseverity_imp) ~ logmigcxcl9 + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(migcxcl9)
confint(migcxcl9, type = "Wald")
exp(coef(migcxcl9))

#MIP1a/CCL3
mip1a = clm(as.factor(whoseverity_imp) ~ logmip1a + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mip1a)
confint(mip1a, type = "Wald")
exp(coef(mip1a))

#MIP1b/CCL4
mip1b = clm(as.factor(whoseverity_imp) ~ logmip1b + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(mip1b)
confint(mip1b, type = "Wald")
exp(coef(mip1b))

#PDGF-AA
pdgfaa = clm(as.factor(whoseverity_imp) ~ logpdgfaa + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(pdgfaa)
confint(pdgfaa, type = "Wald")
exp(coef(pdgfaa))

#PDGF-AB/BB
pdgfabbb = clm(as.factor(whoseverity_imp) ~ logpdgfabbb + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(pdgfabbb)
confint(pdgfabbb, type = "Wald")
exp(coef(pdgfabbb))

#RANTES
rantes = clm(as.factor(whoseverity_imp) ~ lograntes + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(rantes)
confint(rantes, type = "Wald")
exp(coef(rantes))

#TGF-alpha
tgfa = clm(as.factor(whoseverity_imp) ~ logtgfa + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(tgfa)
confint(tgfa, type = "Wald")
exp(coef(tgfa))

#TNF-alpha
tnfa = clm(as.factor(whoseverity_imp) ~ logtnfa + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(tnfa)
confint(tnfa, type = "Wald")
exp(coef(tnfa))

#TNF-beta/LT-alpha
tnfb = clm(as.factor(whoseverity_imp) ~ logtnfb + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(tnfb)
confint(tnfb, type = "Wald")
exp(coef(tnfb))

#VEGF-A
vegfa = clm(as.factor(whoseverity_imp) ~ logvegfa + age + sex + hivrdtresult + phase3 , data = covid, link = "logit")
summary(vegfa)
confint(vegfa, type = "Wald")
exp(coef(vegfa))


