#Figure 1f and S3 - Interactions between mediators and HIV for COVID severity

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients, WHO severity to factor
covid <- subset(combined, pathogencode == 1)
covid$whoseverity_imp <- as.factor(covid$whoseverity_imp)

#Relabel HIV variable
covid <- subset(covid, hivrdtresult==1 | hivrdtresult==0)
covid[covid$hivrdtresult==1, "PLWH"] <- "Yes"
covid[covid$hivrdtresult==0, "PLWH"] <- "No"
covid$PLWH <- as.factor(covid$PLWH)
levels(covid$PLWH)

#Biomarker concentrations across illness severity with HIV*mediator interaction adjusted for age/sex
library(ordinal)
library(sjPlot)
library(sjmisc)
library(ggplot2)

#sCD40L
scd40l = clm(whoseverity_imp ~ logscd40l*PLWH + age + sex, data = covid, link = "logit")
summary(scd40l)
confint(scd40l, type = "Wald")
exp(coef(scd40l))

scd40linteraction <- plot_model(scd40l, type = "int", ci.lvl=0.95)
scd40linteraction <- scd40linteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.188") + xlab("sCD40L log10 pg/ml") + ylab("")
scd40linteraction

#EGF
egf = clm(whoseverity_imp ~ logegf*PLWH + age + sex, data = covid, link = "logit")
summary(egf)
confint(egf, type = "Wald")
exp(coef(egf))

egfinteraction <- plot_model(egf, type = "int", ci.lvl=0.95)
egfinteraction <- egfinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.096") + xlab("EGF log10 pg/ml") + ylab("")
egfinteraction

#Eotaxin
eotaxin = clm(whoseverity_imp ~ logeotaxin*PLWH + age + sex, data = covid, link = "logit")
summary(eotaxin)
confint(eotaxin, type = "Wald")
exp(coef(eotaxin))

eotaxininteraction <- plot_model(eotaxin, type = "int", ci.lvl=0.95)
eotaxininteraction <- eotaxininteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.702") + xlab("Eotaxin log10 pg/ml") + ylab("")
eotaxininteraction

#FGF-2
fgf2 = clm(whoseverity_imp ~ logfgf2*PLWH + age + sex, data = covid, link = "logit")
summary(fgf2)
confint(fgf2, type = "Wald")
exp(coef(fgf2))

fgf2interaction <- plot_model(fgf2, type = "int", ci.lvl=0.95)
fgf2interaction <- fgf2interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.228") + xlab("FGF2 log10 pg/ml") + ylab("")
fgf2interaction

#FLT3L
flt3l = clm(whoseverity_imp ~ logflt3l*PLWH + age + sex, data = covid, link = "logit")
summary(flt3l)
confint(flt3l, type = "Wald")
exp(coef(flt3l))

flt3linteraction <- plot_model(flt3l, type = "int", ci.lvl=0.95)
flt3linteraction <- flt3linteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.174") + xlab("FLT3L log10 pg/ml") + ylab("")
flt3linteraction

#Fractalkine
fractalkine = clm(whoseverity_imp ~ logfractalkine*PLWH + age + sex, data = covid, link = "logit")
summary(fractalkine)
confint(fractalkine, type = "Wald")
exp(coef(fractalkine))

fractalkineinteraction <- plot_model(fractalkine, type = "int", ci.lvl=0.95)
fractalkineinteraction <- fractalkineinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.324") + xlab("Fractalkine log10 pg/ml") + ylab("")
fractalkineinteraction

#G-CSF
gcsf = clm(whoseverity_imp ~ loggcsf*PLWH + age + sex, data = covid, link = "logit")
summary(gcsf)
confint(gcsf, type = "Wald")
exp(coef(gcsf))

gcsfinteraction <- plot_model(gcsf, type = "int", ci.lvl=0.95)
gcsfinteraction <- gcsfinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.149") + xlab("G-CSF log10 pg/ml") + ylab("")
gcsfinteraction

#GM-CSF
gmcsf = clm(whoseverity_imp ~ loggmcsf*PLWH + age + sex, data = covid, link = "logit")
summary(gmcsf)
confint(gmcsf, type = "Wald")
exp(coef(gmcsf))

gmcsfinteraction <- plot_model(gmcsf, type = "int", ci.lvl=0.95)
gmcsfinteraction <- gmcsfinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.672") + xlab("GM-CSF log10 pg/ml") + ylab("")
gmcsfinteraction

#GRO-alpha/CXCL1
groalpha = clm(whoseverity_imp ~ loggroalpha*PLWH + age + sex, data = covid, link = "logit")
summary(groalpha)
confint(groalpha, type = "Wald")
exp(coef(groalpha))

groalphainteraction <- plot_model(groalpha, type = "int", ci.lvl=0.95)
groalphainteraction <- groalphainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.112") + xlab("CXCL1 log10 pg/ml") + ylab("")
groalphainteraction

#IFN-a2
ifna2 = clm(whoseverity_imp ~ logifna2*PLWH + age + sex, data = covid, link = "logit")
summary(ifna2)
confint(ifna2, type = "Wald")
exp(coef(ifna2))

ifna2interaction <- plot_model(ifna2, type = "int", ci.lvl=0.95)
ifna2interaction <- ifna2interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.042") + xlab("IFN-\u03b12 log10 pg/ml") + ylab("")
ifna2interaction

covid$whoseverity_imp <- sjmisc::rec(covid$whoseverity_imp, rec = "1=Asymptomatic;2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ logifna2*PLWH + age + sex, data = covid, link = "logit")
ifna2interaction2 <- plot_model(m2, type = "pred", terms = c("logifna2", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
ifna2interaction2 <- ifna2interaction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.042") + xlab("IFN-\u03b12 log10 pg/ml") + ylab("") + theme(strip.text.x = element_text(size = 12))
ifna2interaction2 

ifna2hiv = clm(whoseverity_imp ~ logifna2 + age + sex, data = subset(covid, PLWH==1), link = "logit")
summary(ifna2hiv)
confint(ifna2hiv, type = "Wald")
exp(coef(ifna2hiv))

ifna2hivneg = clm(whoseverity_imp ~ logifna2 + age + sex, data = subset(covid, PLWH==0), link = "logit")
summary(ifna2hivneg)
confint(ifna2hivneg, type = "Wald")
exp(coef(ifna2hivneg))

#IFN-gamma
ifngamma = clm(whoseverity_imp ~ logifngamma*PLWH + age + sex, data = covid, link = "logit")
summary(ifngamma)
confint(ifngamma, type = "Wald")
exp(coef(ifngamma))

ifngammainteraction <- plot_model(ifngamma, type = "int", ci.lvl=0.95)
ifngammainteraction <- ifngammainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.500") + xlab("IFN-gamma log10 pg/ml") + ylab("")
ifngammainteraction

#IL-1a
il1a = clm(whoseverity_imp ~ logil1a*PLWH + age + sex, data = covid, link = "logit")
summary(il1a)
confint(il1a, type = "Wald")
exp(coef(il1a))

il1ainteraction <- plot_model(il1a, type = "int", ci.lvl=0.95)
il1ainteraction <- il1ainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.056") + xlab("IL-1a log10 pg/ml") + ylab("")
il1ainteraction

covid$whoseverity_imp <- sjmisc::rec(covid$whoseverity_imp, rec = "1=Asymptomatic;2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ logil1a*PLWH + age + sex, data = covid, link = "logit")
il1ainteraction2 <- plot_model(m2, type = "pred", terms = c("logil1a", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
il1ainteraction2 <- il1ainteraction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.056") + xlab("IL-1a log10 pg/ml") + ylab("") + theme(strip.text.x = element_text(size = 12))
il1ainteraction2

il1ahiv = clm(whoseverity_imp ~ logil1a + age + sex, data = subset(covid, PLWH==1), link = "logit")
summary(il1ahiv)
confint(il1ahiv, type = "Wald")
exp(coef(il1ahiv))

il1ahivneg = clm(whoseverity_imp ~ logil1a + age + sex, data = subset(covid, PLWH==0), link = "logit")
summary(il1ahivneg)
confint(il1ahivneg, type = "Wald")
exp(coef(il1ahivneg))

#IL-1b
il1b = clm(whoseverity_imp ~ logil1b*PLWH + age + sex, data = covid, link = "logit")
summary(il1b)
confint(il1b, type = "Wald")
exp(coef(il1b))

il1binteraction <- plot_model(il1b, type = "int", ci.lvl=0.95)
il1binteraction <- il1binteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.177") + xlab("IL-1b log10 pg/ml") + ylab("")
il1binteraction

#IL-1Ra
il1ra = clm(whoseverity_imp ~ logil1ra*PLWH + age + sex, data = covid, link = "logit")
summary(il1ra)
confint(il1ra, type = "Wald")
exp(coef(il1ra))

il1rainteraction <- plot_model(il1ra, type = "int", ci.lvl=0.95)
il1rainteraction <- il1rainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.053") + xlab("IL-1Ra log10 pg/ml") + ylab("")
il1rainteraction

covid$whoseverity_imp <- sjmisc::rec(covid$whoseverity_imp, rec = "1=Asymptomatic;2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ logil1ra*PLWH + age + sex, data = covid, link = "logit")
il1rainteraction2 <- plot_model(m2, type = "pred", terms = c("logil1ra", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
il1rainteraction2 <- il1rainteraction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.053") + xlab("IL-1Ra log10 pg/ml") + ylab("") + theme(strip.text.x = element_text(size = 12))
il1rainteraction2

il1rahiv = clm(whoseverity_imp ~ logil1ra + age + sex, data = subset(covid, PLWH==1), link = "logit")
summary(il1rahiv)
confint(il1rahiv, type = "Wald")
exp(coef(il1rahiv))

il1rahivneg = clm(whoseverity_imp ~ logil1ra + age + sex, data = subset(covid, PLWH==0), link = "logit")
summary(il1rahivneg)
confint(il1rahivneg, type = "Wald")
exp(coef(il1rahivneg))

#IL-2
il2 = clm(whoseverity_imp ~ logil2*PLWH + age + sex, data = covid, link = "logit")
summary(il2)
confint(il2, type = "Wald")
exp(coef(il2))

il2interaction <- plot_model(il2, type = "int", ci.lvl=0.95)
il2interaction <- il2interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.208") + xlab("IL-2 log10 pg/ml") + ylab("")
il2interaction

#IL-3
il3 = clm(whoseverity_imp ~ logil3*PLWH + age + sex, data = covid, link = "logit")
summary(il3)
confint(il3, type = "Wald")
exp(coef(il3))

il3interaction <- plot_model(il3, type = "int", ci.lvl=0.95)
il3interaction <- il3interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.652") + xlab("IL-3 log10 pg/ml") + ylab("")
il3interaction

#IL-4
il4 = clm(whoseverity_imp ~ logil4*PLWH + age + sex, data = covid, link = "logit")
summary(il4)
confint(il4, type = "Wald")
exp(coef(il4))

il4interaction <- plot_model(il4, type = "int", ci.lvl=0.95)
il4interaction <- il4interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.191") + xlab("IL-4 log10 pg/ml") + ylab("")
il4interaction

#IL-5
il5 = clm(whoseverity_imp ~ logil5*PLWH + age + sex, data = covid, link = "logit")
summary(il5)
confint(il5, type = "Wald")
exp(coef(il5))

il5interaction <- plot_model(il5, type = "int", ci.lvl=0.95)
il5interaction <- il5interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.322") + xlab("IL-5 log10 pg/ml") + ylab("")
il5interaction

#IL-6
il6 = clm(whoseverity_imp ~ logil6*PLWH + age + sex, data = covid, link = "logit")
summary(il6)
confint(il6, type = "Wald")
exp(coef(il6))

il6interaction <- plot_model(il6, type = "int", ci.lvl=0.95)
il6interaction <- il6interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.035") + xlab("IL-6 log10 pg/ml") + ylab("")
il6interaction

covid$whoseverity_imp <- sjmisc::rec(covid$whoseverity_imp, rec = "1=Asymptomatic;2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ logil6*PLWH + age + sex, data = covid, link = "logit")
il6interaction2 <- plot_model(m2, type = "pred", terms = c("logil6", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
il6interaction2 <- il6interaction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.035") + xlab("IL-6 log10 pg/ml") + ylab("") + theme(strip.text.x = element_text(size = 12))
il6interaction2

il6hiv = clm(whoseverity_imp ~ logil6 + age + sex, data = subset(covid, PLWH==1), link = "logit")
summary(il6hiv)
confint(il6hiv, type = "Wald")
exp(coef(il6hiv))

il6hivneg = clm(whoseverity_imp ~ logil6 + age + sex, data = subset(covid, PLWH==0), link = "logit")
summary(il6hivneg)
confint(il6hivneg, type = "Wald")
exp(coef(il6hivneg))

#IL-7
il7 = clm(whoseverity_imp ~ logil7*PLWH + age + sex, data = covid, link = "logit")
summary(il7)
confint(il7, type = "Wald")
exp(coef(il7))

il7interaction <- plot_model(il7, type = "int", ci.lvl=0.95)
il7interaction <- il7interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.687") + xlab("IL-7 log10 pg/ml") + ylab("")
il7interaction

#IL-8
il8 = clm(whoseverity_imp ~ logil8*PLWH + age + sex, data = covid, link = "logit")
summary(il8)
confint(il8, type = "Wald")
exp(coef(il8))

il8interaction <- plot_model(il8, type = "int", ci.lvl=0.95)
il8interaction <- il8interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.475") + xlab("IL-1a log10 pg/ml") + ylab("")
il8interaction

#IL-9
il9 = clm(whoseverity_imp ~ logil9*PLWH + age + sex, data = covid, link = "logit")
summary(il9)
confint(il9, type = "Wald")
exp(coef(il9))

il9interaction <- plot_model(il9, type = "int", ci.lvl=0.95)
il9interaction <- il9interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.488") + xlab("IL-9 log10 pg/ml") + ylab("")
il9interaction

#IL-10
il10 = clm(whoseverity_imp ~ logil10*PLWH + age + sex, data = covid, link = "logit")
summary(il10)
confint(il10, type = "Wald")
exp(coef(il10))

il10interaction <- plot_model(il10, type = "int", ci.lvl=0.95)
il10interaction <- il10interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.461") + xlab("IL-10 log10 pg/ml") + ylab("")
il10interaction

#IL-12p40
il12p40 = clm(whoseverity_imp ~ logil12p40*PLWH + age + sex, data = covid, link = "logit")
summary(il12p40)
confint(il12p40, type = "Wald")
exp(coef(il12p40))

il12p40interaction <- plot_model(il12p40, type = "int", ci.lvl=0.95)
il12p40interaction <- il12p40interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.807") + xlab("IL-12p40 log10 pg/ml") + ylab("")
il12p40interaction

#IL-12p70
il12p70 = clm(whoseverity_imp ~ logil12p70*PLWH + age + sex, data = covid, link = "logit")
summary(il12p70)
confint(il12p70, type = "Wald")
exp(coef(il12p70))

il12p70interaction <- plot_model(il12p70, type = "int", ci.lvl=0.95)
il12p70interaction <- il12p70interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.227") + xlab("IL-12p70 log10 pg/ml") + ylab("")
il12p70interaction

#IL-13
il13 = clm(whoseverity_imp ~ logil13*PLWH + age + sex, data = covid, link = "logit")
summary(il13)
confint(il13, type = "Wald")
exp(coef(il13))

il13interaction <- plot_model(il13, type = "int", ci.lvl=0.95)
il13interaction <- il13interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.997") + xlab("IL-13 log10 pg/ml") + ylab("")
il13interaction

#IL-15
il15 = clm(whoseverity_imp ~ logil15*PLWH + age + sex, data = covid, link = "logit")
summary(il15)
confint(il15, type = "Wald")
exp(coef(il15))

il15interaction <- plot_model(il15, type = "int", ci.lvl=0.95)
il15interaction <- il15interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.365") + xlab("IL-15 log10 pg/ml") + ylab("")
il15interaction

#IL-17a
il17a = clm(whoseverity_imp ~ logil17a*PLWH + age + sex, data = covid, link = "logit")
summary(il17a)
confint(il17a, type = "Wald")
exp(coef(il17a))

il17ainteraction <- plot_model(il17a, type = "int", ci.lvl=0.95)
il17ainteraction <- il17ainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.203") + xlab("IL-17a log10 pg/ml") + ylab("")
il17ainteraction

#IL-17E/IL-25
il17eil25 = clm(whoseverity_imp ~ logil17eil25*PLWH + age + sex, data = covid, link = "logit")
summary(il17eil25)
confint(il17eil25, type = "Wald")
exp(coef(il17eil25))

il17eil25interaction <- plot_model(il17eil25, type = "int", ci.lvl=0.95)
il17eil25interaction <- il17eil25interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.281") + xlab("IL-17E/IL-25 log10 pg/ml") + ylab("")
il17eil25interaction

#IL-17F
il17f = clm(whoseverity_imp ~ logil17f*PLWH + age + sex, data = covid, link = "logit")
summary(il17f)
confint(il17f, type = "Wald")
exp(coef(il17f))

il17finteraction <- plot_model(il17f, type = "int", ci.lvl=0.95)
il17finteraction <- il17finteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.665") + xlab("IL-17F log10 pg/ml") + ylab("")
il17finteraction

#IL-18
il18 = clm(whoseverity_imp ~ logil18*PLWH + age + sex, data = covid, link = "logit")
summary(il18)
confint(il18, type = "Wald")
exp(coef(il18))

il18interaction <- plot_model(il18, type = "int", ci.lvl=0.95)
il18interaction <- il18interaction + theme_bw() + theme(plot.title=element_text(size=12))   + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.667") + xlab("IL-18 log10 pg/ml") + ylab("")
il18interaction

#IL-22
il22 = clm(whoseverity_imp ~ logil22*PLWH + age + sex, data = covid, link = "logit")
summary(il22)
confint(il22, type = "Wald")
exp(coef(il22))

il22interaction <- plot_model(il22, type = "int", ci.lvl=0.95)
il22interaction <- il22interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.378") + xlab("IL-22 log10 pg/ml") + ylab("")
il22interaction

#IL-27
il27 = clm(whoseverity_imp ~ logil27*PLWH + age + sex, data = covid, link = "logit")
summary(il27)
confint(il27, type = "Wald")
exp(coef(il27))

il27interaction <- plot_model(il27, type = "int", ci.lvl=0.95)
il27interaction <- il27interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.912") + xlab("IL-27 log10 pg/ml") + ylab("")
il27interaction

#IP-10/CXCL10
ip10 = clm(whoseverity_imp ~ logip10*PLWH + age + sex, data = covid, link = "logit")
summary(ip10)
confint(ip10, type = "Wald")
exp(coef(ip10))

ip10interaction <- plot_model(ip10, type = "int", ci.lvl=0.95)
ip10interaction <- ip10interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.761") + xlab("CXCL10 log10 pg/ml") + ylab("")
ip10interaction

#MCP-1
mcp1 = clm(whoseverity_imp ~ logmcp1*PLWH + age + sex, data = covid, link = "logit")
summary(mcp1)
confint(mcp1, type = "Wald")
exp(coef(mcp1))

mcp1interaction <- plot_model(mcp1, type = "int", ci.lvl=0.95)
mcp1interaction <- mcp1interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.398") + xlab("MCP-1 log10 pg/ml") + ylab("")
mcp1interaction

#MCP-3
mcp3 = clm(whoseverity_imp ~ logmcp3*PLWH + age + sex, data = covid, link = "logit")
summary(mcp3)
confint(mcp3, type = "Wald")
exp(coef(mcp3))

mcp3interaction <- plot_model(mcp3, type = "int", ci.lvl=0.95)
mcp3interaction <- mcp3interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.834") + xlab("MCP-3 log10 pg/ml") + ylab("")
mcp3interaction

#M-CSF
mcsf = clm(whoseverity_imp ~ logmcsf*PLWH + age + sex, data = covid, link = "logit")
summary(mcsf)
confint(mcsf, type = "Wald")
exp(coef(mcsf))

mcsfinteraction <- plot_model(mcsf, type = "int", ci.lvl=0.95)
mcsfinteraction <- mcsfinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.364") + xlab("M-CSF log10 pg/ml") + ylab("")
mcsfinteraction

#MDC
mdc = clm(whoseverity_imp ~ logmdc*PLWH + age + sex, data = covid, link = "logit")
summary(mdc)
confint(mdc, type = "Wald")
exp(coef(mdc))

mdcinteraction <- plot_model(mdc, type = "int", ci.lvl=0.95)
mdcinteraction <- mdcinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.155") + xlab("MDC log10 pg/ml") + ylab("")
mdcinteraction

#MIG/CXCL9
migcxcl9 = clm(whoseverity_imp ~ logmigcxcl9*PLWH + age + sex, data = covid, link = "logit")
summary(migcxcl9)
confint(migcxcl9, type = "Wald")
exp(coef(migcxcl9))

migcxcl9interaction <- plot_model(migcxcl9, type = "int", ci.lvl=0.95)
migcxcl9interaction <- migcxcl9interaction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.735") + xlab("MIG/CXCL9 log10 pg/ml") + ylab("")
migcxcl9interaction

#MIP1a/CCL3
mip1a = clm(whoseverity_imp ~ logmip1a*PLWH + age + sex, data = covid, link = "logit")
summary(mip1a)
confint(mip1a, type = "Wald")
exp(coef(mip1a))

mip1ainteraction <- plot_model(mip1a, type = "int", ci.lvl=0.95)
mip1ainteraction <- mip1ainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.009") + xlab("MIP-1a/CCL3 log10 pg/ml") + ylab("")
mip1ainteraction

covid$whoseverity_imp <- sjmisc::rec(covid$whoseverity_imp, rec = "1=Asymptomatic;2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ logmip1a*PLWH + age + sex, data = covid, link = "logit")
mip1ainteraction2 <- plot_model(m2, type = "pred", terms = c("logmip1a", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
mip1ainteraction2 <- mip1ainteraction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.009") + xlab("CCL3 log10 pg/ml") + ylab("") + theme(strip.text.x = element_text(size = 12))
mip1ainteraction2

mip1ahiv = clm(whoseverity_imp ~ logmip1a + age + sex, data = subset(covid, PLWH==1), link = "logit")
summary(mip1ahiv)
confint(mip1ahiv, type = "Wald")
exp(coef(mip1ahiv))

mip1ahivneg = clm(whoseverity_imp ~ logmip1a + age + sex, data = subset(covid, PLWH==0), link = "logit")
summary(mip1ahivneg)
confint(mip1ahivneg, type = "Wald")
exp(coef(mip1ahivneg))

#MIP1b/CCL4
mip1b = clm(whoseverity_imp ~ logmip1b*PLWH + age + sex, data = covid, link = "logit")
summary(mip1b)
confint(mip1b, type = "Wald")
exp(coef(mip1b))

mip1binteraction <- plot_model(mip1b, type = "int", ci.lvl=0.95)
mip1binteraction <- mip1binteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.068") + xlab("MIP-1b/CCL4 log10 pg/ml") + ylab("")
mip1binteraction

#PDGF-AA
pdgfaa = clm(whoseverity_imp ~ logpdgfaa*PLWH + age + sex, data = covid, link = "logit")
summary(pdgfaa)
confint(pdgfaa, type = "Wald")
exp(coef(pdgfaa))

pdgfaainteraction <- plot_model(pdgfaa, type = "int", ci.lvl=0.95)
pdgfaainteraction <- pdgfaainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.104") + xlab("PDGF-AA log10 pg/ml") + ylab("")
pdgfaainteraction

#PDGF-AB/BB
pdgfabbb = clm(whoseverity_imp ~ logpdgfabbb*PLWH + age + sex, data = covid, link = "logit")
summary(pdgfabbb)
confint(pdgfabbb, type = "Wald")
exp(coef(pdgfabbb))

pdgfabbbinteraction <- plot_model(pdgfabbb, type = "int", ci.lvl=0.95)
pdgfabbbinteraction <- pdgfabbbinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.129") + xlab("PDGF-AB/BB log10 pg/ml") + ylab("")
pdgfabbbinteraction

#RANTES
rantes = clm(whoseverity_imp ~ lograntes*PLWH + age + sex, data = covid, link = "logit")
summary(rantes)
confint(rantes, type = "Wald")
exp(coef(rantes))

rantesinteraction <- plot_model(rantes, type = "int", ci.lvl=0.95)
rantesinteraction <- rantesinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.350") + xlab("RANTES/CCL5 log10 pg/ml") + ylab("")
rantesinteraction

#TGF-alpha
tgfa = clm(whoseverity_imp ~ logtgfa*PLWH + age + sex, data = covid, link = "logit")
summary(tgfa)
confint(tgfa, type = "Wald")
exp(coef(tgfa))

tgfainteraction <- plot_model(tgfa, type = "int", ci.lvl=0.95)
tgfainteraction <- tgfainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.822") + xlab("TGF-a log10 pg/ml") + ylab("")
tgfainteraction

#TNF-alpha
tnfa = clm(whoseverity_imp ~ logtnfa*PLWH + age + sex, data = covid, link = "logit")
summary(tnfa)
confint(tnfa, type = "Wald")
exp(coef(tnfa))

tnfainteraction <- plot_model(tnfa, type = "int", ci.lvl=0.95)
tnfainteraction <- tnfainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.689") + xlab("TNF-a log10 pg/ml") + ylab("")
tnfainteraction

#TNF-beta/LT-alpha
tnfb = clm(whoseverity_imp ~ logtnfb*PLWH + age + sex, data = covid, link = "logit")
summary(tnfb)
confint(tnfb, type = "Wald")
exp(coef(tnfb))

tnfbinteraction <- plot_model(tnfb, type = "int", ci.lvl=0.95)
tnfbinteraction <- tnfbinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.836") + xlab("Lymphotoxin-a log10 pg/ml") + ylab("")
tnfbinteraction

#VEGF-A
vegfa = clm(whoseverity_imp ~ logvegfa*PLWH + age + sex, data = covid, link = "logit")
summary(vegfa)
confint(vegfa, type = "Wald")
exp(coef(vegfa))

vegfainteraction <- plot_model(vegfa, type = "int", ci.lvl=0.95)
vegfainteraction <- vegfainteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.204") + xlab("VEGF-A log10 pg/ml") + ylab("")
vegfainteraction