#Figures 2e and 2f - Correlation matrices - immune cells/clinical variables and immune cells/biomarkers 

#Clear R environment 
rm(list = ls())

#Import master data, PID to fixed row identifier 
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid patients with RNAseq data, exclude asymptomatic patients
covid <- subset(combined, pathogencode==1)
covid <- subset(covid, min_no_sym==0)
covid <- subset(covid, rna_seq==1)

#Correlation plot - immune cell populations and clinical variables
covidsympbiom <- data.matrix(subset(covid)[1:84, c(6,7,286,288, 290:297,301:307, 311,313,318,319,321,322)], rownames.force = NA)
head(covidsympbiom)
colnames(covidsympbiom) <- c("Sex",
                             "Age",
                             "KPS",
                             "Temp.",
                             "HR",
                             "RR",
                             "SBP",
                             "SpO2",
                             "Lactate",
                             "Hgb",
                             "WBC",
                             "Plts",
                             "Naive B cells",
                             "Mem. B cells",
                             "Plasma cells",
                             "CD8+",
                             "Naive CD4+",
                             "Resting CD4+",
                             "Activated CD4+",
                             "Resting NK cells",
                             "Monocytes",
                             "Activated DCs",
                             "Resting mast cells",
                             "Eosinophils",
                             "Neutrophils")
head(covidsympbiom)

#Calculate correlation coefficients 
cormatsevere <- cor(covidsympbiom, method = c("spearman"))

#Unadjusted p-values
library(corrplot)
p3 <- cor.mtest(covidsympbiom)
colnames(cormatsevere)
colnames(p3$p)

#FDR-adjusted P-values 
library(psych)
p3adj <- corr.p(cormatsevere,84,adjust="BH",alpha=.05)
colnames(cormatsevere)
colnames(p3adj$p)

#Plot correlation matrix 
library(corrplot)
corrplot(cormatsevere,title = "", 
         method = "square", 
         outline = T, 
         addgrid.col = "grey40", 
         order="hclust", 
         hclust.method = 'ward.D2',
         mar = c(2,0,2,0), 
         addrect = NULL, 
         type = "lower",
         rect.col = "black", 
         rect.lwd = 5,
         tl.pos = "ld",
         tl.col = "black", 
         tl.srt = 90,
         tl.cex = 1, 
         cl.cex = 1,
         cl.pos = "b",
         bg = "grey96",
         p.mat = p3adj$p, sig.level = 0.05, insig = "blank")

#Correlation plot - immune cell populations and biomarkers
covidsympbiom <- data.matrix(subset(covid)[1:84, c(237:242, 245:251, 253:284, 301:307, 311,313,318,319,321,322)], rownames.force = NA)
head(covidsympbiom)
colnames(covidsympbiom) <- c("sCD40L", 
                             "EGF", 
                             "Eotaxin", 
                             "FGF-2", 
                             "FLT-3L", 
                             "FKN", 
                             "CXCL1", 
                             "IFN-\u03b12", 
                             "IFN-\u03b3", 
                             "IL-1\u03b1", 
                             "IL-1\u03b2", 
                             "IL-1Ra",
                             "IL-2",
                             "IL-4",
                             "IL-5",
                             "IL-6",
                             "IL-7",
                             "IL-8",
                             "IL-9",
                             "IL-10",
                             "IL-12p40",
                             "IL-12p70",
                             "IL-13",
                             "IL-15",
                             "IL-17A",
                             "IL-17E",
                             "IL-17F",
                             "IL-18",
                             "IL-22",
                             "IL-27",
                             "CXCL10",
                             "MCP-1",
                             "MCP-3",
                             "M-CSF",
                             "MDC",
                             "CXCL9",
                             "CCL3",
                             "CCL4",
                             "PDGF-AA",
                             "PDGF-AB/BB",
                             "CCL5",
                             "TGF-\u03b1",
                             "TNF",
                             "LT-\u03b1",
                             "VEGF-A",
                             "Naive B cells",
                             "Mem. B cells",
                             "Plasma cells",
                             "CD8+",
                             "Naive CD4+",
                             "Resting CD4+",
                             "Activated CD4+",
                             "Resting NK cells",
                             "Monocytes",
                             "Activated DCs",
                             "Resting mast cells",
                             "Eosinophils",
                             "Neutrophils")
head(covidsympbiom)

#Calculate correlation coefficients 
cormatsevere <- cor(covidsympbiom, method = c("spearman"))

#Unadjusted p-values
library(corrplot)
p3 <- cor.mtest(covidsympbiom)
colnames(cormatsevere)
colnames(p3$p)

#FDR-adjusted P-values 
library(psych)
p3adj <- corr.p(cormatsevere,84,adjust="BH",alpha=.05)
colnames(cormatsevere)
colnames(p3adj$p)

#Plot correlation matrix 
library(corrplot)
corrplot(cormatsevere,title = "", 
         method = "square", 
         outline = T, 
         addgrid.col = "grey40", 
         order="hclust", 
         hclust.method = 'ward.D2',
         mar = c(2,0,2,0), 
         addrect = NULL, 
         type = "lower",
         rect.col = "black", 
         rect.lwd = 5,
         tl.pos = "ld",
         tl.col = "black", 
         tl.srt = 90,
         tl.cex = 0.8, 
         cl.cex = 0.80,
         cl.pos = "b",
         bg = "grey96",
         p.mat = p3adj$p, sig.level = 0.05, insig = "blank")

