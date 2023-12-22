#Figure 5a-5b - Correlations between clinical variables and mediators 

#Clear R environment 
rm(list = ls())

#Import master data, PID to fixed row identifier 
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid, exclude asymptomatic 
covid <- subset(combined, pathogencode==1 & min_no_sym==0)

#Correlation plot for symptomatic COVID-19 - exclude GM-CSF, G-CSF, and IL-3 given large prop. values < LLD
covidsympbiom <- data.matrix(subset(covid)[1:240, c(6,7,286,288, 290:297,237:242, 245:251, 253:284)], rownames.force = NA)
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
                             "sCD40L", 
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
                             "VEGF-A")
head(covidsympbiom)

#Calculate correlation coefficients 
cormatsevere <- cor(covidsympbiom, method = c("spearman"))

#Unadjusted p-values
library(corrplot)
p3 <- cor.mtest(covidsympbiom)
colnames(cormatsevere)
colnames(p3$p)

#BH-adjusted P-values 
library(psych)
p3adj <- corr.p(cormatsevere,240,adjust="BH",alpha=.05)
colnames(cormatsevere)
colnames(p3adj$p)

#Plot correlation matrix 
#Create optional color palette
colcorr <- colorRampPalette(c("#fcfdbf", "#fc8961", "#b73779", "#51127c", "#000004"))
colcorr2 <- colorRampPalette(c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"))

library(corrplot)
par(mar = c(6, 6, 1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1)) 
symptomaticcorrplot <- corrplot(cormatsevere,title = "", 
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
         tl.cex = 0.9, 
         cl.cex = 0.9,
         cl.pos = "b",
         bg = "grey96",
         p.mat = p3adj$p, sig.level = 0.05, insig = "blank")

#FDR-adjusted network correlation plot
library(qgraph)
covidsympfdr <- FDRnetwork(cormatsevere, cutoff = 0.05, method = 'pval')

qgraphcovidsympfdr <- qgraph(covidsympfdr, 
                               title = "",
                               title.cex = 1.5,
                               graph = "cor",
                               layout = "spring", 
                               repulsion = 0.95,
                               threshold = 0, 
                               directed = FALSE, 
                               edge.labels = FALSE, 
                               esize = 8,
                               vsize = 5,
                               color = "grey96",
                               labels=colnames(covidsympbiom), 
                               label.color = "black",
                               label.prop=0.8,
                               borders = TRUE,
                               edge.width = 1.75,
                               border.color = "black",
                               border.width = 2,
                               posCol = "#2166ac",
                               negCol = "#b2182b",
                               curve = 0.2, curveAll = T)

#Centrality statistics 
centralitystatsseverecovid <- centrality_auto(qgraphcovidsympfdr, 
                                         weighted = TRUE, 
                                         signed = FALSE)
centralityscaledseverecovid <- scale(centralitystatsseverecovid$node.centrality, 
                                center=TRUE, 
                                scale=TRUE)

#Plot these metrics, along with Expected Influence
centralityplot <- centralityPlot(qgraphcovidsympfdr, include = 
                 c("Strength", "ExpectedInfluence"),
               orderBy ="ExpectedInfluence")
centralityplot <- centralityplot + theme(text = element_text(size = 14))
centralityplot

#Repeat analysis including HIV status, exclude patients with indeterminate HIV status
covid <- subset(combined, pathogencode==1)
covid <- subset(covid, min_no_sym==0 & (hivrdtresult==1 | hivrdtresult==0))

#Correlation plot for symptomatic COVID-19
covidsympbiomhiv <- data.matrix(subset(covid)[1:236, c(6,7,286,288, 290:297,237:242, 245:251, 253:284, 130)], rownames.force = NA)
head(covidsympbiomhiv)
colnames(covidsympbiomhiv) <- c("Sex",
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
                             "sCD40L", 
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
                             "HIV")
head(covidsympbiomhiv)

#Calculate correlation coefficients 
cormatseverehiv <- cor(covidsympbiomhiv, method = c("spearman"))

#Unadjusted p-values
library(corrplot)
p3hiv <- cor.mtest(covidsympbiomhiv)
colnames(cormatseverehiv)
colnames(p3hiv$p)

#FDR-adjusted P-values 
library(psych)
p3adjhiv <- corr.p(cormatseverehiv,236,adjust="BH",alpha=.05)
colnames(cormatseverehiv)
colnames(p3adjhiv$p)

#Plot correlation matrix 
#Create optional color palette
colcorr <- colorRampPalette(c("#fcfdbf", "#fc8961", "#b73779", "#51127c", "#000004"))
colcorr2 <- colorRampPalette(c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"))

library(corrplot)
corrplot(cormatseverehiv,title = "", 
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
         cl.cex = 0.8,
         cl.pos = "b",
         bg = "grey96",
         p.mat = p3adjhiv$p, sig.level = 0.05, insig = "blank")

#FDR-adjusted network correlation plot
library(qgraph)
covidsympfdrhiv <- FDRnetwork(cormatseverehiv, cutoff = 0.05, method = 'pval')

qgraphcovidsympfdrhiv <- qgraph(covidsympfdrhiv, 
                             title = "",
                             title.cex = 1.5,
                             graph = "cor",
                             layout = "spring", 
                             repulsion = 0.95,
                             threshold = 0, 
                             directed = FALSE, 
                             edge.labels = FALSE, 
                             esize = 8,
                             vsize = 5,
                             color = "grey96",
                             labels=colnames(covidsympbiomhiv), 
                             label.color = "black",
                             label.prop=0.8,
                             borders = TRUE,
                             edge.width = 1.5,
                             border.color = "black",
                             border.width = 2,
                             posCol = "#2166ac",
                             negCol = "#b2182b",
                             curve = 0.2, curveAll = T)

#Centrality statistics 
centralitystatsseverecovidhiv <- centrality_auto(qgraphcovidsympfdrhiv, 
                                              weighted = TRUE, 
                                              signed = FALSE)
centralityscaledseverecovidhiv <- scale(centralitystatsseverecovidhiv$node.centrality, 
                                     center=TRUE, 
                                     scale=TRUE)

#Lets now plot these metrics, along with Expected Influence
centralityplothiv <- centralityPlot(qgraphcovidsympfdrhiv, include = 
                 c("Strength", "ExpectedInfluence"),
               orderBy ="ExpectedInfluence")
centralityplothiv <- centralityplothiv + theme(text = element_text(size = 14))
centralityplothiv

