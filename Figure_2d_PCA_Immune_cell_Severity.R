#Figure 2d - PCA of immune cell variance stratified by severe vs non-severe COVID

#Clear R environment 
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients with RNAseq data
covid <- subset(combined, pathogencode == 1 & rna_seq==1)

#PCA of CIBERSORTx-inferred immune cell abundance 
covidpcamatrix <- data.matrix(covid[1:100, c(301:307, 311,313,318,319,321,322)], rownames.force = NA)
head(covidpcamatrix)
colnames(covidpcamatrix) <- c("Naive B cells",
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
head(covidpcamatrix)

library("FactoMineR")
res.pca <- PCA(covidpcamatrix, ncp=5, scale.unit = TRUE, graph = TRUE)
print(res.pca)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

#Determination of top variables driving separation in the PCs 
library("corrplot")
library("factoextra")
var <- get_pca_var(res.pca)
head(var$cos2, 13)
head(var$contrib, 13)

fviz_pca_var(res.pca, col.var = "cos2",
             select.var = list(contrib = 13),
             gradient.cols = "aaas", 
             repel = TRUE)

#Corrplot
library(corrplot)
PCAloadmatrix <- var$cos2
PCAloadmatrix
colnames(PCAloadmatrix) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
rownames(PCAloadmatrix) <- c("Naive B cells",
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

par(xpd=TRUE)
pcfactorloadings <- corrplot(PCAloadmatrix[,1:3], 
                             is.corr=FALSE, 
                             tl.col = "black", 
                             tl.cex = 1.25, 
                             tl.srt = 45,
                             cl.cex = 1.25, 
                             cl.ratio = 2, 
                             cl.pos = "r",
                             col = COL2(diverging = c("PRGn"), n = 200),
                             mar = c(0.5, 0.5, 0.5, 1.3),
                             outline=TRUE) 
recordPlot(pcfactorloadings)
pcfactorloadings

#First 2 PCs with density plot
library("mixOmics")
library("PLSDAbatch")

pcacols <- c("#4BB446", "#AF46B4")

covid[covid$whosevere_imp==1, "whosevere_imp2"] <- "Severe"
covid[covid$whosevere_imp==0, "whosevere_imp2"] <- "Non-severe"
covid$whosevere_imp2 <- as.factor(covid$whosevere_imp2)

pcav2 <- pca(covid[1:100, c(301:307, 311,313,318,319,321,322)], ncomp = 3, scale = TRUE)
pcadensity <- Scatter_Density(object = pcav2, 
                batch = covid$whosevere_imp2,
                trt = NULL,
                color.set = pcacols,
                batch.legend.title = "COVID-19 Severity",
                density.lwd = 0.1,
                title = NULL,
                title.cex = 1.1,
                legend.cex = 1.1,
                legend.title.cex = 1.1)
pcadensity < pcadensity + theme(text = element_text(size = 14))
pcadensity

