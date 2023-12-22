#Figure 5b-5l - Cluster-derived COVID Response Signatures 

#Clear R environment 
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select symptomatic COVID patients
covid <- subset(combined, pathogencode == 1 & min_no_sym == 0)

#Select log-transformed biomarker data to analyze - exclude GM-CSF, G-CSF, and IL-3 given large prop. values <LLD
biomarkers <- covid[1:240, c(237:242,245:251,253:284)]
head(biomarkers)
colnames(biomarkers) <- c("sCD40L", 
                             "EGF", 
                             "Eotaxin", 
                             "FGF-2", 
                             "FLT-3L", 
                             "FKN", 
                             "CXCL1", 
                             "IFN-a2", 
                             "IFN-gamma", 
                             "IL-1a", 
                             "IL-1b", 
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
                             "RANTES",
                             "TGF-a",
                             "TNF-a",
                             "LT-a",
                             "VEGF-A")
head(biomarkers)

#Scale/center biomarkers
biomscaled <- scale(biomarkers, center=TRUE, scale=TRUE)
biomscaled <- as.data.frame(biomscaled)

#Explore optimal cluster partition 
library(NbClust)
library(factoextra)

set.seed(12345)
nb <- NbClust(biomscaled, distance = "euclidean", min.nc = 2,
              max.nc = 12, method = "km")

nb$All.index
nb$Best.nc
nb$All.CriticalValues
nb$Best.partition

#Implement k-means consensus clustering  
library(ConsensusClusterPlus)

results = ConsensusClusterPlus(t(biomscaled),
                            maxK=12,
                            reps=1000,
                            pItem=0.8,pFeature=1,
                            title="example3",
                            distance="euclidean",
                            seed=12345,
                            clusterAlg="km")

#Consensus matrices and cluster/item consensuses
results[[2]][["consensusMatrix"]][1:12,1:12] 
icl = calcICL(results,plot=NULL)
icl[["clusterConsensus"]]

#Label patients with their assigned cluster based on 2-cluster partition and extract cluster assignments
hcc <- results[[2]][["consensusClass"]]
hcc

#Merge cluster assignments into dataset 
hcc <- data.frame(hcc)

covidclustassign <- merge(covid,hcc,all=T,by='row.names')

rownames(covidclustassign) = covidclustassign$Row.names
covidclustassign$Row.names = NULL

covidcluster1 <- subset(covidclustassign, hcc==1)
covidcluster2 <- subset(covidclustassign, hcc==2)

#PCA - exclude GM-CSF, G-CSF, and IL-3
covidpcamatrix <- data.matrix(covid[1:240, c(237:242,245:251,253:284)], rownames.force = NA)
head(covidpcamatrix)
colnames(covidpcamatrix) <- c("sCD40L", 
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
                              "RANTES",
                              "TGF-\u03b1",
                              "TNF",
                              "LT-\u03b1",
                              "VEGF-A")
head(covidpcamatrix)

library("FactoMineR")
res.pca <- PCA(covidpcamatrix, ncp=5, scale.unit = TRUE, graph = TRUE)
print(res.pca)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

#Determination of top variables driving separation in the PCs 
library("corrplot")
var <- get_pca_var(res.pca)
head(var$cos2, 48)
head(var$contrib, 48)

fviz_pca_var(res.pca, col.var = "cos2",
             select.var = list(contrib = 23),
             gradient.cols = "aaas", 
             repel = TRUE)

#First 2 PCs stratified by cluster
library("FactoMineR")
res.pca <- PCA(covidclustassign[1:240, c(237:242,245:251,253:284)], scale.unit=TRUE, graph = TRUE)
fviz_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             mean.point = FALSE,
             pointsize = 2.5,
             col.ind = as.factor(covidclustassign$hcc), # color by groups
             palette = c("#00468B99", "#ED000099"),
             addEllipses = TRUE, 
             legend.title = "")

#First 2 PCs by cluster with density plot
library("mixOmics")
library("PLSDAbatch")

pcacols <- c("#00468B99", "#ED000099")

pcav2 <- pca(covidclustassign[1:240, c(237:242,245:251,253:284)], ncomp = 3, scale = TRUE)
pcadensity <- Scatter_Density(object = pcav2, 
                batch = covidclustassign$hcc,
                trt = NULL,
                color.set = pcacols,
                batch.legend.title = "CRS",
                density.lwd = 0.1,
                title = NULL,
                title.cex = 1.5,
                legend.cex = 1.5,
                legend.title.cex = 1.5)
pcadensity < pcadensity + theme(text = element_text(size = 14))
pcadensity

#XGBoost to determine most important mediators in predicting cluster assignment
library(readr)
library(pROC)
library(xgboost)
library(dplyr)
library(caret)
library(ggplot2)
#Create matrix of relevant biomarker data for input to XGBoost - exclude GM-CSF, G-CSF, IL-3 
covidxgbmatrix <- data.matrix(covidclustassign[1:240, c(237:242,245:251,253:284)], rownames.force = NA)
head(covidxgbmatrix)
colnames(covidxgbmatrix) <- c("sCD40L", 
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
                              "RANTES",
                              "TGF-\u03b1",
                              "TNF",
                              "LT-\u03b1",
                              "VEGF-A")
head(covidxgbmatrix)

#Create outcome for variable to be classified/predicted, in this case cluster assignment 
clusterpredict <- covidclustassign[,"hcc"] == "2"

##Performing a grid search for Hyperparameter Tuning
### Tuned Parameters - "nrounds", "eta", and "max_depth"
grid_label <- factor(clusterpredict, labels = c("C1", "C2"))
boost_train_cont = trainControl(method = "cv", number = 10, verboseIter = TRUE, returnData = FALSE, returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
boost_grid <- expand.grid(nrounds = c(10, 100, 200, 400, 800, 1000, 2000),eta = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.001),max_depth = c(2, 3, 4, 6, 8, 10), gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1)
set.seed(12345)
model_train <- train(x=covidxgbmatrix, y=grid_label, trControl = boost_train_cont, tuneGrid = boost_grid, method = "xgbTree")
ggplot(model_train$results, aes(x = as.factor(eta), y = max_depth, size = ROC, color = ROC)) + geom_point() + theme_bw() + scale_size_continuous(guide = "none")
options(max.print = 100000)
model_train$results
model_train$results[order(model_train$results$ROC),] 

#XGBoost to predict cluster assignment based on all 45 biomarkers
library(xgboost)
params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.1, gamma=0, max_depth=2, min_child_weight=1, subsample=1, colsample_bytree=1)
bst <- xgboost(data = covidxgbmatrix, label = clusterpredict, nround = 2000, params = params)    

#Determine variable importance in predicting cluster assignment 
importance <- xgb.importance(feature_names = colnames(covidxgbmatrix), model = bst)
head(importance)

#Plot variable importance - 10 most important variables 
library(Ckmeans.1d.dp)
library(ggplot2)
importanceplot <- xgb.ggplot.importance(importance_matrix = importance, rel_to_first = FALSE, top_n = 20, n_clusters = c(1), xlab = "Split-Gain")

importanceplot <- importanceplot+ 
  scale_fill_manual(values=c("mediumpurple4")) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size=14)) + 
  theme(axis.text.y = element_text(color = "black")) +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "", y = "Split-Gain", x = "", color="black") + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = "blank")) + 
  theme(panel.grid.minor = element_line(size = 0.25, linetype = "blank"))
importanceplot <- importanceplot + theme(axis.title.x = element_text(colour = "black", size = 12), axis.title.y = element_text(colour = "black", size = 12), 
                     axis.text.y = element_text(colour = "black", size = 12), axis.text.x = element_text(colour = "black", size = 12)) 
importanceplot     

#SHAP summary plot
contr <- predict(bst, covidxgbmatrix, predcontrib = TRUE)
xgb.plot.shap(covidxgbmatrix, contr, model = bst, top_n = 20, n_col = 3)
shap <- xgb.ggplot.shap.summary(covidxgbmatrix, contr, model = bst, top_n = 20)
shap <- shap + theme_bw() 
shap <- shap + labs(title = "", y = "SHAP value (impact on XGBoost model)", x = "", color="Mediator\n value") 
shap <- shap + theme(legend.title = element_text(color = "black", size = 12),
                     legend.text = element_text(color = "black", size = 12), legend.key.size = unit(0.5, "cm"))
shap <- shap + theme(axis.title.x = element_text(colour = "black", size = 12), axis.title.y = element_text(colour = "black", size = 12), 
                     axis.text.y = element_text(colour = "black", size = 12), axis.text.x = element_text(colour = "black", size = 12)) 
shap

#Combine importance and SHAP plots
library(ggpubr)
covidclusterboostshap <- ggarrange(importanceplot, shap, ncol = 1, nrow=2, labels = c("", ""))
covidclusterboostshap

#Visualize between-cluster variance in biomarkers using ComplexHeapMap
#Clinical variables to factors for the annotations
#HCC
covidclustassign$hcc <- factor(covidclustassign$hcc)
levels(covidclustassign$hcc)

#Sex
covidclustassign$sex <- factor(covidclustassign$sex)
levels(covidclustassign$sex)

#HIV
covidclustassign$hivrdtresult <- factor(covidclustassign$hivrdtresult)
levels(covidclustassign$hivrdtresult)

#WHO clinical Severity
covidclustassign$whoseverity_imp <- factor(covidclustassign$whoseverity_imp)
levels(covidclustassign$whoseverity_imp)

#Study Phase
covidclustassign$phase3 <- factor(covidclustassign$phase3)
levels(covidclustassign$phase3)

#Create a new biomarker matrix with labels
biomscaledheatmatrix <- data.matrix(biomscaled, rownames.force = NA)
biomscaledheatmatrix
colnames(biomscaledheatmatrix) <- c("sCD40L", 
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
                                    "RANTES",
                                    "TGF-\u03b1",
                                    "TNF",
                                    "LT-\u03b1",
                                    "VEGF-A")
biomscaledheatmatrix

#Plot
library(ComplexHeatmap)

library(RColorBrewer)
display.brewer.pal(9, "Oranges")
brewer.pal(9, "Oranges")

display.brewer.pal(9, "Greens")
brewer.pal(9, "Greens")

brewer.pal(n=5,"Set1")

#Rename cluster assignment variable
covidclustassign[covidclustassign$hcc == 1, "hcc3"] <- "CRS 1"
covidclustassign[covidclustassign$hcc == 2, "hcc3"] <- "CRS 2"

#Make annotation
library(circlize)
col_age = colorRamp2(c(10, 20, 30, 40, 50, 60, 70, 80, 90), c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))
col_duration = colorRamp2(c(0,3,6,9,12,15,18,21,24), c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"))

ha = HeatmapAnnotation("Age" = covidclustassign$age[1:240],
                       "Sex" = covidclustassign$sex[1:240],
                       "Living with HIV" = covidclustassign$hivrdtresult[1:240],
                       "Dominant variant(s)" = covidclustassign$phase3[1:240],
                       "WHO Clinical Severity" = covidclustassign$whoseverity_imp[1:240],
                       show_legend = c("HIV infection" = TRUE, "Age" = TRUE, "Illness Duration" = TRUE, "Sex" = TRUE, "Study Phase" = TRUE, "WHO Clinical Severity" = TRUE),
                       annotation_legend_param = list("Age" = list(title = "Age", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)), "Sex" = list(title = "Sex", at = c("1", "0"), labels = c("Male", "Female"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
                                                      "Living with HIV" = list(title = "Living with HIV", at = c("1", "0"), labels = c("Yes", "No"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
                                                      "Dominant variant(s)" = list(title = "Dominant variant(s)", at = c("1", "2", "3"), labels = c("Varied A/B Lineages", "A.23/A.23.1", "Delta (B.1.617.2/AY)"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
                                                      "WHO Clinical Severity" = list(title = "WHO Clinical Severity", at = c("2", "3", "4"), labels = c("Mild", "Moderate", "Severe"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14))),
                       col = list(
                                  "Sex" = c("1" = "#143d59", "0" = "#f4b41a"),
                                  "Living with HIV" = c("1" = "#E2725A", "0" = "#4A6274"),
                                  "WHO Clinical Severity" = c("2" = "blue4", "3" = "#BC3C29FF", "4" = "#e6ab02"),
                                  "Dominant variant(s)" = c("1" = "#440154", "2" = "#2a788e", "3" = "#7ad151"),
                                  "Age" = col_age),
                                  na_col = "#4A6274",
                       annotation_name_gp = gpar(fontsize = 14))

library(circlize)
col_map = colorRamp2(c(-4, -2, 0, 2, 4), c("#000004", "#ED000099", "#b73779", "#fc8961", "#fcfdbf"))
   
col_map2 = colorRamp2(c(-4, -2, 0, 2, 4), c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"))

col_map3 = colorRamp2(c(-4, -2, 0, 2, 4), c("#053061", "#053061", "#f7f7f7", "#fde725", "#fde725"))

set.seed(1234)
Heatmap(t(biomscaledheatmatrix), 
        name = "Z-score", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_km = 4,
        show_column_names = FALSE, 
        show_column_dend = FALSE, 
        show_row_dend = FALSE, 
        row_names_side = "left", 
        column_split = covidclustassign$hcc3, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
        column_title_gp = gpar(fontsize = 14),
        cluster_column_slices = FALSE,
        row_title = c(),
        row_names_gp = gpar(fontsize = 14),  
        top_annotation = ha,
        border=TRUE)

#Compare biomarker concentrations across clusters
library(tableone)
library(tableHTML)
#Create a variable list 
listbiomarkers <- c("logscd40l", 
                    "logegf",	
                    "logeotaxin",	
                    "logfgf2",	
                    "logflt3l",	
                    "logfractalkine",	
                    "loggcsf",	
                    "loggmcsf",	
                    "loggroalpha",	
                    "logifna2",	
                    "logifngamma",
                    "logil1a",
                    "logil1b",
                    "logil1ra",
                    "logil2",
                    "logil3",
                    "logil4",
                    "logil5",
                    "logil6",
                    "logil7",
                    "logil8",
                    "logil9",
                    "logil10",
                    "logil12p40",	
                    "logil12p70",
                    "logil13",	
                    "logil15",
                    "logil17a",
                    "logil17eil25",
                    "logil17f",
                    "logil18",
                    "logil22",
                    "logil27",
                    "logip10",
                    "logmcp1",
                    "logmcp3",
                    "logmcsf",
                    "logmdc",
                    "logmigcxcl9",
                    "logmip1a",
                    "logmip1b",
                    "logpdgfaa",
                    "logpdgfabbb",
                    "lograntes",
                    "logtgfa",
                    "logtnfa",
                    "logtnfb"	,
                    "logvegfa")

tablebiomarkerscluster <- CreateTableOne(vars = listbiomarkers, 
                                  data = covidclustassign, 
                                  strata = "hcc",
                                  addOverall = TRUE)
tablebiomarkerscluster
summary(tablebiomarkerscluster)

tablebiomarkerscluster <- print(tablebiomarkerscluster, nonnormal = c("logscd40l", 
                                                        "logegf",	
                                                        "logeotaxin",	
                                                        "logfgf2",	
                                                        "logflt3l",	
                                                        "logfractalkine",	
                                                        "loggcsf",	
                                                        "loggmcsf",	
                                                        "loggroalpha",	
                                                        "logifna2",	
                                                        "logifngamma",
                                                        "logil1a",
                                                        "logil1b",
                                                        "logil1ra",
                                                        "logil2",
                                                        "logil3",
                                                        "logil4",
                                                        "logil5",
                                                        "logil6",
                                                        "logil7",
                                                        "logil8",
                                                        "logil9",
                                                        "logil10",
                                                        "logil12p40",	
                                                        "logil12p70",
                                                        "logil13",	
                                                        "logil15",
                                                        "logil17a",
                                                        "logil17eil25",
                                                        "logil17f",
                                                        "logil18",
                                                        "logil22",
                                                        "logil27",
                                                        "logip10",
                                                        "logmcp1",
                                                        "logmcp3",
                                                        "logmcsf",
                                                        "logmdc",
                                                        "logmigcxcl9",
                                                        "logmip1a",
                                                        "logmip1b",
                                                        "logpdgfaa",
                                                        "logpdgfabbb",
                                                        "lograntes",
                                                        "logtgfa",
                                                        "logtnfa",
                                                        "logtnfb"	,
                                                        "logvegfa"))

write_tableHTML(tableHTML(tablebiomarkerscluster), file = 'biomarkersbycluster.html')

#Table comparing clinical characteristics across clusters - TableOne package
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessdurationenroll_imp",
              "phase3",
              "historyoffever", 
              "nightsweats", 
              "headache", 
              "cough", 
              "sorethroat",
              "runnynose",
              "sob",
              "tempmax_imp", 
              "heartrate3_imp", 
              "resprate3_imp", 
              "sbp3_imp",
              "o2sat3", 
              "kpsadmit16plus_imp",
              "lactate_imp",
              "hgb_imp",
              "wbc_imp",
              "platelets_imp",
              "malariardtresult",
              "hivrdtresult", 
              "artprior",
              "priortuberculosis",
              "heartdisease",
              "hypertension",
              "diabetes",
              "supplementaloxygen",
              "supplementaloxygenamount",
              "abx",
              "steroids",
              "hospdeathtransf") 

#Define categorical variables
catVars <- c("sex", 
             "phase3",
             "historyoffever", 
             "nightsweats", 
             "headache", 
             "cough", 
             "sorethroat",
             "runnynose",
             "sob",
             "malariardtresult",
             "hivrdtresult", 
             "artprior",
             "heartdisease",
             "hypertension",
             "diabetes",
             "priortuberculosis",
             "whoseverity_imp",
             "supplementaloxygen",
             "abx",
             "steroids",
             "hospdeathtransf") 

tableclinicalcluster <- CreateTableOne(vars = listVars, 
                                     data = covidclustassign, 
                                     factorVars = catVars,
                                     strata = "hcc",
                                     includeNA = TRUE,
                                     addOverall = TRUE)
tableclinicalcluster
summary(tableclinicalcluster)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tableclinicalcluster <- print(tableclinicalcluster, nonnormal = c("age",
                                                                  "tempmax_imp",
                                                                  "illnessdurationenroll_imp",
                                                                  "heartrate3_imp", 
                                                                  "resprate3_imp", 
                                                                  "sbp3_imp",
                                                                  "o2sat3",
                                                                  "lactate_imp",
                                                                  "hgb_imp",
                                                                  "wbc_imp",
                                                                  "platelets_imp",
                                                                  "kpsadmit16plus_imp",
                                                                  "supplementaloxygenamount"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tableclinicalcluster), file = 'tablecovidclinicalcluster.html')

#Mosaic plots for WHO Clinical Severity groups across clusters
library(ggplot2)
library(dplyr)
mosaicseverity <- covidclustassign %>%
  group_by(hcc, whoseverity_imp) %>%
  summarise(count = n()) %>%
  mutate(hcc.count = sum(count),
         prop = count/sum(count)) %>%
  ungroup()

mosaicseverityplot <- ggplot(mosaicseverity,
       aes(x = hcc, y = prop, width = hcc.count, fill = whoseverity_imp)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(~hcc, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c("blue4", "#BC3C29FF", "#e6ab02"), 
                    name="",
                    breaks=c("2", "3", "4"),
                    labels=c("Mild", "Moderate", "Severe"))
mosaicseverityplot 

#Adjust appearance
mosaicseverityplot<-mosaicseverityplot + 
  theme(axis.title.x = element_text(color = "black", size = 14, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, face = "plain"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  theme(axis.text.x = element_text(color = "black", size = 14, face = "plain"),
        axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  theme(legend.title = element_text(color = "black", size = 14, face = "plain"), legend.text = element_text(color = "black", size = 14, face = "plain"), legend.position = "top") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  theme(strip.text = element_blank())
mosaicseverityplot <- mosaicseverityplot + labs(x = "", y = "Proportion of Patients",) + 
                      scale_x_discrete(labels=c("1" = "CRS 1", "2" = "CRS 2")) + scale_y_continuous(expand = c(0,0))
mosaicseverityplot 

#Mosaic plots for pandemic phase across cluster
library(ggplot2)
library(dplyr)
mosaicphase <- covidclustassign %>%
  group_by(hcc, phase3) %>%
  summarise(count = n()) %>%
  mutate(hcc.count = sum(count),
         prop = count/sum(count)) %>%
  ungroup()

mosaicphaseplot <- ggplot(mosaicphase,
                             aes(x = hcc, y = prop, width = hcc.count, fill = phase3)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(~hcc, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c("#440154", "#2a788e", "#7ad151"), 
                    name="",
                    breaks=c("1", "2", "3"),
                    labels=c("Varied A/B Lineages", "A.23/A.23.1", "Delta (B.1.617.2/AY)"))
mosaicphaseplot 

#Adjust appearance
mosaicphaseplot<-mosaicphaseplot + 
  theme(axis.title.x = element_text(color = "black", size = 14, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, face = "plain"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  theme(axis.text.x = element_text(color = "black", size = 14, face = "plain"),
        axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  theme(legend.title = element_text(color = "black", size = 14, face = "plain"), legend.text = element_text(color = "black", size = 14, face = "plain"), legend.position = "top") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  theme(strip.text = element_blank())
mosaicphaseplot <- mosaicphaseplot + labs(x = "", y = "Proportion of Patients",) + 
  scale_x_discrete(labels=c("1" = "CRS 1", "2" = "CRS 2")) + scale_y_continuous(expand = c(0,0))
mosaicphaseplot 

#Combine mosaic plots 
library(ggpubr)
covidclustermosaicplots <- ggarrange(mosaicseverityplot, mosaicphaseplot, ncol = 1, nrow=2, labels = c("", ""))
covidclustermosaicplots

#Upset plot - format data
library(ggplot2)
library(ggthemes)
library(reshape)
library(gmodels)
library(DescTools)

#Create binary variable for KPS <=50
covidclustassign$kps50less <- NA
covidclustassign[covidclustassign$kpsadmit16plus<=50, "kps50less"] <- 1
covidclustassign[covidclustassign$kpsadmit16plus>=60, "kps50less"] <- 0

#Create binary variable for WHO Severe Covid 
covidclustassign$severecovid_who <- NA
covidclustassign[covidclustassign$whoseverity_imp==4, "severecovid_who"] <- 1
covidclustassign[covidclustassign$whoseverity_imp==2 | covidclustassign$whoseverity_imp==3, "severecovid_who"] <- 0

#Upset plot - ComplexUpset
library(ComplexUpset)

covidupset <- subset(covidclustassign, select=c(severecovid_who, kps50less, assisttowalkoob, supplementaloxygen, hcc))
colnames(covidupset) <- c("Severe COVID-19", "KPS \u2264 50", "Unable to ambulate", "Oxygen Therapy", "CRS")
severity = colnames(covidupset)[1:4]

upsetplot2 <- upset(
  covidupset,
  severity,
  name = "",
  min_degree=1,
  sort_intersections_by='degree',
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=CRS)) + scale_fill_manual(values=c("1" = "#377eb8", "2" = "#df8f44"))
  ),
  width_ratio=0.1,
  themes=upset_default_themes(text=element_text(size=14)),
  set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))
  ))
upsetplot2 

#Cumulative incidence plot of hospital mortality by cluster - first in all symptomatic patients and then in severe patients
library(tidycmprsk)
library(ggsurvfit)

covidclustassign$hospdeathtransf <- as.factor(covidclustassign$hospdeathtransf)
covidclustassign$hcc <- as.factor(covidclustassign$hcc)

mortalityplot2 <- cuminc(Surv(daysspentinhospitalimp, hospdeathtransf) ~ hcc, covidclustassign) %>%
  ggcuminc(linewidth = 1.5) 
mortalityplot2 <- mortalityplot2 + theme_bw() + labs(title = "Symptomatic COVID-19", x ="Days since hospital admission") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.position = "bottom", legend.text=element_text(size=14)) + ylim(0, 0.50) 
mortalityplot2 <- mortalityplot2 + scale_colour_manual(values = brewer.pal(8, "Set1")[4:5])
mortalityplot2

mortalityplotsevere2 <- cuminc(Surv(daysspentinhospitalimp, hospdeathtransf) ~ hcc, subset(covidclustassign, whosevere_imp==1)) %>%
  ggcuminc(linewidth = 1.5) 
mortalityplotsevere2 <- mortalityplotsevere2 + theme_bw() + labs(title = "Severe COVID-19", x ="Days since hospital admission") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.position = "bottom", legend.text=element_text(size=14)) + ylim(0, 0.50) 
mortalityplotsevere2 <- mortalityplotsevere2 + scale_colour_manual(values = brewer.pal(8, "Set1")[4:5])
mortalityplotsevere2

#Combine plots
library(ggpubr)
covidclustermortalityplots2 <- ggarrange(mortalityplot2, mortalityplotsevere2, ncol = 1, nrow=2, labels = c("", ""), legend = "bottom", common.legend = TRUE)
covidclustermortalityplots2

#CRR
crr_symp <- crr(Surv(daysspentinhospitalimp, hospdeathtransf) ~ hcc, covidclustassign)
crr_symp

crr_severe <- crr(Surv(daysspentinhospitalimp, hospdeathtransf) ~ hcc, subset(covidclustassign, whosevere_imp==1))
crr_severe

#Density plot for supplemental O2 levels across clusters
#Make data frame for cluster and supplemental O2 level
covidclustero2levels <- subset(covidclustassign, select=c(supplementaloxygenamount, hcc)) 
covidclustero2levels <- na.omit(covidclustero2levels)

#Calculate median for each cluster and compare medians using Wilcoxon rank-sum test
library(dplyr)
group_by(covidclustero2levels, hcc) %>%
  dplyr::summarize(count = n(),
    median = median(supplementaloxygenamount, na.rm = TRUE),
    IQR = IQR(supplementaloxygenamount, na.rm = TRUE))

wilcox.test(supplementaloxygenamount ~ hcc, data = covidclustero2levels,
                   exact = FALSE)

#Make plot
o2levelsplot <-ggplot(covidclustero2levels, aes(x=supplementaloxygenamount, fill=hcc)) +
  geom_density(alpha=0.4) 
o2levelsplot <- o2levelsplot + scale_fill_manual(labels = c("CRS 1", "CRS 2"), values=c("#767676FF", "#350E20FF"))+ labs(y= "Density", x = "Oxygen Therapy Flow Rate (liters/minute)")
o2levelsplot <- o2levelsplot + annotate("segment", x = 5, xend = 5, y = 0.002, yend = 0.1815, color = "black", linetype = "dashed", size = 1.25, arrow = arrow(ends = "both", angle = 180, length = unit(.2,"cm"))) 
o2levelsplot <- o2levelsplot + annotate("segment", x = 15, xend = 15, y = 0.002, yend = 0.030, color = "black", linetype = "dashed", size = 1.25, arrow = arrow(ends = "both", angle = 180, length = unit(.2,"cm"))) 
o2levelsplot

#Adjust appearance
o2levelsplot<-o2levelsplot + 
  theme(axis.title.x = element_text(color = "black", size = 14, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, face = "plain"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  theme(axis.text.x = element_text(color = "black", size = 14, face = "plain"),
        axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  theme(legend.title = element_blank(), legend.text = element_text(color = "black", size = 14, face = "plain")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")
o2levelsplot 

#Determine if HIV status modifies relationship between cluster assignment and risk of severe COVID-19
#Load packages
library(ordinal)
library(sjPlot)
library(sjmisc)
library(ggplot2)

#Rename HIV variable
covidclustassignhiv <- subset(covidclustassign, hivrdtresult==1 | hivrdtresult==0)
covidclustassignhiv[covidclustassignhiv$hivrdtresult==1, "PLWH"] <- "Yes"
covidclustassignhiv[covidclustassignhiv$hivrdtresult==0, "PLWH"] <- "No"
covidclustassignhiv$PLWH <- as.factor(covidclustassignhiv$PLWH)
levels(covidclustassignhiv$PLWH)

#Interaction between HIV status and CRS assignment for risk of severe COVID-19
severecovidhiv = clm(whoseverity_imp ~ hcc*PLWH, data = covidclustassignhiv, link = "logit")
summary(severecovidhiv)
confint(severecovidhiv, type = "Wald")
exp(coef(severecovidhiv))

severecovidhivinteraction <- plot_model(severecovidhiv, type = "int", ci.lvl=0.95)
severecovidhivinteraction <- severecovidhivinteraction + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12))  + ggtitle("Predicted probabilities of COVID-19 severity\n p-value of interaction term=0.329") + xlab("CRS Assignment") + ylab("")
severecovidhivinteraction <- severecovidhivinteraction + geom_line(size = 2)

#Reformat plot
covidclustassignhiv$whoseverity_imp <- sjmisc::rec(covidclustassignhiv$whoseverity_imp, rec = "2=Mild;3=Moderate;4=Severe", as.num = FALSE)
m2 <- clm(whoseverity_imp ~ hcc*PLWH, data = covidclustassignhiv, link = "logit")
severecovidhivinteraction2 <- plot_model(m2, type = "pred", terms = c("hcc", "PLWH", "whoseverity_imp"), ci.lvl = 0.95)
severecovidhivinteraction2 <- severecovidhivinteraction2 + theme_bw() + theme(plot.title=element_text(size=12)) + theme(text=element_text(size=12)) + ggtitle("p-value of interaction term=0.329") + xlab("CRS Assignment") + ylab("") + theme(strip.text.x = element_text(size = 12))
severecovidhivinteraction2 <- severecovidhivinteraction2 + scale_color_sjplot("system")
severecovidhivinteraction2

#Combine upset plot and interaction plot
library(ggpubr)
covidclusterupsethivinteractionplots <- ggarrange(upsetplot2, severecovidhivinteraction2, ncol = 1, nrow=2, labels = c("", ""))
covidclusterupsethivinteractionplots

#Multivariable logistic regression model for association between Delta phase COVID and CRS2
#Include 4 patients with unknown HIV status with HIV neg
covidclustassign$hivrdtresult[is.na(covidclustassign$hivrdtresult)] = 0

deltaphasecluster <- glm(hcc ~ phasedelta + age + sex + illnessdurationenroll_imp + hivrdtresult + whoseverity_imp, data = covidclustassign, family = binomial)
summary(deltaphasecluster)
exp(deltaphasecluster$coefficients)
exp(confint.default(deltaphasecluster))

