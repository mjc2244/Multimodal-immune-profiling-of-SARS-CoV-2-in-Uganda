#Figures 4b and 4c - XGBoost and SHAP mediator importance for Delta phase COVID

#Clear R environment
rm(list = ls())

#Import master data, PID to fixed row identifier 
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid
covid <- subset(combined, pathogencode==1)

#Load packages
library(readr)
library(pROC)
library(xgboost)
library(dplyr)
library(caret)
library(ggplot2)

#Create matrix of relevant biomarker data for input to XGBoost - exclude GM-CSF, G-CSF, and IL-3 given large prop. values <LLD
covidxgbmatrix <- data.matrix(covid[1:306, c(237:242,245:251,253:284)], rownames.force = NA)
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
                              "TNF-\u03b1",
                              "LT-\u03b1",
                              "VEGF-A")
head(covidxgbmatrix)

#Create outcome for variable to be classified/predicted, in this case Delta phase COVID
predictdelta <- covid[,"phasedelta"] == "1"

## Performing a grid search for Hyperparameter Tuning
### Tuned Parameters were "nrounds", "eta", and "max_depth"
grid_label <- factor(predictdelta, labels = c("nondelta", "delta"))
boost_train_cont = trainControl(method = "cv", number = 10, verboseIter = TRUE, returnData = FALSE, returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
boost_grid <- expand.grid(nrounds = c(10, 100, 200, 400, 800, 1000, 2000),eta = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.001),max_depth = c(2, 3, 4, 6, 8, 10), gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1)
set.seed(12345)
model_train <- train(x=covidxgbmatrix, y=grid_label, trControl = boost_train_cont, tuneGrid = boost_grid, method = "xgbTree")
ggplot(model_train$results, aes(x = as.factor(eta), y = max_depth, size = ROC, color = ROC)) + geom_point() + theme_bw() + scale_size_continuous(guide = "none")
options(max.print = 100000)
model_train$results
model_train$results[order(model_train$results$ROC),] 

#XGBoost to predict Delta COVID-19
library(xgboost)
params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.1, gamma=0, max_depth=3, min_child_weight=1, subsample=1, colsample_bytree=1)
bst <- xgboost(data = covidxgbmatrix, label = predictdelta, nround = 400, params = params)                     

#Variable importance in predicting cluster assignment 
importance <- xgb.importance(feature_names = colnames(covidxgbmatrix), model = bst)
head(importance)

#Plot variable importance - 10 most important variables 
library(Ckmeans.1d.dp)
library(ggplot2)
importanceplot <- xgb.ggplot.importance(importance_matrix = importance, rel_to_first = FALSE, top_n = 20, n_clusters = c(1), xlab = "Split-Gain")

importanceplot <- importanceplot+ 
  scale_fill_manual(values=c("mediumpurple4")) + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_text(colour = "black", size = 12)) + 
  theme(axis.text.y = element_text(color = "black", size = 12)) +
  theme(axis.text.x = element_text(color = "black", size = 12)) +
  labs(title = "", y = "Split-Gain", x = "", color="black") + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = "blank")) + 
  theme(panel.grid.minor = element_line(size = 0.25, linetype = "blank"))
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

#Combine plots
library(ggpubr)
ggarrange(importanceplot, shap, ncol = 1, nrow = 2)



