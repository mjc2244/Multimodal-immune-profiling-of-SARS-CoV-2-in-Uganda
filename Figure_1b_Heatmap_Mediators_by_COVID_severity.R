#Figure 1b - Heatmap - Soluble mediators stratified by COVID-19 severity 

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, PID to fixed row ID
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid patients
covid <- subset(combined, pathogencode==1 )

#Clinical variables to factors for annotations
#Sex
covid$sex <- factor(covid$sex)
levels(covid$sex)

#HIV
covid$hivrdtresult <- factor(covid$hivrdtresult)
levels(covid$hivrdtresult)

#Study Phase
covid$phase3 <- factor(covid$phase3)
levels(covid$phase3)

#WHO clinical severity
covid[covid$whoseverity_imp==1, "whoseverity_imp2"] <- "Asymptomatic"
covid[covid$whoseverity_imp==2, "whoseverity_imp2"] <- "Mild"
covid[covid$whoseverity_imp==3, "whoseverity_imp2"] <- "Mod."
covid[covid$whoseverity_imp==4, "whoseverity_imp2"] <- "Severe"
covid$whoseverity_imp2 <- as.factor(covid$whoseverity_imp2)
levels(covid$whoseverity_imp2)

#Create a new matrix of biomarkers with labels - exclude G-CSF, GM-CSF, IL-3 given large prop. of values < LLD
biomarkers <- covid[1:306, c(c(237:242,245:251,253:284))]
head(biomarkers)

#Scale/center biomarkers, convert to dataframe 
biomscaled <- scale(biomarkers, center=TRUE, scale=TRUE)
biomscaled <- as.data.frame(biomscaled)

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
                                    "CCL5",
                                    "TGF-\u03b1",
                                    "TNF",
                                    "LT-\u03b1",
                                    "VEGF-A")
head(biomscaledheatmatrix)

#Plot
library(ComplexHeatmap)

library(RColorBrewer)
display.brewer.pal(9, "Oranges")
brewer.pal(9, "Oranges")

display.brewer.pal(9, "Greens")
brewer.pal(9, "Greens")

brewer.pal(n=5,"Set1")

#Make annotation
library(circlize)
col_age = colorRamp2(c(10, 20, 30, 40, 50, 60, 70, 80, 90), c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))
col_duration = colorRamp2(c(0,3,6,9,12,15,18,21,24), c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"))

ha = HeatmapAnnotation("Age" = covid$age[1:306],
                       "Sex" = covid$sex[1:306],
                       "Living with HIV" = covid$hivrdtresult[1:306],
                       "Dominant variant(s)" = covid$phase3[1:306],
                       show_legend = c("HIV infection" = TRUE, "Age" = TRUE, "Sex" = TRUE, "Dominant variant(s)" = TRUE),
                       annotation_legend_param = list("Age" = list(title = "Age", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)), "Sex" = list(title = "Sex", at = c("1", "0"), labels = c("Male", "Female"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
                                                      "Living with HIV" = list(title = "Living with HIV", at = c("1", "0"), labels = c("Yes", "No"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
                                                      "Dominant variant(s)" = list(title = "Dominant variant(s)", at = c("1", "2", "3"), labels = c("Varied A/B Lineages", "A.23/A.23.1", "Delta (B.1.617.2/AY)"), labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14))),
                       col = list(
                         "Sex" = c("1" = "#143d59", "0" = "#f4b41a"),
                         "Living with HIV" = c("1" = "#E2725A", "0" = "#4A6274"),
                         "Dominant variant(s)" = c("1" = "#440154", "2" = "#2a788e", "3" = "#7ad151"),
                         "Age" = col_age),
                       na_col = "#4A6274",
                       annotation_name_gp = gpar(fontsize = 14))

library(circlize)
col_map1 = colorRamp2(c(-1.5, -1, 0, 1, 1.5), c("#000004", "#51127c", "#b73779", "#fc8961", "#fcfdbf"))

set.seed(123)
severityheatmap <- Heatmap(t(biomscaledheatmatrix), 
        name = "Z-score", 
        col = col_map1,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_km = 4,
        show_column_names = FALSE, 
        show_column_dend = FALSE, 
        show_row_dend = FALSE, 
        row_names_side = "left", 
        column_split = covid$whoseverity_imp2, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14)),
        column_title_gp = gpar(fontsize = 14),
        cluster_column_slices = FALSE,
        row_title = c(),
        row_names_gp = gpar(fontsize = 14),  
        top_annotation = ha,
        border=TRUE)
severityheatmap

