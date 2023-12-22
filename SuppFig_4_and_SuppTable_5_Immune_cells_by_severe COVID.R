#Figure S4 and Table S6 - Inferred immune cell abundance by severe vs non-severe COVID

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients with RNAseq data
covid <- subset(combined, pathogencode == 1 & rna_seq ==1)

#Immune cell abundance across COVID severity groups
library(ggpubr)
library(rstatix)

#Naive B cells
naivebcells <- ggviolin(covid, x = "whosevere_imp", y = "bcellsnaive", 
                   title = "Naive B cells",
                   color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                   order = c("0", "1"),
                   add = "jitter", draw_quantiles = 0.5,
                   ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 10, label.x = 1.5, size=4.5) 
naivebcells

#Memory B cells
memorybcells <- ggviolin(covid, x = "whosevere_imp", y = "bcellsmemory", 
                        title = "Memory B cells",
                        color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                        order = c("0", "1"),
                        add = "jitter", draw_quantiles = 0.5,
                        ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 2, label.x = 1.5, size=4.5) 
memorybcells

#Plasma cells
plasmacells <- ggviolin(covid, x = "whosevere_imp", y = "plasmacells", 
                        title = "Plasma cells",
                        color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                        order = c("0", "1"),
                        add = "jitter", draw_quantiles = 0.5,
                        ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 7.5, label.x = 1.5, size=4.5) 
plasmacells

#CD8+ T cells
tcellscd8 <- ggviolin(covid, x = "whosevere_imp", y = "tcellscd8", 
                         title = "CD8+ T cells",
                         color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                         order = c("0", "1"),
                         add = "jitter", draw_quantiles = 0.5,
                         ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 10, label.x = 1.5, size=4.5) 
tcellscd8

#CD4+ T cells - Naive 
tcellscd4naive <- ggviolin(covid, x = "whosevere_imp", y = "tcellscd4naive", 
                      title = "Naive CD4+ T cells",
                      color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                      order = c("0", "1"),
                      add = "jitter", draw_quantiles = 0.5,
                      ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 20, label.x = 1.5, size=4.5) 
tcellscd4naive

#CD4+ T cells - resting memory 
tcellscd4resting <- ggviolin(covid, x = "whosevere_imp", y = "tcellscd4memoryresting", 
                           title = "Resting memory CD4+ T cells",
                           color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                           order = c("0", "1"),
                           add = "jitter", draw_quantiles = 0.5,
                           ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 25, label.x = 1.5, size=4.5) 
tcellscd4resting

#CD4+ T cells - activated memory
tcellscd4activated <- ggviolin(covid, x = "whosevere_imp", y = "tcellscd4memoryactivated", 
                           title = "Activated memory CD4+ T cells",
                           color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                           order = c("0", "1"),
                           add = "jitter", draw_quantiles = 0.5,
                           ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, label.x = 1.5, size=4.5) 
tcellscd4activated

#Resting NK cells
nkcellsresting <- ggviolin(covid, x = "whosevere_imp", y = "nkcellsresting", 
                             title = "Resting NK cells",
                             color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                             order = c("0", "1"),
                             add = "jitter", draw_quantiles = 0.5,
                             ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 25, label.x = 1.5, size=4.5) 
nkcellsresting

#Monocytes
monocytes <- ggviolin(covid, x = "whosevere_imp", y = "monocytes", 
                           title = "Monocytes",
                           color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                           order = c("0", "1"),
                           add = "jitter", draw_quantiles = 0.5,
                           ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 22, label.x = 1.5, size=4.5) 
monocytes

#Activated DCs
activateddendritic <- ggviolin(covid, x = "whosevere_imp", y = "dendriticcellsactivated", 
                           title = "Activated dendritic cells",
                           color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                           order = c("0", "1"),
                           add = "jitter", draw_quantiles = 0.5,
                           ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 0.75, label.x = 1.5, size=4.5) 
activateddendritic

#Resting mast cells
restingmastcells <- ggviolin(covid, x = "whosevere_imp", y = "mastcellsresting", 
                               title = "Resting mast cells",
                               color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                               order = c("0", "1"),
                               add = "jitter", draw_quantiles = 0.5,
                               ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, label.x = 1.5, size=4.5) 
restingmastcells

#Eosinophils
eosinophils <- ggviolin(covid, x = "whosevere_imp", y = "eosinophils", 
                             title = "Eosinophils",
                             color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                             order = c("0", "1"),
                             add = "jitter", draw_quantiles = 0.5,
                             ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, label.x = 1.5, size=4.5) 
eosinophils

#Neutrophils
neutrophils <- ggviolin(covid, x = "whosevere_imp", y = "neutrophils", 
                        title = "Neutrophils",
                        color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
                        order = c("0", "1"),
                        add = "jitter", draw_quantiles = 0.5,
                        ylab = "Inferred abundance", xlab = "") +
  scale_x_discrete(labels=c("0" = "Non-severe", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 125, label.x = 1.5, size=4.5) 
neutrophils

#Tables showing these values
library(tableone)
library(tableHTML)
#Create a variable list 
listcells <- c("bcellsnaive",
                    "bcellsmemory",
                    "plasmacells",
                    "tcellscd8",
                    "tcellscd4naive",
                    "tcellscd4memoryresting",
               "tcellscd4memoryactivated",
                    "nkcellsresting",
                    "monocytes",
                    "dendriticcellsactivated",
               "mastcellsresting",
               "eosinophils",
               "neutrophils")

tablecells <- CreateTableOne(vars = listcells, 
                                            data = covid, 
                                            strata = "whosevere_imp",
                                            addOverall = TRUE)
tablecells
summary(tablecells)

tablecells <- print(tablecells, nonnormal = c("bcellsnaive",
                                              "bcellsmemory",
                                              "plasmacells",
                                              "tcellscd8",
                                              "tcellscd4naive",
                                              "tcellscd4memoryresting",
                                              "tcellscd4memoryactivated",
                                              "nkcellsresting",
                                              "monocytes",
                                              "dendriticcellsactivated",
                                              "mastcellsresting",
                                              "eosinophils",
                                              "neutrophils"))

write_tableHTML(tableHTML(tablecells), file = 'immunecellsseverityLM22.html')

#Combine plots
library(ggpubr)

immunecellsbyseverity <- ggarrange(naivebcells, memorybcells, plasmacells, tcellscd8, tcellscd4naive, tcellscd4resting, tcellscd4activated, nkcellsresting, monocytes, activateddendritic, restingmastcells, eosinophils, neutrophils, ncol = 3, nrow =5)
immunecellsbyseverity




