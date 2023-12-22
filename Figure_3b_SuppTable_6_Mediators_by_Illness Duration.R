#Figure 3b and Table S7 - Mediators by severe vs. non-severe COVID enrolled <=7d after illness onset

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, make row ID fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients, symptomatic, enrolled <=7d after illness onset 
covid <- subset(combined, pathogencode == 1 & min_no_sym == 0 & illnessdurationenroll_imp <=7)

#Biomarker concentrations across illness severity groups
library(ggpubr)
library(rstatix)

#sCD40L
scd40l <- ggviolin(covid, x = "whosevere_imp", y = "logscd40l", 
          title = "sCD40L",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "sCD40L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5) 
scd40l

#EGF
egf <- ggviolin(covid, x = "whosevere_imp", y = "logegf", 
          title = "EGF",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "EGF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
egf

#Eotaxin
eotaxin <- ggviolin(covid, x = "whosevere_imp", y = "logeotaxin", 
          title = "Eotaxin",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "Eotaxin (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
eotaxin

#FGF2
fgf2 <- ggviolin(covid, x = "whosevere_imp", y = "logfgf2", 
          title = "FGF-2",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FGF-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
fgf2

#FLT-3L
flt3l <- ggviolin(covid, x = "whosevere_imp", y = "logflt3l", 
          title = "FLT-3L",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FLT-3L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
flt3l

#FKN
fkn <- ggviolin(covid, x = "whosevere_imp", y = "logfractalkine", 
          title = "Fractalkine",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Fractalkine (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
fkn

#G-CSF
gcsf<- ggviolin(covid, x = "whosevere_imp", y = "loggcsf", 
          title = "G-CSF",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "G-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.5, size=4.5)
gcsf

#GM-CSF
gmcsf <- ggviolin(covid, x = "whosevere_imp", y = "loggmcsf", 
          title = "GM-CSF",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "GM-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5)
gmcsf

#GRO-alpha/CXCL1
groalpha <- ggviolin(covid, x = "whosevere_imp", y = "loggroalpha", 
          title = "CXCL1",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
groalpha

#IFN-a2
ifna2 <- ggviolin(covid, x = "whosevere_imp", y = "logifna2", 
          title = "IFN-\u03b12",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
ifna2

#IFN-gamma
ifngamma <- ggviolin(covid, x = "whosevere_imp", y = "logifngamma", 
          title = "IFN-\u03b3",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
ifngamma

#IL-1a
il1a <- ggviolin(covid, x = "whosevere_imp", y = "logil1a", 
          title = "IL-1\u03b1",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.25, size=4.5)
il1a

#IL-1b
il1b <- ggviolin(covid, x = "whosevere_imp", y = "logil1b", 
          title = "IL-1\u03B2",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
il1b

#IL-1Ra
il1ra <- ggviolin(covid, x = "whosevere_imp", y = "logil1ra", 
          title = "IL-1Ra",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1Ra (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
il1ra

#IL-2
il2 <- ggviolin(covid, x = "whosevere_imp", y = "logil2", 
          title = "IL-2",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.75, size=4.5)
il2

#IL-3
il3 <- ggviolin(covid, x = "whosevere_imp", y = "logil3", 
          title = "IL-3",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 2, size=4.5)
il3

#IL-4
il4 <- ggviolin(covid, x = "whosevere_imp", y = "logil4", 
          title = "IL-4",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 2.25, size=4.5)
il4

#IL-5
il5 <- ggviolin(covid, x = "whosevere_imp", y = "logil5", 
          title = "IL-5",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-5 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 2.5, size=4.5)
il5

#IL-6
il6 <- ggviolin(covid, x = "whosevere_imp", y = "logil6", 
          title = "IL-6",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-6 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
il6

#IL-7
il7 <- ggviolin(covid, x = "whosevere_imp", y = "logil7", 
          title = "IL-7",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-7 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 2.5, size=4.5)
il7

#IL-8
il8 <- ggviolin(covid, x = "whosevere_imp", y = "logil8", 
          title = "IL-8",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-8 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.5, size=4.5)
il8

#IL-9
il9 <- ggviolin(covid, x = "whosevere_imp", y = "logil9", 
          title = "IL-9",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
il9

#IL-10
il10 <- ggviolin(covid, x = "whosevere_imp", y = "logil10", 
          title = "IL-10",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
il10

#IL-12p40
il12p40 <- ggviolin(covid, x = "whosevere_imp", y = "logil12p40", 
          title = "IL-12p40",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p40 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
il12p40

#IL-12p70
il12p70 <- ggviolin(covid, x = "whosevere_imp", y = "logil12p70", 
          title = "IL-12p70",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p70 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
il12p70

#IL-13
il13 <- ggviolin(covid, x = "whosevere_imp", y = "logil13", 
          title = "IL-13",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-13 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
il13

#IL-15
il15 <- ggviolin(covid, x = "whosevere_imp", y = "logil15", 
          title = "IL-15",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-15 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
il15

#IL-17A
il17a <- ggviolin(covid, x = "whosevere_imp", y = "logil17a", 
          title = "IL-17A",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
il17a

#IL-17E
il17e <- ggviolin(covid, x = "whosevere_imp", y = "logil17eil25", 
          title = "IL-17E",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17E (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5)
il17e

#IL-17F
il17f <- ggviolin(covid, x = "whosevere_imp", y = "logil17f", 
          title = "IL-17F",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17F (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5)
il17f

#IL-18
il18 <- ggviolin(covid, x = "whosevere_imp", y = "logil18", 
          title = "IL-18",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-18 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
il18

#IL-22
il22 <- ggviolin(covid, x = "whosevere_imp", y = "logil22", 
          title = "IL-22",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-22 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
il22

#IL-27
il27 <- ggviolin(covid, x = "whosevere_imp", y = "logil27", 
          title = "IL-27",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-27 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5)
il27

#IP-10/CXCL10
ip10 <- ggviolin(covid, x = "whosevere_imp", y = "logip10", 
          title = "CXCL10",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 6, size=4.5)
ip10

#MCP-1
mcp1 <- ggviolin(covid, x = "whosevere_imp", y = "logmcp1", 
          title = "MCP-1",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
mcp1

#MCP-3
mcp3 <- ggviolin(covid, x = "whosevere_imp", y = "logmcp3", 
          title = "MCP-3",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3, size=4.5)
mcp3

#M-CSF
mcsf <- ggviolin(covid, x = "whosevere_imp", y = "logmcsf", 
          title = "M-CSF",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "M-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
mcsf

#MDC
mdc <- ggviolin(covid, x = "whosevere_imp", y = "logmdc", 
          title = "MDC",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MDC (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
mdc

#MIG/CXCL9
migcxcl9 <- ggviolin(covid, x = "whosevere_imp", y = "logmigcxcl9", 
          title = "CXCL9",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5, size=4.5)
migcxcl9

#MIP-1a/CCL3
mip1a <- ggviolin(covid, x = "whosevere_imp", y = "logmip1a", 
          title = "CCL3",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
mip1a

#MIP-1b/CCL4
mip1b <- ggviolin(covid, x = "whosevere_imp", y = "logmip1b", 
          title = "CCL4",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.25, size=4.5)
mip1b

#PDGF-AA
pdgfaa <- ggviolin(covid, x = "whosevere_imp", y = "logpdgfaa", 
          title = "PDGF-AA",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AA (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.5, size=4.5)
pdgfaa

#PDGF-AB/BB
pdgfabbb <- ggviolin(covid, x = "whosevere_imp", y = "logpdgfabbb", 
          title = "PDGF-AB/BB",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 5.35, size=4.5)
pdgfabbb

#RANTES
rantes <- ggviolin(covid, x = "whosevere_imp", y = "lograntes", 
          title = "RANTES",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "RANTES (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.5, size=4.5)
rantes

#TGF-alpha
tgfa <- ggviolin(covid, x = "whosevere_imp", y = "logtgfa", 
          title = "TGF-\u03b1",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4.25, size=4.5)
tgfa

#TNF-alpha
tnfa <- ggviolin(covid, x = "whosevere_imp", y = "logtnfa", 
          title = "TNF",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TNF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.5, size=4.5)
tnfa

#Lymphotoxin-alpha/TNF-beta
lta <- ggviolin(covid, x = "whosevere_imp", y = "logtnfb", 
          title = "Lymphotoxin-\u03b1",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 3.25, size=4.5)
lta

#VEGF-A
vegfa <- ggviolin(covid, x = "whosevere_imp", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "whosevere_imp", alpha = 0.5, palette =  c("#F39B7FFF","#8491B4FF"),
          order = c("0", "1"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("0" = "Mild-Moderate", "1" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means(label.y = 4, size=4.5)
vegfa

#Tables showing these values
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

tablebiomarkers <- CreateTableOne(vars = listbiomarkers, 
                                            data = covid, 
                                            strata = "whosevere_imp",
                                            addOverall = TRUE)
tablebiomarkers
summary(tablebiomarkers)

tablebiomarkers <- print(tablebiomarkers, nonnormal = c("logscd40l", 
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

write_tableHTML(tableHTML(tablebiomarkers), file = 'biomarkersseverity7d.html')




