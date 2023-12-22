#Figure 1c and Table S1 - Mediators by COVID severity_KW

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients
covid <- subset(combined, pathogencode == 1)

#Biomarker concentrations across illness severity groups, including asymptomatic controls
library(ggpubr)
library(rstatix)

#sCD40L
ggviolin(covid, x = "whoseverity_imp", y = "logscd40l", 
          title = "sCD40L",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "sCD40L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalscd40l <- covid %>% kruskal_test(logscd40l ~ whoseverity_imp)
res.kruskalscd40l

pwcscd40l <- covid %>% dunn_test(logscd40l ~ whoseverity_imp, p.adjust.method = "BH") 
pwcscd40l

pwcscd40l <- pwcscd40l %>% add_xy_position(x = "whoseverity_imp")
pwcscd40lplot <- ggviolin(covid, x = "whoseverity_imp", y = "logscd40l", 
         title = "sCD40L",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "sCD40L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcscd40l, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalscd40l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcscd40lplot

#EGF
ggviolin(covid, x = "whoseverity_imp", y = "logegf", 
          title = "EGF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "EGF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalegf <- covid %>% kruskal_test(logegf ~ whoseverity_imp)
res.kruskalegf

pwcegf <- covid %>% dunn_test(logegf ~ whoseverity_imp, p.adjust.method = "BH") 
pwcegf

pwcegf <- pwcegf %>% add_xy_position(x = "whoseverity_imp")
pwcegfplot <- ggviolin(covid, x = "whoseverity_imp", y = "logegf", 
         title = "EGF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "EGF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcegf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalegf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcegfplot

#Eotaxin
ggviolin(covid, x = "whoseverity_imp", y = "logeotaxin", 
          title = "Eotaxin",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "Eotaxin (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskaleotaxin <- covid %>% kruskal_test(logeotaxin ~ whoseverity_imp)
res.kruskaleotaxin

pwceotaxin <- covid %>% dunn_test(logeotaxin ~ whoseverity_imp, p.adjust.method = "BH") 
pwceotaxin

pwceotaxin <- pwceotaxin %>% add_xy_position(x = "whoseverity_imp")
pwceotaxinplot <- ggviolin(covid, x = "whoseverity_imp", y = "logeotaxin", 
         title = "Eotaxin",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Eotaxin (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwceotaxin, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaleotaxin, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwceotaxinplot

#FGF2
ggviolin(covid, x = "whoseverity_imp", y = "logfgf2", 
          title = "FGF-2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FGF-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalfgf2 <- covid %>% kruskal_test(logfgf2 ~ whoseverity_imp)
res.kruskalfgf2

pwcfgf2 <- covid %>% dunn_test(logfgf2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcfgf2

pwcfgf2 <- pwcfgf2 %>% add_xy_position(x = "whoseverity_imp")
pwcfgf2plot <- ggviolin(covid, x = "whoseverity_imp", y = "logfgf2", 
         title = "FGF-2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FGF-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcfgf2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfgf2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcfgf2plot

#FLT-3L
ggviolin(covid, x = "whoseverity_imp", y = "logflt3l", 
          title = "FLT-3L",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FLT-3L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalflt3l <- covid %>% kruskal_test(logflt3l ~ whoseverity_imp)
res.kruskalflt3l

pwcflt3l <- covid %>% dunn_test(logflt3l ~ whoseverity_imp, p.adjust.method = "BH") 
pwcflt3l

pwcflt3l <- pwcflt3l %>% add_xy_position(x = "whoseverity_imp")
pwcflt3lplot <- ggviolin(covid, x = "whoseverity_imp", y = "logflt3l", 
         title = "FLT-3L",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FLT-3L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcflt3l, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalflt3l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcflt3lplot

#FKN
ggviolin(covid, x = "whoseverity_imp", y = "logfractalkine", 
          title = "Fractalkine",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Fractalkine (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalfractalkine <- covid %>% kruskal_test(logfractalkine ~ whoseverity_imp)
res.kruskalfractalkine

pwcfractalkine <- covid %>% dunn_test(logfractalkine ~ whoseverity_imp, p.adjust.method = "BH") 
pwcfractalkine

pwcfractalkine <- pwcfractalkine %>% add_xy_position(x = "whoseverity_imp")
pwcfractalkineplot <- ggviolin(covid, x = "whoseverity_imp", y = "logfractalkine", 
         title = "Fractalkine",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Fractalkine (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcfractalkine, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfractalkine, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcfractalkineplot

#G-CSF
ggviolin(covid, x = "whoseverity_imp", y = "loggcsf", 
          title = "G-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "G-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgcsf <- covid %>% kruskal_test(loggcsf ~ whoseverity_imp)
res.kruskalgcsf

pwcgcsf <- covid %>% dunn_test(loggcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwcgcsf

pwcgcsf <- pwcgcsf %>% add_xy_position(x = "whoseverity_imp")
pwcgcsfplot <- ggviolin(covid, x = "whoseverity_imp", y = "loggcsf", 
         title = "G-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "G-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcgcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcgcsfplot

#GM-CSF
ggviolin(covid, x = "whoseverity_imp", y = "loggmcsf", 
          title = "GM-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "GM-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgmcsf <- covid %>% kruskal_test(loggmcsf ~ whoseverity_imp)
res.kruskalgmcsf

pwcgmcsf <- covid %>% dunn_test(loggmcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwcgmcsf

pwcgmcsf <- pwcgmcsf %>% add_xy_position(x = "whoseverity_imp")
pwcgmcsfplot <- ggviolin(covid, x = "whoseverity_imp", y = "loggmcsf", 
         title = "GM-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "GM-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcgmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcgmcsfplot

#GRO-alpha/CXCL1
ggviolin(covid, x = "whoseverity_imp", y = "loggroalpha", 
          title = "CXCL1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgroalpha <- covid %>% kruskal_test(loggroalpha ~ whoseverity_imp)
res.kruskalgroalpha

pwcgroalpha <- covid %>% dunn_test(loggroalpha ~ whoseverity_imp, p.adjust.method = "BH") 
pwcgroalpha

pwcgroalpha <- pwcgroalpha %>% add_xy_position(x = "whoseverity_imp")
pwcgroalphaplot <- ggviolin(covid, x = "whoseverity_imp", y = "loggroalpha", 
         title = "CXCL1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcgroalpha, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgroalpha, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcgroalphaplot

#IFN-a2
ggviolin(covid, x = "whoseverity_imp", y = "logifna2", 
          title = "IFN-\u03b12",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalifna2 <- covid %>% kruskal_test(logifna2 ~ whoseverity_imp)
res.kruskalifna2

pwcifna2 <- covid %>% dunn_test(logifna2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcifna2

pwcifna2 <- pwcifna2 %>% add_xy_position(x = "whoseverity_imp")
pwcifna2plot <- ggviolin(covid, x = "whoseverity_imp", y = "logifna2", 
         title = "IFN-\u03b12",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcifna2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifna2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcifna2plot

#IFN-gamma
ggviolin(covid, x = "whoseverity_imp", y = "logifngamma", 
          title = "IFN-\u03b3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalifngamma <- covid %>% kruskal_test(logifngamma ~ whoseverity_imp)
res.kruskalifngamma

pwcifngamma <- covid %>% dunn_test(logifngamma ~ whoseverity_imp, p.adjust.method = "BH") 
pwcifngamma

pwcifngamma <- pwcifngamma %>% add_xy_position(x = "whoseverity_imp")
pwcifngammaplot <- ggviolin(covid, x = "whoseverity_imp", y = "logifngamma", 
         title = "IFN-\u03b3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcifngamma, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifngamma, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcifngammaplot

#IL-1a
ggviolin(covid, x = "whoseverity_imp", y = "logil1a", 
          title = "IL-1\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1a <- covid %>% kruskal_test(logil1a ~ whoseverity_imp)
res.kruskalil1a

pwcil1a <- covid %>% dunn_test(logil1a ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil1a

pwcil1a <- pwcil1a %>% add_xy_position(x = "whoseverity_imp")
pwcil1aplot <- ggviolin(covid, x = "whoseverity_imp", y = "logil1a", 
         title = "IL-1\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil1aplot

#IL-1b
ggviolin(covid, x = "whoseverity_imp", y = "logil1b", 
          title = "IL-1\u03B2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1b <- covid %>% kruskal_test(logil1b ~ whoseverity_imp)
res.kruskalil1b

pwcil1b <- covid %>% dunn_test(logil1b ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil1b

pwcil1b <- pwcil1b %>% add_xy_position(x = "whoseverity_imp")
pwcil1bplot <- ggviolin(covid, x = "whoseverity_imp", y = "logil1b", 
         title = "IL-1\u03B2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil1bplot

#IL-1Ra
ggviolin(covid, x = "whoseverity_imp", y = "logil1ra", 
          title = "IL-1Ra",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1Ra (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1ra <- covid %>% kruskal_test(logil1ra ~ whoseverity_imp)
res.kruskalil1ra

pwcil1ra <- covid %>% dunn_test(logil1ra ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil1ra

pwcil1ra <- pwcil1ra %>% add_xy_position(x = "whoseverity_imp")
pwcil1raplot <- ggviolin(covid, x = "whoseverity_imp", y = "logil1ra", 
         title = "IL-1Ra",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1Ra (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil1ra, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1ra, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil1raplot

#IL-2
ggviolin(covid, x = "whoseverity_imp", y = "logil2", 
          title = "IL-2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil2 <- covid %>% kruskal_test(logil2 ~ whoseverity_imp)
res.kruskalil2

pwcil2 <- covid %>% dunn_test(logil2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil2

pwcil2 <- pwcil2 %>% add_xy_position(x = "whoseverity_imp")
pwcil2plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil2", 
         title = "IL-2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil2plot 

#IL-3
ggviolin(covid, x = "whoseverity_imp", y = "logil3", 
          title = "IL-3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil3 <- covid %>% kruskal_test(logil3 ~ whoseverity_imp)
res.kruskalil3

pwcil3 <- covid %>% dunn_test(logil3 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil3

pwcil3 <- pwcil3 %>% add_xy_position(x = "whoseverity_imp")
pwcil3plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil3", 
         title = "IL-3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil3, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil3plot

#IL-4
ggviolin(covid, x = "whoseverity_imp", y = "logil4", 
          title = "IL-4",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil4 <- covid %>% kruskal_test(logil4 ~ whoseverity_imp)
res.kruskalil4

pwcil4 <- covid %>% dunn_test(logil4 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil4

pwcil4 <- pwcil4 %>% add_xy_position(x = "whoseverity_imp")
pwcil4plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil4", 
         title = "IL-4",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil4, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil4, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil4plot

#IL-5
ggviolin(covid, x = "whoseverity_imp", y = "logil5", 
          title = "IL-5",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-5 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil5 <- covid %>% kruskal_test(logil5 ~ whoseverity_imp)
res.kruskalil5

pwcil5 <- covid %>% dunn_test(logil5 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil5

pwcil5 <- pwcil5 %>% add_xy_position(x = "whoseverity_imp")
pwcil5plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil5", 
         title = "IL-5",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-5 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil5, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0) + 
  labs(
    subtitle = get_test_label(res.kruskalil5, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil5plot

#IL-6
ggviolin(covid, x = "whoseverity_imp", y = "logil6", 
          title = "IL-6",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-6 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil6 <- covid %>% kruskal_test(logil6 ~ whoseverity_imp)
res.kruskalil6

pwcil6 <- covid %>% dunn_test(logil6 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil6

pwcil6 <- pwcil6 %>% add_xy_position(x = "whoseverity_imp")
pwcil6plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil6", 
         title = "IL-6",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-6 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil6, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.8, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil6, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil6plot

#IL-7
ggviolin(covid, x = "whoseverity_imp", y = "logil7", 
          title = "IL-7",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-7 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil7 <- covid %>% kruskal_test(logil7 ~ whoseverity_imp)
res.kruskalil7

pwcil7 <- covid %>% dunn_test(logil7 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil7

pwcil7 <- pwcil7 %>% add_xy_position(x = "whoseverity_imp")
pwcil7plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil7", 
         title = "IL-7",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-7 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil7, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil7, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil7plot

#IL-8
ggviolin(covid, x = "whoseverity_imp", y = "logil8", 
          title = "IL-8",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-8 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil8 <- covid %>% kruskal_test(logil8 ~ whoseverity_imp)
res.kruskalil8

pwcil8 <- covid %>% dunn_test(logil8 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil8

pwcil8 <- pwcil8 %>% add_xy_position(x = "whoseverity_imp")
pwcil8plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil8", 
         title = "IL-8",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-8 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil8, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil8, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil8plot

#IL-9
ggviolin(covid, x = "whoseverity_imp", y = "logil9", 
          title = "IL-9",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil9 <- covid %>% kruskal_test(logil9 ~ whoseverity_imp)
res.kruskalil9

pwcil9 <- covid %>% dunn_test(logil9 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil9

pwcil9 <- pwcil9 %>% add_xy_position(x = "whoseverity_imp")
pwcil9plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil9", 
         title = "IL-9",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil9, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil9plot

#IL-10
ggviolin(covid, x = "whoseverity_imp", y = "logil10", 
          title = "IL-10",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil10 <- covid %>% kruskal_test(logil10 ~ whoseverity_imp)
res.kruskalil10

pwcil10 <- covid %>% dunn_test(logil10 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil10

pwcil10 <- pwcil10 %>% add_xy_position(x = "whoseverity_imp")
pwcil10plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil10", 
         title = "IL-10",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil10, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil10plot

#IL-12p40
ggviolin(covid, x = "whoseverity_imp", y = "logil12p40", 
          title = "IL-12p40",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p40 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil12p40 <- covid %>% kruskal_test(logil12p40 ~ whoseverity_imp)
res.kruskalil12p40

pwcil12p40 <- covid %>% dunn_test(logil12p40 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil12p40

pwcil12p40 <- pwcil12p40 %>% add_xy_position(x = "whoseverity_imp")
pwcil12p40plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil12p40", 
         title = "IL-12p40",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p40 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil12p40, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p40, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil12p40plot

#IL-12p70
ggviolin(covid, x = "whoseverity_imp", y = "logil12p70", 
          title = "IL-12p70",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p70 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil12p70 <- covid %>% kruskal_test(logil12p70 ~ whoseverity_imp)
res.kruskalil12p70

pwcil12p70 <- covid %>% dunn_test(logil12p70 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil12p70

pwcil12p70 <- pwcil12p70 %>% add_xy_position(x = "whoseverity_imp")
pwcil12p70plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil12p70", 
         title = "IL-12p70",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p70 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil12p70, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p70, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil12p70plot

#IL-13
ggviolin(covid, x = "whoseverity_imp", y = "logil13", 
          title = "IL-13",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-13 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil13 <- covid %>% kruskal_test(logil13 ~ whoseverity_imp)
res.kruskalil13

pwcil13 <- covid %>% dunn_test(logil13 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil13

pwcil13 <- pwcil13 %>% add_xy_position(x = "whoseverity_imp")
pwcil13plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil13", 
         title = "IL-13",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-13 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil13, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil13, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil13plot

#IL-15
ggviolin(covid, x = "whoseverity_imp", y = "logil15", 
          title = "IL-15",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-15 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil15 <- covid %>% kruskal_test(logil15 ~ whoseverity_imp)
res.kruskalil15

pwcil15 <- covid %>% dunn_test(logil15 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil15

pwcil15 <- pwcil15 %>% add_xy_position(x = "whoseverity_imp")
pwcil15plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil15", 
         title = "IL-15",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-15 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil15, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil15, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil15plot

#IL-17A
ggviolin(covid, x = "whoseverity_imp", y = "logil17a", 
          title = "IL-17A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17a <- covid %>% kruskal_test(logil17a ~ whoseverity_imp)
res.kruskalil17a

pwcil17a <- covid %>% dunn_test(logil17a ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil17a

pwcil17a <- pwcil17a %>% add_xy_position(x = "whoseverity_imp")
pwcil17aplot <- ggviolin(covid, x = "whoseverity_imp", y = "logil17a", 
         title = "IL-17A",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil17a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil17aplot

#IL-17E
ggviolin(covid, x = "whoseverity_imp", y = "logil17eil25", 
          title = "IL-17E",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17E (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17eil25 <- covid %>% kruskal_test(logil17eil25 ~ whoseverity_imp)
res.kruskalil17eil25

pwcil17eil25 <- covid %>% dunn_test(logil17eil25 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil17eil25

pwcil17eil25 <- pwcil17eil25 %>% add_xy_position(x = "whoseverity_imp")
pwcil17eil25plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil17eil25", 
         title = "IL-17E",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17E (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil17eil25, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17eil25, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil17eil25plot

#IL-17F
ggviolin(covid, x = "whoseverity_imp", y = "logil17f", 
          title = "IL-17F",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17F (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17f <- covid %>% kruskal_test(logil17f ~ whoseverity_imp)
res.kruskalil17f

pwcil17f <- covid %>% dunn_test(logil17f ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil17f

pwcil17f <- pwcil17f %>% add_xy_position(x = "whoseverity_imp")
pwcil17fplot <- ggviolin(covid, x = "whoseverity_imp", y = "logil17f", 
         title = "IL-17F",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17F (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil17f, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17f, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil17fplot

#IL-18
ggviolin(covid, x = "whoseverity_imp", y = "logil18", 
          title = "IL-18",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-18 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil18 <- covid %>% kruskal_test(logil18 ~ whoseverity_imp)
res.kruskalil18

pwcil18 <- covid %>% dunn_test(logil18 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil18

pwcil18 <- pwcil18 %>% add_xy_position(x = "whoseverity_imp")
pwcil18plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil18", 
         title = "IL-18",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-18 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil18, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil18, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil18plot

#IL-22
ggviolin(covid, x = "whoseverity_imp", y = "logil22", 
          title = "IL-22",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-22 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil22 <- covid %>% kruskal_test(logil22 ~ whoseverity_imp)
res.kruskalil22

pwcil22 <- covid %>% dunn_test(logil22 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil22

pwcil22 <- pwcil22 %>% add_xy_position(x = "whoseverity_imp")
pwcil22plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil22", 
         title = "IL-22",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-22 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil22, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil22, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil22plot

#IL-27
ggviolin(covid, x = "whoseverity_imp", y = "logil27", 
          title = "IL-27",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-27 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil27 <- covid %>% kruskal_test(logil27 ~ whoseverity_imp)
res.kruskalil27

pwcil27 <- covid %>% dunn_test(logil27 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcil27

pwcil27 <- pwcil27 %>% add_xy_position(x = "whoseverity_imp")
pwcil27plot <- ggviolin(covid, x = "whoseverity_imp", y = "logil27", 
         title = "IL-27",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-27 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcil27, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil27, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcil27plot

#IP-10/CXCL10
ggviolin(covid, x = "whoseverity_imp", y = "logip10", 
          title = "CXCL10",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalip10 <- covid %>% kruskal_test(logip10 ~ whoseverity_imp)
res.kruskalip10

pwcip10 <- covid %>% dunn_test(logip10 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcip10

pwcip10 <- pwcip10 %>% add_xy_position(x = "whoseverity_imp")
pwcip10plot <- ggviolin(covid, x = "whoseverity_imp", y = "logip10", 
         title = "CXCL10",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcip10, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.25, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalip10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcip10plot

#MCP-1
ggviolin(covid, x = "whoseverity_imp", y = "logmcp1", 
          title = "MCP-1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcp1 <- covid %>% kruskal_test(logmcp1 ~ whoseverity_imp)
res.kruskalmcp1

pwcmcp1 <- covid %>% dunn_test(logmcp1 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmcp1

pwcmcp1 <- pwcmcp1 %>% add_xy_position(x = "whoseverity_imp")
pwcmcp1plot <- ggviolin(covid, x = "whoseverity_imp", y = "logmcp1", 
         title = "MCP-1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmcp1, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp1, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmcp1plot

#MCP-3
ggviolin(covid, x = "whoseverity_imp", y = "logmcp3", 
          title = "MCP-3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcp3 <- covid %>% kruskal_test(logmcp3 ~ whoseverity_imp)
res.kruskalmcp3

pwcmcp3 <- covid %>% dunn_test(logmcp3 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmcp3

pwcmcp3 <- pwcmcp3 %>% add_xy_position(x = "whoseverity_imp")
pwcmcp3plot <- ggviolin(covid, x = "whoseverity_imp", y = "logmcp3", 
         title = "MCP-3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmcp3, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmcp3plot

#M-CSF
ggviolin(covid, x = "whoseverity_imp", y = "logmcsf", 
          title = "M-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "M-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcsf <- covid %>% kruskal_test(logmcsf ~ whoseverity_imp)
res.kruskalmcsf

pwcmcsf <- covid %>% dunn_test(logmcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmcsf

pwcmcsf <- pwcmcsf %>% add_xy_position(x = "whoseverity_imp")
pwcmcsfplot <- ggviolin(covid, x = "whoseverity_imp", y = "logmcsf", 
         title = "M-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "M-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmcsfplot

#MDC
ggviolin(covid, x = "whoseverity_imp", y = "logmdc", 
          title = "MDC",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MDC (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmdc <- covid %>% kruskal_test(logmdc ~ whoseverity_imp)
res.kruskalmdc

pwcmdc <- covid %>% dunn_test(logmdc ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmdc

pwcmdc <- pwcmdc %>% add_xy_position(x = "whoseverity_imp")
pwcmdcplot <- ggviolin(covid, x = "whoseverity_imp", y = "logmdc", 
         title = "MDC",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MDC (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmdc, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.25, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmdc, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmdcplot

#MIG/CXCL9
ggviolin(covid, x = "whoseverity_imp", y = "logmigcxcl9", 
          title = "CXCL9",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmigcxcl9 <- covid %>% kruskal_test(logmigcxcl9 ~ whoseverity_imp)
res.kruskalmigcxcl9

pwcmigcxcl9 <- covid %>% dunn_test(logmigcxcl9 ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmigcxcl9

pwcmigcxcl9 <- pwcmigcxcl9 %>% add_xy_position(x = "whoseverity_imp")
pwcmigcxcl9plot <- ggviolin(covid, x = "whoseverity_imp", y = "logmigcxcl9", 
         title = "CXCL9",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmigcxcl9, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmigcxcl9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmigcxcl9plot

#MIP-1a/CCL3
ggviolin(covid, x = "whoseverity_imp", y = "logmip1a", 
          title = "CCL3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmip1a <- covid %>% kruskal_test(logmip1a ~ whoseverity_imp)
res.kruskalmip1a

pwcmip1a <- covid %>% dunn_test(logmip1a ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmip1a

pwcmip1a <- pwcmip1a %>% add_xy_position(x = "whoseverity_imp")
pwcmip1aplot <- ggviolin(covid, x = "whoseverity_imp", y = "logmip1a", 
         title = "CCL3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmip1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmip1aplot

#MIP-1b/CCL4
ggviolin(covid, x = "whoseverity_imp", y = "logmip1b", 
          title = "CCL4",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmip1b <- covid %>% kruskal_test(logmip1b ~ whoseverity_imp)
res.kruskalmip1b

pwcmip1b <- covid %>% dunn_test(logmip1b ~ whoseverity_imp, p.adjust.method = "BH") 
pwcmip1b

pwcmip1b <- pwcmip1b %>% add_xy_position(x = "whoseverity_imp")
pwcmip1bplot <- ggviolin(covid, x = "whoseverity_imp", y = "logmip1b", 
         title = "CCL4",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcmip1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcmip1bplot

#PDGF-AA
ggviolin(covid, x = "whoseverity_imp", y = "logpdgfaa", 
          title = "PDGF-AA",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AA (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalpdgfaa <- covid %>% kruskal_test(logpdgfaa ~ whoseverity_imp)
res.kruskalpdgfaa

pwcpdgfaa <- covid %>% dunn_test(logpdgfaa ~ whoseverity_imp, p.adjust.method = "BH") 
pwcpdgfaa

pwcpdgfaa <- pwcpdgfaa %>% add_xy_position(x = "whoseverity_imp")
pwcpdgfaaplot <- ggviolin(covid, x = "whoseverity_imp", y = "logpdgfaa", 
         title = "PDGF-AA",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AA (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcpdgfaa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfaa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcpdgfaaplot

#PDGF-AB/BB
ggviolin(covid, x = "whoseverity_imp", y = "logpdgfabbb", 
          title = "PDGF-AB/BB",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalpdgfabbb <- covid %>% kruskal_test(logpdgfabbb ~ whoseverity_imp)
res.kruskalpdgfabbb

pwcpdgfabbb <- covid %>% dunn_test(logpdgfabbb ~ whoseverity_imp, p.adjust.method = "BH") 
pwcpdgfabbb

pwcpdgfabbb <- pwcpdgfabbb %>% add_xy_position(x = "whoseverity_imp")
pwcpdgfabbbplot <- ggviolin(covid, x = "whoseverity_imp", y = "logpdgfabbb", 
         title = "PDGF-AB/BB",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcpdgfabbb, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfabbb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcpdgfabbbplot

#RANTES
ggviolin(covid, x = "whoseverity_imp", y = "lograntes", 
          title = "RANTES",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "RANTES (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalrantes <- covid %>% kruskal_test(lograntes ~ whoseverity_imp)
res.kruskalrantes

pwcrantes <- covid %>% dunn_test(lograntes ~ whoseverity_imp, p.adjust.method = "BH") 
pwcrantes

pwcrantes <- pwcrantes %>% add_xy_position(x = "whoseverity_imp")
pwcrantesplot <- ggviolin(covid, x = "whoseverity_imp", y = "lograntes", 
         title = "RANTES",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "RANTES (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcrantes, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalrantes, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcrantesplot

#TGF-alpha
ggviolin(covid, x = "whoseverity_imp", y = "logtgfa", 
          title = "TGF-\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltgfa <- covid %>% kruskal_test(logtgfa ~ whoseverity_imp)
res.kruskaltgfa

pwctgfa <- covid %>% dunn_test(logtgfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwctgfa

pwctgfa <- pwctgfa %>% add_xy_position(x = "whoseverity_imp")
pwctgfaplot <- ggviolin(covid, x = "whoseverity_imp", y = "logtgfa", 
         title = "TGF-\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwctgfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltgfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwctgfaplot

#TNF-alpha
ggviolin(covid, x = "whoseverity_imp", y = "logtnfa", 
          title = "TNF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TNF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltnfa <- covid %>% kruskal_test(logtnfa ~ whoseverity_imp)
res.kruskaltnfa

pwctnfa <- covid %>% dunn_test(logtnfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwctnfa

pwctnfa <- pwctnfa %>% add_xy_position(x = "whoseverity_imp")
pwctnfaplot <- ggviolin(covid, x = "whoseverity_imp", y = "logtnfa", 
         title = "TNF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TNF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwctnfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwctnfaplot

#Lymphotoxin-alpha/TNF-beta
ggviolin(covid, x = "whoseverity_imp", y = "logtnfb", 
          title = "Lymphotoxin-\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltnfb <- covid %>% kruskal_test(logtnfb ~ whoseverity_imp)
res.kruskaltnfb

pwctnfb <- covid %>% dunn_test(logtnfb ~ whoseverity_imp, p.adjust.method = "BH") 
pwctnfb

pwctnfb <- pwctnfb %>% add_xy_position(x = "whoseverity_imp")
pwctnfbplot <- ggviolin(covid, x = "whoseverity_imp", y = "logtnfb", 
         title = "Lymphotoxin-\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwctnfb, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwctnfbplot

#VEGF-A
ggviolin(covid, x = "whoseverity_imp", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalvegfa <- covid %>% kruskal_test(logvegfa ~ whoseverity_imp)
res.kruskalvegfa

pwcvegfa <- covid %>% dunn_test(logvegfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwcvegfa

pwcvegfa <- pwcvegfa %>% add_xy_position(x = "whoseverity_imp")
pwcvegfaplot <- ggviolin(covid, x = "whoseverity_imp", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwcvegfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalvegfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwcvegfaplot

#Table showing these values
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
                                            strata = "whoseverity_imp",
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

write_tableHTML(tableHTML(tablebiomarkers), file = 'biomarkersseverity.html')




