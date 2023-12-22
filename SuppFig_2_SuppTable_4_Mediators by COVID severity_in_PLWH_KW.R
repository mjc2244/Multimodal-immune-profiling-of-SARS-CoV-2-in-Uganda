#Figure S2 and Table S4 - Mediators by COVID severity in PLWH

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients
covid <- subset(combined, pathogencode == 1)

#Biomarker concentrations across illness severity groups in PLWH, including asymptomatic controls
library(ggpubr)
library(rstatix)

#sCD40L
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logscd40l", 
          title = "sCD40L",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "sCD40L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalscd40l <- subset(covid, hivrdtresult==1) %>% kruskal_test(logscd40l ~ whoseverity_imp)
res.kruskalscd40l

pwchivscd40l <- subset(covid, hivrdtresult==1) %>% dunn_test(logscd40l ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivscd40l

pwchivscd40l <- pwchivscd40l %>% add_xy_position(x = "whoseverity_imp")
pwchivscd40lplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logscd40l", 
         title = "sCD40L",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "sCD40L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivscd40l, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalscd40l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#EGF
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logegf", 
          title = "EGF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "EGF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalegf <- subset(covid, hivrdtresult==1) %>% kruskal_test(logegf ~ whoseverity_imp)
res.kruskalegf

pwchivegf <- subset(covid, hivrdtresult==1) %>% dunn_test(logegf ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivegf

pwchivegf <- pwchivegf %>% add_xy_position(x = "whoseverity_imp")
pwchivegfplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logegf", 
         title = "EGF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "EGF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivegf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalegf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#Eotaxin
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logeotaxin", 
          title = "Eotaxin",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "Eotaxin (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskaleotaxin <- subset(covid, hivrdtresult==1) %>% kruskal_test(logeotaxin ~ whoseverity_imp)
res.kruskaleotaxin

pwchiveotaxin <- subset(covid, hivrdtresult==1) %>% dunn_test(logeotaxin ~ whoseverity_imp, p.adjust.method = "BH") 
pwchiveotaxin

pwchiveotaxin <- pwchiveotaxin %>% add_xy_position(x = "whoseverity_imp")
pwchiveotaxinplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logeotaxin", 
         title = "Eotaxin",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Eotaxin (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchiveotaxin, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaleotaxin, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#FGF2
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logfgf2", 
          title = "FGF-2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FGF-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalfgf2 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logfgf2 ~ whoseverity_imp)
res.kruskalfgf2

pwchivfgf2 <- subset(covid, hivrdtresult==1) %>% dunn_test(logfgf2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivfgf2

pwchivfgf2 <- pwchivfgf2 %>% add_xy_position(x = "whoseverity_imp")
pwchivfgf2plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logfgf2", 
         title = "FGF-2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FGF-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivfgf2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfgf2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#FLT-3L
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logflt3l", 
          title = "FLT-3L",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FLT-3L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalflt3l <- subset(covid, hivrdtresult==1) %>% kruskal_test(logflt3l ~ whoseverity_imp)
res.kruskalflt3l

pwchivflt3l <- subset(covid, hivrdtresult==1) %>% dunn_test(logflt3l ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivflt3l

pwchivflt3l <- pwchivflt3l %>% add_xy_position(x = "whoseverity_imp")
pwchivflt3lplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logflt3l", 
         title = "FLT-3L",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FLT-3L (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivflt3l, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalflt3l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#FKN
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logfractalkine", 
          title = "Fractalkine",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Fractalkine (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalfractalkine <- subset(covid, hivrdtresult==1) %>% kruskal_test(logfractalkine ~ whoseverity_imp)
res.kruskalfractalkine

pwchivfractalkine <- subset(covid, hivrdtresult==1) %>% dunn_test(logfractalkine ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivfractalkine

pwchivfractalkine <- pwchivfractalkine %>% add_xy_position(x = "whoseverity_imp")
pwchivfractalkineplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logfractalkine", 
         title = "Fractalkine",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Fractalkine (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivfractalkine, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfractalkine, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#G-CSF
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggcsf", 
          title = "G-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "G-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgcsf <- subset(covid, hivrdtresult==1) %>% kruskal_test(loggcsf ~ whoseverity_imp)
res.kruskalgcsf

pwchivgcsf <- subset(covid, hivrdtresult==1) %>% dunn_test(loggcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivgcsf

pwchivgcsf <- pwchivgcsf %>% add_xy_position(x = "whoseverity_imp")
pwchivgcsfplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggcsf", 
         title = "G-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "G-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivgcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#GM-CSF
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggmcsf", 
          title = "GM-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "GM-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgmcsf <- subset(covid, hivrdtresult==1) %>% kruskal_test(loggmcsf ~ whoseverity_imp)
res.kruskalgmcsf

pwchivgmcsf <- subset(covid, hivrdtresult==1) %>% dunn_test(loggmcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivgmcsf

pwchivgmcsf <- pwchivgmcsf %>% add_xy_position(x = "whoseverity_imp")
pwchivgmcsfplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggmcsf", 
         title = "GM-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "GM-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivgmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#GRO-alpha/CXCL1
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggroalpha", 
          title = "CXCL1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalgroalpha <- subset(covid, hivrdtresult==1) %>% kruskal_test(loggroalpha ~ whoseverity_imp)
res.kruskalgroalpha

pwchivgroalpha <- subset(covid, hivrdtresult==1) %>% dunn_test(loggroalpha ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivgroalpha

pwchivgroalpha <- pwchivgroalpha %>% add_xy_position(x = "whoseverity_imp")
pwchivgroalphaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "loggroalpha", 
         title = "CXCL1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivgroalpha, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgroalpha, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IFN-a2
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logifna2", 
          title = "IFN-\u03b12",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalifna2 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logifna2 ~ whoseverity_imp)
res.kruskalifna2

pwchivifna2 <- subset(covid, hivrdtresult==1) %>% dunn_test(logifna2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivifna2

pwchivifna2 <- pwchivifna2 %>% add_xy_position(x = "whoseverity_imp")
pwchivifna2plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logifna2", 
         title = "IFN-\u03b12",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivifna2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifna2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IFN-gamma
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logifngamma", 
          title = "IFN-\u03b3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalifngamma <- subset(covid, hivrdtresult==1) %>% kruskal_test(logifngamma ~ whoseverity_imp)
res.kruskalifngamma

pwchivifngamma <- subset(covid, hivrdtresult==1) %>% dunn_test(logifngamma ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivifngamma

pwchivifngamma <- pwchivifngamma %>% add_xy_position(x = "whoseverity_imp")
pwchivifngammaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logifngamma", 
         title = "IFN-\u03b3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivifngamma, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifngamma, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-1a
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1a", 
          title = "IL-1\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1a <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil1a ~ whoseverity_imp)
res.kruskalil1a

pwchivil1a <- subset(covid, hivrdtresult==1) %>% dunn_test(logil1a ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil1a

pwchivil1a <- pwchivil1a %>% add_xy_position(x = "whoseverity_imp")
pwchivil1aplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1a", 
         title = "IL-1\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-1b
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1b", 
          title = "IL-1\u03B2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1b <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil1b ~ whoseverity_imp)
res.kruskalil1b

pwchivil1b <- subset(covid, hivrdtresult==1) %>% dunn_test(logil1b ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil1b

pwchivil1b <- pwchivil1b %>% add_xy_position(x = "whoseverity_imp")
pwchivil1bplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1b", 
         title = "IL-1\u03B2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-1Ra
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1ra", 
          title = "IL-1Ra",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1Ra (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil1ra <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil1ra ~ whoseverity_imp)
res.kruskalil1ra

pwchivil1ra <- subset(covid, hivrdtresult==1) %>% dunn_test(logil1ra ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil1ra

pwchivil1ra <- pwchivil1ra %>% add_xy_position(x = "whoseverity_imp")
pwchivil1raplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil1ra", 
         title = "IL-1Ra",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1Ra (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil1ra, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1ra, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-2
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil2", 
          title = "IL-2",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil2 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil2 ~ whoseverity_imp)
res.kruskalil2

pwchivil2 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil2 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil2

pwchivil2 <- pwchivil2 %>% add_xy_position(x = "whoseverity_imp")
pwchivil2plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil2", 
         title = "IL-2",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-2 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil2, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-3
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil3", 
          title = "IL-3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil3 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil3 ~ whoseverity_imp)
res.kruskalil3

pwchivil3 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil3 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil3

pwchivil3 <- pwchivil3 %>% add_xy_position(x = "whoseverity_imp")
pwchivil3plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil3", 
         title = "IL-3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil3, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-4
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil4", 
          title = "IL-4",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil4 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil4 ~ whoseverity_imp)
res.kruskalil4

pwchivil4 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil4 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil4

pwchivil4 <- pwchivil4 %>% add_xy_position(x = "whoseverity_imp")
pwchivil4plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil4", 
         title = "IL-4",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil4, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil4, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-5
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil5", 
          title = "IL-5",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-5 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil5 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil5 ~ whoseverity_imp)
res.kruskalil5

pwchivil5 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil5 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil5

pwchivil5 <- pwchivil5 %>% add_xy_position(x = "whoseverity_imp")
pwchivil5plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil5", 
         title = "IL-5",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-5 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil5, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil5, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-6
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil6", 
          title = "IL-6",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-6 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil6 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil6 ~ whoseverity_imp)
res.kruskalil6

pwchivil6 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil6 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil6

pwchivil6 <- pwchivil6 %>% add_xy_position(x = "whoseverity_imp")
pwchivil6plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil6", 
         title = "IL-6",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-6 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil6, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.8, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil6, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-7
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil7", 
          title = "IL-7",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-7 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil7 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil7 ~ whoseverity_imp)
res.kruskalil7

pwchivil7 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil7 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil7

pwchivil7 <- pwchivil7 %>% add_xy_position(x = "whoseverity_imp")
pwchivil7plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil7", 
         title = "IL-7",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-7 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil7, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil7, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-8
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil8", 
          title = "IL-8",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-8 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil8 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil8 ~ whoseverity_imp)
res.kruskalil8

pwchivil8 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil8 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil8

pwchivil8 <- pwchivil8 %>% add_xy_position(x = "whoseverity_imp")
pwchivil8plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil8", 
         title = "IL-8",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-8 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil8, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil8, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-9
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil9", 
          title = "IL-9",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil9 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil9 ~ whoseverity_imp)
res.kruskalil9

pwchivil9 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil9 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil9

pwchivil9 <- pwchivil9 %>% add_xy_position(x = "whoseverity_imp")
pwchivil9plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil9", 
         title = "IL-9",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil9, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-10
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil10", 
          title = "IL-10",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil10 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil10 ~ whoseverity_imp)
res.kruskalil10

pwchivil10 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil10 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil10

pwchivil10 <- pwchivil10 %>% add_xy_position(x = "whoseverity_imp")
pwchivil10plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil10", 
         title = "IL-10",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil10, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-12p40
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil12p40", 
          title = "IL-12p40",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p40 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil12p40 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil12p40 ~ whoseverity_imp)
res.kruskalil12p40

pwchivil12p40 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil12p40 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil12p40

pwchivil12p40 <- pwchivil12p40 %>% add_xy_position(x = "whoseverity_imp")
pwchivil12p40plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil12p40", 
         title = "IL-12p40",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p40 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil12p40, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p40, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-12p70
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil12p70", 
          title = "IL-12p70",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p70 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means() 

res.kruskalil12p70 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil12p70 ~ whoseverity_imp)
res.kruskalil12p70

pwchivil12p70 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil12p70 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil12p70

pwchivil12p70 <- pwchivil12p70 %>% add_xy_position(x = "whoseverity_imp")
pwchivil12p70plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil12p70", 
         title = "IL-12p70",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p70 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil12p70, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p70, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-13
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil13", 
          title = "IL-13",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-13 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil13 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil13 ~ whoseverity_imp)
res.kruskalil13

pwchivil13 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil13 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil13

pwchivil13 <- pwchivil13 %>% add_xy_position(x = "whoseverity_imp")
pwchivil13plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil13", 
         title = "IL-13",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-13 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil13, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil13, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-15
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil15", 
          title = "IL-15",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-15 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil15 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil15 ~ whoseverity_imp)
res.kruskalil15

pwchivil15 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil15 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil15

pwchivil15 <- pwchivil15 %>% add_xy_position(x = "whoseverity_imp")
pwchivil15plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil15", 
         title = "IL-15",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-15 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil15, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil15, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-17A
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17a", 
          title = "IL-17A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17a <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil17a ~ whoseverity_imp)
res.kruskalil17a

pwchivil17a <- subset(covid, hivrdtresult==1) %>% dunn_test(logil17a ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil17a

pwchivil17a <- pwchivil17a %>% add_xy_position(x = "whoseverity_imp")
pwchivil17aplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17a", 
         title = "IL-17A",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil17a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-17E
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17eil25", 
          title = "IL-17E",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17E (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17eil25 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil17eil25 ~ whoseverity_imp)
res.kruskalil17eil25

pwchivil17eil25 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil17eil25 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil17eil25

pwchivil17eil25 <- pwchivil17eil25 %>% add_xy_position(x = "whoseverity_imp")
pwchivil17eil25plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17eil25", 
         title = "IL-17E",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17E (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil17eil25, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17eil25, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-17F
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17f", 
          title = "IL-17F",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17F (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil17f <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil17f ~ whoseverity_imp)
res.kruskalil17f

pwchivil17f <- subset(covid, hivrdtresult==1) %>% dunn_test(logil17f ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil17f

pwchivil17f <- pwchivil17f %>% add_xy_position(x = "whoseverity_imp")
pwchivil17fplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil17f", 
         title = "IL-17F",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17F (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil17f, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17f, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-18
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil18", 
          title = "IL-18",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-18 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil18 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil18 ~ whoseverity_imp)
res.kruskalil18

pwchivil18 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil18 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil18

pwchivil18 <- pwchivil18 %>% add_xy_position(x = "whoseverity_imp")
pwchivil18plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil18", 
         title = "IL-18",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-18 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil18, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil18, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-22
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil22", 
          title = "IL-22",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-22 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil22 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil22 ~ whoseverity_imp)
res.kruskalil22

pwchivil22 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil22 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil22

pwchivil22 <- pwchivil22 %>% add_xy_position(x = "whoseverity_imp")
pwchivil22plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil22", 
         title = "IL-22",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-22 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil22, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil22, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IL-27
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil27", 
          title = "IL-27",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-27 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalil27 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logil27 ~ whoseverity_imp)
res.kruskalil27

pwchivil27 <- subset(covid, hivrdtresult==1) %>% dunn_test(logil27 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivil27

pwchivil27 <- pwchivil27 %>% add_xy_position(x = "whoseverity_imp")
pwchivil27plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logil27", 
         title = "IL-27",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-27 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivil27, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil27, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#IP-10/CXCL10
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logip10", 
          title = "CXCL10",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalip10 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logip10 ~ whoseverity_imp)
res.kruskalip10

pwchivip10 <- subset(covid, hivrdtresult==1) %>% dunn_test(logip10 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivip10

pwchivip10 <- pwchivip10 %>% add_xy_position(x = "whoseverity_imp")
pwchivip10plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logip10", 
         title = "CXCL10",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL10 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivip10, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalip10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MCP-1
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcp1", 
          title = "MCP-1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcp1 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmcp1 ~ whoseverity_imp)
res.kruskalmcp1

pwchivmcp1 <- subset(covid, hivrdtresult==1) %>% dunn_test(logmcp1 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmcp1

pwchivmcp1 <- pwchivmcp1 %>% add_xy_position(x = "whoseverity_imp")
pwchivmcp1plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcp1", 
         title = "MCP-1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmcp1, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp1, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MCP-3
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcp3", 
          title = "MCP-3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcp3 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmcp3 ~ whoseverity_imp)
res.kruskalmcp3

pwchivmcp3 <- subset(covid, hivrdtresult==1) %>% dunn_test(logmcp3 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmcp3

pwchivmcp3 <- pwchivmcp3 %>% add_xy_position(x = "whoseverity_imp")
pwchivmcp3plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcp3", 
         title = "MCP-3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmcp3, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#M-CSF
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcsf", 
          title = "M-CSF",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "M-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmcsf <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmcsf ~ whoseverity_imp)
res.kruskalmcsf

pwchivmcsf <- subset(covid, hivrdtresult==1) %>% dunn_test(logmcsf ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmcsf

pwchivmcsf <- pwchivmcsf %>% add_xy_position(x = "whoseverity_imp")
pwchivmcsfplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmcsf", 
         title = "M-CSF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "M-CSF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MDC
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmdc", 
          title = "MDC",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MDC (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmdc <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmdc ~ whoseverity_imp)
res.kruskalmdc

pwchivmdc <- subset(covid, hivrdtresult==1) %>% dunn_test(logmdc ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmdc

pwchivmdc <- pwchivmdc %>% add_xy_position(x = "whoseverity_imp")
pwchivmdcplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmdc", 
         title = "MDC",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MDC (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmdc, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmdc, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MIG/CXCL9
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmigcxcl9", 
          title = "CXCL9",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmigcxcl9 <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmigcxcl9 ~ whoseverity_imp)
res.kruskalmigcxcl9

pwchivmigcxcl9 <- subset(covid, hivrdtresult==1) %>% dunn_test(logmigcxcl9 ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmigcxcl9

pwchivmigcxcl9 <- pwchivmigcxcl9 %>% add_xy_position(x = "whoseverity_imp")
pwchivmigcxcl9plot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmigcxcl9", 
         title = "CXCL9",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL9 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmigcxcl9, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmigcxcl9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MIP-1a/CCL3
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmip1a", 
          title = "CCL3",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmip1a <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmip1a ~ whoseverity_imp)
res.kruskalmip1a

pwchivmip1a <- subset(covid, hivrdtresult==1) %>% dunn_test(logmip1a ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmip1a

pwchivmip1a <- pwchivmip1a %>% add_xy_position(x = "whoseverity_imp")
pwchivmip1aplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmip1a", 
         title = "CCL3",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL3 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmip1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#MIP-1b/CCL4
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmip1b", 
          title = "CCL4",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalmip1b <- subset(covid, hivrdtresult==1) %>% kruskal_test(logmip1b ~ whoseverity_imp)
res.kruskalmip1b

pwchivmip1b <- subset(covid, hivrdtresult==1) %>% dunn_test(logmip1b ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivmip1b

pwchivmip1b <- pwchivmip1b %>% add_xy_position(x = "whoseverity_imp")
pwchivmip1bplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logmip1b", 
         title = "CCL4",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL4 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivmip1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#PDGF-AA
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logpdgfaa", 
          title = "PDGF-AA",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AA (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalpdgfaa <- subset(covid, hivrdtresult==1) %>% kruskal_test(logpdgfaa ~ whoseverity_imp)
res.kruskalpdgfaa

pwchivpdgfaa <- subset(covid, hivrdtresult==1) %>% dunn_test(logpdgfaa ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivpdgfaa

pwchivpdgfaa <- pwchivpdgfaa %>% add_xy_position(x = "whoseverity_imp")
pwchivpdgfaaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logpdgfaa", 
         title = "PDGF-AA",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AA (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivpdgfaa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfaa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#PDGF-AB/BB
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logpdgfabbb", 
          title = "PDGF-AB/BB",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalpdgfabbb <- subset(covid, hivrdtresult==1) %>% kruskal_test(logpdgfabbb ~ whoseverity_imp)
res.kruskalpdgfabbb

pwchivpdgfabbb <- subset(covid, hivrdtresult==1) %>% dunn_test(logpdgfabbb ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivpdgfabbb

pwchivpdgfabbb <- pwchivpdgfabbb %>% add_xy_position(x = "whoseverity_imp")
pwchivpdgfabbbplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logpdgfabbb", 
         title = "PDGF-AB/BB",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivpdgfabbb, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfabbb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#RANTES
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "lograntes", 
          title = "RANTES",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "RANTES (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalrantes <- subset(covid, hivrdtresult==1) %>% kruskal_test(lograntes ~ whoseverity_imp)
res.kruskalrantes

pwchivrantes <- subset(covid, hivrdtresult==1) %>% dunn_test(lograntes ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivrantes

pwchivrantes <- pwchivrantes %>% add_xy_position(x = "whoseverity_imp")
pwchivrantesplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "lograntes", 
         title = "RANTES",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "RANTES (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivrantes, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalrantes, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#TGF-alpha
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtgfa", 
          title = "TGF-\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltgfa <- subset(covid, hivrdtresult==1) %>% kruskal_test(logtgfa ~ whoseverity_imp)
res.kruskaltgfa

pwchivtgfa <- subset(covid, hivrdtresult==1) %>% dunn_test(logtgfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivtgfa

pwchivtgfa <- pwchivtgfa %>% add_xy_position(x = "whoseverity_imp")
pwchivtgfaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtgfa", 
         title = "TGF-\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivtgfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltgfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#TNF-alpha
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtnfa", 
          title = "TNF-\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TNF-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltnfa <- subset(covid, hivrdtresult==1) %>% kruskal_test(logtnfa ~ whoseverity_imp)
res.kruskaltnfa

pwchivtnfa <- subset(covid, hivrdtresult==1) %>% dunn_test(logtnfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivtnfa

pwchivtnfa <- pwchivtnfa %>% add_xy_position(x = "whoseverity_imp")
pwchivtnfaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtnfa", 
         title = "TNF",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TNF (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivtnfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#Lymphotoxin-alpha/TNF-beta
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtnfb", 
          title = "Lymphotoxin-\u03b1",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskaltnfb <- subset(covid, hivrdtresult==1) %>% kruskal_test(logtnfb ~ whoseverity_imp)
res.kruskaltnfb

pwchivtnfb <- subset(covid, hivrdtresult==1) %>% dunn_test(logtnfb ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivtnfb

pwchivtnfb <- pwchivtnfb %>% add_xy_position(x = "whoseverity_imp")
pwchivtnfbplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logtnfb", 
         title = "Lymphotoxin-\u03b1",
         color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
         order = c("1", "2", "3", "4"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivtnfb, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))

#VEGF-A
ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "WHO Clinical Severity") +
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_compare_means()

res.kruskalvegfa <- subset(covid, hivrdtresult==1) %>% kruskal_test(logvegfa ~ whoseverity_imp)
res.kruskalvegfa

pwchivvegfa <- subset(covid, hivrdtresult==1) %>% dunn_test(logvegfa ~ whoseverity_imp, p.adjust.method = "BH") 
pwchivvegfa

pwchivvegfa <- pwchivvegfa %>% add_xy_position(x = "whoseverity_imp")
pwchivvegfaplot <- ggviolin(subset(covid, hivrdtresult==1), x = "whoseverity_imp", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "whoseverity_imp", alpha = 0.5, palette =  c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),
          order = c("1", "2", "3", "4"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "WHO Clinical Severity") + 
  scale_x_discrete(labels=c("1" = "Asymptomatic", "2" = "Mild", "3" = "Moderate", "4" = "Severe")) + 
  stat_pvalue_manual(pwchivvegfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 4, bracket.nudge.y = 0.75, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalvegfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 12))
pwchivvegfaplot

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

tablebiomarkershiv <- CreateTableOne(vars = listbiomarkers, 
                                            data = subset(covid, hivrdtresult==1), 
                                            strata = "whoseverity_imp",
                                            addOverall = TRUE)
tablebiomarkershiv
summary(tablebiomarkershiv)

tablebiomarkershiv <- print(tablebiomarkershiv, nonnormal = c("logscd40l", 
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

write_tableHTML(tableHTML(tablebiomarkershiv), file = 'biomarkersseverityhiv.html')



