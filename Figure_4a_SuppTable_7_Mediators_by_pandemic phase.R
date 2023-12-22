#Figure 4a and Table S9 - Mediator concentrations stratified by SARS-CoV-2 variant-driven pandemic phase 

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select COVID patients
covid <- subset(combined, pathogencode == 1)

#Biomarker concentrations across pandemic phases, including asymptomatic controls
library(ggpubr)
library(rstatix)

#sCD40L
ggviolin(covid, x = "phase3", y = "logscd40l", 
          title = "sCD40L",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "sCD40L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalscd40l <- covid %>% kruskal_test(logscd40l ~ phase3)
res.kruskalscd40l

pwcscd40l <- covid %>% dunn_test(logscd40l ~ phase3, p.adjust.method = "BH") 
pwcscd40l

pwcscd40l <- pwcscd40l %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logscd40l", 
         title = "sCD40L",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "sCD40L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcscd40l, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalscd40l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#EGF
ggviolin(covid, x = "phase3", y = "logegf", 
          title = "EGF",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "EGF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalegf <- covid %>% kruskal_test(logegf ~ phase3)
res.kruskalegf

pwcegf <- covid %>% dunn_test(logegf ~ phase3, p.adjust.method = "BH") 
pwcegf

pwcegf <- pwcegf %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logegf", 
         title = "EGF",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "EGF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcegf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalegf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#Eotaxin
ggviolin(covid, x = "phase3", y = "logeotaxin", 
          title = "Eotaxin",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "Eotaxin (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskaleotaxin <- covid %>% kruskal_test(logeotaxin ~ phase3)
res.kruskaleotaxin

pwceotaxin <- covid %>% dunn_test(logeotaxin ~ phase3, p.adjust.method = "BH") 
pwceotaxin

pwceotaxin <- pwceotaxin %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logeotaxin", 
         title = "Eotaxin",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Eotaxin (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwceotaxin, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaleotaxin, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#FGF2
ggviolin(covid, x = "phase3", y = "logfgf2", 
          title = "FGF-2",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FGF-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalfgf2 <- covid %>% kruskal_test(logfgf2 ~ phase3)
res.kruskalfgf2

pwcfgf2 <- covid %>% dunn_test(logfgf2 ~ phase3, p.adjust.method = "BH") 
pwcfgf2

pwcfgf2 <- pwcfgf2 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logfgf2", 
         title = "FGF-2",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FGF-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcfgf2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfgf2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#FLT-3L
ggviolin(covid, x = "phase3", y = "logflt3l", 
          title = "FLT-3L",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FLT-3L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalflt3l <- covid %>% kruskal_test(logflt3l ~ phase3)
res.kruskalflt3l

pwcflt3l <- covid %>% dunn_test(logflt3l ~ phase3, p.adjust.method = "BH") 
pwcflt3l

pwcflt3l <- pwcflt3l %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logflt3l", 
         title = "FLT-3L",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FLT-3L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcflt3l, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalflt3l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#FKN
ggviolin(covid, x = "phase3", y = "logfractalkine", 
          title = "Fractalkine",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Fractalkine (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalfractalkine <- covid %>% kruskal_test(logfractalkine ~ phase3)
res.kruskalfractalkine

pwcfractalkine <- covid %>% dunn_test(logfractalkine ~ phase3, p.adjust.method = "BH") 
pwcfractalkine

pwcfractalkine <- pwcfractalkine %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logfractalkine", 
         title = "Fractalkine",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Fractalkine (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcfractalkine, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfractalkine, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#G-CSF
ggviolin(covid, x = "phase3", y = "loggcsf", 
          title = "G-CSF",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "G-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalgcsf <- covid %>% kruskal_test(loggcsf ~ phase3)
res.kruskalgcsf

pwcgcsf <- covid %>% dunn_test(loggcsf ~ phase3, p.adjust.method = "BH") 
pwcgcsf

pwcgcsf <- pwcgcsf %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "loggcsf", 
         title = "G-CSF",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "G-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcgcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#GM-CSF
ggviolin(covid, x = "phase3", y = "loggmcsf", 
          title = "GM-CSF",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "GM-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalgmcsf <- covid %>% kruskal_test(loggmcsf ~ phase3)
res.kruskalgmcsf

pwcgmcsf <- covid %>% dunn_test(loggmcsf ~ phase3, p.adjust.method = "BH") 
pwcgmcsf

pwcgmcsf <- pwcgmcsf %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "loggmcsf", 
         title = "GM-CSF",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "GM-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcgmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#GRO-alpha/CXCL1
ggviolin(covid, x = "phase3", y = "loggroalpha", 
          title = "CXCL1",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalgroalpha <- covid %>% kruskal_test(loggroalpha ~ phase3)
res.kruskalgroalpha

pwcgroalpha <- covid %>% dunn_test(loggroalpha ~ phase3, p.adjust.method = "BH") 
pwcgroalpha

pwcgroalpha <- pwcgroalpha %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "loggroalpha", 
         title = "CXCL1",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcgroalpha, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgroalpha, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IFN-a2
ggviolin(covid, x = "phase3", y = "logifna2", 
          title = "IFN-\u03b12",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalifna2 <- covid %>% kruskal_test(logifna2 ~ phase3)
res.kruskalifna2

pwcifna2 <- covid %>% dunn_test(logifna2 ~ phase3, p.adjust.method = "BH") 
pwcifna2

pwcifna2 <- pwcifna2 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logifna2", 
         title = "IFN-\u03b12",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcifna2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifna2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IFN-gamma
ggviolin(covid, x = "phase3", y = "logifngamma", 
          title = "IFN-\u03b3",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalifngamma <- covid %>% kruskal_test(logifngamma ~ phase3)
res.kruskalifngamma

pwcifngamma <- covid %>% dunn_test(logifngamma ~ phase3, p.adjust.method = "BH") 
pwcifngamma

pwcifngamma <- pwcifngamma %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logifngamma", 
         title = "IFN-\u03b3",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcifngamma, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifngamma, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-1a
ggviolin(covid, x = "phase3", y = "logil1a", 
          title = "IL-1\u03b1",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil1a <- covid %>% kruskal_test(logil1a ~ phase3)
res.kruskalil1a

pwcil1a <- covid %>% dunn_test(logil1a ~ phase3, p.adjust.method = "BH") 
pwcil1a

pwcil1a <- pwcil1a %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil1a", 
         title = "IL-1\u03b1",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-1b
ggviolin(covid, x = "phase3", y = "logil1b", 
          title = "IL-1\u03B2",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil1b <- covid %>% kruskal_test(logil1b ~ phase3)
res.kruskalil1b

pwcil1b <- covid %>% dunn_test(logil1b ~ phase3, p.adjust.method = "BH") 
pwcil1b

pwcil1b <- pwcil1b %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil1b", 
         title = "IL-1\u03B2",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-1Ra
ggviolin(covid, x = "phase3", y = "logil1ra", 
          title = "IL-1Ra",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1Ra (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil1ra <- covid %>% kruskal_test(logil1ra ~ phase3)
res.kruskalil1ra

pwcil1ra <- covid %>% dunn_test(logil1ra ~ phase3, p.adjust.method = "BH") 
pwcil1ra

pwcil1ra <- pwcil1ra %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil1ra", 
         title = "IL-1Ra",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1Ra (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil1ra, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1ra, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-2
ggviolin(covid, x = "phase3", y = "logil2", 
          title = "IL-2",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil2 <- covid %>% kruskal_test(logil2 ~ phase3)
res.kruskalil2

pwcil2 <- covid %>% dunn_test(logil2 ~ phase3, p.adjust.method = "BH") 
pwcil2

pwcil2 <- pwcil2 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil2", 
         title = "IL-2",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-3
ggviolin(covid, x = "phase3", y = "logil3", 
          title = "IL-3",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil3 <- covid %>% kruskal_test(logil3 ~ phase3)
res.kruskalil3

pwcil3 <- covid %>% dunn_test(logil3 ~ phase3, p.adjust.method = "BH") 
pwcil3

pwcil3 <- pwcil3 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil3", 
         title = "IL-3",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil3, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-4
ggviolin(covid, x = "phase3", y = "logil4", 
          title = "IL-4",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil4 <- covid %>% kruskal_test(logil4 ~ phase3)
res.kruskalil4

pwcil4 <- covid %>% dunn_test(logil4 ~ phase3, p.adjust.method = "BH") 
pwcil4

pwcil4 <- pwcil4 %>% add_xy_position(x = "phase3")
pwcil4plot <- ggviolin(covid, x = "phase3", y = "logil4", 
         title = "IL-4",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil4, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil4, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-5
ggviolin(covid, x = "phase3", y = "logil5", 
          title = "IL-5",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-5 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil5 <- covid %>% kruskal_test(logil5 ~ phase3)
res.kruskalil5

pwcil5 <- covid %>% dunn_test(logil5 ~ phase3, p.adjust.method = "BH") 
pwcil5

pwcil5 <- pwcil5 %>% add_xy_position(x = "phase3")
pwcil5plot <- ggviolin(covid, x = "phase3", y = "logil5", 
         title = "IL-5",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-5 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil5, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0) + 
  labs(
    subtitle = get_test_label(res.kruskalil5, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-6
ggviolin(covid, x = "phase3", y = "logil6", 
          title = "IL-6",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-6 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil6 <- covid %>% kruskal_test(logil6 ~ phase3)
res.kruskalil6

pwcil6 <- covid %>% dunn_test(logil6 ~ phase3, p.adjust.method = "BH") 
pwcil6

pwcil6 <- pwcil6 %>% add_xy_position(x = "phase3")
pwcil6plot <- ggviolin(covid, x = "phase3", y = "logil6", 
         title = "IL-6",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-6 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil6, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.8, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil6, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-7
ggviolin(covid, x = "phase3", y = "logil7", 
          title = "IL-7",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-7 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil7 <- covid %>% kruskal_test(logil7 ~ phase3)
res.kruskalil7

pwcil7 <- covid %>% dunn_test(logil7 ~ phase3, p.adjust.method = "BH") 
pwcil7

pwcil7 <- pwcil7 %>% add_xy_position(x = "phase3")
pwcil7plot <- ggviolin(covid, x = "phase3", y = "logil7", 
         title = "IL-7",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-7 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil7, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil7, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-8
ggviolin(covid, x = "phase3", y = "logil8", 
          title = "IL-8",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-8 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil8 <- covid %>% kruskal_test(logil8 ~ phase3)
res.kruskalil8

pwcil8 <- covid %>% dunn_test(logil8 ~ phase3, p.adjust.method = "BH") 
pwcil8

pwcil8 <- pwcil8 %>% add_xy_position(x = "phase3")
pwcil8plot <- ggviolin(covid, x = "phase3", y = "logil8", 
         title = "IL-8",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-8 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil8, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil8, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-9
ggviolin(covid, x = "phase3", y = "logil9", 
          title = "IL-9",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil9 <- covid %>% kruskal_test(logil9 ~ phase3)
res.kruskalil9

pwcil9 <- covid %>% dunn_test(logil9 ~ phase3, p.adjust.method = "BH") 
pwcil9

pwcil9 <- pwcil9 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil9", 
         title = "IL-9",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil9, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-10
ggviolin(covid, x = "phase3", y = "logil10", 
          title = "IL-10",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil10 <- covid %>% kruskal_test(logil10 ~ phase3)
res.kruskalil10

pwcil10 <- covid %>% dunn_test(logil10 ~ phase3, p.adjust.method = "BH") 
pwcil10

pwcil10 <- pwcil10 %>% add_xy_position(x = "phase3")
pwcil10plot <- ggviolin(covid, x = "phase3", y = "logil10", 
         title = "IL-10",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil10, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-12p40
ggviolin(covid, x = "phase3", y = "logil12p40", 
          title = "IL-12p40",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p40 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil12p40 <- covid %>% kruskal_test(logil12p40 ~ phase3)
res.kruskalil12p40

pwcil12p40 <- covid %>% dunn_test(logil12p40 ~ phase3, p.adjust.method = "BH") 
pwcil12p40

pwcil12p40 <- pwcil12p40 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil12p40", 
         title = "IL-12p40",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p40 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil12p40, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p40, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-12p70
ggviolin(covid, x = "phase3", y = "logil12p70", 
          title = "IL-12p70",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p70 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means() 

res.kruskalil12p70 <- covid %>% kruskal_test(logil12p70 ~ phase3)
res.kruskalil12p70

pwcil12p70 <- covid %>% dunn_test(logil12p70 ~ phase3, p.adjust.method = "BH") 
pwcil12p70

pwcil12p70 <- pwcil12p70 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil12p70", 
         title = "IL-12p70",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p70 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil12p70, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p70, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-13
ggviolin(covid, x = "phase3", y = "logil13", 
          title = "IL-13",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-13 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil13 <- covid %>% kruskal_test(logil13 ~ phase3)
res.kruskalil13

pwcil13 <- covid %>% dunn_test(logil13 ~ phase3, p.adjust.method = "BH") 
pwcil13

pwcil13 <- pwcil13 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil13", 
         title = "IL-13",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-13 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil13, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil13, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-15
ggviolin(covid, x = "phase3", y = "logil15", 
          title = "IL-15",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-15 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil15 <- covid %>% kruskal_test(logil15 ~ phase3)
res.kruskalil15

pwcil15 <- covid %>% dunn_test(logil15 ~ phase3, p.adjust.method = "BH") 
pwcil15

pwcil15 <- pwcil15 %>% add_xy_position(x = "phase3")
pwcil15plot <- ggviolin(covid, x = "phase3", y = "logil15", 
         title = "IL-15",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-15 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil15, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil15, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-17A
ggviolin(covid, x = "phase3", y = "logil17a", 
          title = "IL-17A",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil17a <- covid %>% kruskal_test(logil17a ~ phase3)
res.kruskalil17a

pwcil17a <- covid %>% dunn_test(logil17a ~ phase3, p.adjust.method = "BH") 
pwcil17a

pwcil17a <- pwcil17a %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil17a", 
         title = "IL-17A",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil17a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-17E
ggviolin(covid, x = "phase3", y = "logil17eil25", 
          title = "IL-17E",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17E (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil17eil25 <- covid %>% kruskal_test(logil17eil25 ~ phase3)
res.kruskalil17eil25

pwcil17eil25 <- covid %>% dunn_test(logil17eil25 ~ phase3, p.adjust.method = "BH") 
pwcil17eil25

pwcil17eil25 <- pwcil17eil25 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil17eil25", 
         title = "IL-17E",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17E (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil17eil25, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17eil25, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-17F
ggviolin(covid, x = "phase3", y = "logil17f", 
          title = "IL-17F",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17F (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil17f <- covid %>% kruskal_test(logil17f ~ phase3)
res.kruskalil17f

pwcil17f <- covid %>% dunn_test(logil17f ~ phase3, p.adjust.method = "BH") 
pwcil17f

pwcil17f <- pwcil17f %>% add_xy_position(x = "phase3")
pwcil17fplot <- ggviolin(covid, x = "phase3", y = "logil17f", 
         title = "IL-17F",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17F (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil17f, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17f, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-18
ggviolin(covid, x = "phase3", y = "logil18", 
          title = "IL-18",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-18 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil18 <- covid %>% kruskal_test(logil18 ~ phase3)
res.kruskalil18

pwcil18 <- covid %>% dunn_test(logil18 ~ phase3, p.adjust.method = "BH") 
pwcil18

pwcil18 <- pwcil18 %>% add_xy_position(x = "phase3")
pwcil18plot <- ggviolin(covid, x = "phase3", y = "logil18", 
         title = "IL-18",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-18 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil18, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil18, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-22
ggviolin(covid, x = "phase3", y = "logil22", 
          title = "IL-22",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-22 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil22 <- covid %>% kruskal_test(logil22 ~ phase3)
res.kruskalil22

pwcil22 <- covid %>% dunn_test(logil22 ~ phase3, p.adjust.method = "BH") 
pwcil22

pwcil22 <- pwcil22 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logil22", 
         title = "IL-22",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-22 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil22, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil22, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IL-27
ggviolin(covid, x = "phase3", y = "logil27", 
          title = "IL-27",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-27 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalil27 <- covid %>% kruskal_test(logil27 ~ phase3)
res.kruskalil27

pwcil27 <- covid %>% dunn_test(logil27 ~ phase3, p.adjust.method = "BH") 
pwcil27

pwcil27 <- pwcil27 %>% add_xy_position(x = "phase3")
pwcil27plot <- ggviolin(covid, x = "phase3", y = "logil27", 
         title = "IL-27",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-27 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcil27, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil27, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#IP-10/CXCL10
ggviolin(covid, x = "phase3", y = "logip10", 
          title = "CXCL10",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalip10 <- covid %>% kruskal_test(logip10 ~ phase3)
res.kruskalip10

pwcip10 <- covid %>% dunn_test(logip10 ~ phase3, p.adjust.method = "BH") 
pwcip10

pwcip10 <- pwcip10 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logip10", 
         title = "CXCL10",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcip10, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.25, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalip10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MCP-1
ggviolin(covid, x = "phase3", y = "logmcp1", 
          title = "MCP-1",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmcp1 <- covid %>% kruskal_test(logmcp1 ~ phase3)
res.kruskalmcp1

pwcmcp1 <- covid %>% dunn_test(logmcp1 ~ phase3, p.adjust.method = "BH") 
pwcmcp1

pwcmcp1 <- pwcmcp1 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logmcp1", 
         title = "MCP-1",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmcp1, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp1, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MCP-3
ggviolin(covid, x = "phase3", y = "logmcp3", 
          title = "MCP-3",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmcp3 <- covid %>% kruskal_test(logmcp3 ~ phase3)
res.kruskalmcp3

pwcmcp3 <- covid %>% dunn_test(logmcp3 ~ phase3, p.adjust.method = "BH") 
pwcmcp3

pwcmcp3 <- pwcmcp3 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logmcp3", 
         title = "MCP-3",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmcp3, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#M-CSF
ggviolin(covid, x = "phase3", y = "logmcsf", 
          title = "M-CSF",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "M-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmcsf <- covid %>% kruskal_test(logmcsf ~ phase3)
res.kruskalmcsf

pwcmcsf <- covid %>% dunn_test(logmcsf ~ phase3, p.adjust.method = "BH") 
pwcmcsf

pwcmcsf <- pwcmcsf %>% add_xy_position(x = "phase3")
pwcmcsfplot <- ggviolin(covid, x = "phase3", y = "logmcsf", 
         title = "M-CSF",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "M-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MDC
ggviolin(covid, x = "phase3", y = "logmdc", 
          title = "MDC",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MDC (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmdc <- covid %>% kruskal_test(logmdc ~ phase3)
res.kruskalmdc

pwcmdc <- covid %>% dunn_test(logmdc ~ phase3, p.adjust.method = "BH") 
pwcmdc

pwcmdc <- pwcmdc %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logmdc", 
         title = "MDC",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MDC (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmdc, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.25, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmdc, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MIG/CXCL9
ggviolin(covid, x = "phase3", y = "logmigcxcl9", 
          title = "CXCL9",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmigcxcl9 <- covid %>% kruskal_test(logmigcxcl9 ~ phase3)
res.kruskalmigcxcl9

pwcmigcxcl9 <- covid %>% dunn_test(logmigcxcl9 ~ phase3, p.adjust.method = "BH") 
pwcmigcxcl9

pwcmigcxcl9 <- pwcmigcxcl9 %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logmigcxcl9", 
         title = "CXCL9",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmigcxcl9, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmigcxcl9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MIP-1a/CCL3
ggviolin(covid, x = "phase3", y = "logmip1a", 
          title = "CCL3",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmip1a <- covid %>% kruskal_test(logmip1a ~ phase3)
res.kruskalmip1a

pwcmip1a <- covid %>% dunn_test(logmip1a ~ phase3, p.adjust.method = "BH") 
pwcmip1a

pwcmip1a <- pwcmip1a %>% add_xy_position(x = "phase3")
pwcmip1aplot <- ggviolin(covid, x = "phase3", y = "logmip1a", 
         title = "CCL3",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmip1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#MIP-1b/CCL4
ggviolin(covid, x = "phase3", y = "logmip1b", 
          title = "CCL4",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalmip1b <- covid %>% kruskal_test(logmip1b ~ phase3)
res.kruskalmip1b

pwcmip1b <- covid %>% dunn_test(logmip1b ~ phase3, p.adjust.method = "BH") 
pwcmip1b

pwcmip1b <- pwcmip1b %>% add_xy_position(x = "phase3")
pwcmip1bplot <- ggviolin(covid, x = "phase3", y = "logmip1b", 
         title = "CCL4",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcmip1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#PDGF-AA
ggviolin(covid, x = "phase3", y = "logpdgfaa", 
          title = "PDGF-AA",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AA (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalpdgfaa <- covid %>% kruskal_test(logpdgfaa ~ phase3)
res.kruskalpdgfaa

pwcpdgfaa <- covid %>% dunn_test(logpdgfaa ~ phase3, p.adjust.method = "BH") 
pwcpdgfaa

pwcpdgfaa <- pwcpdgfaa %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logpdgfaa", 
         title = "PDGF-AA",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AA (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcpdgfaa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfaa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#PDGF-AB/BB
ggviolin(covid, x = "phase3", y = "logpdgfabbb", 
          title = "PDGF-AB/BB",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalpdgfabbb <- covid %>% kruskal_test(logpdgfabbb ~ phase3)
res.kruskalpdgfabbb

pwcpdgfabbb <- covid %>% dunn_test(logpdgfabbb ~ phase3, p.adjust.method = "BH") 
pwcpdgfabbb

pwcpdgfabbb <- pwcpdgfabbb %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "logpdgfabbb", 
         title = "PDGF-AB/BB",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcpdgfabbb, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfabbb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#RANTES
ggviolin(covid, x = "phase3", y = "lograntes", 
          title = "RANTES",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "RANTES (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalrantes <- covid %>% kruskal_test(lograntes ~ phase3)
res.kruskalrantes

pwcrantes <- covid %>% dunn_test(lograntes ~ phase3, p.adjust.method = "BH") 
pwcrantes

pwcrantes <- pwcrantes %>% add_xy_position(x = "phase3")
ggviolin(covid, x = "phase3", y = "lograntes", 
         title = "RANTES",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "RANTES (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcrantes, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalrantes, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#TGF-alpha
ggviolin(covid, x = "phase3", y = "logtgfa", 
          title = "TGF-\u03b1",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskaltgfa <- covid %>% kruskal_test(logtgfa ~ phase3)
res.kruskaltgfa

pwctgfa <- covid %>% dunn_test(logtgfa ~ phase3, p.adjust.method = "BH") 
pwctgfa

pwctgfa <- pwctgfa %>% add_xy_position(x = "phase3")
pwctgfaplot <- ggviolin(covid, x = "phase3", y = "logtgfa", 
         title = "TGF-\u03b1",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwctgfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltgfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#TNF-alpha
ggviolin(covid, x = "phase3", y = "logtnfa", 
          title = "TNF",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TNF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskaltnfa <- covid %>% kruskal_test(logtnfa ~ phase3)
res.kruskaltnfa

pwctnfa <- covid %>% dunn_test(logtnfa ~ phase3, p.adjust.method = "BH") 
pwctnfa

pwctnfa <- pwctnfa %>% add_xy_position(x = "phase3")
pwctnfaplot <- ggviolin(covid, x = "phase3", y = "logtnfa", 
         title = "TNF",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TNF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwctnfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#Lymphotoxin-alpha/TNF-beta
ggviolin(covid, x = "phase3", y = "logtnfb", 
          title = "Lymphotoxin-\u03b1",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskaltnfb <- covid %>% kruskal_test(logtnfb ~ phase3)
res.kruskaltnfb

pwctnfb <- covid %>% dunn_test(logtnfb ~ phase3, p.adjust.method = "BH") 
pwctnfb

pwctnfb <- pwctnfb %>% add_xy_position(x = "phase3")
pwctnfbplot <- ggviolin(covid, x = "phase3", y = "logtnfb", 
         title = "Lymphotoxin-\u03b1",
         color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwctnfb, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

#VEGF-A
ggviolin(covid, x = "phase3", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  theme(legend.position = "none", text = element_text(size = 9)) +
  stat_compare_means()

res.kruskalvegfa <- covid %>% kruskal_test(logvegfa ~ phase3)
res.kruskalvegfa

pwcvegfa <- covid %>% dunn_test(logvegfa ~ phase3, p.adjust.method = "BH") 
pwcvegfa

pwcvegfa <- pwcvegfa %>% add_xy_position(x = "phase3")
pwcvegfaplot <- ggviolin(covid, x = "phase3", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "phase3", alpha = 0.5, palette =  c("#440154", "#2a788e", "#7ad151"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "Var. A/B", "2" = "A.23/A.23.1", "3" = "Delta")) + 
  stat_pvalue_manual(pwcvegfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalvegfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 9))

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
                                            strata = "phase3",
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

write_tableHTML(tableHTML(tablebiomarkers), file = 'biomarkersphase.html')





