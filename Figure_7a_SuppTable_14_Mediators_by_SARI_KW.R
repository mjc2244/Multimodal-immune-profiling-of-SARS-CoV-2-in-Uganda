#Figure 7a and Table S17 - Mediators by SARI etiology

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Drop asymptomatic covid patients
covidflusari <- subset(combined, min_no_sym==0 | (is.na(min_no_sym)))

#Rename pathogen code 
covidflusari[covidflusari$pathogencode == 1, "pathogencode2"] <- 1
covidflusari[covidflusari$pathogencode == 3, "pathogencode2"] <- 2
covidflusari[covidflusari$pathogencode == 2 | covidflusari$pathogencode == 4, "pathogencode2"] <- 3

#Mediator concentrations across SARI groups
library(ggpubr)
library(rstatix)

#sCD40L
ggviolin(covidflusari, x = "pathogencode2", y = "logscd40l", 
          title = "sCD40L",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "sCD40L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalscd40l <- covidflusari %>% kruskal_test(logscd40l ~ pathogencode2)
res.kruskalscd40l

pwcscd40l <- covidflusari %>% dunn_test(logscd40l ~ pathogencode2, p.adjust.method = "BH") 
pwcscd40l

pwcscd40l <- pwcscd40l %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logscd40l", 
         title = "sCD40L",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "sCD40L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcscd40l, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalscd40l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#EGF
ggviolin(covidflusari, x = "pathogencode2", y = "logegf", 
          title = "EGF",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "EGF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalegf <- covidflusari %>% kruskal_test(logegf ~ pathogencode2)
res.kruskalegf

pwcegf <- covidflusari %>% dunn_test(logegf ~ pathogencode2, p.adjust.method = "BH") 
pwcegf

pwcegf <- pwcegf %>% add_xy_position(x = "pathogencode2")
pwcegfplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logegf", 
         title = "EGF",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "EGF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcegf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalegf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))
pwcegfplot

#Eotaxin
ggviolin(covidflusari, x = "pathogencode2", y = "logeotaxin", 
          title = "Eotaxin",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "Eotaxin (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskaleotaxin <- covidflusari %>% kruskal_test(logeotaxin ~ pathogencode2)
res.kruskaleotaxin

pwceotaxin <- covidflusari %>% dunn_test(logeotaxin ~ pathogencode2, p.adjust.method = "BH") 
pwceotaxin

pwceotaxin <- pwceotaxin %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logeotaxin", 
         title = "Eotaxin",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Eotaxin (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwceotaxin, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaleotaxin, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#FGF2
ggviolin(covidflusari, x = "pathogencode2", y = "logfgf2", 
          title = "FGF-2",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FGF-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalfgf2 <- covidflusari %>% kruskal_test(logfgf2 ~ pathogencode2)
res.kruskalfgf2

pwcfgf2 <- covidflusari %>% dunn_test(logfgf2 ~ pathogencode2, p.adjust.method = "BH") 
pwcfgf2

pwcfgf2 <- pwcfgf2 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logfgf2", 
         title = "FGF-2",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FGF-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcfgf2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfgf2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#FLT-3L
ggviolin(covidflusari, x = "pathogencode2", y = "logflt3l", 
          title = "FLT-3L",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "FLT-3L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalflt3l <- covidflusari %>% kruskal_test(logflt3l ~ pathogencode2)
res.kruskalflt3l

pwcflt3l <- covidflusari %>% dunn_test(logflt3l ~ pathogencode2, p.adjust.method = "BH") 
pwcflt3l

pwcflt3l <- pwcflt3l %>% add_xy_position(x = "pathogencode2")
pwcflt3lplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logflt3l", 
         title = "FLT-3L",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "FLT-3L (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcflt3l, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalflt3l, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#FKN
ggviolin(covidflusari, x = "pathogencode2", y = "logfractalkine", 
          title = "Fractalkine",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Fractalkine (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalfractalkine <- covidflusari %>% kruskal_test(logfractalkine ~ pathogencode2)
res.kruskalfractalkine

pwcfractalkine <- covidflusari %>% dunn_test(logfractalkine ~ pathogencode2, p.adjust.method = "BH") 
pwcfractalkine

pwcfractalkine <- pwcfractalkine %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logfractalkine", 
         title = "Fractalkine",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Fractalkine (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcfractalkine, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalfractalkine, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#G-CSF
ggviolin(covidflusari, x = "pathogencode2", y = "loggcsf", 
          title = "G-CSF",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "G-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalgcsf <- covidflusari %>% kruskal_test(loggcsf ~ pathogencode2)
res.kruskalgcsf

pwcgcsf <- covidflusari %>% dunn_test(loggcsf ~ pathogencode2, p.adjust.method = "BH") 
pwcgcsf

pwcgcsf <- pwcgcsf %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "loggcsf", 
         title = "G-CSF",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "G-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcgcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#GM-CSF
ggviolin(covidflusari, x = "pathogencode2", y = "loggmcsf", 
          title = "GM-CSF",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "GM-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalgmcsf <- covidflusari %>% kruskal_test(loggmcsf ~ pathogencode2)
res.kruskalgmcsf

pwcgmcsf <- covidflusari %>% dunn_test(loggmcsf ~ pathogencode2, p.adjust.method = "BH") 
pwcgmcsf

pwcgmcsf <- pwcgmcsf %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "loggmcsf", 
         title = "GM-CSF",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "GM-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcgmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#GRO-alpha/CXCL1
ggviolin(covidflusari, x = "pathogencode2", y = "loggroalpha", 
          title = "CXCL1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalgroalpha <- covidflusari %>% kruskal_test(loggroalpha ~ pathogencode2)
res.kruskalgroalpha

pwcgroalpha <- covidflusari %>% dunn_test(loggroalpha ~ pathogencode2, p.adjust.method = "BH") 
pwcgroalpha

pwcgroalpha <- pwcgroalpha %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "loggroalpha", 
         title = "CXCL1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcgroalpha, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalgroalpha, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IFN-a2
ggviolin(covidflusari, x = "pathogencode2", y = "logifna2", 
          title = "IFN-\u03b12",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalifna2 <- covidflusari %>% kruskal_test(logifna2 ~ pathogencode2)
res.kruskalifna2

pwcifna2 <- covidflusari %>% dunn_test(logifna2 ~ pathogencode2, p.adjust.method = "BH") 
pwcifna2

pwcifna2 <- pwcifna2 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logifna2", 
         title = "IFN-\u03b12",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b12 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcifna2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifna2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IFN-gamma
ggviolin(covidflusari, x = "pathogencode2", y = "logifngamma", 
          title = "IFN-\u03b3",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalifngamma <- covidflusari %>% kruskal_test(logifngamma ~ pathogencode2)
res.kruskalifngamma

pwcifngamma <- covidflusari %>% dunn_test(logifngamma ~ pathogencode2, p.adjust.method = "BH") 
pwcifngamma

pwcifngamma <- pwcifngamma %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logifngamma", 
         title = "IFN-\u03b3",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IFN-\u03b3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcifngamma, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalifngamma, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-1a
ggviolin(covidflusari, x = "pathogencode2", y = "logil1a", 
          title = "IL-1\u03b1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil1a <- covidflusari %>% kruskal_test(logil1a ~ pathogencode2)
res.kruskalil1a

pwcil1a <- covidflusari %>% dunn_test(logil1a ~ pathogencode2, p.adjust.method = "BH") 
pwcil1a

pwcil1a <- pwcil1a %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil1a", 
         title = "IL-1\u03b1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-1b
ggviolin(covidflusari, x = "pathogencode2", y = "logil1b", 
          title = "IL-1\u03B2",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil1b <- covidflusari %>% kruskal_test(logil1b ~ pathogencode2)
res.kruskalil1b

pwcil1b <- covidflusari %>% dunn_test(logil1b ~ pathogencode2, p.adjust.method = "BH") 
pwcil1b

pwcil1b <- pwcil1b %>% add_xy_position(x = "pathogencode2")
pwcil1bplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil1b", 
         title = "IL-1\u03B2",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1\u03B2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-1Ra
ggviolin(covidflusari, x = "pathogencode2", y = "logil1ra", 
          title = "IL-1Ra",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5,
          ylab = "IL-1Ra (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil1ra <- covidflusari %>% kruskal_test(logil1ra ~ pathogencode2)
res.kruskalil1ra

pwcil1ra <- covidflusari %>% dunn_test(logil1ra ~ pathogencode2, p.adjust.method = "BH") 
pwcil1ra

pwcil1ra <- pwcil1ra %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil1ra", 
         title = "IL-1Ra",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-1Ra (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil1ra, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil1ra, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-2
ggviolin(covidflusari, x = "pathogencode2", y = "logil2", 
          title = "IL-2",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil2 <- covidflusari %>% kruskal_test(logil2 ~ pathogencode2)
res.kruskalil2

pwcil2 <- covidflusari %>% dunn_test(logil2 ~ pathogencode2, p.adjust.method = "BH") 
pwcil2

pwcil2 <- pwcil2 %>% add_xy_position(x = "pathogencode2")
pwcil2plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil2", 
         title = "IL-2",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-2 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil2, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil2, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-3
ggviolin(covidflusari, x = "pathogencode2", y = "logil3", 
          title = "IL-3",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil3 <- covidflusari %>% kruskal_test(logil3 ~ pathogencode2)
res.kruskalil3

pwcil3 <- covidflusari %>% dunn_test(logil3 ~ pathogencode2, p.adjust.method = "BH") 
pwcil3

pwcil3 <- pwcil3 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil3", 
         title = "IL-3",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil3, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-4
ggviolin(covidflusari, x = "pathogencode2", y = "logil4", 
          title = "IL-4",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil4 <- covidflusari %>% kruskal_test(logil4 ~ pathogencode2)
res.kruskalil4

pwcil4 <- covidflusari %>% dunn_test(logil4 ~ pathogencode2, p.adjust.method = "BH") 
pwcil4

pwcil4 <- pwcil4 %>% add_xy_position(x = "pathogencode2")
pwcil4plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil4", 
         title = "IL-4",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil4, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil4, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-5
ggviolin(covidflusari, x = "pathogencode2", y = "logil5", 
          title = "IL-5",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-5 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil5 <- covidflusari %>% kruskal_test(logil5 ~ pathogencode2)
res.kruskalil5

pwcil5 <- covidflusari %>% dunn_test(logil5 ~ pathogencode2, p.adjust.method = "BH") 
pwcil5

pwcil5 <- pwcil5 %>% add_xy_position(x = "pathogencode2")
pwcil5plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil5", 
         title = "IL-5",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-5 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil5, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0) + 
  labs(
    subtitle = get_test_label(res.kruskalil5, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-6
ggviolin(covidflusari, x = "pathogencode2", y = "logil6", 
          title = "IL-6",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-6 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil6 <- covidflusari %>% kruskal_test(logil6 ~ pathogencode2)
res.kruskalil6

pwcil6 <- covidflusari %>% dunn_test(logil6 ~ pathogencode2, p.adjust.method = "BH") 
pwcil6

pwcil6 <- pwcil6 %>% add_xy_position(x = "pathogencode2")
pwcil6plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil6", 
         title = "IL-6",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-6 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil6, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.8, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil6, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-7
ggviolin(covidflusari, x = "pathogencode2", y = "logil7", 
          title = "IL-7",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-7 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil7 <- covidflusari %>% kruskal_test(logil7 ~ pathogencode2)
res.kruskalil7

pwcil7 <- covidflusari %>% dunn_test(logil7 ~ pathogencode2, p.adjust.method = "BH") 
pwcil7

pwcil7 <- pwcil7 %>% add_xy_position(x = "pathogencode2")
pwcil7plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil7", 
         title = "IL-7",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-7 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil7, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil7, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-8
ggviolin(covidflusari, x = "pathogencode2", y = "logil8", 
          title = "IL-8",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-8 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil8 <- covidflusari %>% kruskal_test(logil8 ~ pathogencode2)
res.kruskalil8

pwcil8 <- covidflusari %>% dunn_test(logil8 ~ pathogencode2, p.adjust.method = "BH") 
pwcil8

pwcil8 <- pwcil8 %>% add_xy_position(x = "pathogencode2")
pwcil8plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil8", 
         title = "IL-8",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-8 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil8, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil8, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-9
ggviolin(covidflusari, x = "pathogencode2", y = "logil9", 
          title = "IL-9",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil9 <- covidflusari %>% kruskal_test(logil9 ~ pathogencode2)
res.kruskalil9

pwcil9 <- covidflusari %>% dunn_test(logil9 ~ pathogencode2, p.adjust.method = "BH") 
pwcil9

pwcil9 <- pwcil9 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil9", 
         title = "IL-9",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil9, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-10
ggviolin(covidflusari, x = "pathogencode2", y = "logil10", 
          title = "IL-10",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil10 <- covidflusari %>% kruskal_test(logil10 ~ pathogencode2)
res.kruskalil10

pwcil10 <- covidflusari %>% dunn_test(logil10 ~ pathogencode2, p.adjust.method = "BH") 
pwcil10

pwcil10 <- pwcil10 %>% add_xy_position(x = "pathogencode2")
pwcil10plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil10", 
         title = "IL-10",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil10, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-12p40
ggviolin(covidflusari, x = "pathogencode2", y = "logil12p40", 
          title = "IL-12p40",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p40 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil12p40 <- covidflusari %>% kruskal_test(logil12p40 ~ pathogencode2)
res.kruskalil12p40

pwcil12p40 <- covidflusari %>% dunn_test(logil12p40 ~ pathogencode2, p.adjust.method = "BH") 
pwcil12p40

pwcil12p40 <- pwcil12p40 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil12p40", 
         title = "IL-12p40",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p40 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil12p40, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p40, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-12p70
ggviolin(covidflusari, x = "pathogencode2", y = "logil12p70", 
          title = "IL-12p70",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-12p70 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means() 

res.kruskalil12p70 <- covidflusari %>% kruskal_test(logil12p70 ~ pathogencode2)
res.kruskalil12p70

pwcil12p70 <- covidflusari %>% dunn_test(logil12p70 ~ pathogencode2, p.adjust.method = "BH") 
pwcil12p70

pwcil12p70 <- pwcil12p70 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil12p70", 
         title = "IL-12p70",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-12p70 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil12p70, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil12p70, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-13
ggviolin(covidflusari, x = "pathogencode2", y = "logil13", 
          title = "IL-13",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-13 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil13 <- covidflusari %>% kruskal_test(logil13 ~ pathogencode2)
res.kruskalil13

pwcil13 <- covidflusari %>% dunn_test(logil13 ~ pathogencode2, p.adjust.method = "BH") 
pwcil13

pwcil13 <- pwcil13 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil13", 
         title = "IL-13",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-13 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil13, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil13, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-15
ggviolin(covidflusari, x = "pathogencode2", y = "logil15", 
          title = "IL-15",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-15 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil15 <- covidflusari %>% kruskal_test(logil15 ~ pathogencode2)
res.kruskalil15

pwcil15 <- covidflusari %>% dunn_test(logil15 ~ pathogencode2, p.adjust.method = "BH") 
pwcil15

pwcil15 <- pwcil15 %>% add_xy_position(x = "pathogencode2")
pwcil15plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil15", 
         title = "IL-15",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-15 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil15, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil15, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-17A
ggviolin(covidflusari, x = "pathogencode2", y = "logil17a", 
          title = "IL-17A",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil17a <- covidflusari %>% kruskal_test(logil17a ~ pathogencode2)
res.kruskalil17a

pwcil17a <- covidflusari %>% dunn_test(logil17a ~ pathogencode2, p.adjust.method = "BH") 
pwcil17a

pwcil17a <- pwcil17a %>% add_xy_position(x = "pathogencode2")
pwcil17aplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil17a", 
         title = "IL-17A",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil17a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-17E
ggviolin(covidflusari, x = "pathogencode2", y = "logil17eil25", 
          title = "IL-17E",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17E (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil17eil25 <- covidflusari %>% kruskal_test(logil17eil25 ~ pathogencode2)
res.kruskalil17eil25

pwcil17eil25 <- covidflusari %>% dunn_test(logil17eil25 ~ pathogencode2, p.adjust.method = "BH") 
pwcil17eil25

pwcil17eil25 <- pwcil17eil25 %>% add_xy_position(x = "pathogencode2")
pwcil17eil25plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil17eil25", 
         title = "IL-17E",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17E (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil17eil25, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17eil25, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))
pwcil17eil25plot

#IL-17F
ggviolin(covidflusari, x = "pathogencode2", y = "logil17f", 
          title = "IL-17F",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-17F (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil17f <- covidflusari %>% kruskal_test(logil17f ~ pathogencode2)
res.kruskalil17f

pwcil17f <- covidflusari %>% dunn_test(logil17f ~ pathogencode2, p.adjust.method = "BH") 
pwcil17f

pwcil17f <- pwcil17f %>% add_xy_position(x = "pathogencode2")
pwcil17fplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil17f", 
         title = "IL-17F",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-17F (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil17f, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil17f, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-18
ggviolin(covidflusari, x = "pathogencode2", y = "logil18", 
          title = "IL-18",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-18 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil18 <- covidflusari %>% kruskal_test(logil18 ~ pathogencode2)
res.kruskalil18

pwcil18 <- covidflusari %>% dunn_test(logil18 ~ pathogencode2, p.adjust.method = "BH") 
pwcil18

pwcil18 <- pwcil18 %>% add_xy_position(x = "pathogencode2")
pwcil18plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil18", 
         title = "IL-18",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-18 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil18, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil18, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-22
ggviolin(covidflusari, x = "pathogencode2", y = "logil22", 
          title = "IL-22",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-22 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil22 <- covidflusari %>% kruskal_test(logil22 ~ pathogencode2)
res.kruskalil22

pwcil22 <- covidflusari %>% dunn_test(logil22 ~ pathogencode2, p.adjust.method = "BH") 
pwcil22

pwcil22 <- pwcil22 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logil22", 
         title = "IL-22",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-22 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil22, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil22, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IL-27
ggviolin(covidflusari, x = "pathogencode2", y = "logil27", 
          title = "IL-27",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "IL-27 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalil27 <- covidflusari %>% kruskal_test(logil27 ~ pathogencode2)
res.kruskalil27

pwcil27 <- covidflusari %>% dunn_test(logil27 ~ pathogencode2, p.adjust.method = "BH") 
pwcil27

pwcil27 <- pwcil27 %>% add_xy_position(x = "pathogencode2")
pwcil27plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logil27", 
         title = "IL-27",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "IL-27 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcil27, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalil27, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#IP-10/CXCL10
ggviolin(covidflusari, x = "pathogencode2", y = "logip10", 
          title = "CXCL10",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalip10 <- covidflusari %>% kruskal_test(logip10 ~ pathogencode2)
res.kruskalip10

pwcip10 <- covidflusari %>% dunn_test(logip10 ~ pathogencode2, p.adjust.method = "BH") 
pwcip10

pwcip10 <- pwcip10 %>% add_xy_position(x = "pathogencode2")
pwcip10plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logip10", 
         title = "CXCL10",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL10 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcip10, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.85, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalip10, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))
pwcip10plot

#MCP-1
ggviolin(covidflusari, x = "pathogencode2", y = "logmcp1", 
          title = "MCP-1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmcp1 <- covidflusari %>% kruskal_test(logmcp1 ~ pathogencode2)
res.kruskalmcp1

pwcmcp1 <- covidflusari %>% dunn_test(logmcp1 ~ pathogencode2, p.adjust.method = "BH") 
pwcmcp1

pwcmcp1 <- pwcmcp1 %>% add_xy_position(x = "pathogencode2")
pwcmcp1plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logmcp1", 
         title = "MCP-1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmcp1, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp1, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#MCP-3
ggviolin(covidflusari, x = "pathogencode2", y = "logmcp3", 
          title = "MCP-3",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MCP-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmcp3 <- covidflusari %>% kruskal_test(logmcp3 ~ pathogencode2)
res.kruskalmcp3

pwcmcp3 <- covidflusari %>% dunn_test(logmcp3 ~ pathogencode2, p.adjust.method = "BH") 
pwcmcp3

pwcmcp3 <- pwcmcp3 %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logmcp3", 
         title = "MCP-3",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MCP-3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmcp3, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcp3, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#M-CSF
ggviolin(covidflusari, x = "pathogencode2", y = "logmcsf", 
          title = "M-CSF",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "M-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmcsf <- covidflusari %>% kruskal_test(logmcsf ~ pathogencode2)
res.kruskalmcsf

pwcmcsf <- covidflusari %>% dunn_test(logmcsf ~ pathogencode2, p.adjust.method = "BH") 
pwcmcsf

pwcmcsf <- pwcmcsf %>% add_xy_position(x = "pathogencode2")
pwcmcsfplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logmcsf", 
         title = "M-CSF",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "M-CSF (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmcsf, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmcsf, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#MDC
ggviolin(covidflusari, x = "pathogencode2", y = "logmdc", 
          title = "MDC",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "MDC (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmdc <- covidflusari %>% kruskal_test(logmdc ~ pathogencode2)
res.kruskalmdc

pwcmdc <- covidflusari %>% dunn_test(logmdc ~ pathogencode2, p.adjust.method = "BH") 
pwcmdc

pwcmdc <- pwcmdc %>% add_xy_position(x = "pathogencode2")
ggviolin(covidflusari, x = "pathogencode2", y = "logmdc", 
         title = "MDC",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "MDC (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmdc, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.25, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmdc, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#MIG/CXCL9
ggviolin(covidflusari, x = "pathogencode2", y = "logmigcxcl9", 
          title = "CXCL9",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CXCL9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmigcxcl9 <- covidflusari %>% kruskal_test(logmigcxcl9 ~ pathogencode2)
res.kruskalmigcxcl9

pwcmigcxcl9 <- covidflusari %>% dunn_test(logmigcxcl9 ~ pathogencode2, p.adjust.method = "BH") 
pwcmigcxcl9

pwcmigcxcl9 <- pwcmigcxcl9 %>% add_xy_position(x = "pathogencode2")
pwcmigcxcl9plot <- ggviolin(covidflusari, x = "pathogencode2", y = "logmigcxcl9", 
         title = "CXCL9",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CXCL9 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmigcxcl9, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmigcxcl9, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#MIP-1a/CCL3
ggviolin(covidflusari, x = "pathogencode2", y = "logmip1a", 
          title = "CCL3",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmip1a <- covidflusari %>% kruskal_test(logmip1a ~ pathogencode2)
res.kruskalmip1a

pwcmip1a <- covidflusari %>% dunn_test(logmip1a ~ pathogencode2, p.adjust.method = "BH") 
pwcmip1a

pwcmip1a <- pwcmip1a %>% add_xy_position(x = "pathogencode2")
pwcmip1aplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logmip1a", 
         title = "CCL3",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL3 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmip1a, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1a, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#MIP-1b/CCL4
ggviolin(covidflusari, x = "pathogencode2", y = "logmip1b", 
          title = "CCL4",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "CCL4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalmip1b <- covidflusari %>% kruskal_test(logmip1b ~ pathogencode2)
res.kruskalmip1b

pwcmip1b <- covidflusari %>% dunn_test(logmip1b ~ pathogencode2, p.adjust.method = "BH") 
pwcmip1b

pwcmip1b <- pwcmip1b %>% add_xy_position(x = "pathogencode2")
pwcmip1bplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logmip1b", 
         title = "CCL4",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "CCL4 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcmip1b, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalmip1b, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#PDGF-AA
ggviolin(covidflusari, x = "pathogencode2", y = "logpdgfaa", 
          title = "PDGF-AA",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AA (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalpdgfaa <- covidflusari %>% kruskal_test(logpdgfaa ~ pathogencode2)
res.kruskalpdgfaa

pwcpdgfaa <- covidflusari %>% dunn_test(logpdgfaa ~ pathogencode2, p.adjust.method = "BH") 
pwcpdgfaa

pwcpdgfaa <- pwcpdgfaa %>% add_xy_position(x = "pathogencode2")
pwcpdgfaaplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logpdgfaa", 
         title = "PDGF-AA",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AA (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcpdgfaa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfaa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#PDGF-AB/BB
ggviolin(covidflusari, x = "pathogencode2", y = "logpdgfabbb", 
          title = "PDGF-AB/BB",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalpdgfabbb <- covidflusari %>% kruskal_test(logpdgfabbb ~ pathogencode2)
res.kruskalpdgfabbb

pwcpdgfabbb <- covidflusari %>% dunn_test(logpdgfabbb ~ pathogencode2, p.adjust.method = "BH") 
pwcpdgfabbb

pwcpdgfabbb <- pwcpdgfabbb %>% add_xy_position(x = "pathogencode2")
pwcpdgfabbbplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logpdgfabbb", 
         title = "PDGF-AB/BB",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "PDGF-AB/BB (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcpdgfabbb, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalpdgfabbb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#RANTES
ggviolin(covidflusari, x = "pathogencode2", y = "lograntes", 
          title = "RANTES",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "RANTES (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalrantes <- covidflusari %>% kruskal_test(lograntes ~ pathogencode2)
res.kruskalrantes

pwcrantes <- covidflusari %>% dunn_test(lograntes ~ pathogencode2, p.adjust.method = "BH") 
pwcrantes

pwcrantes <- pwcrantes %>% add_xy_position(x = "pathogencode2")
pwcrantesplot <- ggviolin(covidflusari, x = "pathogencode2", y = "lograntes", 
         title = "RANTES",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "RANTES (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcrantes, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalrantes, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#TGF-alpha
ggviolin(covidflusari, x = "pathogencode2", y = "logtgfa", 
          title = "TGF-\u03b1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskaltgfa <- covidflusari %>% kruskal_test(logtgfa ~ pathogencode2)
res.kruskaltgfa

pwctgfa <- covidflusari %>% dunn_test(logtgfa ~ pathogencode2, p.adjust.method = "BH") 
pwctgfa

pwctgfa <- pwctgfa %>% add_xy_position(x = "pathogencode2")
pwctgfaplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logtgfa", 
         title = "TGF-\u03b1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TGF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwctgfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltgfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#TNF-alpha
ggviolin(covidflusari, x = "pathogencode2", y = "logtnfa", 
          title = "TNF-\u03b1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "TNF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskaltnfa <- covidflusari %>% kruskal_test(logtnfa ~ pathogencode2)
res.kruskaltnfa

pwctnfa <- covidflusari %>% dunn_test(logtnfa ~ pathogencode2, p.adjust.method = "BH") 
pwctnfa

pwctnfa <- pwctnfa %>% add_xy_position(x = "pathogencode2")
pwctnfaplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logtnfa", 
         title = "TNF-\u03b1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "TNF-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwctnfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#Lymphotoxin-alpha/TNF-beta
ggviolin(covidflusari, x = "pathogencode2", y = "logtnfb", 
          title = "Lymphotoxin-\u03b1",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskaltnfb <- covidflusari %>% kruskal_test(logtnfb ~ pathogencode2)
res.kruskaltnfb

pwctnfb <- covidflusari %>% dunn_test(logtnfb ~ pathogencode2, p.adjust.method = "BH") 
pwctnfb

pwctnfb <- pwctnfb %>% add_xy_position(x = "pathogencode2")
pwctnfbplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logtnfb", 
         title = "Lymphotoxin-\u03b1",
         color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
         order = c("1", "2", "3"),
         add = "jitter", draw_quantiles = 0.5, size = 1,
         ylab = "Lymphotoxin-\u03b1 (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwctnfb, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.02, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskaltnfb, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

#VEGF-A
ggviolin(covidflusari, x = "pathogencode2", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  theme(legend.position = "none", text = element_text(size = 10)) +
  stat_compare_means()

res.kruskalvegfa <- covidflusari %>% kruskal_test(logvegfa ~ pathogencode2)
res.kruskalvegfa

pwcvegfa <- covidflusari %>% dunn_test(logvegfa ~ pathogencode2, p.adjust.method = "BH") 
pwcvegfa

pwcvegfa <- pwcvegfa %>% add_xy_position(x = "pathogencode2")
pwcvegfaplot <- ggviolin(covidflusari, x = "pathogencode2", y = "logvegfa", 
          title = "VEGF-A",
          color = "black", fill = "pathogencode2", alpha = 0.5, palette =  c("#7570b3", "#d95f02", "#1b9e77"),
          order = c("1", "2", "3"),
          add = "jitter", draw_quantiles = 0.5, size = 1,
          ylab = "VEGF-A (log10 pg/ml)", xlab = "") +
  scale_x_discrete(labels=c("1" = "COVID-19", "2" = "Influenza", "3" = "Non-Influenza")) + 
  stat_pvalue_manual(pwcvegfa, hide.ns = TRUE, bracket.size = 0.5, label.size = 5, bracket.nudge.y = 0.5, tip.length = 0.03, step.increase = 0.03) + 
  labs(
    subtitle = get_test_label(res.kruskalvegfa, detailed = FALSE)) + 
  theme(legend.position = "none", text = element_text(size = 10))

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

tablebiomarkerspathogen <- CreateTableOne(vars = listbiomarkers, 
                                            data = covidflusari, 
                                            strata = "pathogencode2",
                                            addOverall = TRUE)
tablebiomarkerspathogen
summary(tablebiomarkerspathogen)

tablebiomarkerspathogen <- print(tablebiomarkerspathogen, nonnormal = c("logscd40l", 
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

write_tableHTML(tableHTML(tablebiomarkerspathogen), file = 'biomarkerspathogen.html')





