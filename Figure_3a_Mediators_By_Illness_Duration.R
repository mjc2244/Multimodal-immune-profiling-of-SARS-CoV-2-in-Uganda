#Figure 3a - Mediators by Illness Duration - Robust regression

#Clear R environment
rm(list = ls())

#Import master data, PID to fixed row identifier 
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid, exclude asymptomatic 
covid <- subset(combined, pathogencode==1)
covid <- subset(covid, min_no_sym==0)

#Combine mild-moderate COVID to compare with severe 
covid[covid$whoseverity_imp==2 | covid$whoseverity_imp==3, "whoseverity_imp_bin"] <- 0
covid[covid$whoseverity_imp==4, "whoseverity_imp_bin"] <- 1

covid$whoseverity_imp_bin <- as.factor(covid$whoseverity_imp_bin)
levels(covid$whoseverity_imp_bin)

#Determine outliers in pre-enrollment illness duration variable
library(ggplot2)
library(ggsci)
library(MASS)
ggplot(covid, aes(x=illnessdurationenroll_imp, y=logil6)) + geom_point() + geom_text(label=rownames(covid))

#PIDs CL221, CL458, CL363 are outliers - remove from this analysis
coviddurationout <- covid[row.names(covid) != "CL221",]
coviddurationout <- coviddurationout[row.names(coviddurationout) != "CL458",]
coviddurationout <- coviddurationout[row.names(coviddurationout) != "CL363",]

#Robust regression plots stratified by illness severity
#sCD40L
scd40llong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logscd40l, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="sCD40L", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
scd40llong

#EGF
egflong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logegf, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="EGF", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
egflong

#Eotaxin
eotaxinlong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logeotaxin, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="Eotaxin", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
eotaxinlong

#FGF-2
fgf2long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logfgf2, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="FGF-2", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
fgf2long

#FLT-3L
flt3llong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logflt3l, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="FLT-3L", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
flt3llong

#Fractalkine
fractalkinelong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logfractalkine, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="Fractalkine", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
fractalkinelong

#G-CSF
gcsflong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=loggcsf, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="G-CSF", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
gcsflong

#GM-CSF
gmcsflong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=loggmcsf, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="GM-CSF", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
gmcsflong

#GRO-alpha/CXCL1
groalphalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=loggroalpha, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="CXCL1", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
groalphalong

#IFN-a2
ifna2long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logifna2, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IFN-\u03b12", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
ifna2long

#IFN-gamma
ifngammalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logifngamma, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IFN-\u03b3", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
ifngammalong

#IL-1a
il1along <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil1a, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-1\u03b1", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il1along

#IL-1b
il1blong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil1b, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-1\u03B2", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il1blong

#IL-1Ra
il1ralong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil1ra, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-1Ra", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il1ralong

#IL-2
il2long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil2, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-2", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il2long

#IL-3
il3long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil3, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-3", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il3long

#IL-4
il4long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil4, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-4", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il4long

#IL-5
il5long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil5, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-5", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il5long

#IL-6
il6long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil6, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-6", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il6long

#IL-7
il7long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil7, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-7", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il7long

#IL-8
il8long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil8, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-8", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il8long

#IL-9
il9long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil9, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-9", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il9long

#IL-10
il10long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil10, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-10", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il10long

#IL-12p40
il12p40long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil12p40, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-12p40", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il12p40long

#IL-12p70
il12p70long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil12p70, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-12p70", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il12p70long

#IL-13
il13long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil13, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-13", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il13long

#IL-15
il15long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil15, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-15", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il15long

#IL-17A
il17along <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil17a, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-17A", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il17along

#IL-17E
il17eil25long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil17eil25, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-17E", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il17eil25long

#IL-17F
il17flong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil17f, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-17F", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il17flong

#IL-18
il18long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil18, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-18", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il18long

#IL-22
il22long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil22, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-22", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il22long

#IL-27
il27long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logil27, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="IL-27", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
il27long

#CXCL10/IP-10
ip10long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logip10, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="CXCL10", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
ip10long

#MCP-1
mcp1long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmcp1, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="MCP-1", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mcp1long

#MCP-3
mcp3long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmcp3, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="MCP-3", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mcp3long

#M-CSF
mcsflong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmcsf, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="M-CSF", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mcsflong

#MDC
mdclong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmdc, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="MDC", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mdclong

#MIG/CXCL9
migcxcl9long <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmigcxcl9, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="CXCL9", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
migcxcl9long

#MIP-1a/CCL3
mip1along <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmip1a, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="CCL3", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mip1along

#MIP-1b/CCL4
mip1blong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logmip1b, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="CCL4", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
mip1blong


#PDGF-AA
pdgfaalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logpdgfaa, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="PDGF-AA", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
pdgfaalong

#PDGF-AB/BB
pdgfabbblong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logpdgfabbb, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="PDGF-AB/BB", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
pdgfabbblong

#RANTES
ranteslong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=lograntes, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="RANTES", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
ranteslong

#TGF-alpha 
tgfalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logtgfa, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="TGF-\u03b1", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
tgfalong

#TNF-alpha 
tnfalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logtnfa, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="TNF", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
tnfalong

#TNF-beta/LT-alpha
tnfblong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logtnfb, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="Lymphotoxin-\u03b1", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
tnfblong

#VEGF-A
vegfalong <- ggplot(coviddurationout, aes(x=illnessdurationenroll_imp, y=logvegfa, color=whoseverity_imp_bin)) + 
  geom_point() + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe"))+
  geom_smooth(method=rlm, fullrange=TRUE, aes(fill=whoseverity_imp_bin)) + 
  scale_color_manual(values=c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"), labels = c("Mild-Moderate", "Severe")) +
  labs(title="VEGF-A", x="Days since symptom onset", y = "Concentration (log10 pg/ml)") + 
  labs(color = "Illness Severity", fill = "Illness Severity") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.5, linetype = "blank"),
        panel.grid.minor = element_line(size = 0.25, linetype = "blank"),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom" ) + 
  theme(legend.text = element_text(color = "black", size = 12))
vegfalong
