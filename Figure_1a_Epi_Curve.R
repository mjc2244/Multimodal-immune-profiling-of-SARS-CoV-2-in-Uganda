#Figure 1a - Epi Curve 

#Clear R environment
rm(list = ls())

#Import COVID data
coviddates <- read.csv(file.choose(), header=TRUE)

#PID to fixed row ID
rownames(coviddates) = coviddates$pid
coviddates$pid = NULL

#Load packages 
library(incidence)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)
library(lubridate)
library(dplyr)

#Weekly curve - ungrouped  
i.7 <- incidence(coviddates$dateofhospitalization, 
                 interval = "1 week")
plot(i.7)

class(coviddates$dateofhospitalization)


#Weekly curve - grouped by WHO clinical severity 
i.7.severity <- incidence(coviddates$dateofhospitalization, 
                           interval = "1 sunday week", 
                           groups = coviddates$whoseverity_imp)

i.7.severity.plot <- plot(i.7.severity, 
     stack = TRUE, 
     border = "black", 
     ylab = "No. enrolled participants", 
     labels_week = FALSE,
     n_breaks = 12) + 
  scale_fill_manual(values = c("#01665e", "blue4", "#BC3C29FF", "#e6ab02"),  
                    name = "",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Asymptomatic", "Mild", "Moderate", "Severe")) +
  theme(axis.text.y = element_text(color = "black", size = 14, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, face = "plain")) + 
  theme(legend.title = element_text(color = "black", size = 14, face = "plain"), 
        legend.text = element_text(color = "black", size = 14, face = "plain"))
i.7.severity.plot <- i.7.severity.plot + cowplot::theme_cowplot() + theme(legend.position = "top") + theme(legend.direction = "horizontal") + theme(legend.text=element_text(size=14)) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14, face = "plain")) + theme(axis.text.y = element_text(size = 14, face = "plain"))
i.7.severity.plot

#Add lines for pandemic phases
i.7.severity.plot <- i.7.severity.plot + geom_rect(xmin = as.Date("2020-03-22"), ymin = 0, xmax = as.Date("2020-09-04"), ymax = 30,
          fill = NA, col = "grey", linetype = 1, size = 1)
i.7.severity.plot <- i.7.severity.plot + geom_rect(xmin = as.Date("2020-09-04"), ymin = 0, xmax = as.Date("2021-04-28"), ymax = 30,
                                                     fill = NA, col = "grey", linetype = 1, size = 1)
i.7.severity.plot <- i.7.severity.plot + geom_rect(xmin = as.Date("2021-04-28"), ymin = 0, xmax = as.Date("2021-07-19"), ymax = 30,
                                                     fill = NA, col = "grey", linetype = 1, size = 1)
i.7.severity.plot

i.7.severity.plot <- i.7.severity.plot + annotate("text", x = as.Date("2020-06-06"), y = 30.6, label = "Varied A/B Lineages", size = 5.5)
i.7.severity.plot <- i.7.severity.plot + annotate("text", x = as.Date("2020-12-25"), y = 30.6, label = "A.23/A.23.1", size = 5.5)
i.7.severity.plot <- i.7.severity.plot + annotate("text", x = as.Date("2021-06-01"), y = 30.6, label = "Delta*", size = 5.5)
i.7.severity.plot

#Plot for national cases
#Import national data
coviddateswho <- read.csv(file.choose(), header=TRUE)
coviddateswho$Date <- as.Date(coviddateswho$Date)

#Generate weekly dates without using repeating format 
coviddateswhoweekly <- coviddateswho %>% group_by(week=floor_date(coviddateswho$Date, "week"))

#Plot national cases using parameters from weekly study curve using geom_point
i.plot.7.who <-ggplot(data=coviddateswhoweekly, aes(week, New_cases)) +
  stat_summary(fun = sum, geom = "line", size = 0.75) + 
  incidence::scale_x_incidence(i.7.severity, n_breaks = nrow(i.7.severity)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_y_continuous(position = "right") +             
  labs(x = "",
       y = "No. national cases") 
i.plot.7.who <- i.plot.7.who + theme(axis.text.y = element_text(size = 14, face = "plain"))
i.plot.7.who <- i.plot.7.who + theme_cowplot()
i.plot.7.who

#Combine study, national, and variant plots
aligned_plots <- cowplot::align_plots(i.plot.7.who + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 14, face = "plain"), axis.ticks.x = element_blank(), axis.title.x = element_blank()), i.7.severity.plot, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
aligned_plots





