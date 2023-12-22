#Figure 2b - GSEA - Metabolic, Thrombotic, Endothelial Pathways Severe vs. non-severe COVID

#Clear R environment
rm(list = ls())

#Import enriched metabolic, thrombotic, endothelial pathways 
metabolic <- read.csv(file.choose(), header=TRUE)

#Barplot
library(ggplot2)
metabolicbarplot <- ggplot(metabolic, aes(x = process, y = nes, color = process)) + 
  geom_point(aes(reorder(process,nes)), size = 5.5, position = position_dodge(preserve = "single")) 
metabolicbarplot <- metabolicbarplot + coord_flip() + theme_bw() 
metabolicbarplot <- metabolicbarplot + theme(legend.position = "none") + ylim(-2, 2) +
  scale_color_manual("pathway", values = c("Pos. regulation of ROS metabolic process" = "#1874d2",
                                          "Cellular polysaccharide catabolic process" = "#1874d2",
                                          "Platelet aggregation" = "#1874d2",
                                          "Respiratory burst" = "#1874d2",
                                          "NOS biosynthetic process" = "#1874d2",
                                          "Platelet activation" = "#1874d2",
                                          "Cellular oxidant detoxification" = "#1874d2",
                                          "Collagen catabolic process" = "#1874d2",
                                          "Temperature homeostasis" = "#1874d2",
                                          "Chondroitin sulf./proteoglycan biosynthesis" = "#1874d2",
                                          "Cellular response to oxygen levels" = "#1874d2",
                                          "Carbohydrate catabolic process" = "#1874d2",
                                          "Glycolipid catabolic process" = "#1874d2",
                                          "Adaptive thermogenesis" = "#1874d2",
                                          "Response to angiotensin" = "#1874d2",
                                          "Cellular response to VEGF stimulus" = "#1874d2",
                                          "Regulation of coagulation" = "#1874d2",
                                          "Heparan sulf./proteoglycan biosynthesis" = "#1874d2",
                                          "Icosanoid transport" = "#1874d2",
                                          "Leukotriene metabolic process" = "#1874d2",
                                          "Ceramide catabolic process" = "#1874d2",
                                          "Glycoprotein catabolic process" = "#1874d2",
                                          "Cellular response to corticosteroid stimulus" = "#1874d2",
                                          "Superoxide metabolic process" = "#1874d2")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(color = "black", size = 14, face = "plain"), axis.text.x = element_text(color = "black", size = 14, face = "plain"), axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  ylab("Normalized enrichment score, severe vs. non-severe COVID-19") 
metabolicbarplot <- metabolicbarplot + ggtitle("Metabolic, thrombotic, and endothelial function") + theme(plot.title = element_text(size=12, face = "bold"))
metabolicbarplot <- metabolicbarplot + geom_hline(yintercept = 0, size=1) + geom_point(pch=21,size=5.5,colour="black")
metabolicbarplot





