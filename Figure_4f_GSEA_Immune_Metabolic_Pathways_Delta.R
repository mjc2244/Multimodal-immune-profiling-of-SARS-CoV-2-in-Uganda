#Figure 4f - GSEA Immune and metabolic pathways - Delta vs. non-Delta COVID

#Clear R environment
rm(list = ls())

#Import enriched Delta vs. non-Delta pathways 
immunometabolicdelta <- read.csv(file.choose(), header=TRUE)

#Barplot
library(ggplot2)
immunometabolicdeltabarplot <- ggplot(immunometabolicdelta, aes(x = process, y = nes, color = process)) + 
  geom_point(aes(reorder(process,nes)), size = 5.5, position = position_dodge(preserve = "single")) 
immunometabolicdeltabarplot <- immunometabolicdeltabarplot + coord_flip() + theme_bw() 
immunometabolicdeltabarplot <- immunometabolicdeltabarplot + theme(legend.position = "none") + ylim(-2.5, 2.5) +
  scale_color_manual("pathway", values = c("Inflammasome complex assembly" = "#008B00",
                                           "Negative regulation of coagulation" = "#008B00",
                                           "Positive regulation of IL-1 production" = "#008B00",
                                           "Phagocytosis" = "#008B00",
                                          "Cytokine production in inflammatory response" = "#008B00",
                                          "Complement activation" = "#008B00",
                                          "Leukotriene metabolic process" = "#008B00",
                                          "Superoxide metabolic process" = "#008B00",
                                          "Neg. reg. of extrinsic apoptotic signaling via DRR"= "#008B00",
                                          "Receptor-mediated endocytosis" = "#008B00",
                                          "Myeloid cell activation in immune response" = "#008B00",
                                          "Myeloid leukocyte-mediated immunity" = "#008B00",
                                          "Reg. of NLRP3 inflammasome complex assembly" = "#008B00",
                                          "Autophagosome organization" = "#008B00",
                                          "Regulation of leukocyte degranulation" = "#008B00",
                                          "Cytoplasmic translation" = "#8B00B8",
                                          "Ribosome biogenesis" = "#8B00B8",
                                          "tRNA modification" = "#8B00B8",
                                          "Mitochondrial gene expression" = "#8B00B8",
                                          "RNA methylation" = "#8B00B8",
                                          "Ribosomal small subunit assembly" = "#8B00B8",
                                          "Ribosomal large subunit assembly" = "#8B00B8",
                                          "Positive T cell selection" = "#8B00B8",
                                          "Peroxisome organization" = "#8B00B8",
                                          "Ribonucleoprotein complex biogenesis" = "#8B00B8",
                                          "Pos. reg. of syncytium formation by PM fusion" = "#008B00",
                                          "Histone H3 K14 acetylation" = "#8B00B8",
                                          "Regulation of translational fidelity" = "#8B00B8")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "plain"), axis.text.x = element_text(color = "black", size = 12, face = "plain"), axis.text.y = element_text(color = "black", size = 12, face = "plain")) + 
  ylab("Normalized enrichment score, Delta vs. non-Delta phase COVID-19") 
immunometabolicdeltabarplot <- immunometabolicdeltabarplot + ggtitle("Immune pathways") + theme(plot.title = element_text(size=12, face = "bold"))
immunometabolicdeltabarplot <- immunometabolicdeltabarplot + geom_hline(yintercept = 0, size=1) + geom_point(pch=21,size=5.5,colour="black")
immunometabolicdeltabarplot





