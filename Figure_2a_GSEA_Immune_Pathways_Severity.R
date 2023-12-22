#Figure 2a - GSEA - Immune Pathways Severe vs. non-severe COVID

#Clear R environment
rm(list = ls())

#Import enriched immune pathways 
immune <- read.csv(file.choose(), header=TRUE)

#Barplot
library(ggplot2)
immunebarplot <- ggplot(immune, aes(x = process, y = nes, color = process)) + 
  geom_point(aes(reorder(process,nes)), size = 5.5, position = position_dodge(preserve = "single")) 
immunebarplot <- immunebarplot + coord_flip() + theme_bw() 
immunebarplot <- immunebarplot + theme(legend.position = "none") + ylim(-2, 2) +
  scale_color_manual("pathway", values = c("Phagocytosis" = "#1874d2",
                                          "Pos. regulation of PRR signaling" = "#1874d2",
                                          "Regulation of TLR4 signaling" = "#1874d2",
                                          "Neutrophil activation in immune response" = "#1874d2",
                                          "Pos. regulation of IL-6 production" = "#1874d2",
                                          "Neutrophil migration" = "#1874d2",
                                          "MYD88-dependent TLR signaling" = "#1874d2",
                                          "Innate immune respose in mucosa" = "#1874d2",
                                          "Complement activation" = "#1874d2",
                                          "Humoral immune response via circulating Ig" = "#1874d2",
                                          "Regulation of Th1-type immune response" = "#1874d2",
                                          "Macrophage cytokine production" = "#1874d2",
                                          "Th1 cell differentiation" = "#1874d2",
                                          "TNF superfamily cytokine production" = "#1874d2",
                                          "IL-1-mediated signaling" = "#1874d2",
                                          "Granulocyte activation" = "#1874d2",
                                          "Pos. regulation of IL-8 production" = "#1874d2",
                                          "T cell activation in immune response" = "#1874d2",
                                          "Pos. regulation of NF-kB activity" = "#1874d2",
                                          "Regulation of CD4+ T cell differentiation" = "#1874d2",
                                          "Inflammasome complex assembly" = "#1874d2",
                                          "FcR-mediated stimulatory signaling" = "#1874d2",
                                          "Receptor signaling pathway via STAT" = "#1874d2",
                                          "Phagocytosis" = "#1874d2",
                                          "IFN-mediated signaling" = "#1874d2",
                                          "Integrated stress response signaling" = "#1874d2",
                                          "Myeloid leukocyte cytokine production" = "#1874d2",
                                          "Dendritic cell differentiation" = "#1874d2")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(color = "black", size = 14, face = "plain"), axis.text.x = element_text(color = "black", size = 14, face = "plain"), axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  ylab("Normalized enrichment score, severe vs. non-severe COVID-19") 
immunebarplot <- immunebarplot + ggtitle("Innate and adaptive immune functions") + theme(plot.title = element_text(size=12, face = "bold"))
immunebarplot <- immunebarplot + geom_hline(yintercept = 0, size=1) + geom_point(pch=21,size=5.5,colour="black")
immunebarplot





