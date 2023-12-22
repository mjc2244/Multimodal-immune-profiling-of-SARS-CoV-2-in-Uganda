#Figure 6b - GSEA Immune and metabolic pathways - PLWH vs. those without HIV

#Clear R environment
rm(list = ls())

#Important pathways 
immunometabolichiv <- read.csv(file.choose(), header=TRUE)

#Barplot
library(ggplot2)
immunometabolichivbarplot <- ggplot(immunometabolichiv, aes(x = process, y = nes, color = process)) + 
  geom_point(aes(reorder(process,nes)), size = 5.5, position = position_dodge(preserve = "single")) 
immunometabolichivbarplot <- immunometabolichivbarplot + coord_flip() + theme_bw() 
immunometabolichivbarplot <- immunometabolichivbarplot + theme(legend.position = "none") + ylim(-2.5, 2.5) +
  scale_color_manual("pathway", values = c("*Macrophage activation in immune response" = "lightsteelblue4",
                                          "*Acute inflammatory response to antigenic stimulus" = "lightsteelblue4",
                                          "*Protein O-linked mannosylation" = "lightsteelblue4",
                                          "*FcR-mediated stimulatory signaling pathway" = "lightsteelblue4",
                                          "CD40 signaling pathway" = "lightsteelblue4",
                                          "Neutrophil-mediated immunity" = "lightsteelblue4",
                                          "Inflammasome complex assembly" = "lightsteelblue4",
                                          "Complement activation" = "lightsteelblue4",
                                          "Pyroptosis" = "lightsteelblue4",
                                          "B-cell mediated immunity" = "lightsteelblue4",
                                          "TLR4 signaling pathway" = "lightsteelblue4",
                                          "Pos. reg. of NO metabolic process" = "lightsteelblue4",
                                          "Pos. reg. of ox. stress-induced cell death" = "lightsteelblue4",
                                          "Response to IFN-alpha" = "lightsteelblue4",
                                          "Response to IFN-gamma" = "lightsteelblue4",
                                          "IL-1 production" = "lightsteelblue4",
                                          "PRR signaling pathway" = "lightsteelblue4",
                                          "Response to IFN-beta" = "lightsteelblue4",
                                          "Neg. reg. of viral genome replication" = "lightsteelblue4",
                                          "Protein N-linked mannosylation" = "#00008B",
                                          "Lymphocyte chemotaxis" = "#00008B",
                                          "Hydrogen peroxide metabolic process" = "#00008B",
                                          "PDGF receptor signaling pathway" = "#00008B",
                                          "Activated T cell proliferation" = "#00008B",
                                          "Cellular response to oxidative stress" = "#00008B",
                                          "Ubiquitin-dependent ERAD pathway" = "#00008B",
                                          "Reg. of CD8+ alpha-beta T cell activation" = "#00008B",
                                          "TLR9 signaling pathway" = "#00008B",
                                          "Maintenance of protein localization in organelle" = "#00008B",
                                          "Monocyte differentiation" = "#00008B",
                                          "PERK-mediated unfolded protein response" = "#00008B",
                                          "Pos. reg. of myeloid cell differentiation" = "#00008B")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "plain"), axis.text.x = element_text(color = "black", size = 12, face = "plain"), axis.text.y = element_text(color = "black", size = 12, face = "plain")) + 
  ylab("Normalized enrichment score, PLWH vs. those without HIV") 
immunometabolichivbarplot <- immunometabolichivbarplot + ggtitle("Immune pathways") + theme(plot.title = element_text(size=12, face = "bold"))
immunometabolichivbarplot <- immunometabolichivbarplot + geom_hline(yintercept = 0, size=1) + geom_point(pch=21,size=5.5,colour="black")
immunometabolichivbarplot





