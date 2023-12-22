#Figure 2c - GSEA - Protein Synthesis, Cell Transport, Proteolysis Pathways Severe vs. non-severe COVID

#Clear R environment
rm(list = ls())

#Import enriched Protein Synthesis, Cell Transport, Proteolysis Pathways
synth <- read.csv(file.choose(), header=TRUE)

#Barplot
library(ggplot2)
synthbarplot <- ggplot(synth, aes(x = process, y = nes, color = process)) + 
  geom_point(aes(reorder(process,nes)), size = 5.5, position = position_dodge(preserve = "single")) 
synthbarplot <- synthbarplot + coord_flip() + theme_bw() 
synthbarplot <- synthbarplot + theme(legend.position = "none") + ylim(-2.5, 2.5) +
  scale_color_manual("pathway", values = c("Regulation of receptor-mediated endocytosis" = "#1874d2",
                                          "Vacuolar acidification" = "#1874d2",
                                          "Ubiquitin-dependent ERAD pathway" = "#1874d2",
                                          "Extrinsic apoptotic signaling via DDR" = "#1874d2",
                                          "Pos. regulation of syncytium formation by PM fusion" = "#1874d2",
                                          "Cell adhesion mediated by integrin" = "#1874d2",
                                          "Ribosomal small subunit biogenesis" = "#ed7953",
                                          "tRNA modification" = "#ed7953",
                                          "Ribonucleoprotein complex biogenesis" = "#ed7953",
                                          "Mitochondrial translation" = "#ed7953",
                                          "tRNA metabolic process" = "#ed7953",
                                          "Ribosomal large subunit biogenesis" = "#ed7953",
                                          "RNA methylation" = "#ed7953",
                                          "Autophagosome organization" = "#1874d2",
                                          "ER to cytosol transport" = "#1874d2",
                                          "Pos. regulation of proteolysis" = "#1874d2",
                                          "Regulation of lysosomal lumen pH" = "#1874d2",
                                          "Viral release from host cell" = "#1874d2",
                                          "Transcription by RNA polymerase III" = "#ed7953",
                                          "Cytoplasmic translational initiation" = "#ed7953",
                                          "Ribosomal subunit export from nucleus" = "#ed7953",
                                          "Ribosomal large subunit assembly" = "#ed7953",
                                          "Pos. regulation of autophagy" = "#1874d2",
                                          "Virion assembly" = "#1874d2",
                                          "rRNA metabolic process" = "#ed7953",
                                          "Mitochondrial RNA metabolic process" = "#ed7953")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(color = "black", size = 14, face = "plain"), axis.text.x = element_text(color = "black", size = 14, face = "plain"), axis.text.y = element_text(color = "black", size = 14, face = "plain")) + 
  ylab("Normalized enrichment score, severe vs. non-severe COVID-19") 
synthbarplot <- synthbarplot + ggtitle("Cell entry, transport, proteolysis, and protein synthesis") + theme(plot.title = element_text(size=12, face = "bold"))
synthbarplot <- synthbarplot + geom_hline(yintercept = 0, size=1) + geom_point(pch=21,size=5.5,colour="black")
synthbarplot





