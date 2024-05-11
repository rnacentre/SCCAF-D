cairo_pdf(file = 'integrated4_0727_sccaf.pdf',width = 15,height = 12,family = 'Arial')#original15，12
p1 <- DimPlot(integrated4_0727_intersect, reduction = "umap", label=F,group.by = "cellType",cols = cols,repel = T,label.size = 20,pt.size = 1)+theme_dr()+theme(panel.grid=element_blank(),
                                                                                                                     axis.title.x = element_text(size=40),
                                                                                                                     axis.title.y = element_text(size=40),axis.text.x = element_text(colour = 'black',size = 40),
                                                                                                                     axis.text.y = element_text(colour = 'black',size = 40),
                                                                                                                     legend.text = element_text(size=40),legend.title = element_text(size=40),
                                                                                                                     plot.title = element_text(colour = 'black',size = 40,hjust = 0.5),legend.position = 'bottom',legend.spacing.y = unit(1,'cm'))+labs(title = 'Cell type')+guides(color=guide_legend(nrow = 7,override.aes = list(size = 8)))
p2 <- DimPlot(integrated4_0727_counts0910_umapdim55_0.15_2, reduction = "umap", label=T,group.by = "cellType",cols = cols,repel = T,label.size = 12)+theme_dr()+theme(panel.grid=element_blank(),
                                                                                                                             axis.title.x = element_text(size=30),
                                                                                                                             axis.title.y = element_text(size=30),axis.text.x = element_text(colour = 'black',size = 30),
                                                                                                                             axis.text.y = element_text(colour = 'black',size = 30),
                                                                                                                             legend.text = element_text(size=30),legend.title = element_text(size=30),
                                                                                                                             plot.title = element_text(colour = 'black',size = 30,hjust = 0.5))+labs(title = 'Cell type')+guides(color = guide_legend(override.aes = list(size = 10)))
p3 <- DimPlot(ad, reduction = "umap", label=T,group.by = "Cell.type",cols = cols,repel = T,label.size = 12)+theme_dr()+theme(panel.grid=element_blank(),
                                                                                                                             axis.title.x = element_text(size=30),
                                                                                                                             axis.title.y = element_text(size=30),axis.text.x = element_text(colour = 'black',size = 30),
                                                                                                                             axis.text.y = element_text(colour = 'black',size = 30),
                                                                                                                             legend.text = element_text(size=30),legend.title = element_text(size=30),
                                                                                                                             plot.title = element_text(colour = 'black',size = 30,hjust = 0.5))+labs(title = 'Cell type')
p1
p2
p3
dev.off()
#guides(color=guide_legend(nrow = 4)
integrated4_0727_counts0910_umapdim55_0.15_2$cellType=factor(integrated4_0727_counts0910_umapdim55_0.15_2$cellType,levels = c('Hepatocyte','Cholangiocyte','HSC','B','Plasma','T','NK','Macrophage','Monocyte','Neutrophil','pDC','cDC','Endothelial','Fibroblast'))