
#heatmap 
library(ComplexHeatmap)

col_fun <- circlize::colorRamp2(
  breaks = c( 0, 1), 
  colors = c('white',"#2166AC")
)
###modify name
# colnames(dat_nnls) <- gsub('Marquina','Marquina \n -Sanchez',colnames(dat_nnls))
# rownames(dat_nnls) <- gsub('Marquina','Marquina \n -Sanchez',rownames(dat_nnls))
rownames(dat_dwls) <- gsub('Integrated','SCCAF-D',rownames(dat_dwls))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_dwls <- dat_dwls[new_order,,drop=FALSE]
###
rownames(dat_epidish) <- gsub('Integrated','SCCAF-D',rownames(dat_epidish))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_epidish <- dat_epidish[new_order,,drop=FALSE]
###
rownames(dat_fardeep) <- gsub('Integrated','SCCAF-D',rownames(dat_fardeep))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_fardeep <- dat_fardeep[new_order,,drop=FALSE]
###
rownames(dat_music) <- gsub('Integrated','SCCAF-D',rownames(dat_music))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_music <- dat_music[new_order,,drop=FALSE]
###
rownames(dat_nnls) <- gsub('Integrated','SCCAF-D',rownames(dat_nnls))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_nnls <- dat_nnls[new_order,,drop=FALSE]
###
rownames(dat_rlr) <- gsub('Integrated','SCCAF-D',rownames(dat_rlr))
new_order <- c('SCCAF-D',"Banovich","Krasnow","Lafyatis","Meyer")
dat_rlr <- dat_rlr[new_order,,drop=FALSE]
####
dat_fardeep <- dat_fardeep[-1,]
dat_epidish <- dat_epidish[-1,]
dat_music <- dat_music[-1,]
dat_nnls <- dat_nnls[-1,]
dat_rlr <- dat_rlr[-1,]
###
cairo_pdf(filename = 'lung_heatmap_all.pdf',width = 26,height = 6,family = 'Arial')
p1 <- Heatmap(t(as.matrix(dat_dwls)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_split = c(rep('a',1),rep('b',4)),column_gap = unit(8, "mm"),border = TRUE,
              column_title = 'DWLS',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_dwls)))*unit(1.15, "cm"), 
              height = nrow(t(as.matrix(dat_dwls)))*unit(1, "cm"))


p3 <- Heatmap(t(as.matrix(dat_fardeep)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_gap = unit(5, "mm"),border = TRUE,column_title = 'FARDEEP',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_fardeep)))*unit(1, "cm"), 
              height = nrow(t(as.matrix(dat_fardeep)))*unit(1, "cm"))

p2 <- Heatmap(t(as.matrix(dat_epidish)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_gap = unit(5, "mm"),border = TRUE,column_title = 'EpiDISH',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_epidish)))*unit(1, "cm"), 
              height = nrow(t(as.matrix(dat_epidish)))*unit(1, "cm"))

p4 <- Heatmap(t(as.matrix(dat_music)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_gap = unit(5, "mm"),border = TRUE,column_title = 'MuSiC',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_music)))*unit(1, "cm"), 
              height = nrow(t(as.matrix(dat_music)))*unit(1, "cm"))

p5 <- Heatmap(t(as.matrix(dat_nnls)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_gap = unit(5, "mm"),border = TRUE,column_title = 'NNLS',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_nnls)))*unit(1, "cm"), 
              height = nrow(t(as.matrix(dat_nnls)))*unit(1, "cm"))

p6 <- Heatmap(t(as.matrix(dat_rlr)),col=col_fun,name = 'Pearson',cluster_rows = F,cluster_columns = F,column_gap = unit(5, "mm"),border = TRUE,column_title = 'RLR',row_names_side = 'left',row_names_gp = gpar(fontsize=20),
              column_names_gp =gpar(fontsize=20),rect_gp = gpar(col='white',lwd=1),heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20)),column_title_gp = gpar(fontsize = 20),width = ncol(t(as.matrix(dat_rlr)))*unit(1, "cm"), 
              height = nrow(t(as.matrix(dat_rlr)))*unit(1, "cm"))

ht_list_pancreas = p1+p2+p3+p4+p5+p6
draw(ht_list, ht_gap = unit(2, "cm"))

dev.off()
####
ggplot(data)+geom_bar(aes(x=pseudobulk,y=1/RMSE,fill=group),stat = "identity",position = "dodge",width = 0.7)+scale_fill_manual(values = pal_npg()(5),#modify legend title text
                                                                                                                                name="Reference")+theme_bw()+theme(text=element_text(family='Arial'),axis.title.x = element_text(family='Arial',size=25,margin = margin(0.3,1,0.3,2,'cm'),colour = 'black'),
                                                                                                                                                                   axis.title.y=element_text(family='Arial',size=25,colour = 'black'),axis.text.x=element_text(family = 'Arial',size=25,angle = 90,colour = 'black'),axis.text.y=element_text(family='Arial',size=25,colour = 'black'),strip.text.x = element_text(family = 'Arial',size=25,colour = 'black'),
                                                                                                                                                                   legend.title=element_text(family='Arial',size=25,colour = 'black'),legend.text=element_text(family='Arial',size=25,colour = 'black'),panel.spacing.x = unit(2, "cm"),strip.background = element_rect(colour = 'white',fill = 'white'))#+facet_wrap(.~Methods,nrow = 1)

cairo_pdf(filename = 'pancreas_fardeep_dwls_rmse.pdf',width = 10,height = 6,family = 'Arial')
###

cairo_pdf(filename = 'pancreas_dwls_rmse.pdf',width = 7,height = 6,family = 'Arial')
ggplot(dwls)+geom_bar(aes(x=pseudobulk,y=1/RMSE,fill=group),stat = "identity",position = "dodge",width = 0.7,)+scale_fill_manual(values = pal_npg()(5),#modify legend title text
                                                                                                                                 name="Reference")+theme_bw()+theme(axis.title.x = element_text(size=20,margin = margin(0.2,1,0,1,'cm')),legend.text = element_text(size=20,colour = 'black'),legend.title = element_text(size=20,colour = 'black'),axis.title.y = element_text(size=20,colour = 'black'),
                                                                                                                                                                    strip.text.x = element_text(size=20,colour = 'black'),#调整分页标题
                                                                                                                                                                    axis.text.x = element_text(size=20,angle = 90,hjust = 1,vjust = 0.5,colour = 'black'),strip.background = element_rect(colour = 'black',fill = 'grey90'),
                                                                                                                                                                    axis.text.y=element_text(size=20,colour = 'black'),plot.title = element_text(hjust = 0.5,size=20,colour = 'black'))+labs(x='Pseudobulk',title = 'DWLS')#+facet_wrap(.~Methods)
dev.off()

