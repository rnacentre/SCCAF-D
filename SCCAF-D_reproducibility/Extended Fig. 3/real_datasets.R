###real_bulk
library(ggpubr)
library(ggplot2)
###
mydata$group <- factor(mydata$group,levels = c('Arunachalam','Lee','Schulte','Wilk','Reference'))
cairo_pdf(file = 'real_bulk_rmse.pdf',width = 17,height = 6,family = 'Arial')
ggplot(mydata1)+geom_bar(aes(x=Methods,y=1/RMSE,fill=group),stat = "identity",position = "dodge",width = 0.7)+scale_fill_manual(values = cols,#modify legend title text
                                                                                                                                  name="Dataset")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.title.x = element_text(family='Arial',size=18,margin = margin(0.5,1,0,1,'cm')),legend.text = element_text(family = 'Arial',size=18),legend.title = element_text(family='Arial',size=18),axis.title.y = element_text(family = 'Arial',size=18),
                                                                                                                                                                          strip.text.x = element_text(family = 'Arial',size=18,colour = 'black'),#调整分页标题
                                                                                                                                                                          axis.text.x = element_text(family='Arial',size=18,angle = 90,hjust = 1,colour = 'black',vjust = 0.5),strip.background = element_rect(colour = 'black',fill = 'grey90'),
                                                                                                                                                                          axis.text.y=element_text(family='Arial',size=18,colour = 'black'))+facet_wrap(.~study,nrow = 1)+labs(x='Method')
dev.off()

###
cairo_pdf(file = 'all_tissue_bulk_rmse.pdf',width = 17,height = 18,family = 'Arial')
cols <- c("#DD6091" ,"#006091","#6C00BF","#53A39D",  "#FFB307")
