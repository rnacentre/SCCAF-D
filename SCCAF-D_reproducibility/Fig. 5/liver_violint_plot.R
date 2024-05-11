#violin_plot_liver
options(repr.plot.width=7,repr.plot.height=7)
library(ggpubr)
my_comparisons <- list(c('1','0'),c('2','0'),c('3','0'),c('4','0'))
###
plot_violint <- function(data,n=9,ncols=3){
  parray <- ggplot()
  length(parray) <- n
  i=1
  for(ts in unique(fib$CT)){
    y <- fib[fib$CT==ts,]
    parray[[i]] <- ggviolin(y,x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),add=c('jitter'),ylim=c(0,0.12),
             xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.08,0.09,0.1,0.116),
                                                     aes(label = after_stat(p.signif)))+
      theme(text=element_text(family='Arial'),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 20,family = 'Arial'),axis.title.y = element_text(size=20,family = 'Arial'),axis.text.x = element_text(size = 12,family='Arial'),axis.text.y = element_text(size = 12,family='Arial'),legend.text = element_text(size=12,family = 'Arial'),
            legend.title = element_text(size=12,family = 'Arial'))+labs(fill='Fibrosis stage')+
      ggtitle(ts)
    i=i+1
  }
  g<-DeconRNASeq::multiplot(plotlist = parray, cols=ncols)
  return(g)
}

#####
ggviolin(fib[fib$CT=='cDC',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),add=c('jitter'),ylim=c(0,0.02),
         xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.001,0.003,0.005,0.015),
                                                 aes(label = after_stat(p.signif)))+
  theme(text=element_text(family='Arial'),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 20,family = 'Arial'),axis.title.y = element_text(size=20,family = 'Arial'),axis.text.x = element_text(size = 12,family='Arial'),axis.text.y = element_text(size = 12,family='Arial'),legend.text = element_text(size=12,family = 'Arial'),
        legend.title = element_text(size=12,family = 'Arial'))+labs(fill='Fibrosis stage')+
  ggtitle('cDC')




#####
B           cDC Cholangiocyte   endothelial    fibroblast    Hepatocyte           HSC    macrophage 
87            87            87            87            87            87            87            87 
monocyte    neutrophil            NK           pDC        Plasma             T 
87            87            87            87            87            87 
###
jpeg(filename = 'liver_T.jpeg',width = 30,height = 30,units = 'cm',res = 300)
pdf(file = 'liver_T.pdf',width = 30,height = 30)
####CDC
ggviolin(fib[fib$CT=='T',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
         add=c('jitter'),ylim=c(0,0.15),
         xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.09,0.1,0.11,0.12),
                                                 aes(label = after_stat(p.signif)))+
  theme(text=element_text(family='Arial',size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 20,family = 'Arial'),
        axis.title.y = element_text(size=20,family = 'Arial'),axis.text.x = element_text(size = 30,family='Arial'),axis.text.y = element_text(size = 30,family='Arial'),legend.text = element_text(size=30,family = 'Arial'),
        legend.title = element_text(size=30,family = 'Arial'))+labs(fill='Fibrosis stages')+
  ggtitle('T')
dev.off()

#####pdf
cairo_pdf(filename = 'liver_T.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='T',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
         add=c('jitter'),ylim=c(0,0.045),
         xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.01,0.02,0.03,0.04),
                                                 aes(label = after_stat(p.signif)))+
  theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
        legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
  ggtitle('T')
dev.off()



####
cairo_pdf(filename = 'liver_Plasma.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='Plasma',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.045),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.01,0.02,0.03,0.04),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Plasma')
dev.off()


###
cairo_pdf(filename = 'liver_pDC.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='pDC',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.045),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.01,0.02,0.03,0.04),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('pDC')
dev.off()
###
cairo_pdf(filename = 'liver_HSC.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='HSC',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.45),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.1,0.2,0.3,0.4),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('HSC')
dev.off()
###
cairo_pdf(filename = 'liver_Hepatocyte.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='Hepatocyte',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,2.2),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(1.0,1.3,1.6,1.9),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Hepatocyte')
dev.off()
###
cairo_pdf(filename = 'liver_fibroblast.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='fibroblast',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.06),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.02,0.03,0.04,0.05),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Fibroblast')
dev.off()
###
cairo_pdf(filename = 'liver_endothelial.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='endothelial',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.2),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.09,0.12,0.15,0.18),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Endothelial')
dev.off()
###
cairo_pdf(filename = 'liver_Cholangiocyte.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='Cholangiocyte',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.4),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.15,0.2,0.30,0.3),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Cholangiocyte')
dev.off()

##
cairo_pdf(filename = 'liver_pDC.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='pDC',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.02),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.01,0.013,0.016,0.019),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('pDC')
dev.off()

###
cairo_pdf(filename = 'liver_Cholangiocyte.pdf',width = 7,height = 6,family = 'Arial')
ggviolin(fib[fib$CT=='Cholangiocyte',],x="fib",y="observed_values",fill="fib",palette =cols,add.params = list(fill="white"),
add=c('jitter'),ylim=c(0,0.4),
xlab = '',ylab = '')+stat_compare_means(comparisons = my_comparisons,label.y = c(0.15,0.2,0.25,0.3),
aes(label = after_stat(p.signif)))+
theme(text=element_text(size=30),plot.margin = margin(3,15,0,0),axis.title.x = element_text(size = 30),
axis.title.y = element_text(size=30),axis.text.x = element_text(size = 30,colour = 'black'),axis.text.y = element_text(size = 30,colour = 'black'),legend.text = element_text(size=30),
legend.title = element_text(size=30))+labs(fill='Fibrosis stage')+ylab('Cell proportion')+
ggtitle('Cholangiocyte')
dev.off()

