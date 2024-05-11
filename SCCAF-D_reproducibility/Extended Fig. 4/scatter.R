scaleFUN <- function(x) sprintf("%.3f", x)



plotTopResults<-function(x, n=9, ncols=3){
  parray <- ggplot()
  length(parray) <- n
  i = 1
  for(ts in head(unique(x$CT), n=n)){
    y = x[x$CT==ts,]
    parray[[i]]<-ggplot(y, aes(x=observed_values, y=expected_values, color=CT)) + 
      geom_point(size=7)+
      geom_smooth(method = lm,colour='black')+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            text = element_text(size = 35,family = 'Arial'),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
            axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
            axis.title.y = element_text(family = 'Arial',size = 30),plot.margin = margin(3,15,0,0),axis.line.x =element_line(size = 1),
            axis.line.y = element_line(size = 1),
            legend.position = "none")+xlab(label = '')+ylab(label = '')+
      scale_color_manual(values='#1f77b4')+labs(title = ts)+scale_x_continuous(labels=scaleFUN)
    i=i+1
  }
  g<-DeconRNASeq::multiplot(plotlist = parray, cols=ncols)
  return(g)
}
###
jpeg(filename = 'finotello_arun.jpeg',width = 65,height = 60,res = 300,units = 'cm')

###] B          CD4T       CD8T       Monocyte   Neutrophil NK    
parray <- ggplot()
p6<- ggplot(newman_lee[newman_lee$CT=='NK',], aes(x=observed_values, y=expected_values, color=CT)) + 
  geom_point(size=3)+
  geom_smooth(method = lm,colour='black')+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 30,family = 'Arial'),axis.text.x = element_text(family = 'Arial',size = 25),
        axis.text.y = element_text(family = 'Arial',size = 25),axis.title.x = element_text(family = 'Arial',size = 30),
        axis.title.y = element_text(family = 'Arial',size = 30),plot.margin = margin(1,20,0,0),
        legend.position = "none")+xlab(label = '')+ylab(label = '')+
  scale_color_manual(values='#1f77b4')+labs(title = 'NK')#+scale_x_continuous()

####
cairo_pdf(filename = 'finotello_lee_r.pdf',width = 17,height = 16,family = 'Arial')#finotello==12
#plotTopResults(finotello_arun)

####
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(finotello_lee, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
          ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
      axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                                               axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))


p1+facetted_pos_scales(y=list(CT == "B" ~ scale_y_continuous(limits=c(0.03,0.09),breaks=seq(0.03,0.09,0.03)),
                              CT == "CD4T" ~ scale_y_continuous(limits=c(0, 0.5),breaks=seq(0,0.5,0.2)),
                              CT == "CD8T" ~ scale_y_continuous(limits=c(0.15,0.5),breaks=seq(0.15,0.5,0.15)),
                              CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
                              CT == "Monocyte" ~ scale_y_continuous(limits=c(0.15,0.42),breaks=seq(0.15,0.42,0.1)),
                              CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.02,0.07),breaks=seq(0.02,0.07,0.02))
                              # CT == "NK" ~ scale_y_continuous(limits=c(0.03,0.17),breaks=seq(0.03,0.17,0.06))),
),
  x = list(
           # CT == "B" ~ scale_x_continuous(limits=c(0.03,0.082),breaks=seq(0.03,0.082,0.02)),
           # CT == "CD4T" ~ scale_x_continuous(limits=c(0.13, 0.38),breaks=seq(0.13,0.38,0.1)),
           CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.18),breaks=seq(0,0.18,0.09)),
           CT == "cDC" ~ scale_x_continuous(limits=c(0,0.025),breaks=seq(0,0.025,0.01)),
           CT == "Monocyte" ~ scale_x_continuous(limits=c(0.25,0.42),breaks=seq(0.25,0.42,0.1)),
           CT == "Neutrophil" ~ scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.025)),
           # CT == "NK" ~ scale_x_continuous(limits=c(0.04,0.23),breaks=seq(0.04,0.24,0.1))
           ))
# p1
dev.off()
####全部细胞的线性图
cairo_pdf(filename = 'monaco_array_plotall.pdf',width = 17,height = 10,family = 'Arial')
ggscatter(mydata, x = "observed_values", y = "expected_values",size = 5,ylab='',xlab = '',color = 'CT',palette = cols,
          ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=9,)+facet_wrap(.~reference,ncol=3,scales = 'free')+
  theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
        axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
        axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=30,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'),legend.text = element_text(size=30),legend.title = element_text(size=30))



dev.off()

###
monaco_rnaseq_arun$reference <- rep('Arunachalam',72)
monaco_rnaseq_lee$reference <- rep('Lee',72)
monaco_rnaseq_schulte$reference <- rep('Schulte',72)
monaco_rnaseq_wilk$reference <- rep('Wilk',72)
monaco_rnaseq_integrated0908mean5$reference <- rep('Reference',72)

mydata$reference <- factor(mydata$reference,levels = c('Arunachalam','Lee','Schulte','Wilk','Reference'))
###
mydata <- rbind(rbind(rbind(rbind(finotello_arun,finotello_integrated0908mean5),finotello_lee),finotello_schulte),finotello_wilk)
mydata <- rbind(rbind(rbind(rbind(monaco_array_arun,monaco_array_integrated0908mean5),monaco_array_lee),monaco_array_schulte),monaco_array_wilk)
mydata <- rbind(rbind(rbind(rbind(monaco_rnaseq_arun,monaco_rnaseq_integrated0908mean5),monaco_rnaseq_lee),monaco_rnaseq_schulte),monaco_rnaseq_wilk)
mydata <- rbind(rbind(rbind(rbind(newman_arun,newman_integrated0908mean5),newman_lee),newman_schulte),newman_wilk)


###test with r
ggscatter(newman_integrated0908mean5,x = "observed_values", y = "expected_values",size = 10,ylab='',xlab = '',color = "#006091",
          ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=6)+facet_wrap(.~CT, scales = "free",nrow=3)+
  theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 25,colour = 'black'),
        axis.text.y = element_text(family = 'Arial',size = 25,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
        axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))





#####original
cairo_pdf(filename = 'monaco_array_plotall.pdf',width = 17,height = 10,family = 'Arial')
p1 <- ggscatter(mydata, x = "observed_values", y = "expected_values",size = 5,ylab='',xlab = '',color = 'CT',palette = cols,
ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=9,)+facet_wrap(.~reference,ncol=3,scales = 'free')+
theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=30,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'),legend.text = element_text(size=30),legend.title = element_text(size=30))
##adjust axis
p1+facetted_pos_scales(y=list(reference=='Arunachalam'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Lee'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Schulte'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Wilk'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Reference'~scale_y_continuous(breaks = seq(0,0.5,0.2))
),
x=list(reference=='Arunachalam'~scale_x_continuous(breaks = seq(0,0.6,0.2)),
reference=='Lee'~scale_x_continuous(breaks = seq(0,0.6,0.3)),
reference=='Schulte'~scale_x_continuous(breaks = seq(0,0.5,0.25)),
reference=='Wilk'~scale_x_continuous(breaks = seq(0,0.6,0.3)),
reference=='Reference'~scale_x_continuous(breaks = seq(0,0.5,0.25))
))
dev.off()
####
cairo_pdf(filename = 'newman_plotall.pdf',width = 17,height = 10,family = 'Arial')
p1 <- ggscatter(mydata, x = "observed_values", y = "expected_values",size = 5,ylab='',xlab = '',color = 'CT',palette = cols,
ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=9,)+facet_wrap(.~reference,ncol=3,scales = 'free')+
theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=30,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'),legend.text = element_text(size=30),legend.title = element_text(size=30))
##adjust axis
p1+facetted_pos_scales(y=list(reference=='Arunachalam'~scale_y_continuous(breaks = seq(0,0.9,0.4)),
reference=='Lee'~scale_y_continuous(breaks = seq(0,0.8,0.4)),
reference=='Schulte'~scale_y_continuous(breaks = seq(0,0.8,0.4)),
reference=='Wilk'~scale_y_continuous(breaks = seq(0,0.8,0.4)),
reference=='Reference'~scale_y_continuous(breaks = seq(0,0.8,0.4))
),
x=list(reference=='Arunachalam'~scale_x_continuous(breaks = seq(0,0.8,0.4)),
reference=='Lee'~scale_x_continuous(breaks = seq(0,0.8,0.4)),
reference=='Schulte'~scale_x_continuous(breaks = seq(0,0.7,0.3)),
reference=='Wilk'~scale_x_continuous(breaks = seq(0,0.8,0.4)),
reference=='Reference'~scale_x_continuous(breaks = seq(0,0.7,0.3))
))
dev.off()
###
cairo_pdf(filename = 'finotello_plotall.pdf',width = 17,height = 10,family = 'Arial')
p1 <- ggscatter(mydata, x = "observed_values", y = "expected_values",size = 5,ylab='',xlab = '',color = 'CT',palette = cols,
ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=9,)+facet_wrap(.~reference,ncol=3,scales = 'free')+
theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=30,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'),legend.text = element_text(size=30),legend.title = element_text(size=30))
##adjust axis
p1+facetted_pos_scales(y=list(reference=='Arunachalam'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Lee'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Schulte'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Wilk'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Reference'~scale_y_continuous(breaks = seq(0,0.5,0.2))
),
x=list(reference=='Arunachalam'~scale_x_continuous(breaks = seq(0,0.6,0.25)),
reference=='Lee'~scale_x_continuous(breaks = seq(0,0.5,0.25)),
reference=='Schulte'~scale_x_continuous(breaks = seq(0,0.4,0.2)),
reference=='Wilk'~scale_x_continuous(breaks = seq(0,0.5,0.2)),
reference=='Reference'~scale_x_continuous(breaks = seq(0,0.5,0.2))
))
dev.off()
####
cairo_pdf(filename = 'monaco_rnaseq_plotall.pdf',width = 17,height = 10,family = 'Arial')
p1 <- ggscatter(mydata, x = "observed_values", y = "expected_values",size = 5,ylab='',xlab = '',color = 'CT',palette = cols,
ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=9,)+facet_wrap(.~reference,ncol=3,scales = 'free')+
theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 30,colour = 'black'),
axis.text.y = element_text(family = 'Arial',size = 30,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=30,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'),legend.text = element_text(size=30),legend.title = element_text(size=30))
##adjust axis
p1+facetted_pos_scales(y=list(reference=='Arunachalam'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Lee'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Schulte'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Wilk'~scale_y_continuous(breaks = seq(0,0.5,0.2)),
reference=='Reference'~scale_y_continuous(breaks = seq(0,0.5,0.2))
),
x=list(reference=='Arunachalam'~scale_x_continuous(breaks = seq(0,0.8,0.3)),
reference=='Lee'~scale_x_continuous(breaks = seq(0,0.6,0.25)),
reference=='Schulte'~scale_x_continuous(breaks = seq(0,0.6,0.3)),
reference=='Wilk'~scale_x_continuous(breaks = seq(0,0.5,0.25)),
reference=='Reference'~scale_x_continuous(breaks = seq(0,0.7,0.3))
))
dev.off()

