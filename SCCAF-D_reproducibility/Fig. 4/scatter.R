######
cairo_pdf(filename = 'newman_arun_r.pdf',width = 17,height = 10.67,family = 'Arial')#finotello==16,10.67
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(newman_arun, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
                ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
                                                 axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                 axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))
p1+facetted_pos_scales(y=list(CT == "B" ~ scale_y_continuous(limits=c(-0.02,0.12),breaks=seq(-0.02,0.12,0.07)),
                              CT == "CD4T" ~ scale_y_continuous(limits=c(0, 0.22),breaks=seq(0,0.22,0.1)),
                              CT == "CD8T" ~ scale_y_continuous(limits=c(0,0.18),breaks=seq(0,0.18,0.07)),
                              # CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
                              CT == "Monocyte" ~ scale_y_continuous(limits=c(0.04,0.18),breaks=seq(0.04,0.18,0.06)),
                              # CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.02,0.07),breaks=seq(0.02,0.07,0.02))
                              CT == "NK" ~ scale_y_continuous(limits=c(0.01,0.12),breaks=seq(0.01,0.12,0.05))
),
x = list(
  CT == "B" ~ scale_x_continuous(limits=c(0.02,0.24),breaks=seq(0.02,0.24,0.1)),
  CT == "CD4T" ~ scale_x_continuous(limits=c(0.1, 0.47),breaks=seq(0.1,0.47,0.15)),
  CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.11),breaks=seq(0,0.11,0.05)),
  # CT == "cDC" ~ scale_x_continuous(limits=c(0,0.025),breaks=seq(0,0.025,0.01)),
  # CT == "Monocyte" ~ scale_x_continuous(limits=c(0.25,0.42),breaks=seq(0.25,0.42,0.1)),
  CT == "Neutrophil" ~ scale_x_continuous(limits=c(0,0.12),breaks=seq(0,0.12,0.06)),
  CT == "NK" ~ scale_x_continuous(limits=c(0,0.0345),breaks=seq(0,0.0345,0.015))
))
# p1
dev.off()
####
cairo_pdf(filename = 'newman_lee_r.pdf',width = 17,height = 10.67,family = 'Arial')#finotello==12
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(newman_lee, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
                ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
                                                 axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                 axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))
p1+facetted_pos_scales(y=list(
  CT == "B" ~ scale_y_continuous(limits=c(0,0.125),breaks=seq(0,0.125,0.06)),
  CT == "CD4T" ~ scale_y_continuous(limits=c(0.03, 0.25),breaks=seq(0.03,0.25,0.1)),
  CT == "CD8T" ~ scale_y_continuous(limits=c(0,0.25),breaks=seq(0,0.25,0.1)),
  # CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
  CT == "Monocyte" ~ scale_y_continuous(limits=c(0.04,0.22),breaks=seq(0.04,0.22,0.08)),
  CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.4,0.9),breaks=seq(0.4,0.9,0.2)),
  CT == "NK" ~ scale_y_continuous(limits=c(0.01,0.13),breaks=seq(0.01,0.13,0.05))
),
x = list(
  CT == "B" ~ scale_x_continuous(limits=c(0.02,0.3),breaks=seq(0.02,0.3,0.13)),
  # CT == "CD4T" ~ scale_x_continuous(limits=c(0.1, 0.47),breaks=seq(0.1,0.47,0.15)),
  CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.008),breaks=seq(0,0.008,0.004)),
  # CT == "cDC" ~ scale_x_continuous(limits=c(0,0.025),breaks=seq(0,0.025,0.01)),
  CT == "Monocyte" ~ scale_x_continuous(limits=c(0,0.33),breaks=seq(0,0.33,0.15)),
  CT == "Neutrophil" ~ scale_x_continuous(limits=c(0.3,0.9),breaks=seq(0.3,0.9,0.3)),
  CT == "NK" ~ scale_x_continuous(limits=c(0,0.045),breaks=seq(0,0.045,0.02))
))
# p1
dev.off()
###newman_schulte
cairo_pdf(filename = 'newman_schulte_r.pdf',width = 17,height = 10.67,family = 'Arial')#finotello==12
####
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(newman_schulte, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
                ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
                                                 axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                 axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))
p1+facetted_pos_scales(y=list(
  # CT == "B" ~ scale_y_continuous(limits=c(0,0.125),breaks=seq(0,0.125,0.06)),
  CT == "CD4T" ~ scale_y_continuous(limits=c(0.03, 0.25),breaks=seq(0.03,0.25,0.1)),
  CT == "CD8T" ~ scale_y_continuous(limits=c(0,0.15),breaks=seq(0,0.15,0.06)),
  # CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
  # CT == "Monocyte" ~ scale_y_continuous(limits=c(0.04,0.22),breaks=seq(0.04,0.22,0.06)),
  CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.4,0.9),breaks=seq(0.4,0.9,0.2)),
  CT == "NK" ~ scale_y_continuous(limits=c(0.01,0.13),breaks=seq(0.01,0.13,0.05))
),
x = list(
  CT == "B" ~ scale_x_continuous(limits=c(0,0.75),breaks=seq(0,0.75,0.3)),
  CT == "CD4T" ~ scale_x_continuous(limits=c(0.06, 0.3),breaks=seq(0.06,0.3,0.1)),
  CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.1),breaks=seq(0,0.1,0.05)),
  # CT == "cDC" ~ scale_x_continuous(limits=c(0,0.025),breaks=seq(0,0.025,0.01)),
  CT == "Monocyte" ~ scale_x_continuous(limits=c(0.03,0.5),breaks=seq(0.03,0.5,0.2)),
  CT == "Neutrophil" ~ scale_x_continuous(limits=c(0.15,0.75),breaks=seq(0.15,0.75,0.25)),
  CT == "NK" ~ scale_x_continuous(limits=c(0,0.047),breaks=seq(0,0.047,0.02))
))
# p1
dev.off()
###
cairo_pdf(filename = 'newman_wilk_r.pdf',width = 17,height = 10.67,family = 'Arial')#finotello==12
####
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(newman_wilk, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
                ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
                                                 axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                 axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))
p1+facetted_pos_scales(y=list(
  CT == "B" ~ scale_y_continuous(limits=c(0.01,0.125),breaks=seq(0.01,0.125,0.05)),
  CT == "CD4T" ~ scale_y_continuous(limits=c(0.02, 0.25),breaks=seq(0.02,0.25,0.1)),
  CT == "CD8T" ~ scale_y_continuous(limits=c(0,0.16),breaks=seq(0,0.16,0.08)),
  # CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
  # CT == "Monocyte" ~ scale_y_continuous(limits=c(0.15,0.42),breaks=seq(0.15,0.42,0.1)),
  CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.3,0.9),breaks=seq(0.3,0.85,0.25)),
  CT == "NK" ~ scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,0.06))),
  x = list(
    CT == "B" ~ scale_x_continuous(limits=c(0,0.8),breaks=seq(0,0.8,0.4)),
    # CT == "CD4T" ~ scale_x_continuous(limits=c(0.06, 0.28),breaks=seq(0.06,0.28,0.1)),
    CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.06),breaks=seq(0,0.06,0.03)),
    # CT == "cDC" ~ scale_x_continuous(limits=c(0,0.025),breaks=seq(0,0.025,0.01)),
    CT == "Monocyte" ~ scale_x_continuous(limits=c(0,0.02),breaks=seq(0,0.02,0.01)),
    CT == "Neutrophil" ~ scale_x_continuous(limits=c(0.2,0.9),breaks=seq(0.2,0.9,0.3)),
    CT == "NK" ~ scale_x_continuous(limits=c(0,0.022),breaks=seq(0,0.022,0.01))
  ))
# p1
dev.off()
###
cairo_pdf(filename = 'newman_integrated0908mean5_r.pdf',width = 17,height = 10.67,family = 'Arial')#finotello==12
####
#facet_grid(drv ~ cyl, scales = "free")
p1 <- ggscatter(newman_integrated0908mean5, x = "observed_values", y = "expected_values",size = 7,ylab='',xlab = '',color = "#006091",
                ggtheme =  theme_classic())+geom_smooth(method = lm,colour='black')+ggpubr::stat_cor(color="black",p.accuracy = 0.001,size=12)+
  facet_wrap(.~CT, scales = "free",ncol=3)+theme(axis.line.x = element_line(colour = "black",size = 1),axis.line.y = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(family = 'Arial',size = 35,colour = 'black'),
                                                 axis.text.y = element_text(family = 'Arial',size = 35,colour = 'black'),axis.title.x = element_text(family = 'Arial',size = 30),
                                                 axis.title.y = element_text(family = 'Arial',size = 30),strip.text.x = element_text(family = 'Arial',size=40,colour = 'black'),strip.background = element_rect(colour = 'white',fill = 'white'))
p1+facetted_pos_scales(y=list(CT == "B" ~ scale_y_continuous(limits=c(0,0.125),breaks=seq(0,0.125,0.06)),
                              CT == "CD4T" ~ scale_y_continuous(limits=c(0, 0.25),breaks=seq(0.,0.25,0.1)),
                              CT == "CD8T" ~ scale_y_continuous(limits=c(0,0.2),breaks=seq(0,0.2,0.1)),
                              # CT == "cDC" ~ scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0.01,0.06,0.02)),
                              CT == "Monocyte" ~ scale_y_continuous(limits=c(0.04,0.2),breaks=seq(0.04,0.2,0.08)),
                              CT == "Neutrophil" ~ scale_y_continuous(limits=c(0.3,0.85),breaks=seq(0.3,0.85,0.2)),
                              CT == "NK" ~ scale_y_continuous(limits=c(0.015,0.125),breaks=seq(0.015,0.125,0.05))
),
x = list(
  CT == "B" ~ scale_x_continuous(limits=c(0.01,0.13),breaks=seq(0.01,0.13,0.05)),
  CT == "CD4T" ~ scale_x_continuous(limits=c(0.08, 0.36),breaks=seq(0.08,0.36,0.12)),
  CT == "CD8T" ~ scale_x_continuous(limits=c(0, 0.066),breaks=seq(0,0.066,0.03)),
  # CT == "cDC" ~ scale_x_continuous(limits=c(0.018,0.036),breaks=seq(0.018,0.036,0.009)),
  CT == "Monocyte" ~ scale_x_continuous(limits=c(0.1,0.31),breaks=seq(0.1,0.31,0.1)),
  CT == "Neutrophil" ~ scale_x_continuous(limits=c(0.3,0.72),breaks=seq(0.3,0.72,0.2)),
  CT == "NK" ~ scale_x_continuous(limits=c(0.005,0.082),breaks=seq(0.005,0.082,0.03))
))
# p1
dev.off()