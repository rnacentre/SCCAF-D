
####different donor
P <- readRDS("/Users/fsh/Documents/Fig_deconvolution_final/figure_10-10/pseudobulk_sample1_P.rds")
P = P[gtools::mixedsort(rownames(P)),]
P$CT = rownames(P)
P = data.table::melt(P, id.vars="CT")
colnames(P) <-c("CT","tissue","expected_values")
s = sum(P$expected_values)/1000
P$expected_values = P$expected_values/s
data <- readRDS('/Users/fsh/Documents/Fig_deconvolution_final/figure_10-10/sample4_DWLS_sc_TMM_cross.rds')
data <- arrange(data,data$tissue)
data <- data[,-4]
data <- merge(data,P)
###different study
data <- data <- readRDS('/Users/fsh/Documents/Fig_deconvolution_final/figure_10-10/figure-2024-1-15/results_marquina_muraro_dwls_SC_TMM.rds')
data$expected_values <- data$expected_values/10000
###
cairo_pdf(filename = 'marquina-muraro_sc_dwls-2024-3-22.pdf',width = 6,height = 4,family = 'Arial')# 6,4 and 8,5
p1 <- plotAllResults(data,title = 'Marquina-Sanchez vs Muraro')
# p1+scale_y_continuous(breaks = seq(0,0.45,0.2))+scale_x_continuous(breaks = seq(0,0.45,0.2))
# p1+scale_y_continuous(breaks = seq(0,0.35,0.15))+scale_x_continuous(breaks = seq(0,0.35,0.15))#sample1 vs sample4
# p1+scale_y_continuous(breaks = seq(0,0.45,0.2))+scale_x_continuous(breaks = seq(0,0.35,0.15))#marquina vs muraro
p1+scale_y_continuous(breaks = seq(0,0.4,0.2))+scale_x_continuous(breaks = seq(0,0.3,0.15))#Muraro vs Marquina-Sanchez
dev.off()
###
# plotAllResults<-function(x, n=9, ncols=3, title=''){
#   p<-ggplot(x, aes(x=observed_values, y=expected_values, color=CT)) + 
#     geom_point(alpha=.8,size=8)+
#     geom_smooth(method=lm, color='red')+ggpubr::stat_cor(method = "pearson",color="black",p.accuracy = 0.001,size=6)+
#     ggtitle(title)+ 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           text=element_text(size=25),
#           axis.title.x = element_text(size=25,colour = 'black'),axis.title.y = element_text(size=25,colour = 'black'),
#           axis.text.x = element_text(size=25,colour = 'black'),axis.text.y = element_text(size=25,colour = 'black'),
#           legend.text = element_text(size=25,colour = 'black'),legend.title = element_text(size=25,colour = 'black'),
#           legend.key = element_rect(colour = "transparent", fill = "white"))+
#     scale_color_manual(values=cols)
#   return(p)
# }


###original codes
cairo_pdf(filename = 'cross_reference Muraro vs Marquina-Sanchez.pdf',width = 6,height =4,family = 'Arial')

##
plotAllResults<-function(x, n=9, ncols=3, title=''){
p<-ggplot(x, aes(x=observed_values, y=expected_values, color=CT)) +
geom_point(alpha=.9,size=6)+
geom_smooth(method=lm, color='red')+ggpubr::stat_cor(method = "pearson",color="black",p.accuracy = 0.001,size=6)+
ggtitle(title)+
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
axis.line = element_line(colour = "black"),
text=element_text(size = 18,colour = 'black'),legend.text = element_text(size=18,colour = 'black'),axis.text.x= element_text(size=18,colour = 'black'),axis.text.y= element_text(size=18,colour = 'black'),legend.title = element_text(size=15,colour = 'black'),plot.title = element_text(size=18,colour = 'black'),
legend.key = element_rect(colour = "transparent", fill = "white"))+guides(color=guide_legend(ncol=1))+
scale_color_manual(values=cols)
return(p)
}

#####color
cols <- c("#DD6091" ,"#00B8C3", "#006091", "#FF5F1F","#FFB307" , "#AFE1AF" ,"#005C07" ,"#B4C424" ,
"#74A0FF","#0BDA51",
"#9900B3",
"#FFB6C1" ,"#FFEA00" ,"#7779BF" , "#475D4B", "#802600" ,"#9A0F03",
"#966919" ,"#D2B48C" ,"#F0E68C", "#708090" ,"#A9A9A9", "#231F20")

##
cairo_pdf(filename = 'cross_reference Marquina-Sanchez vs Muraro_rlr.pdf',width = 6,height =4,family = 'Arial')
data <- data <- readRDS('/Users/fsh/Documents/Fig_deconvolution_final/figure_10-10/figure-2024-1-15/File_decon_2024-1-18/marquina_muraro_RLR_bulk_TMM_all_cross_2024-3-18.rds')
data$expected_values <- data$expected_values/10000
p1 <- plotAllResults(data,title = 'Marquina-Sanchez vs Muraro')
p1+scale_y_continuous(breaks = seq(0,0.5,0.2))+scale_x_continuous(breaks = seq(0,0.5,0.25))
dev.off()
##
