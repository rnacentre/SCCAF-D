
  for (i in 1:nrow(data1@meta.data)){
    data1@meta.data$Donor[i] <- strsplit(data1@meta.data$Donor[i],split = "_")[[1]][2]
  }

###
data1$disease=ifelse(str_detect(data1$Donor,'cirrhotic'),'Diseased','Healthy')

###

  

###
cluster=c(
  'B'='B',
  'cDC'='cDC',
  'Cholangiocyte'='Cholangiocyte',
  'endothelial'='Endothelial',
  'fibroblast'='Fibroblast',
  'Hepatocyte'='Hepatocyte',
  'HSC'='HSC',
  'macrophage'='Macrophage',
  'monocyte'='Monocyte',
  'neutrophil'='Neutrophil',
  'NK'='NK',
  'pDC'='pDC',
  'Plasma'='Plasma',
  'T'='T'
)

integrated4_0727_intersect <- RenameIdents(integrated4_0727_intersect, cluster) 
integrated4_0727_intersect[['cellType']] = unname(cluster[integrated4_0727_intersect@meta.data$cellType])

jjDotPlot(data,gene = marker,id='cellType',ytree = FALSE,cluster.order = rev(levels(data$cellType)),legend.position = 'top')+ylab('cell type')+xlab('marker gene')+theme(text = element_text(size=18,family = 'Arial'))+ theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(vjust = 0.5,size = 18,family = 'Arial'),axis.title.x = element_text(size=18,family = 'Arial'),axis.title.y = element_text(size=18,family = 'Arial'),axis.text.y = element_text(size=18,family = 'Arial'))
cairo_pdf(filename = 'marker_liver1.pdf',width = 12,height = 10)
dev.off()

data$cellType=factor(data$cellType,levels = c('Hepatocyte','Cholangiocyte','HSC','B','Plasma','T','NK','Macrophage','Monocyte','Neutrophil','pDC','cDC','Endothelial','Fibroblast'))
####sccaf-d marker
Hepatocyte：'TF','APOC3'
Cholangiocyte:'DEFB1','ANXA4','ELF3',
HSC:'DCN','TIMP1','BGN'
B:'CD79A','CD79B'
Plasma:'MZB1','JCHAIN',
T:'CD3D','CD3E'
  NK:'NKG7','KLRF1','KLRD1'
  macrophage:'CD163','MS4A7'
  monocyte:'CD68','FCN1','APOBEC3A','CD14'
  neutrophil:'S100A8','LTF','LYZ'
  pDC:'IRF7','ITM2C','GZMB'
  cDC:'CLEC9A','IRF8'
  endothelial:'FCN3','LDB2','VWF'
  fibroblast:'COL3A1','LUM'

  marker=c('TF','APOC3','DEFB1','ANXA4','ELF3','DCN','TIMP1','BGN','CD79A','CD79B','MZB1','JCHAIN','CD3D','CD3E','NKG7','KLRF1','KLRD1','CD163','MS4A7','CD68','FCN1','APOBEC3A','CD14','S100A8','LTF','LYZ','IRF7','ITM2C','GZMB','CLEC9A','IRF8','FCN3','LDB2','VWF','COL3A1','LUM')
  