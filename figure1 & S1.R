fdir = ""
## I. load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(CellChat)
library(pheatmap)

## II. load data
escc = read_rds("obj/escc_merge.rds")
sub_meta = escc@sub_meta.data

## Fig1B ## 
escc$subtype<-factor(escc$subtype,levels = c(
  'Epithelia',
  'Endothelial',
  'Pericyte',
  'Fibroblast',
  'Bcell',
  'Plasma',
  'CD8+T',
  'Cycling',
  'CD4+T',
  'Mast',
  'Macrophage',
  'Monocyte',
  'Neutrophil',
  'DC'
))
escc@active.ident=escc$subtype

subtype_color=c(
  '#73c5d3',
  '#a56e40',
  '#2da94d',
  '#a9d591',
  '#f4bbbb',
  '#da463e',
  '#60a098',
  '#ec8843',
  '#d099c3',
  '#e1cc77',
  '#bda6cc',
  '#4584b3',
  '#a1c6d9',
  '#b29773'
)

DimPlot(escc,group.by = 'subtype',label = T,raster=FALSE,reduction = 'tsne')+
scale_color_manual(values = c(subtype_color))

## FigS1A ##

escc$Tissue=factor(escc$Tissue,levels = c('N','T_NM','T_M'))
DimPlot(escc,group.by = 'Tissue',label = F,raster=F,reduction = 'tsne')+
scale_color_manual(values = c('#a9d591','#ec8843','#73c5d3'))

sample_color=c(
  '#1f77b4', 
  '#ff7f0e', 
  '#279e68', 
  '#d62728', 
  '#aa40fc', 
  '#8c564b', 
  '#e377c2', 
  '#b5bd61', 
  '#17becf', 
  '#9edae5',
  '#ff9896',
  '#ffbb78', 
  '#98df8a', 
  '#c5b0d5', 
  '#c49c94', 
  '#aec7e8',  
  '#ad494a'  
)
DimPlot(escc,group.by = 'Sample',label = F,raster=FALSE,reduction = 'tsne')+
scale_color_manual(values =sample_color))

## FigS1B ##

list_genes=list(  Epithelial=c('EPCAM','KRT5', 'KRT15', 'KRT18'), # epithelial cells
                  Endothelial=c('PECAM1', 'EGFL7','VWF'), # endothelial cells
                  Pericyte=c('ACTA2','MCAM', 'RGS5', 'CSPG4'), # pericytes
                  Fibroblast=c('DCN', 'COL1A1','COL3A1'), # fibroblast cells
                  Bcell=c('CD79A', 'MS4A1', 'CD19'), # B cells
                  Plasma=c('IGHG1', 'MZB1','DERL3'),# plasma
                  CD8T=c('CD3D', 'CD8A','CCL5'),# CD8T cells
                  Cycling=c('MKI67','CDK1' ), # cycling cells
                  CD4T=c('FOXP3','IL7R','ICOS' ), # CD4T cells
                  Mast=c('TPSAB1', 'CPA3' , 'KIT'), # Mast
                  Macrophage=c('CD163','C1QB','LGMN'),# Macrophage
                  Monocyte=c('IL1B','EREG','LYZ'), # Monocyte
                  Neutrophil=c('FCGR3B','G0S2','S100A9'),# Neu 
                  DC=c('CD1C','CD1E','S100B')# DC
)

DotPlot(escc,
           features=list_genes,
           #cols = c("grey", "red"),
           cluster.idents = F)+
  RotatedAxis()+
  theme(
    panel.border = element_rect(color="black"), 
    panel.spacing = unit(1, "mm"), 
    #strip.background = element_rect(color="red"),
    strip.text = element_text(margin=margin(b=1, unit="mm")),
    strip.placement = 'outlet', #
    axis.line = element_blank(),
  )+labs(x="", y="")+
  scale_color_gradientn(colors = c("#41b6e6","white","#de4307"))

## FigS1C ## 

as.data.frame(sub_meta) %>%
  group_by(Sample) %>%
  summarise(n=n()) ->week.cell

as.data.frame(sub_meta) %>%
  group_by(Sample,subtype) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("Sample"="Sample")) %>%
  dplyr::mutate(ratio = cls_n/n) ->data.plot

ggplot(data.plot,aes(Sample,cls_n,fill=subtype))+
  geom_bar(stat="identity",position="fill")+theme_pubr()+
  theme(legend.position = "right",legend.title = element_blank())+
  scale_fill_manual(values=subtype_color)+
  xlab("")+
  ylab("Cell Compartment")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

## Fig1C

myloid=subset(escc,subtype %in% c(
'Mast',
'Macrophage',
'Monocyte',
'Neutrophil',
'DC'))
stromal=subset(escc,subtype %in% c(
'Endothelial',
'Pericyte',
'Fibroblast'))
lymphocyte=subset(escc,subtype %in% c(
'Bcell',
'Plasma',
'CD8+T',
'CD4+T'))
 
sub_meta=stromal@sub_meta.data
# the following two steps are performed separately
# sub_meta=lymphocyte@sub_meta.data
# sub_meta=myeloid@sub_meta.data

sub_meta$Tissue<-factor(sub_meta$Tissue,levels = c('M_N','P_T','M_T'))
as.data.frame(sub_meta) %>%
  group_by(Sample) %>%
  summarise(n=n()) ->week.cell

as.data.frame(sub_meta) %>%
  group_by(Sample,subtype) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("Sample"="Sample")) %>%
  dplyr::mutate(ratio = cls_n/n) ->ratio

sub=unique(sub_meta[,c('Sample','Tissue')])
ratio=left_join(ratio,sub,by='Sample')

ggboxplot(ratio, x = 'Tissue', y = 'ratio',size = 0.3,add.params=list(size=0.5),
             fill = 'Tissue', add = "jitter",alpha=0.8)+
  facet_wrap(~subtype,nrow = 1)+
  theme(axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_fill_manual(values = c('#94D0FF','#F7E967','#F89D4A'))+ylim(-0.01,1)
  
  
## Fig1D

cellchat.MN <- readRDS('cellchat_N.RDS')
cellchat.PT <- readRDS('cellchat_T_NM.RDS')
cellchat.MT <- readRDS('cellchat_T_M.RDS')
object.list <- list(MN=cellchat.MN,PT = cellchat.PT, MT = cellchat.MT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
res=rbind(gg[[1]]$data,gg[[2]]$data,gg[[3]]$data)
res$tissue=rep(c('N','T_NM','T_M'),each=13)
# N
N_strength= data.frame(celltype= rep(c("Bcell","Plasma","CD4+T","CD8+T",'Epithelia',"Endothelial","Fibroblast",'Pericyte',
                         'DC',"Macrophage","Mast","Monocyte","Neutrophil"),2),
                         Interaction_strength = c(-(res[1:13,1]),res[1:13,2]),
                         type= rep(c('outgoing','incoming'),each=13))
N_strength$type=factor(N_strength$type,levels = c('outgoing','incoming'))
ggplot(N_strength,aes(x=celltype,y=Interaction_strength,fill=type))+
  geom_bar(stat = "identity",width = 0.65)+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
  theme_classic()+
  coord_flip()
# T_NM
TNM_strength= data.frame(celltype= rep(c("Bcell","Plasma","CD4+T","CD8+T",'Epithelia',"Endothelial","Fibroblast",'Pericyte',
                         'DC',"Macrophage","Mast","Monocyte","Neutrophil"),2),
                         Interaction_strength = c(-(res[14:26,1]),res[14:26,2]),
                         type= rep(c('outgoing','incoming'),each=13))
TNM_strength$type=factor(TNM_strength$type,levels = c('outgoing','incoming'))
ggplot(TNM_strength,aes(x=celltype,y=Interaction_strength,fill=type))+
  geom_bar(stat = "identity",width = 0.65)+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
  theme_classic()+
  coord_flip()
# T_M 
TM_strength= data.frame(celltype= rep(c("Bcell","Plasma","CD4+T","CD8+T",'Epithelia',"Endothelial","Fibroblast",'Pericyte',
                         'DC',"Macrophage","Mast","Monocyte","Neutrophil"),2),
                         Interaction_strength = c(-(res[27:39,1]),res[27:39,2]),
                         type= rep(c('outgoing','incoming'),each=13))
TM_strength$type=factor(TM_strength$type,levels = c('outgoing','incoming'))
ggplot(TM_strength,aes(x=celltype,y=Interaction_strength,fill=type))+
  geom_bar(stat = "identity",width = 0.65)+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
  theme_classic()+
  coord_flip()

# heatmap of comparation between tissues
out_t_n=res[27:39,1]-res[1:13,1]
out_mt_pt=res[27:39,1]-res[14:26,1]
in_t_n=res[27:39,2]-res[1:13,2]
in_mt_pt=res[27:39,2]-res[14:26,2]

mat=matrix(c(out_t_n,out_mt_pt,in_t_n,in_mt_pt),13,4,byrow = F)
rownames(mat)=rownames(res)[1:13]
library(pheatmap)
pheatmap(mat,cluster_rows = F,cluster_cols = F,border_color = 'white',scale = 'column')

  
## Fig1E
netVisual_heatmap(cellchat, measure = "weight",comparison = c(2,3))


## FigS1D ## 
groupSize <- table(cellchat_PT@idents) %>% as.numeric()
mat <- matrix(0, 
              nrow = nrow(cellchat_PT@net$weight),
              ncol = ncol(cellchat_PT@net$weight),
              dimnames = dimnames(cellchat_PT@net$weight))
mat[8, ] <- cellchat_PT@net$weight[8, ]
netVisual_circle(mat, 
                 vertex.weight = groupSize,
                 weight.scale = T,label.edge = T,
                 edge.weight.max = max(cellchat_P_T@net$weight),vertex.label.cex=0.8,
                 title.name = 'T_NM')


groupSize <- table(cellchat_MT@idents) %>% as.numeric()
mat <- matrix(0, 
              nrow = nrow(cellchat_MT@net$weight),
              ncol = ncol(cellchat_MT@net$weight),
              dimnames = dimnames(cellchat_MT@net$weight))
mat[8, ] <- cellchat_MT@net$weight[8, ]
netVisual_circle(mat, 
                 vertex.weight = groupSize,
                 weight.scale = T,label.edge = T,
                 edge.weight.max = max(cellchat_MT@net$weight),vertex.label.cex=0.8,
                 title.name = 'T_NM')

## Fig1H

res_scale=read.table('malignant_around100um_existRatio_celltype_splitsample3_scale_result.txt',header = T,sep = '\t')
res_scale$Group=factor(res_scale$Group,levels = c('NM','M'))
res_scale$Celltype=factor(res_scale$Celltype,levels = c('Pericytes','Myeloid','Endothelial','Fibroblasts','T.NK','Bcells'))
ggplot(res_scale,aes(x=Celltype,y=Ratio,color=Group,size=Sample,shape=Sample))+
geom_point()+theme_bw()+
  scale_size_manual(values = c(2,2,2,2,2,4,4))+
  scale_shape_manual(values = c(20,20,20,20,20,17,17))+
  scale_color_manual(values = c('#66aff9','#ff9000'))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))