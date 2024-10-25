fdir = ""
## I. load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(harmony)
library(RColorBrewer)

## II. load data
escc = read_rds("obj/escc_merge.rds")
peri <- readRDS("obj/Pericyte.RDS")

## Fig S8A  ## 

CD8T=subset(escc,subtype_plot=='CD8+T')

dim(CD8T) 
CD8T$Tissue=factor(CD8T$Tissue,levels = c('N','T_NM','T_M'))

CD8T <- NormalizeData(CD8T)
CD8T <- FindVariableFeatures(CD8T, selection.method = "vst", nfeatures = 2000)
CD8T <- ScaleData(CD8T)
CD8T <- RunPCA(CD8T)
library(harmony)
CD8T<-RunHarmony(CD8T,"Sample")
gc()
CD8T <- FindNeighbors(CD8T, dims = 1:10,reduction = "harmony")
CD8T <- FindClusters(CD8T, resolution = 0.4) 
CD8T<- RunUMAP(CD8T,reduction = "harmony", dims=1:10)

DimPlot(CD8T,label = T)

mk=FindAllMarkers(CD8T)

CD8T$anno=factor(CD8T$anno,levels = c('CD8T_CD69','CD8T_CTLA4','CD8T_STMN1','CD8T_IL7R','CD8T_IL32','CD8T_FGFBP2','NK_FCER1G'))
CD8T@active.ident=CD8T$anno


DimPlot(CD8T,label = T)+scale_color_manual(values  = c('#73c5d3',
  '#d099c3',
  '#a9d591',
  '#f4bbbb',
  '#ec8843',
  '#e1cc77',
  '#b29773',
  '#e53e91',
  '#4584b3',
  '#a1c6d9',
  '#bda6cc',
  '#7a519e'))


## Fig S8B  ## 

DotPlot(CD8T,features = c(
  'CD69',
  'TNFSF9',
  'IFNG', 
  'FOS',
  'CCL4L2',
  'NR4A1',
  'FOSB',
  'CTLA4',
  'CXCL13',
  'HAVCR2',
  'LAG3',
  'TIGIT',
  'PDCD1',
  'LAYN',
  'TOX',
  'STMN1',
  'MCM5',
  'TYMS',
  'DUT',
  'GZMB',
  'IL7R',
  'CXCR4',
  'EZR',
  'CD44',
  'TGFB1',
  'IL32',
  'CCL5',
  'CD52', 
  'EVL',
  'RNF213',
  'GIMAP7',
  'IKZF3',
  'NKG7',
  'PRF1',
  'FGFBP2',
  'FCGR3A',
  'KLF2',
  'KLRF1',
  'TNF',
  'GZMH',
  'KLRD1',
  'FCER1G',
  'TYROBP',
  'AREG',
  'GNLY',
  'KLRC1'
))+scale_color_distiller(palette = "Spectral")+
theme(axis.text.x = element_text(angle = 60,hjust = 1))+
coord_flip()

## Fig S8C  ## 

as.data.frame(CD8T@meta.data) %>%
  group_by(Sample) %>%
  summarise(n=n()) ->week.cell

CD8T$Sample=factor(CD8T$Sample,levels = c(
  "P03_T", "P06_T", "P07_T", "P21_T", "P22_T", "P25_T", 
  "P04_T", "P08_T", "P14_T", "P26_T", "P27_T", "P28_T",
  "P04_N", "P08_N", "P26_N", "P27_N"))

as.data.frame(CD8T@meta.data) %>%
  group_by(Sample,anno) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("Sample"="Sample")) %>%
  dplyr::mutate(ratio = cls_n/n) ->data.plot

ggplot(data.plot,aes(Sample,cls_n,fill=anno))+
  geom_bar(stat="identity",position="fill")+theme_pubr()+
  theme(legend.position = "right",legend.title = element_blank())+
  scale_fill_manual(values=cluster_color)+
  xlab("")+
  ylab("Cell Compartment")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))


## Fig S8E  ## 

CD4T=subset(escc,subtype=='CD4T') 
dim(CD4T) 
CD4T$Tissue=factor(CD4T$Tissue,levels = c('N','T_NM','T_M'))

CD4T <- NormalizeData(CD4T)
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 2000)
CD4T <- ScaleData(CD4T)
CD4T <- RunPCA(CD4T)

CD4T<-RunHarmony(CD4T,"Sample")
gc()
CD4T <- FindNeighbors(CD4T, dims = 1:20,reduction = "harmony")
CD4T <- FindClusters(CD4T, resolution = 0.4) 
CD4T<- RunUMAP(CD4T,reduction = "harmony", dims=1:20)

DimPlot(CD4T,label = T)

mk=FindAllMarkers(CD4T)

CD4T$anno=factor(CD4T$anno,levels = c('CD4T_CCR7','CD4T_FOXP3','CD4T_MALAT1','CD4T_KLRB1','CD4T_CXCL13','CD4T_STMN1'))
CD4T@active.ident=CD4T$anno

DimPlot(CD4T,label = T)+scale_color_manual(values  = c('#73c5d3',
  '#d099c3',
  '#a9d591',
  '#f4bbbb',
  '#ec8843',
  '#e1cc77',
  '#b29773',
  '#e53e91',
  '#4584b3',
  '#a1c6d9',
  '#bda6cc',
  '#7a519e'))

## Fig S8F  ## 

DotPlot(CD4T,features = c(
'CCR7',
'CD55',
'KLF2',
'ANXA1',
'IL7R',
'FOXP3',
'BATF',
'IKZF2',
'TNFRSF4',
'IL2RA',
'MALAT1',
'NEAT1',
'CCL5',
'KLRB1',
'GZMA',
'GZMB',
'CCL20',
'CAPG',
'IFNG',
'CXCL13',
'PTPN13',
'TOX2',
'PDCD1',
'COTL1',
'TOX',
'STMN1',
'PCLAF',
'HMGN2',
'TUBB'
))+scale_color_distiller(palette = "Spectral")+
theme(axis.text.x = element_text(angle = 60,hjust = 1))+
coord_flip()


## Fig S8G  ## 

as.data.frame(CD4T@meta.data) %>%
  group_by(Sample) %>%
  summarise(n=n()) ->week.cell

CD4T$Sample=factor(CD4T$Sample,levels = c(
  "P03_T", "P06_T", "P07_T", "P21_T", "P22_T", "P25_T", 
  "P04_T", "P08_T", "P14_T", "P26_T", "P27_T", "P28_T",
  "P04_N", "P08_N", "P26_N", "P27_N"))

as.data.frame(CD4T@meta.data) %>%
  group_by(Sample,anno) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("Sample"="Sample")) %>%
  dplyr::mutate(ratio = cls_n/n) ->data.plot

ggplot(data.plot,aes(Sample,cls_n,fill=anno))+
  geom_bar(stat="identity",position="fill")+theme_pubr()+
  theme(legend.position = "right",legend.title = element_blank())+
  scale_fill_manual(values=cluster_color)+
  xlab("")+
  ylab("Cell Compartment")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

## Fig S8D & S8H  ## 

ratio=read.table('ratio_subcluster_patient_tumor.txt',sep='\t',header=T)

cor_pearson=matrix(8,dim(ratio)[2],dim(ratio)[2])
p_pearson=matrix(8,dim(ratio)[2],dim(ratio)[2])
for (i in 1:dim(ratio)[2]){
  for (j in 1:dim(ratio)[2]){
    cor_pearson[i,j]=cor.test(ratio[,i],ratio[,j],method = 'pearson')$estimate
    p_pearson[i,j]=cor.test(ratio[,i],ratio[,j],method = 'pearson')$p.value
  }
}
rownames(cor_pearson)=colnames(ratio)
colnames(cor_pearson)=colnames(ratio)
rownames(p_pearson)=colnames(ratio)
colnames(p_pearson)=colnames(ratio)

cor_spearman=matrix(8,dim(ratio)[2],dim(ratio)[2])
p_spearman=matrix(8,dim(ratio)[2],dim(ratio)[2])
for (i in 1:dim(ratio)[2]){
  for (j in 1:dim(ratio)[2]){
    cor_spearman[i,j]=cor.test(ratio[,i],ratio[,j],method = 'spearman')$estimate
    p_spearman[i,j]=cor.test(ratio[,i],ratio[,j],method = 'spearman')$p.value
  }
}
rownames(cor_spearman)=colnames(ratio)
colnames(cor_spearman)=colnames(ratio)
rownames(p_spearman)=colnames(ratio)
colnames(p_spearman)=colnames(ratio)

# 13 - GPR116+ PC
pearson1= cbind(p_pearson[,13],cor_pearson[,13])
spearman1= cbind(p_spearman[,13],cor_spearman[,13])

# 14 - MUSTN+ PC
pearson2= cbind(p_pearson[,14],cor_pearson[,14])
spearman2= cbind(p_spearman[,14],cor_spearman[,14])

# 16 - STEAP4+ PC
pearson3= cbind(p_pearson[,16],cor_pearson[,16])
spearman3= cbind(p_spearman[,16],cor_spearman[,16])


data1=as.data.frame((ratio[,'GPR116+PC']))
colnames(data1)[1]='ratio_GPR116_PC'

data2=as.data.frame((ratio[,'CD8T_CTLA4']))
colnames(data2)[1]='ratio_CD8T_CTLA4'

data3=as.data.frame((ratio[,'CD4T_FOXP3']))
colnames(data3)[1]='ratio_CD4T_FOXP3'

data=cbind(data1,data2)
data$patient=rownames(data)
ggplot(data,aes(x=ratio_GPR116_PC,y=ratio_CD8T_CTLA4))+
geom_smooth(method = "lm",fill='#A7D8F2',color='#4A90E2')+
geom_point(aes(color=patient))+
theme_pubr()+
stat_cor(method = "pearson")

data=cbind(data1,data3)
data$patient=rownames(data)
ggplot(data,aes(x=ratio_GPR116_PC,y=ratio_CD4T_FOXP3))+
geom_smooth(method = "lm",fill='#A7D8F2',color='#4A90E2')+
geom_point(aes(color=patient))+
theme_pubr()+
stat_cor(method = "pearson")


## Fig S8J & S8L  ## 

treg=read.table('corr_116PC_CD4Treg_GEPIA2.txt',header = T)

PD1=read.table('corr_116PC_CD8TPD1_GEPIA2.txt',header = T)

treg$sig=sapply(treg$p_value,function(x){
  if(x<0.05){y='yes'}
  else{y='no'}
  y
})
treg=treg[order(treg$Spearman,decreasing = T),]
treg$Cancer_type=factor(treg$Cancer_type,levels = treg$Cancer_type)
ggplot(treg,aes(x=Spearman,y=Cancer_type))+
geom_bar(stat='identity',width = 0.8,aes(fill=sig))+
  scale_fill_manual(values = c('Grey','Orange'))+
  ylab(label = NULL)+xlab(label = 'Spearman_corr')+
  theme (axis.text.x = element_text (color="black"), axis.text.y = element_text ( color="black", size=10.5))+
  geom_vline(xintercept=c(0.3), linetype="dotted")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  coord_flip()
  



PD1$sig=sapply(PD1$p_value,function(x){
  if(x<0.05){y='yes'}
  else{y='no'}
  y
})
PD1=PD1[order(PD1$Spearman,decreasing = T),]
PD1$Cancer_type=factor(PD1$Cancer_type,levels = PD1$Cancer_type)
ggplot(PD1,aes(x=Spearman,y=Cancer_type))+
geom_bar(stat='identity',width = 0.8,aes(fill=sig))+theme_bw()+
  scale_fill_manual(values = c('Grey','#4EB9F9'))+
  ylab(label = NULL)+xlab(label = 'Spearman_corr')+
  theme (axis.text.x = element_text (color="black"), axis.text.y = element_text ( color="black", size=10.5))+
  geom_vline(xintercept=c(0.3), linetype="dotted")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  coord_flip()
