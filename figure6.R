fdir = ""
## I. load packages
library("Seurat")
library('ggplot2')
library('ggpubr')
library('dplyr')
library('UCell')
library('irGSEA')

## II. load data
escc = read_rds("obj/escc_merge.rds")
tumor <- readRDS("obj/Tumor_cells.RDS")

## Fig 6A ##
cellchat_EMT <- readRDS("cellchat_peri_Epi_subset.RDS")
table(cellchat_EMT@idents)
# MUSTN1+PC STEAP4+PC GPR116+PC  nEpi  nonEMT    lowEMT   highEMT 

pairs=read.table('LRpair_related_to_EGFL6.txt',
                header = T,sep = '\t')

netVisual_bubble(cellchat_EMT, sources.use = c(1:3), targets.use = c(4:7), 
                 remove.isolate = FALSE, 
                 pairLR.use = pair)
				 
				
## Fig 6B ##

EGFL6_receptors=c('EGFR',
'HDGF',
'ALG13',
'LDLR',
'SYVN1',
'ITGB1',
'ITGB2',
'ITGB3',
'ITGB4',
'ITGB5',
'ITGB6',
'ITGB7',
'ITGB8')

a=as.data.frame(AverageExpression(tumor,group.by = 'EMT_group',
features = EGFL6_receptors,
assays = 'RNA',slot = 'data')$RNA)

a=a[order(a$high_EMT,decreasing = T),]
pheatmap(a,cluster_rows = F,cluster_cols = F,border_color = 'black',scale = 'column')

b=as.data.frame(AverageExpression(tumor,group.by = 'Tissue',
features = EGFL6_receptors,
assays = 'RNA',slot = 'data')$RNA)

b=b[order(b$M_T,decreasing = T),]
pheatmap(b,cluster_rows = F,cluster_cols = F,border_color = 'black',scale = 'column')


## Fig 6D ##

VlnPlot(tumor,features='ITGB1',group.by = 'Tissue',pt.size = 0)+
  geom_boxplot(width=0.3,color='black',position=position_dodge(1))+
  stat_compare_means(label="p.signif")
  
VlnPlot(tumor,features='ITGB1',group.by = 'EMT_group',pt.size = 0)+
  geom_boxplot(width=0.3,color='black',position=position_dodge(1))+
  stat_compare_means(label="p.signif")
  
## Fig 6F ##

score.msidb_H_ssgsea <- irGSEA.score(object = tumor, assay = "RNA",
                              slot = "data", seeds = 123, ncores = 1,
                              min.cells = 3, min.feature = 0,
                              custom = F, 
                              geneset = NULL, # change to user geneset
                              msigdb = T,
                              species = "Homo sapiens", 
                              category = "H",  # msigdbr::msigdbr_collections() to view all available collections gene sets
                              subcategory = NULL, geneid = "symbol",
                              method = c("UCell"),
                              #"ssgsea","singscore", "UCell",'AUCell'),
                              aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                              kcdf = 'Gaussian')

# 返回一个Seurat对象，富集分数矩阵存放在RNA外的assay中
Seurat::Assays(score.msidb_H)

meta<-as.data.frame(score.msidb_H_UCell@assays$UCell@data)
meta<-as.data.frame(t(meta))

sum(rownames(meta)==colnames(tumor))
meta$anno<-tumor$EMT_group
meta$Tissue<-tumor$Tissue

# 使用UCell数据绘制箱式图
ggplot(meta,aes(Tissue,`HALLMARK-TNFA-SIGNALING-VIA-NFKB`))+
  geom_violin(aes(fill=Tissue))+
  geom_boxplot(width=0.3,fill='white')+
  facet_wrap(~anno)+
  stat_compare_means(method = 't.test')+
  theme_pubr()


