fdir = ""
## I. load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

## II. load data
peri <- readRDS("obj/Pericyte.RDS")

## Fig 2A ##
DimPlot(peri,group.by = 'celltype_plot',label = T)+
  scale_color_manual(values = c('#ec8843','#a9d591','#4584b3'))


## Fig S2A ##
peri$celltype_plot=factor(peri$celltype_plot,
  levels = c('GPR116+PC','MUSTN1+PC','STEAP4+PC'))
peri@active.ident=factor(peri$celltype_plot)
marker<-FindAllMarkers(peri,only.pos = T)
marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
  
set.seed(521)
obj<-subset(peri,cells= c(sample(colnames(peri)[peri$celltype_plot!='STEAP4+PC'],200),
                          sample(colnames(peri)[peri$celltype_plot=='STEAP4+PC'],60)))
obj<-ScaleData(obj,features = top10$gene)
DoHeatmap(obj, features = c(top10$gene),#disp.min=-2.5, disp.max=2.5,
          size = 3,draw.lines = T)+
  scale_fill_gradientn(colors = c("#41b6e6","white","#de4307"))


## Fig 2B ##
DotPlot(peri,features=c(
'ADGRF5',
'EGFL6',
'THY1',
'SORBS2',
'RERGL',
'MUSTN1',
'STEAP4',
'CCL21',
'CCL19',
'CXCL12'))+
  theme(axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_color_gradientn(colors = c("#41b6e6","white","#de4307"))
  
  
## Fig S2B ##
data<-cbind(Embeddings(object = peri[['umap',]]),FetchData(peri,'Tissue'))

ggplot(data,aes(x=UMAP_1,y=UMAP_2))+stat_density_2d(aes(fill =..level..),geom = "polygon", colour= "white")+
  facet_wrap(~Tissue,ncol = 5)+
  scale_fill_distiller(palette=4,direction=1)+theme_bw()


## Fig S2C ##
meta=peri@meta.data
mm=table(meta$Tissue,meta$celltype_plot)
mm <- matrix(mm, ncol=ncol(mm), dimnames=dimnames(mm))
OR_mm=matrix(rep(-1,3*3),3,3)
rownames(OR_mm)=rownames(mm)
colnames(OR_mm)=colnames(mm)

p_mm=matrix(rep(-1,3*3),3,3)
rownames(p_mm)=rownames(mm)
colnames(p_mm)=colnames(mm)

for (i in 1:3){
  for (j in 1:3){
    cell=rownames(mm)[i]
    tissue=colnames(mm)[j]
    
    n1=mm[i,j]        # celltype i in Tissue j
    
    if(j==1){n2=sum(mm[i,c(2:3)])}
    if(j==3){n2=sum(mm[i,c(1:2)])}
    if(j %in% c(2)){n2=sum(mm[i,c(1,3)])}   # celltype i not in Tissue j
    
    if(i==1){n3=sum(mm[c(2:3),j])}
    if(i==3){n3=sum(mm[c(1:2),j])}
    if(i %in% c(2)){n3=sum(mm[c(1,3),j])}   # not celltype i in Tissue j
    
    n4= sum(mm)-mm[i,j] # not celltype i not in Tissue j
    
    data <- matrix(c(n1, n2, n3, n4), byrow = T,nrow = 2) 
    colnames(data) <- c("celltype 1", "celltype 2") 
    rownames(data) <- c("tissue 1", "tissue 2")
    
    test=fisher.test(data)
    
    OR_mm[i,j]=test$estimate
    
    p_mm[i,j]=test$p.value
    
    #test=chisq.test(data)
    
    #test$observed/test$expected
    
    #test$p.value
  }
}

# significance
adjustp<-matrix(p.adjust(p_mm,"BH"),byrow = F,3,3)
pmt=adjustp
if (!is.null(pmt)){
  
  sssmt <- pmt< 0.001
  pmt[sssmt] <-'***'
  ssmt <- pmt >=0.001& pmt <0.01
  pmt[ssmt] <-'**'
  smt <- pmt >=0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt&!sssmt]<- ''
} else {
  pmt <- F
}

pheatmap(OR_mm,display_numbers = pmt,cluster_rows = F,cluster_cols = F,scale = 'column',border_color = 'white')


## Fig S2E ##

# merge 5 data sets

# Datasets Filtering Criteria
# 1. 读取 clean后data
# 2. our data 不要N
# 3. 筛选n>50样本
# 4. 进行聚类分群

# 1. our escc
data1 <- readRDS("obj/pericyte.RDS")
data1$batch=rep('Our',dim(data1)[2])

# 2. Liu_2023
data2 <- readRDS("public_scRNA/peri_liu2023_clean.RDS")
data2$batch=rep('Liu_2023',dim(data2)[2])

# 3. Zhang_2021
data3 <- readRDS("public_scRNA/peri_zhang2021_clean.RDS")
data3$batch=rep('Zhang_2021',dim(data3)[2])

# 4. Pan_2022
data4 <- readRDS("public_scRNA/peri_pan2022_clean.RDS")
data4$batch=rep('Pan_2022',dim(data4)[2])

# 5. Dinh_2021
data5 <- readRDS("public_scRNA/peri_dinh2021_clean.RDS")
data5$batch=rep('Dinh_2021',dim(data5)[2])

final=merge(data5,y=c(data4,data3,data2,data1)) 

final <- NormalizeData(final)
final <- FindVariableFeatures(final, selection.method = "vst", nfeatures = 2000)
final <- ScaleData(final)
final <- RunPCA(final)
final<-RunHarmony(final,"batch")
gc()
final <- FindNeighbors(final, dims = 1:10,reduction = "harmony")
final <- FindClusters(final, resolution = 0.2)
final<- RunUMAP(final,reduction = "harmony", dims=1:10)
# DimPlot(final, reduction = "umap",label = T)
final<- RunTSNE(final,reduction = "harmony", dims=1:10)
DimPlot(final, reduction = "tsne",label = T,pt.size = 0.1)+
  scale_color_manual(values = c('#e1cc77','#c5b0d5','#9edae5','#b29773','#ffbb78'))

DimPlot(final, reduction = "tsne",label =F,group.by = 'batch',pt.size = 0.01)+
  scale_color_manual(values = c('#dbdb8d','#aa40fc','#d62728','#279e68','#ff7f0e'))


## Fig S2F ##

final=AddModuleScore(final,features = list(c(
'ADGRF5',
'EGFL6',
'SPARC',
'COL1A1',
'COL3A1',
'THY1',
'IFI27',
'TIMP1',
'PRSS23',
'TNFAIP6')),
   ctrl = 100,name='GPR116+PC_score')
final=AddModuleScore(final,features = list(c(
'SORBS2',
'MUSTN1',
'MYH11',
'ADIRF',
'RERGL',
'PLN',
'SNCG',
'TAGLN',
'GADD45B',
'MT1M')),
  ctrl = 100,name='SORBS2+PC_score')
final=AddModuleScore(final,features = list(c(
"STEAP4",
"CCL21",
"CCL19",
"CFD",
"CD36",
"FABP4",
"CFH",
"CXCL12",
"CCDC80",
"APOE")),
  ctrl = 100,name='STEAP4+PC_score')
  
## Fig S2G ##

av=AverageExpression(final,assays = 'RNA')
av=av[[1]]
# Select the top 200 genes with the biggest standard deviation
cg=names(tail(sort(apply(av, 1, sd)),200))
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),border_color = '#FFF8DC')

## Fig S2H ##

final$anno=sapply(final$seurat_clusters,function(x){
  if(x ==1){y='STEAP4+PC'}
  else if(x==2){y='SORBS2+PC'}
  else{y='GPR116+PC'}
  y
})

DimPlot(final, reduction = "tsne",label =F,group.by = 'anno',pt.size = 0.01)+
  scale_color_manual(values = c())
  
## Fig S2I ##

DotPlot(final,features = c(
'EGFL6',
'THY1',
'ADGRF5',
'MUSTN1',
'MYH11',
'RERGL',
'SORBS2',
'STEAP4',
'CCL21',
'CCL19',
'CXCL12'),group.by = 'anno')+
  theme(axis.text.x = element_text(angle = 35,hjust = 1))+scale_color_gradientn(colors = c("#41b6e6","white","#de4307"))

## Fig S2J ##

test=subset(final,batch=='Our')

a=table(test$celltype_plot,test$anno)

b= a/rowSums(a)

pheatmap::pheatmap(b,cluster_rows = F,cluster_cols = F,scale='none',display_numbers = T,border_color = 'white')

## Fig S2K & Fig 2C ##

final$group=factor(final$group,levels = c('T_NM','T_M'))
meta=final@meta.data

as.data.frame(meta) %>%
  group_by(Sample) %>%
  summarise(n=n()) ->week.cell

as.data.frame(meta) %>%
  group_by(Sample,anno) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("Sample"="Sample")) %>%
  dplyr::mutate(ratio = cls_n/n) ->ratio

sub=unique(meta[,c('Sample','group')])

ratio=left_join(ratio,sub,by='Sample')

ratio=ratio[ratio$n>50,]

p<-ggboxplot(ratio, x = 'group', y = 'ratio',
             color = 'group', add = "jitter")+
  facet_wrap(~anno,nrow = 1)+
  theme(axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_color_manual(values = c('#66aff9','#ff9000'))

p+stat_compare_means(method = 't.test',label.y = 1.1)+ylim(0,1.2)

## Fig 2M ##

data=read.table('GPR116_PC-TABLE-DATAANALYZE.txt',header=T,sep='\t')

plot_data=data[!(is.na(data$GPR116PC_percent)),]
plot_data$Survival_state <- as.numeric(plot_data$Status)
plot_data$Survival_Month <- as.numeric(plot_data$survival.month)

surv_data=plot_data
surv_data$Survival_state <- as.numeric(surv_data$Survival_state)
surv_data$Survival_Month <- as.numeric(surv_data$Survival_Month)
#########################################################
# segregate exp by median

# add new col "high"
surv_data <- transform(surv_data, GPR116PC_percent_group = "high")

# change the value based exp median
surv_data[which(surv_data[,'GPR116PC_percent'] <= median(surv_data[,'GPR116PC_percent'])),'GPR116PC_percent_group'] <- "low"

# plot
library("survival")
library("survminer")

fit <- survfit(Surv(Survival_Month, Survival_state) ~ surv_data[,c('GPR116PC_percent_group')], data = surv_data)
print(fit)
# n events median 0.95LCL 0.95UCL
# surv_data[, c("GPR116PC_percent_group")]=high 100     85   14.5      11      19
# surv_data[, c("GPR116PC_percent_group")]=low  100     69   26.0      18      42

ggsurvplot(fit,
           pval = TRUE, #conf.int = TRUE,# confidence interval
           #cumevents=T,
           risk.table = TRUE, # add risk table
           risk.table.col = "strata", # 
           #lineTissue = "strata", # 
           surv.median.line = "hv", # median survival
           ggtheme = theme_bw(), # theme
           #palette = c("#35a12b", "#0d79b4","#f39999"),
           #fun = "event" #
           legend.labs=c("high","low"),# legend labels
           #legend.titel ='',
           xlab="Time(Months)" )