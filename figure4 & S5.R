fdir = ""
## I. load packages
library("infercnv")
library("Seurat")
library('ggplot2')
library('ggpubr')
library('CellChat')
library('harmony')

## II. load data
escc = read_rds("obj/escc_merge.rds")
Pericyte <- readRDS("obj/Pericyte.RDS")

## Fig S5A ##
epi=subset(escc,celltype=='Epithelial')
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi)

# run harmony merge data
epi$Sample<-factor(epi$Sample,levels = c('P03_T_0','P06_T_0','P07_T_0','P21_T_0','P22_T_0','P25_T_0',
                                           'P04_T_0','P08_T_0','P14_T_0','P26_T_0','P27_T_0','P28_T_0'))
epi<-RunHarmony(epi,"Sample")
gc()
epi <- FindNeighbors(epi, dims = 1:20,reduction = "harmony")
epi <- FindClusters(epi, resolution = 0.3)
epi<- RunUMAP(epi,reduction = "harmony", dims=1:20)
epi<- RunTSNE(epi,reduction = "harmony", dims=1:20)

DimPlot(epi,reduction = 'tsne')

sample_color=c(
  '#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#9edae5', 
  '#ff9896','#ffbb78', '#98df8a', '#c5b0d5', '#c49c94', '#aec7e8',  '#ad494a')
DimPlot(epi,group.by = 'Patient',reduction = 'tsne',cols = sample_color)


## Fig S5B,C ##
# run infercnv -----------------------------------------------------------------
# scRNA data
setwd("/home/lijy/data/infercnv/20240520")
seurat_obj_sub <- readRDS("sub_dataset_for_computeCNV.RDS") # celltype %in% c("Epithelial", "Bcells", "T_NK")

# DimPlot(seurat_obj_sub, reduction = "umap", label = T, repel = T, group.by = "celltype")

# cell annotation file
ann <- data.frame(seurat_obj_sub@meta.data)
ann <- cbind(row.names(ann), ann)
colnames(ann)[1] <- "barcode"
ann <- ann[, c("barcode", "celltype")]
write.table(ann, file = "ann.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# exp matrix
count <- GetAssayData(seurat_obj_sub, slot = 'counts')
# count_dat <- as.matrix(count)
write.table(count, file = "exp.txt", sep = "\t", quote = F)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count,
                                    annotations_file = "ann.txt",
                                    delim = "\t",
                                    gene_order_file = "~/data/infercnv/genome_position/human_genes_pos.txt",
                                    ref_group_names = c("Bcells", "T_NK"))

# perform infercnv operations to reveal cnv signal
infercnv_obj_re = infercnv::run(infercnv_obj,
                                cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = "./result/update",  # dir is auto-created for storing outputs
                                cluster_by_groups = F,   # cluster
                                analysis_mode = "subclusters",
                                denoise = T,
                                HMM = T,
                                num_threads = 1
)

cnv_score=read.table('inferCNV_cnv_obs_score.txt',header = T)
epi$CNV_score=cnv_score
DimPlot(epi,group.by = 'TN_type',reduction = 'tsne')
DimPlot(epi,group.by = 'TN_type',reduction = 'tsne')

## Fig 4A ##

name3.0=colnames(Pericyte)[which(Pericyte$cluster=='ADGRF5+PC')]
name3.1=colnames(Pericyte)[which(Pericyte$cluster=='MUSTN1+PC')]
name3.2=colnames(Pericyte)[which(Pericyte$cluster=='STEAP4+PC')]

escc$edit_subclass=apply(escc@meta.data,1,FUN = function(x){
  # x[24]:cellID
  if (x[24] %in% name3.0) {y='GPR116+PC'}
  else if (x[24] %in% name3.1) {y='MUSTN1+PC'}
  else if (x[24] %in% name3.2) {y='STEAP4+PC'}
  else if (!(x[24] %in% c(name3.0,name3.1,name3.2))){y=x[17] }#x[17]:celltype
  y
})

table(escc$edit_subclass)

escc=NormalizeData(escc)
escc=FindVariableFeatures(escc, selection.method = "vst", nfeatures = 2000)
escc=ScaleData(escc)

# cellchat need normalized data  
meta = escc@meta.data[,c('edit_subclass','Tissue')] # a dataframe with rownames containing cell mata data
data.input = escc@assays$RNA@data # normalized data matrix

cellchat <- createCellChat(object = data.input, meta = meta,group.by = 'edit_subclass')

cellchat <- setIdent(cellchat, ident.use = 'edit_subclass') # set "anno_big_edit" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
edit_intereaction=read.table('interaction_addEGFL6.txt',header = T,sep = '\t')
rownames(edit_intereaction)=edit_intereaction$interaction_name
CellChatDB$interaction=edit_intereaction

showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

options(future.globals.maxSize= 1000*1024^2)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 50)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_scatter(cellchat)

# ## Fig 4B ##

escc = read_rds("obj/escc_merge.rds")
T_NM<-subset(count,Tissue=='T_NM')
T_M<-subset(count,Tissue=='T_M')
# For T_NM and T_M, construct cellchat object separately
# ident.use = 'edit_subclass'

cellchat.TNM <- readRDS('/home/peixy/20230525_ESCC6v6_new/only_NT/Cellchat/res_Tissue/cellchat_P_T.RDS')
cellchat.TM <- readRDS('/home/peixy/20230525_ESCC6v6_new/only_NT/Cellchat/res_Tissue/cellchat_M_T.RDS')

a=cellchat.TNM@net$weight
b=cellchat.TM@net$weight

res= data.frame(edit_subtype=colnames(a),
strength_difference=b['GPR116+PC',]-a['GPR116+PC',])
ggdotchart(res, x = "edit_subtype", y = "strength_difference",
           color = "edit_subtype",                                
           sorting = "descending",                      
           add = "segments",                             
           add.params = list(color = "lightgray", size = 1.5),
           rotate = TRUE,                                
           #group = "cyl",                                
           dot.size = 6,                                 
           label = round(res$strength_difference,2),                      
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               
           ggtheme = theme_pubr()                       
)

# ## Fig S5D ##
cellchat.NM <- readRDS("cellchat_peri_Epi_split_NM.RDS")
cellchat.M <- readRDS("cellchat_peri_Epi_split_M.RDS")

object.list <- list(NM = cellchat.NM, M = cellchat.M)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# ## Fig S5E ##
epi_mk<-c(
  'CDH1',
  'CLDN3',
  'CLDN4',
  'CLDN7',
  'OCLN',
  'TJP1',
  'KRT8',
  'KRT9',
  'KRT18',
  'KRT19',
  'KRT20',
  'DSP',
  'SDC1',
  'MUC1'
  
)

mes_mk<-c(
  #'CDH2',
  'VIM',
  'FN1',
  'S100A4',
  'ACTA2',
  'COL1A1',
  'COL1A2',
  'COL3A1'
)

ent_reg<-c(
  'SNAI1',
  'SNAI2',
  'TWIST1',
  'ZEB1',
  'ZEB2',
  'EPCAM',
  'MMP2',
  'MMP3',
  'MMP9'
)
gene<-list(intersect(epi_mk,rownames(tumor)))
tumor<-AddModuleScore(
  object = tumor,
  features = gene,
  ctrl = 100,
  name = ('score_Epithelial')
)

gene<-list(intersect(mes_mk,rownames(tumor)))
tumor<-AddModuleScore(
  object = tumor,
  features = gene,
  ctrl = 100,
  name = ('score_Mesenchymal')
)

tumor$EMT_score<-tumor$score_Mesenchymal1-tumor$score_Epithelial1

tumor$EMT_group=sapply(tumor$EMT_score,FUN = function(x){
  if (x<0){y='non_EMT'}
  if (x>=0 & x<1){y='low_EMT'}
  if (x>=1){y='high_EMT'}
  y
})
table(tumor$Tissue,tumor$EMT_group)                          

tumor$EMT_group<-factor(tumor$EMT_group,levels = c('non_EMT','low_EMT','high_EMT'))

VlnPlot(tumor,'EMT_score',group.by = 'EMT_group',pt.size = 0)+scale_fill_manual(values = c('#B2B2B2','#74A0A1','#D69C9B'))+
  geom_boxplot(width=0.2,fill='white')+
  geom_hline(aes(yintercept=0),linetype="dashed") + geom_hline(aes(yintercept=1),linetype="dashed")
  
 
# ## Fig S5G ##
tumor$CellType_edit=paste(tumor$EMT_group,tumor$Tissue,sep = '_')

tumor$Tissue=factor(tumor$Tissue,levels = c('P_T','M_T'))

tumor$EMT_group=factor(tumor$EMT_group,levels = c('non_EMT','low_EMT','high_EMT'))

tumor$CellType_edit=factor(tumor$CellType_edit,levels = c('non_EMT_P_T','non_EMT_M_T',
'low_EMT_P_T','low_EMT_M_T','high_EMT_P_T','high_EMT_M_T')

a=as.data.frame(AverageExpression(tumor,group.by = 'CellType_edit',features = c(epi_mk,mes_mk,ent_reg),
                                  assays = 'RNA',slot = 'data')$RNA)
pheatmap(a,cluster_rows = F,cluster_cols = F,border_color = 'black',scale = 'row',gaps_col = c(2,4),gaps_row = c(length(epi_mk),length(epi_mk)+length(mes_mk)))

# ## Fig 4C##

cellchat.NM_EMT <- readRDS("cellchat_peri_Epi_Group_split_NM.RDS")
cellchat.M_EMT <- readRDS("cellchat_peri_Epi_Group_split_M.RDS")

object.list <- list(NM = cellchat.NM_EMT, M = cellchat.M_EMT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#> Do heatmap based on a merged object
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# ## Fig S5H ##

obj <- readRDS("bin50_G100M20K20R0.6D30_nosd_Cluster.rds")
DefaultAssay(obj) <- "Spatial";obj[["SCT"]] <- NULL
obj_mali <- readRDS("T_0_bin50Maligant.rds")
cell_all <- colnames(obj)
cell_mali <- colnames(obj_mali)
cell_other <- cell_all[!(cell_all %in% cell_mali)]
obj_other <- subset(obj,cells = cell_other)
cd=obj_other@images$slice1@coordinates
head(cd);dim(cd)
obj_other@meta.data$x <- cd$col
obj_other@meta.data$y <- cd$row
head(obj_other)
obj;obj_mali;obj_other
EMT_score_marker <- list(EMT_score = c('ABI3BP', 'ACTA2', 'ADAM12', 'ANPEP', 'APLP1', 'AREG', 'BASP1', 'BDNF', 'BGN', 'BMP1', 'CADM1', 'CALD1', 'CALU', 'CAP2', 'CAPG', 'CCN1', 'CCN2', 'CD44', 'CD59', 'CDH11', 'CDH2', 'CDH6', 'COL11A1', 'COL12A1', 'COL16A1', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A2', 'COL6A3', 'COL7A1', 'COL8A2', 'COLGALT1', 'COMP', 'COPA', 'CRLF1', 'CTHRC1', 'CXCL1', 'CXCL12', 'CXCL6', 'CXCL8', 'DAB2', 'DCN', 'DKK1', 'DPYSL3', 'DST', 'ECM1', 'ECM2', 'EDIL3', 'EFEMP2', 'ELN', 'EMP3', 'ENO2', 'FAP', 'FAS', 'FBLN1', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT2', 'FGF2', 'FLNA', 'FMOD', 'FN1', 'FOXC2', 'FSTL1', 'FSTL3', 'FUCA1', 'FZD8', 'GADD45A', 'GADD45B', 'GAS1', 'GEM', 'GJA1', 'GLIPR1', 'GPC1', 'GPX7', 'GREM1', 'HTRA1', 'ID2', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IL15', 'IL32', 'IL6', 'INHBA', 'ITGA2', 'ITGA5', 'ITGAV', 'ITGB1', 'ITGB3', 'ITGB5', 'JUN', 'LAMA1', 'LAMA2', 'LAMA3', 'LAMC1', 'LAMC2', 'LGALS1', 'LOX', 'LOXL1', 'LOXL2', 'LRP1', 'LRRC15', 'LUM', 'MAGEE1', 'MATN2', 'MATN3', 'MCM7', 'MEST', 'MFAP5', 'MGP', 'MMP1', 'MMP14', 'MMP2', 'MMP3', 'MSX1', 'MXRA5', 'MYL9', 'MYLK', 'NID2', 'NNMT', 'NOTCH2', 'NT5E', 'NTM', 'OXTR', 'P3H1', 'PCOLCE', 'PCOLCE2', 'PDGFRB', 'PDLIM4', 'PFN2', 'PLAUR', 'PLOD1', 'PLOD2', 'PLOD3', 'PMEPA1', 'PMP22', 'POSTN', 'PPIB', 'PRRX1', 'PRSS2', 'PTHLH', 'PTX3', 'PVR', 'QSOX1', 'RGS4', 'RHOB', 'SAT1', 'SCG2', 'SDC1', 'SDC4', 'SERPINE1', 'SERPINE2', 'SERPINH1', 'SFRP1', 'SFRP4', 'SGCB', 'SGCD', 'SGCG', 'SLC6A8', 'SLIT2', 'SLIT3', 'SNAI2', 'SNTB1', 'SPARC', 'SPOCK1', 'SPP1', 'TAGLN', 'TFPI2', 'TGFB1', 'TGFBI', 'TGFBR3', 'TGM2', 'THBS1', 'THBS2', 'THY1', 'TIMP1', 'TIMP3', 'TNC', 'TNFAIP3', 'TNFRSF11B', 'TNFRSF12A', 'TPM1', 'TPM2', 'TPM4', 'VCAM1', 'VCAN', 'VEGFA', 'VEGFC', 'VIM', 'WIPF1', 'WNT5A'))
obj_mali <- AddModuleScore(object = obj_mali,features = EMT_score_marker,nbin = 24)
colnames(obj_mali@meta.data)[(ncol(obj_mali@meta.data)-length(EMT_score_marker)+1):ncol(obj_mali@meta.data)] <- names(EMT_score_marker)
head(obj_mali)
peri_score_marker <- list(peri_score =  c('MCAM','RGS5',"ACTA2",'ADGRF5','EGFL6','SPARC','COL1A1','THY1','IFI27','ISG15','IFI6','COL3A1','TIMP1'))
obj_other <- AddModuleScore(object = obj_other,features = peri_score_marker,nbin = 24)
colnames(obj_other@meta.data)[(ncol(obj_other@meta.data)-length(peri_score_marker)+1):ncol(obj_other@meta.data)] <- names(peri_score_marker)
head(obj_other)
SpatialFeaturePlot(obj_mali, features = "EMT_score",pt.size.factor = 1,stroke=NA) + ggtitle(paste0(pre ,"_mali", 'EMT_gene_set_score'))
SpatialFeaturePlot(obj_other, features = "peri_score",pt.size.factor = 1,stroke=NA) + ggtitle(paste0(pre ,"_non_mali", "peri_score"))

plot(density(obj_mali@meta.data$EMT_score))
plot(density(obj_other@meta.data$peri_score))

# ## Fig 4D ##

cd_mali = obj_mali@images$slice1@coordinates
cd_mali$id=paste0(cd_mali$imagerow,'-',cd_mali$imagecol)
head(cd_mali)
cd <- cd_mali[cd_mali$row > 10 & cd_mali$row < 60 & cd_mali$col > 175 & cd_mali$col < 245, ]
obj_mali_1 <- subset(obj_mali,cells = interset)
obj_mali_1
head(obj_mali_1)

cd_peri = obj_other@images$slice1@coordinates
cd_peri$id=paste0(cd_peri$imagerow,'-',cd_peri$imagecol)
head(cd_peri)
cd <- cd_peri[cd_peri$row > 10 & cd_peri$row < 60 & cd_peri$col > 175 & cd_peri$col < 245, ]
obj_peri_1 <- subset(obj_other,cells = interset) 
obj_peri_1
head(obj_peri_1)

EMT_score <- obj_mali_1@meta.data$EMT_score
df_EMT <- data.frame(EMT_score)
row_names <- rownames(obj_mali_1@meta.data)
rownames(df_EMT) <- row_names
df_EMT$x <- obj_mali_1@meta.data$x
df_EMT$y <- obj_mali_1@meta.data$y
df_EMT$celltype <- "mali"
colnames(df_EMT)[colnames(df_EMT) == "EMT_score"] <- "score"
head(df_EMT);dim(df_EMT)

peri_score <- obj_peri_1@meta.data$peri_score
df_peri <- data.frame(peri_score)
row_names <- rownames(obj_peri_1@meta.data)
rownames(df_peri) <- row_names
df_peri$x <- obj_peri_1@meta.data$x
df_peri$y <- obj_peri_1@meta.data$y
df_peri$celltype <- "other"
colnames(df_peri)[colnames(df_peri) == "peri_score"] <- "score"
head(df_peri);dim(df_peri)

df_all <- rbind(df_peri,df_EMT)
df_all$id=paste0(df_all$x,'-',df_all$y)
head(df_all);dim(df_all)
table(df_all$celltype)

p4 <- ggplot() +
  geom_point(data=df_EMT,aes(x = x, y = -y, color = EMT_score),size=2) +  scale_color_continuous(low="yellow",high="red")+
  new_scale_color() +
  geom_point(data=df_peri,aes(x = x, y = -y, color = peri_score),size=2) +scale_color_continuous(low="green",high="blue")
p4

# dismat是一个矩阵，非肿瘤细胞与边界出肿瘤细胞的距离
dismat=readRDS('distance_malignantcells_normalcells.RDS')

p1 <- c()
for (i in 1:ncol(dismat)) {
  if (any(dismat[, i] >= 0 & dismat[, i] <= 5)){
    p1 <- c(p1, colnames(dismat)[i])
  }
}
p1 <- as.list(p1)
length(p1);head(p1);class(p1)

p2 <- c()
for (i in 1:ncol(dismat)) {
  if (any(dismat[, i] > 5 & dismat[, i] <= 10) && !any(dismat[, i] >= 0 & dismat[, i] <= 5)) {
    p2 <- c(p2, colnames(dismat)[i])
  }
}
p2 <- as.list(p2)
length(p2);head(p2);class(p2)


p3 <- c()
for (i in 1:ncol(dismat)) {
  if (any(dismat[, i] >10 & dismat[, i] <= 15) && !any(dismat[, i] >= 0 & dismat[, i] <= 10)) {
    p3 <- c(p3, colnames(dismat)[i])
  }
}
p3 <- as.list(p3)
length(p3);head(p3);class(p3)

p4 <- c()
for (i in 1:ncol(dismat)) {
  if (any(dismat[, i] >15) && !any(dismat[, i] <= 15)) {
    p4 <- c(p4, colnames(dismat)[i])
  }
}

p4 <- as.list(p4)
length(p4);head(p4);class(p4)

rownames(df_peri) <- sub("^....", "", rownames(df_peri)) 

df_peri_p1 <- df_peri[rownames(df_peri) %in% p1, ]
df_peri_p2 <- df_peri[rownames(df_peri) %in% p2, ]
df_peri_p3 <- df_peri[rownames(df_peri) %in% p3, ]
df_peri_p4 <- df_peri[rownames(df_peri) %in% p4, ]
new_df_1$level <- "b0" 
# new_df_1$level <- NULL
df_peri_p1$level <- "level_1"
df_peri_p2$level <- "level_2"
df_peri_p3$level <- "level_3"
df_peri_p4$level <- "level_4"
dim(new_df_1) 
dim(df_peri) 
dim(df_peri_p1)
dim(df_peri_p2)
dim(df_peri_p3)
dim(df_peri_p4)

df <- rbind(new_df_1,df_peri_p1,df_peri_p2,df_peri_p3,df_peri_p4) 
head(df)
table(df$level)

p6 <- ggplot(data = df, aes(x = x, y = -y, color = level)) +  
  geom_point(size = 2) +  
  scale_color_manual(values = c(  
    "b0" = "blue",  
    "level_1" = "red",  
    "level_2" = "orange",  
    "level_3" = "yellow",  
    "level_4" = "green",  
    "level_5" = "purple"  
  )) +  
  labs(color = "Level") +  
  theme_minimal()  
p6

obj_1 <- subset(obj,cells = rownames(df));obj_1;head(obj_1) 
df$score <- obj_1@meta.data$score
df_1 <- df[df$level != "b0", ]
head(df_1);table(df_1$level)
df_avg <- aggregate(score ~ level, data = df_1, FUN = mean)

ggplot(data = df_avg, mapping = aes(x = level, y = score, group = 1)) +
  geom_line() +
  geom_point() +
  labs(x = "Group", y = "Average Score", title = "Average Score by Group")
write.table(df_1,file = paste0("data_for_",pre,"_line_graph.txt"),sep = "\t",quote = F)

score_NM <- read.delim("data_for_NM_bin50_line_graph.txt",sep="\t",header=TRUE,check.names=FALSE)
score_M <- read.delim("data_for_M_bin50_line_graph.txt",sep="\t",header=TRUE,check.names=FALSE)
score_NM$score_scaled <- 2 * (score_NM$score - min(score_NM$score)) / (max(score_NM$score) - min(score_NM$score)) - 1 
score_M$score_scaled <- 2 * (score_M$score - min(score_M$score)) / (max(score_M$score) - min(score_M$score)) - 1 
df <- rbind(score_NM,score_M)
head(df)

score_summary <- score_NM %>%  
  group_by(level) %>%  
  summarise(  
    mean_score = mean(score_scaled, na.rm = TRUE),  
    se_score = sd(score_scaled, na.rm = TRUE) / sqrt(n())  
  )  

score_M_summary <- score_M %>%  
  group_by(level) %>%  
  summarise(  
    mean_score = mean(score_scaled, na.rm = TRUE),  
    se_score = sd(score_scaled, na.rm = TRUE) / sqrt(n())  
  )  
  
score_summary$level <- factor(score_summary$level, levels = c("level_4", "level_3", "level_2", "level_1"))  
score_M_summary$level <- factor(score_M_summary$level, levels = c("level_4", "level_3", "level_2", "level_1"))  

combined_data <- rbind(score_summary, score_M_summary)  

p <- ggplot(combined_data, aes(x = level, y = mean_score, color = type)) +    
  geom_point() +  
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), width = 0.2) +  
  
  geom_line(aes(group = type)) +   
  
  labs(title = "non_mali GPR116 score",  
       x = "Level",  
       y = "scaled_score",  
       color = "Pacient") +
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5,size = 16), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p <- p+theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())
p
