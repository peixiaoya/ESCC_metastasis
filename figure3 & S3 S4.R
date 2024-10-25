fdir = ""
## I. load packages
library(SCENIC)
library(AUCell)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

## II. load data
peri <- readRDS("obj/Pericyte.RDS")
peri_merge <- readRDS("obj/Pericyte_merge.RDS")

## Fig 3A ##

regulonAUC<- readRDS("SCENIC_res/int/3.4_regulonAUC.Rds")
regulonAUC<- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# regulon numbers
dim(regulonAUC@assays@data@listData$AUC)
# 179 regulons

cellInfo<-readRDS("SCENIC_res/int/cellinfo.Rds")
cellInfo<-cellInfo[colnames(regulonAUC),]

c0<-colnames(regulonAUC)[which(cellInfo$celltype_plot=='SORBS2+PC')]
c1<-colnames(regulonAUC)[which(cellInfo$celltype_plot=='STEAP4+PC')]
c2<-colnames(regulonAUC)[which(cellInfo$celltype_plot=='GPR116+PC')]

getAUC<-getAUC(regulonAUC)
colnames(getAUC)<-colnames(regulonAUC)

AUC_mean=data.frame()
for(x in rownames(regulonAUC@assays@data@listData$AUC)){
  auc_c0<-mean(getAUC[x,c0]) 
  auc_c1<-mean(getAUC[x,c1])
  auc_c2<-mean(getAUC[x,c2])
  a=c(auc_c0,auc_c1,auc_c2)
  AUC_mean=rbind(AUC_mean,a) 
}
colnames(AUC_mean)=c('SORBS2+PC', 'STEAP4+PC','GPR116+PC' )
rownames(AUC_mean)<-rownames(getAUC)

pheatmap(AUC_mean,scale = 'row',cluster_cols = F,show_rownames = F)


## Fig 3B ##

# intersection
# 1. SCENIC result, regluons activated in GPR116+PC, filtering: p.adjust<0.05 & log2FC>0
# 2. pericyte DEGs, up DEGs in pericyte of metastatic patients, filtering: p.adjust<0.05 & log2FC>log2(1.5)
# 3. pericyte DEGs, up DEGs in GPR116+PC, filtering: p.adjust<0.05 & log2FC>log2(1.5)

# 1. 
p<-c()
logFC=c()
getAUC<-getAUC(regulonAUC)
colnames(getAUC)<-colnames(regulonAUC)

for(x in rownames(regulonAUC@assays@data@listData$AUC)){
  GPR116Neg<-getAUC[x,c(c0,c1)]
  GPR116Pos<-getAUC[x,c2]
  C1QNCM_p<-t.test(GPR116Neg,GPR116Pos)
  p<-c(p,C1QNCM_p$p.value)
  logFC= c(logFC,log2( mean(GPR116Pos) / mean(GPR116Neg) ))
}
p_adjust<-fdr<-p.adjust(p,'BH')

frame<-data.frame(regulon=rownames(regulonAUC@assays@data@listData$AUC),
                  p.value=p,p_adjust=p_adjust,log2FC=logFC
)

res1=frame$regulon[frame$p_adjust<0.05 & frame$log2FC>0]
res1=sapply(res1,function(x){
  if(length(grep('_',x))==0){y=strsplit(x," ")[[1]][1]}
  else{y=strsplit(x,"_")[[1]][1]}
  y
}) # 33 genes

# 2. 
mk=readRDS("peri/markers_tissue.Rds") # T-NM high
res2=mk$gene[mk$p_val_adj<0.05 & mk$avg_log2FC>log2(1.5)]  # 66 genes 

# 3.
mk2=readRDS("peri/markers_celltype.Rds") # GPR116+PC high
res3=mk2$gene[mk2$p_val_adj<0.05 & mk2$avg_log2FC>log2(1.5)]  # 282 genes 

intersect(intersect(res1,res2),res3) # # "SOX4"  "PRRX1" "HIF1A"


## Fig 3C ##

exprset=peri@assays$RNA@data
exprset=as.data.frame(t(exprset))

y=as.numeric(exprset[,'PRRX1'])
y2=as.numeric(exprset[,'HIF1A'])
y3=as.numeric(exprset[,'SOX4'])

colnames=colnames(exprset)
cor_data_df=data.frame(colnames)
for(i in 1:length(colnames)){
  test=cor.test(as.numeric(exprset[,i]),y,type='spearman')
  cor_data_df[i,2]=test$estimate
  cor_data_df[i,3]=test$p.value
}
names(cor_data_df)=c('symbol','correlation','pvalue')


cor_data_df2=data.frame(colnames)
for(i in 1:length(colnames)){
  test=cor.test(as.numeric(exprset[,i]),y2,type='spearman')
  cor_data_df2[i,2]=test$estimate
  cor_data_df2[i,3]=test$p.value
}
names(cor_data_df2)=c('symbol','correlation','pvalue')


cor_data_df3=data.frame(colnames)
for(i in 1:length(colnames)){
  test=cor.test(as.numeric(exprset[,i]),y3,type='spearman')
  cor_data_df3[i,2]=test$estimate
  cor_data_df3[i,3]=test$p.value
}
names(cor_data_df3)=c('symbol','correlation','pvalue')

cor_data_df$TF=rep('PRRX1','36601')
cor_data_df2$TF=rep('HIF1A','36601')
cor_data_df3$TF=rep('SOX4','36601')

genes=c('ADGRF5','EGFL6','SPARC','COL1A1','THY1','COL3A1','IFI27','TIMP1','TDO2')

bind=rbind(cor_data_df[cor_data_df$symbol %in% genes,],
           cor_data_df2[cor_data_df2$symbol %in% genes,],
           cor_data_df3[cor_data_df3$symbol %in% genes,])

bind_cor=bind[,-3]

mat_bind=matrix(bind_cor$correlation,3,length(genes),byrow = T)
rownames(mat_bind)=c('PRRX1','HIF1A','SOX4')
colnames(mat_bind)=bind_cor$symbol[1:length(genes)]
mat_bind=mat_bind[,genes]
pheatmap::pheatmap(mat_bind,cluster_rows = F,cluster_cols = F,border_color = 'white',scale = 'column')


## Fig S4 ##

addArchRThreads(threads = 12)
addArchRGenome("hg38")

# create Arrow Files
inputfiles=c(
  paste0(path,samples,'/atac_fragments.tsv.gz'
  ))
names(inputfiles)=samples
ArrowFiles <- createArrowFiles(
  inputFiles = inputfiles,
  sampleNames = names(inputfiles),
  outputNames= names(inputfiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Inferring Doublets with ArchR
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# Creating An ArchRProject
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "DownStream_analysis",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1


# Fig S4F 
# Plotting Sample Fragment Size Distribution
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

# Fig S4F 
# Plotting TSS enrichment profiles
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

# subset Pericyte

idxPass <- which(projHeme2$celltype == 'Pericyte')
cellsPass <- projHeme2$cellNames [idxPass ] 
projHeme4 <- projHeme2[cellsPass, ] 

projHeme4 <- addIterativeLSI(
  ArchRProj = projHeme4,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

projHeme4 <- addHarmony(
  ArchRProj = projHeme4,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

projHeme4 <- addClusters(
  input = projHeme4,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)
table(projHeme4$Clusters) 

projHeme4 <- addTSNE(
  ArchRProj = projHeme4, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

markersGS <- getMarkerFeatures(
  ArchRProj = projHeme4, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
A=getMarkers(markersGS)
C1=data.frame(A$C1$name,A$C1$Log2FC,A$C1$FDR,A$C1$MeanDiff)
C2=data.frame(A$C2$name,A$C2$Log2FC,A$C2$FDR,A$C2$MeanDiff)

genes=c('ADGRF5','HIGD1B','SPON2','MDK','NAT14',
        'SPEG','HRH2','MOB2','ECHS1','C1QTNF2','TUBB4B','MFAP4','BCAM'
        )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.8", 
  labelMarkers = genes,   
  transpose = TRUE,
  plotLog2FC = TRUE
) # Must use plotLog2FC = TRUE when ncol(seMarker) <= 2!
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")


# Fig 3H & S4H
# ArchRBrowser tracks

# 'ADGRF5','EGFL6','THY1','MUSTN1'
p=plotBrowserTrack(
  ArchRProj = projHeme4, 
  groupBy = "Clusters", 
  geneSymbol = c('ADGRF5'), # 'EGFL6','THY1','MUSTN1'
  upstream = 15000,
  downstream = 100000
)
grid::grid.newpage()
grid::grid.draw(p$ADGRF5)
grid::grid.draw(p$EGFL6)
grid::grid.draw(p$THY1)
grid::grid.draw(p$MUSTN1) 

# Fig 3J

# Motif Deviations
if("Motif" %ni% names(projHeme5@peakAnnotation)){
  projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}

projHeme5 <- addBgdPeaks(projHeme5,force=T)

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
data=   getVarDeviations(projHeme5, name = "MotifMatrix",plot = F)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

motifs <- c("MYB", "TCF7", "FOXP1", "NFATC1", "TRPS1", "BACH2",'PRRX1')
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

arkerMotifs <- grep("z:",markerMotifs, value = TRUE)
markerMotifs <- markerMotifs [markerMotifs %in% c("z:PRRX1_440")]
markerMotifs

p <- plotGroups(ArchRProj = projHeme5, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projHeme5)
)
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 8, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
