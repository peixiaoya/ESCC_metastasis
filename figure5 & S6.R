fdir = ""
## I. load packages
library("Seurat")
library('ggplot2')
library('ggpubr')
library('ggrepel')
library('dplyr')

## II. load data
escc = read_rds("obj/escc_merge.rds")
Pericyte <- readRDS("obj/Pericyte.RDS")

## Fig 5A ##

mk=FindMarkers(pericyte,ident.1 = 'GPR116+PC')

mk$difference <- mk$pct.1 - mk$pct.2
mk_sig <- mk[which(mk$p_val_adj<10^(-10) & abs(mk$avg_log2FC) >2 & (mk$difference>0.5 | mk$difference<(-0.3))),]
mk_sig$label <- rownames(mk_sig)

ggplot(mk, aes(x=difference, y=avg_log2FC)) + xlim(-1,1)+
  geom_point(size=2, color="grey60") + 
  geom_point(data=mk[which(mk$p_val_adj<0.05 & mk$avg_log2FC>1 & mk$difference>0.1),],
             aes(x=difference, y=avg_log2FC,size=avg_log2FC),
             color="#ffb549")+
  geom_point(data=mk[which(mk$p_val_adj<0.05 & mk$avg_log2FC< -1 & mk$difference< -0.1),],
             aes(x=difference, y=avg_log2FC,size=-(avg_log2FC)),
             color="#41b6e6")+
  geom_text_repel(data=mk_sig, aes(label=label),
                  color="black",fontface="italic",size=3.5,max.overlaps=300)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 12),
        axis.line = element_line(color = 'black', size = 0.6))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
  
## Fig 5B ##
VlnPlot(pericyte,features = 'EGFL6',pt.size = 0,group.by='celltype_plot')+geom_boxplot(width=0.3,fill='white')

## Fig 5D ##
VlnPlot(pericyte,features = 'EGFL6',pt.size = 0,group.by='Tissue')+geom_boxplot(width=0.3,fill='white')

## Fig 5I ##
data=read.table('Chip-EGFL6-GPR116PC.txt',header = F,sep = '\t')
ggplot(data,aes(x=GPR116pc,y=egfl6))+geom_smooth(method = "lm",fill='#efcab1',color='#fe654c')+
  geom_point(size=2,color='#fe654c')+
  theme_pubr()+stat_cor(method = "pearson")

## Fig 5K, S6G, S6I ##
Elisa_res_ESCC <- read.table("ESCC_EGFL6_ELISA.txt", sep = "\t", header = T)
Elisa_res_ESCC <- read.table("ESCC_EGFL6_ELISA.txt", sep = "\t", header = T)
Elisa_res_ESCC <- read.table("ESCC_EGFL6_ELISA.txt", sep = "\t", header = T)

# For each result
data_roc <- roc(Elisa_res$Group, Elisa_res$EGFL6)

p <- ggroc(data_roc, legacy.axes = TRUE, color = "#e64b35")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "darkgrey", linetype = 4)+
  theme_bw()+#ggtitle('')+
  # ggsci::scale_color_aaas()+
  # scale_colour_manual(values = c("red", "blue", "black", "green"))+
  # theme(legend.position="none")+
  annotate("text", x = 0.2, y = 0.7, hjust = "left", size = 5, label=paste("AUC = ", round(data_roc$auc, 3)))+
  theme(legend.title = element_blank())+
  theme(panel.grid = element_blank())
p

## Fig S6A ##
Dotplot(escc,group='subtype',features='EGFL6')

## Fig S6B ##

pheno<-read.csv('TCGA_ESCA_clinical.csv',row.names = 1,sep = ',',header = T)
colnames(pheno)[grep('pathologic_',colnames(pheno))]
colnames(pheno)[grep('clinical_',colnames(pheno))]
rownames(pheno)=pheno$submitter_id

pheno=pheno[-(which(is.na(pheno$ajcc_pathologic_n))),]
pheno=pheno[which(pheno$primary_diagnosis=='Squamous cell carcinoma, NOS'),]

expr_tpm<-read.table('E:/Need_data/ESCC_Public data/ESCA_TCGA/TCGA_ESCA_gene_expression_Tpm.tsv',header = T,sep = '\t')
expr=expr_tpm
rownames(expr)<-expr$ID

old=colnames(expr)
new=sapply(colnames(expr), FUN=function(x){
  a=strsplit(x,'[.]')[[1]]
  b=paste0(a[1],'-',a[2],"-",a[3])
  b
  })
colnames(expr)=new

state=sapply(old, FUN=function(x){
  a=strsplit(x,'[.]')[[1]]
  b=a[4]
  b
})

expr<-expr[,which(state=='01A')]

pheno$N_stage_merge<-sapply(pheno$ajcc_pathologic_n,FUN = function(x){
  if(x=='N0'){A='N0'}
  if(x %in% c('N1','N3','N2','NX')){A='N1+N2+N3'}
  if(x %in% c('')){A='NULL'}
  A
})
table(pheno$N_stage_merge)

gene='EGFL6'
exp=expr[gene,]
exp=as.data.frame(t(exp))
exp$N_stage_merge=pheno[rownames(exp),which(colnames(pheno)=='N_stage_merge')]
colnames(exp)[1]='Gene'
exp=exp[-(which(is.na(exp$N_stage_merge))),]
  
ggplot(exp,aes(x=N_stage_merge,y=Gene))+geom_boxplot(aes(fill=N_stage_merge),outlier.colour = "red")+theme_pubr()+
    stat_compare_means()+ggtitle(gene)+ geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5)


a=exp[which(exp$N_stage_merge=='N0'),1]
IQR_A=quantile(a)[4]-quantile(a)[2]
b=exp[which(exp$N_stage_merge=='N1+N2+N3'),1]
IQR_B=quantile(b)[4]-quantile(b)[2]

exp=exp[-which(exp$N_stage_merge=='N0' & (exp$Gene<quantile(a)[2]-1.5*IQR_A | exp$Gene>quantile(a)[4]+1.5*IQR_A)),]
#exp=exp[-which(exp$N_stage_merge=='N1+N2+N3' & (exp$Gene<quantile(b)[2]-1.5*IQR_B | exp$Gene>quantile(b)[4]+1.5*IQR_B)),]

ggplot(exp,aes(x=N_stage_merge,y=Gene))+geom_boxplot(aes(fill=N_stage_merge),outlier.colour = "red")+theme_pubr()+
  stat_compare_means()+ggtitle(gene)+ geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5)