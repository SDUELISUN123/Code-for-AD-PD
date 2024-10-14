
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
#library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(monocle)
library(ggpubr)


cors<-ggsci::pal_igv("default", alpha =0.25)(51) #加了透明度


plan("multicore", workers = 8) ###set the compute core
#options(future.globals.maxSize = 14000 * 1024^2)






############################ 1 读取数据 #######################################
library(data.table)
library(tidyverse)
sqy<-fread("IPDCO_hg_midbrain_UMI.tsv",header = T)
ens<-read.table("ens.txt",header = T)
symbol<-read.table("symbol.txt",header = T)
ens_symbol<-left_join(ens,symbol,by="row")
ens_symbol_intersect<-intersect(ens_symbol$row,rownames(sqy))
row_name<-subset(ens_symbol,ens_symbol$row%in%ens_symbol_intersect)
merge<-cbind(row_name$gene.y,sqy)
colnames(merge)[1]<-"ID"
head(merge)[1:5,1:5]
merge<-subset(merge,!duplicated(merge$ID))
merge<-subset(merge,!merge$ID%in%NA)
merge_data_frame<-as.data.frame(merge)
head(merge_data_frame)[1:5,1:5]
rownames(merge_data_frame)<-merge_data_frame$ID
merge_data_frame<-merge_data_frame[,-1]
saveRDS(merge_data_frame,"raw_exp_mtx.rds")
rm(list=ls())
gc()

sqy<-readRDS("raw_exp_mtx.rds")
head(sqy)[1:5,1:5]
clinical<-read.table("IPDCO_hg_midbrain_cell.tsv",header = T,sep="\t")
patients<-clinical$patient
celltype_previous<-clinical$cell_ontology


############################ 2 标准化与PCA ################################


sqy=CreateSeuratObject(counts = sqy,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
gc()
sqy$patient<-patients
sqy$celltype_previous<-celltype_previous

library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
gc()

#sqy@meta.data[1:5,]

#使用PercentageFeatureSet函数计算线粒体基因的百分比
sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^MT-")#鼠源的换成mt


#质控
pdf("qc.pdf",width = 8,height = 4)
VlnPlot(object = sqy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "patient",cols = cors)  #ncol是指每行放几张图
dev.off()
sqy=subset(x = sqy, subset = nFeature_RNA > 50 & percent.mt < 10)    #对数据进行过滤

#测序深度的相关性图
plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5,cols = cors,group.by = "patient")+
  RotatedAxis()
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5,cols = cors,group.by = "patient")+
  RotatedAxis()
pdf("测序深度.pdf",width = 8,height = 4)
CombinePlots(plots = list(plot1, plot2))
dev.off()
rm(plot1,plot2)

#看细胞属于哪一期，并加到矩阵里
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sqy <- CellCycleScoring(sqy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sqy@meta.data[1:5,]
######标准化
sqy<-NormalizeData(sqy,verbose = T)   #标准化

####PCA
sqy<-FindVariableFeatures(sqy,selection.method = "vst", nfeatures = 4000)   #找前2000个差异显著的基因
sqy<-ScaleData(sqy,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = T) #去除线粒体基因和分裂期的影响
sqy<-RunPCA(sqy,verbose = T,npcs = 70)  #pca降维
ElbowPlot(sqy,ndims = 70)  #看拐点
#pc_select=25



############################ 3 UMAP降维 ################################



#UMAP/tsne
sqy <- FindNeighbors(object = sqy, dims = 1:25)       #计算邻接距离
sqy <- FindClusters(object = sqy, resolution = 0.8)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
#sqy <- RunUMAP(object = sqy, dims = 1:25) 
sqy<-RunTSNE(object = sqy, dims = 1:25)

sqy@active.ident<-sqy$seurat_clusters
pdf(file = "tsne_cluster.pdf",width=4,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors)+
  theme_bw()+
  NoLegend()
dev.off()

sqy@active.ident<-factor(sqy$celltype_previous)
pdf(file = "tsne_cluster_annotation_previous.pdf",width=4,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors,repel=T)+
  theme_bw()+
  NoLegend()
dev.off()

sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  "0"="Oligodendrocytes",
                  "1"="Oligodendrocytes",
                  "2"="Oligodendrocytes",
                  "3"="Oligodendrocytes",
                  "4"="Neurons",
                  "5"="OPCs",
                  "6"="Astrocytes",
                  "7"="Microglia",
                  "8"="Oligodendrocytes",
                  "9"="ECs",
                  "10"="Astrocytes",
                  "11"="Neurons",
                  "12"="Pericytes",
                  "13"="Microglia",
                  "14"="Astrocytes",
                  "15"="Neurons",
                  "16"="Ependymal",
                  "17"="Neurons",
                  "18"="Neurons",
                  "19"="Neurons",
                  "20"="Microglia",
                  "21"="Pericytes",
                  "22"="Microglia")
sqy$celltype<-sqy@active.ident
table(sqy$celltype)

saveRDS(sqy,"注释完毕.rds")

sqy@active.ident<-factor(sqy$celltype)
pdf(file = "tsne_cluster_annotation.pdf",width=4,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors,repel=T)+
  theme_bw()+
  NoLegend()
dev.off()

sqy@active.ident<-factor(sqy$patient)
sqy<-RenameIdents(sqy,
                  "C1"="Normal",
                  "C2"="Normal",
                  "C3"="Normal",
                  "C4"="Normal",
                  "C5"="Normal",
                  "C6"="Normal",
                  "PD1"="PD",
                  "PD2"="PD",
                  "PD3"="PD",
                  "PD4"="PD",
                  "PD5"="PD")
sqy$group<-sqy@active.ident

sqy@active.ident<-factor(sqy$celltype)
pdf(file = "tsne_cluster_annotation_group.pdf",width=8,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors,repel=T,split.by="group")+
  theme_bw()+
  NoLegend()
dev.off()

sqy@active.ident<-factor(sqy$Phase)
pdf(file = "tsne_cluster_Phase.pdf",width=5,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors,repel=T)+
  theme_bw()
dev.off()







############################ 5 看每种细胞在每个样本中的比例 ################################

##### 1 堆砌柱状图
library(ggplot2)

sqy@active.ident<-factor(sqy$celltype)
tab<-table(Idents(sqy),sqy$patient)
tab#每个样本中各种细胞数目
write.csv(tab,file = "每个样本各细胞数目.csv")
tab<-as.data.frame(tab)#改成绘图所需要的格式
tab


# 1.2绘制百分比的柱状图
tab<-table(Idents(sqy),sqy$patient)
tab
tab<-prop.table(tab,2)*100
tab
tab<-as.data.frame(tab)
p<-ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='group',y='Cell Proportion (%)')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ scale_fill_manual(values=cors)

pdf("百分比堆砌柱状图，每个样本.pdf",width = 6,height = 4)
p
dev.off()

tab<-table(Idents(sqy),sqy$group)
tab
tab<-prop.table(tab,2)*100
tab
tab<-as.data.frame(tab)
p<-ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='group',y='Cell Proportion (%)')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ scale_fill_manual(values=cors)

pdf("百分比堆砌柱状图，con和pd.pdf",width = 6,height = 4)
p
dev.off()




##### 6 提取神经元 ####
pdf("feat_GFRA2_exp.pdf",width = 4,height = 4)
FeaturePlot(sqy,features = "GFRA2",reduction = "tsne",cols = c("lightgrey","red"))+
  theme_bw()+
  NoLegend()
dev.off()

pdf("Vln_GFRA2_exp.pdf",width = 4,height = 4)
VlnPlot(sqy,"GFRA2",cols = cors,group.by = "celltype",pt.size = 0)+
  theme_bw()+
  geom_boxplot(outlier.size = -1)+
  NoLegend()+
  stat_compare_means(label.x = 1.5,label.y = 2.5)+
  RotatedAxis()
dev.off()

Neurons<-sqy[,sqy$celltype%in%"Neurons"]

#UMAP/tsne
Neurons <- FindNeighbors(object = Neurons, dims = 1:25)       #计算邻接距离
Neurons <- FindClusters(object = Neurons, resolution = 2)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
#sqy <- RunUMAP(object = sqy, dims = 1:25) 
Neurons<-RunTSNE(object = Neurons, dims = 1:25)

Neurons@active.ident<-factor(Neurons$seurat_clusters)
pdf(file = "tsne_cluster_neurons.pdf",width=4,height = 4)
TSNEPlot(object = Neurons, label = TRUE,cols=cors)+
  theme_bw()+
  NoLegend()
dev.off()

Neurons@active.ident<-Neurons$seurat_clusters
Neurons<-RenameIdents(Neurons,
                      "0"="Excitatory",
                      "1"="Inhibitory",
                      "2"="Excitatory",
                      "3"="Inhibitory",
                      "4"="Inhibitory",
                      "5"="Excitatory",
                      "6"="Excitatory",
                      "7"="Inhibitory",
                      "8"="Inhibitory",
                      "9"="Excitatory",
                      "10"="GABA",
                      "11"="Inhibitory",
                      "12"="Excitatory",
                      "13"='Inhibitory',
                      '14'='GABA',
                      '15'='Excitatory',
                      '16'='Excitatory',
                      '17'='Inhibitory',
                      '18'='GABA',
                      '19'='Glutamatergic',
                      '20'='Excitatory',
                      '21'='CADPS2+',
                      '22'='DaNs',
                      '23'='Excitatory',
                      '24'='Inhibitory',
                      '25'='Excitatory')

Neurons$neuron_type<-Neurons@active.ident
saveRDS(Neurons,"注释完毕_神经元.rds")

Neurons@active.ident<-Neurons$neuron_type
pdf(file = "tsne_cluster_neurons_annotation.pdf",width=4,height = 4)
TSNEPlot(object = Neurons, label = TRUE,cols=cors)+
  theme_bw()+
  NoLegend()
dev.off()

pdf("dot_不同神经元的marker表达.pdf",width = 6,height = 4)
DotPlot(Neurons,
        features = c("MAP2","SCN2A",'SLC17A6','GAD2','GRIK1','SLC17A7','CADPS2','TH','SLC6A3'),
        group.by = "neuron_type",
        cols = c("white","black"))+
  RotatedAxis()
dev.off()

pdf("fea_GFRA2在神经元表达.pdf",width = 4,height = 4)
FeaturePlot(Neurons,features = "GFRA2",reduction = "tsne",label = T)+
  theme_bw()+
  NoLegend()
dev.off()

pdf("Vln_GFRA2在神经元表达.pdf",width = 4,height = 4)
VlnPlot(Neurons,"GFRA2",cols = cors,group.by = "neuron_type",pt.size = 0)+
  geom_boxplot(outlier.size = 0.1)+
  NoLegend()+
  stat_compare_means(label.x = 1.5,label.y = 2.3)
dev.off()

##### 7 反卷积 #####
dir.create("cibersort")
setwd("cibersort")
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library('e1071')
library(future)
library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

plan("multicore", workers = 5)
#options(future.globals.maxSize = 14000 * 1024^2)

#导入单细胞数据,形成系数表
Idents(sqy) <- sqy$celltype#细胞类型是哪一列
X <- AverageExpression(sqy)[[1]]
Y<-CreateSeuratObject(counts = X,project = "seurat", min.cells=1, min.features=1)
Y<- FindVariableFeatures(object = Y, selection.method = "vst", nfeatures = 600)#调整反卷积要用的高变基因数目
Y<-ScaleData(Y) 
features<-Y@commands[["ScaleData.RNA"]]@params[["features"]]
X<-X[features,]
symbol<-rownames(X)
X<-cbind("ID"=symbol,X)
write.table(X,"sig.txt",sep = "\t",col.names = T,row.names = F,quote = F)

# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
results <- cibersort(sig_matrix = "sig.txt", mixture_file = "merge.normalize.txt")
results<-results[,1:(ncol(results)-3)]
write.csv(results,"result.csv",quote = F,row.names = T)

#读取免疫细胞浸润文件
rt=read.csv("result.csv", header=T, sep=",", check.names=F, row.names=1)
rt<-rt[,1:(ncol(rt)-3)]

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

#绘制柱状图
pdf(file="barplot.pdf", width=13, height=7.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.5)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="#0088FF")
text(a1[conNum]/2,-0.035,"Normal",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#FF5555")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"PD",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1)
dev.off()

###绘制箱线图
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#绘制箱线图
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=c("#0088FF", "#FF5555"))+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()

#看GFRA2和细胞的关系
GFRA2_exp<-read.table("merge.normalize.txt",header = T,sep = "\t",row.names = 1)
GFRA2_exp<-GFRA2_exp["GFRA2",]
cibersort_result<-read.csv("result.csv",sep = ",",header = T,row.names = 1)
cibersort_result<-t(cibersort_result)
mydata<-rbind(GFRA2_exp,cibersort_result)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(mydata))
mydata2<-rbind("Type"=Type,mydata)
mydata2<-t(mydata2)
mydata2<-as.data.frame(mydata2)
mydata2<-subset(mydata2,mydata2$Type=="Treat")
mydata3<-as.data.frame(t(mydata2))
mydata3<-mydata3[-1,]
rm(mydata,mydata2,Type)
write.csv(mydata3,"cor_input.csv",row.names = T)

#install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)#加载包
mydata3<-read.csv("cor_input.csv",sep = ",",header = T,row.names = 1)
mydata3<-t(mydata3)
head(mydata3)#列名为基因名，行名为样本名
pdf("GFRA2和细胞丰度相关性.pdf",width = 5,height = 5)
chart.Correlation(mydata3, histogram=TRUE, pch=19)
dev.off()

#### 8 看病人GFRA2的高低,然后分析他们的神经元的激活功能 ####
setwd("~/sqy/24.3.6 pd_macaca/GSE157783")

pdf("Vln_GFRA2_exp_每个患者之间.pdf",width = 4,height = 3)
VlnPlot(Neurons[,Neurons$group%in%"PD"],"GFRA2",cols = cors,group.by = "patient",sort=T,pt.size = 0)+
  geom_boxplot(outlier.size = 0.1)+
  NoLegend()+
  stat_compare_means(label.x = 1.2,label.y = 1.75)
dev.off()

Neurons@active.ident<-factor(Neurons$patient)
Neurons<-RenameIdents(Neurons,
                      "PD1"="Low",
                      "PD2"="High",
                      "PD3"="Low",
                      "PD4"="High",
                      "PD5"="High")
Neurons$GFRA2_group<-Neurons@active.ident

Neurons@active.ident<-Neurons$neuron_type
pdf("tsne_高低GFRA2_患者之间的split.pdf",width = 6,height = 4)
TSNEPlot(object = Neurons[,Neurons$group%in%"PD"], label = TRUE,cols=cors,split.by="GFRA2_group")+
  theme_bw()+
  NoLegend()
dev.off()

#分析两种患者神经元之间的差异
dir.create("分析高低GFRA2患者之间神经元差异基因")
setwd("分析高低GFRA2患者之间神经元差异基因")
Neurons@active.ident<-Neurons$GFRA2_group
logFCfilter=0
adjPvalFilter=0.05
sqy.markers <- FindAllMarkers(object = Neurons[,Neurons$group%in%"PD"],
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="cluster_markers.csv")
top10 <- sqy.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "cluster_markers_top10.csv")
write.csv(sqy.markers,"all_markers.csv")

head(sig.markers)


#vol
#三列分别是“基因名”，“logFC”，“pvalue”
cut_off_pvalue = 0.05
cut_off_logFC = 0.5

vol_input<-subset(sqy.markers,sqy.markers$cluster%in%"High")
vol_input<-vol_input[,c("gene","avg_log2FC","p_val_adj")]
colnames(vol_input)<-c("Gene","logFC","Pvalue")
vol_input$logFC<-as.numeric(vol_input$logFC)
vol_input$Pvalue<-as.numeric(vol_input$Pvalue)


vol_input$change = ifelse(vol_input$Pvalue < cut_off_pvalue & abs(vol_input$logFC) >= cut_off_logFC, 
                          ifelse(vol_input$logFC> cut_off_logFC ,'Up','Down'),
                          'Not.sig')
vol_input$nega_logP=-log10(as.numeric(vol_input$Pvalue))
vol_input<-subset(vol_input,vol_input$change!="NA")

pdf(file = "vol_input.pdf",height = 3,width = 3)
ggplot(
  #设置数据
  vol_input, 
  aes(x = logFC, 
      y = nega_logP, 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-as.numeric(cut_off_logFC),as.numeric(cut_off_logFC)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2FC",
       y="-log10 (adj.Pvalue)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+
  NoLegend()
dev.off()

up_down_in_high=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>0.5 & as.numeric(as.vector(sqy.markers$p_val_adj))<0.05),]
up_down_in_high<-subset(up_down_in_high,up_down_in_high$cluster%in%"High")
up_in_high<-subset(up_down_in_high,up_down_in_high$avg_log2FC>0)
down_in_high<-subset(up_down_in_high,up_down_in_high$avg_log2FC<0)
up_in_high<-up_in_high$gene
down_in_high<-down_in_high$gene

up_down_in_high_genes<-data.frame("gene"=c(up_in_high,down_in_high),"Type"=c(rep("Up_in_high",length(up_in_high)),rep("Down_in_high",length(down_in_high))))
write.csv(up_down_in_high_genes,"up_down_in_high_genes.csv",quote=F,row.names = F)
#save.image("24.3.7.RData")


###富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #p值的过滤条件
adjPvalFilter=1        #矫正后p值的过滤条件

#定义颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
  colorSel="pvalue"
}


#+++++++++++++++++++++++++++up_in_high_go+++++++++++++++++++++++
dir.create("up_in_high_go")
setwd("up_in_high_go")

rt=up_in_high     #读取差异分析的结果文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(rt)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt[,2]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<adjPvalFilter),]
#输出显著富集分析的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#柱状图
pdf(file="barplot.pdf", width=8, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


###绘制GO圈图
ontology.col=c("#00CC33FF", "#FFC20AFF", "#CC33FFFF")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="GO.circlize.pdf", width=10, height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制富集显著性pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

#+++++++++++++++++++++++++++down_in_high_go+++++++++++++++++++++++
setwd("~/sqy/24.3.6 pd_macaca/GSE157783/分析高低GFRA2患者之间神经元差异基因")
dir.create("down_in_high_go")
setwd("down_in_high_go")

rt=down_in_high     #读取差异分析的结果文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(rt)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt[,2]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<adjPvalFilter),]
#输出显著富集分析的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#柱状图
pdf(file="barplot.pdf", width=8, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


###绘制GO圈图
ontology.col=c("#00CC33FF", "#FFC20AFF", "#CC33FFFF")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="GO.circlize.pdf", width=10, height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制富集显著性pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


#+++++++++++++++++++++++++++up_in_high_kegg+++++++++++++++++++++++
dir.create("up_in_high_kegg")
setwd("up_in_high_kegg")

rt=up_in_high  #读取基因列表文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(rt)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
colnames(rt)[1]="gene"
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt[,2]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt[1][match(strsplit(x,"/")[[1]],as.character(rt[,2]))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
#输出显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#设置展示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=9, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#+++++++++++++++++++++++++++down_in_high_kegg+++++++++++++++++++++++
setwd("~/sqy/24.3.6 pd_macaca/GSE157783/分析高低GFRA2患者之间神经元差异基因")
dir.create("down_in_high_kegg")
setwd("down_in_high_kegg")

rt=down_in_high  #读取基因列表文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(rt)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
colnames(rt)[1]="gene"
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt[,2]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt[1][match(strsplit(x,"/")[[1]],as.character(rt[,2]))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
#输出显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#设置展示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=9, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

setwd("~/sqy/24.3.6 pd_macaca/GSE157783")
save.image("24.3.7.RData")


################# 9 GSEA #############3
# 1 提取分析的通路
library(devtools)
#install_github("arc85/singleseqgset")##https://arc85.github.io/singleseqgset/articles/singleseqgset.html
library(singleseqgset)
library(dplyr)
#?msigdbr
msigdbr_collections() #查找能用的基因集的category和subcategory,在下一句话中更改gmt文件

h.human <- msigdbr(species="Homo sapiens",category="C2",subcategory="CP:KEGG")  #可更改gmt文件
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
#h.names  #查看所有选中的通路
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])}
#循环计算
#h.sets #查看每个通路下的基因


# 2 计算每类细胞的通路得分，并绘制热图
expr.mat=Neurons[,Neurons$group%in%"PD"]@assays$RNA@data
cluster.ids=Neurons[,Neurons$group%in%"PD"]$GFRA2_group
#rm(sqy)
logfc.data <- logFC(cluster.ids=cluster.ids,expr.mat=expr.mat)  #大矩阵可能报错，内存不够
#logfc.data <- logFC(cluster.ids=sqy$Seurat_harmony[1:1000],expr.mat=expr.mat[,c(1:1000)])
#names(logfc.data)
gse.res <- wmw_gsea(expr.mat=expr.mat,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
#rm(logfc.data)

names(gse.res)
res.stats <- gse.res[["GSEA_statistics"]]
res.stats<-res.stats[,1:2]
res.stats[1:5,1:2]  #前五行，前五列

ac<-data.frame(cluster=factor(colnames(res.stats))) #变成作图需要的矩阵格式
#ac
rownames(ac)=colnames(res.stats)

ann_colors = list(
  cluster = c("Low"="blue","High"="red")
  #GeneClass = c("CD4+" = "#7570B3", "CD8+" = "#E7298A","NK" = "#66A61E")
) #更改颜色，复制前面的即可

#cors

#dev.off()
library(stringr)
rownames(res.stats)<-str_replace_all(rownames(res.stats),pattern = "KEGG_",replacement ="") #把"HALLMARK_"，"KEGG_"等前缀变成空格
rownames(res.stats)<-str_replace_all(rownames(res.stats),pattern = "_",replacement =" ") #把"_"变成空格
res.stats.sig<-subset(res.stats,abs(res.stats$Low)>1&abs(res.stats$High)>1)


pdf("GSEA_enrichment.pdf",width = 5.5,height = 5)
p<-pheatmap::pheatmap(res.stats.sig,
                      fontsize_row = 8,
                      annotation_col = ac,
                      annotation_legend = F,
                      annotation_colors = ann_colors,
                      cluster_rows = T,
                      cluster_cols = F,
                      scale = "column")
print(p)
dev.off()
#write.csv(res.stats,file="gsea.result.NES.keratinocyte.csv")
#sb的r报错的话直接导出自己画热图，爷nb爷就是这么任性
write.csv(gse.res$GSEA_statistics,"gsea.result.NES.csv")
write.csv(gse.res$GSEA_p_values,"gsea.result.pval.csv")
save.image("24.3.7.RData")


###### 9 寻找PD相关神经元 #####
dir.create("寻找pd相关神经元")
setwd("寻找pd相关神经元/")

Neurons@active.ident<-Neurons$seurat_cluster
pdf("tsne_neurons_group_group.pdf",width = 4,height = 3)
TSNEPlot(Neurons,group.by="group",cols=cors,label=F)+
  theme_bw()
dev.off()

pdf("fea_XIST_exp.pdf",width = 3,height = 3)
FeaturePlot(Neurons,features = "XIST",reduction = "tsne",cols = c("lightgrey","red"))+
  theme_bw()+
  NoLegend()
dev.off()

pdf("fea_GFRA2_exp.pdf",width = 3,height = 3)
FeaturePlot(Neurons,features = "GFRA2",reduction = "tsne",cols = c("lightgrey","red"))+
  theme_bw()+
  NoLegend()
dev.off()

V1<-VlnPlot(Neurons[,Neurons$neuron_type%in%c("Excitatory","DaNs")],
            features = "GFRA2",cols = cors,pt.size = -1,
            group.by = "seurat_clusters",sort=F)+
  geom_boxplot(outlier.size = -1)+
  stat_compare_means(label.x = 2,label.y = 2.25)+
  ggtitle(label = "Excitatory Neurons and DaNs")+
  xlab("")+
  ylab("GFRA2")+
  theme_bw()+
  NoLegend()
V2<-VlnPlot(Neurons[,Neurons$neuron_type%in%c("Excitatory","DaNs")],
            features = "XIST",cols = cors,pt.size = -1,
            group.by = "seurat_clusters",sort=F)+
  geom_boxplot(outlier.size = -1)+
  stat_compare_means(label.x = 2,label.y = 3.5)+
  ggtitle("")+
  ylab("XIST")+
  theme_bw()+
  NoLegend()
pdf("vln_兴奋神经元中GFRA2,XIST表达情况_重新聚类.pdf",width = 4.5,height =6)
V1/V2
dev.off()

##查找每个聚类的差异基因
logFCfilter=0.5
adjPvalFilter=0.05
Neurons@active.ident<-Neurons$seurat_clusters
sqy.markers <- FindAllMarkers(object = Neurons[,Neurons$neuron_type%in%c("Excitatory","DaNs")],
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="cluster_markers.csv")

top10 <- sqy.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "cluster_markers_top10.csv")

##### 7 反卷积 #####
dir.create("cibersort")
setwd("cibersort")
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library('e1071')
library(future)
library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(reshape)

plan("multicore", workers = 5)
#options(future.globals.maxSize = 14000 * 1024^2)

#导入单细胞数据,形成系数表
Idents(sqy) <- Neurons$seurat_clusters#细胞类型是哪一列
X <- AverageExpression(sqy)[[1]]
Y<-CreateSeuratObject(counts = X,project = "seurat", min.cells=1, min.features=1)
Y<- FindVariableFeatures(object = Y, selection.method = "vst", nfeatures = 600)#调整反卷积要用的高变基因数目
Y<-ScaleData(Y) 
features<-Y@commands[["ScaleData.RNA"]]@params[["features"]]
X<-X[features,]
symbol<-rownames(X)
X<-cbind("ID"=symbol,X)
write.table(X,"sig.txt",sep = "\t",col.names = T,row.names = F,quote = F)

# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
results <- cibersort(sig_matrix = "sig.txt", mixture_file = "merge.normalize.txt")
results<-results[,1:(ncol(results)-3)]
write.csv(results,"result.csv",quote = F,row.names = T)

#读取免疫细胞浸润文件
rt=read.csv("result.csv", header=T, sep=",", check.names=F, row.names=1)
rt<-rt[,1:(ncol(rt)-3)]

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

###绘制箱线图
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#绘制箱线图
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=c("#0088FF", "#FF5555"))+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()




######## 10 拟时序（神经元，失败） #######

dir.create("monocle")
setwd("monocle")
mmt<-Neurons

#把TAM，monocyte等sub细胞分类信息给R
mmt$cell_type_val<-Idents(mmt)
mmt@meta.data[1:5,]

#寻找高变基因，作为拟时序降维的基因
mmt<-FindVariableFeatures(mmt)
#>这一步，不光是可以使用FindVariableFeatures寻找多变基因
#>也可以用前面的用于pca降维的2000个多变基因作为降维基因
#>二者都要试试，看哪个更符合生物学进程和假说，就用哪个
#>生信需要反复调整参数，而不是去跑一个流程就行
#>?FindVariableFeatures

#提取数据
matrix<-as.matrix(mmt@assays$RNA@counts)  #提取最原始的count数据
dim(matrix) #看多少行多少列
matrix[1:500,1:10]


#基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))

#细胞注释
#oct[["cell_group"]]<-Idents(oct)
mmt@meta.data[1:5,]
sample_ann <- mmt@meta.data

#建立monocle研究对象
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)  #建立monocle研究对象，并进行归一化和质控
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)

#质控
#?detectGenes：某基因必须在多少个细胞里有表达，才算有表达，并删掉没表达的，继续缩小数据
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 10)) #必须至少在10个细胞里（num_cells_expressed）有表达，才纳入该基因,可以自己调
fData(sc_cds_2)[1:5,]
#saveRDS(sc_cds_2,file = "9.06sc_cds_2.rds") #save一下sc_cds_2防止后边报错，闪退等
#sc_cds_2<-readRDS(file = "sc_cds_2.rds")

#设置高变基因
###以monocle下的FindVariableFeatures的高变基因为ordering
ordering_genes<-mmt@assays$RNA@var.features
###>或者
###>以组之间差异基因为ordering gene
#diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
#                                      fullModelFormulaStr = "~cell_type_val",cores = 4,
#                                     verbose = T) #+num_genes_expressed+orig.ident
#ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2)) #这个p值需要进行摸索，摸索到基因数在1000-2000
###>或者选logfc的top200作为高变基因
###>ordering_genes<-markers %>% group_by(cluster) %>% top_n(n=200,wt=avg_logFC)

#拟时序降维
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
#plot_ordering_genes(sc_cds2)  #看降维基因表达水平图？
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")   #降到2维
sc_cds2 <- orderCells(sc_cds2)  #把cell的顺序排出来、
saveRDS(sc_cds2,"monocle_for_plot.rds")
#beepr::beep(1)

