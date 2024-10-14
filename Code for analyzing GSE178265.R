

###### 1 导入原来的数据并重新处理 #####
setwd("~/sqy/24.3.6 pd/GSE178265_MACACA")

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


cors<-ggsci::pal_igv(alpha = 0.5)(51)
cors2<-ggsci::pal_simpsons(alpha = 0.5)(11)


load("~/sqy/24.3.6 pd/GSE178265_MACACA/11.18.RData")
saveRDS(sqy,"Macaca_以前的数据.rds")

rm(list=ls())
sqy<-readRDS("Macaca_以前的数据.rds")

sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^MT-")#鼠源的换成mt

#质控
VlnPlot(object = sqy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")  #ncol是指每行放几张图
sqy=subset(x = sqy, subset = nFeature_RNA > 50 & percent.mt < 10)    #对数据进行过滤

#测序深度的相关性图
plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
rm(plot1,plot2)

#看细胞属于哪一期，并加到矩阵里
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sqy <- CellCycleScoring(sqy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sqy@meta.data[1:5,]
######标准化
sqy<-NormalizeData(sqy,verbose = T)   #标准化

####PCA
sqy<-FindVariableFeatures(sqy,selection.method = "vst", nfeatures = 2000)   #找前2000个差异显著的基因
sqy<-ScaleData(sqy,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = T) #去除线粒体基因和分裂期的影响
sqy<-RunPCA(sqy,verbose = T,npcs = 70)  #pca降维
ElbowPlot(sqy,ndims = 70)  #看拐点
#pc=20

#矫正前的降维图和vln
p1 <- DimPlot(object = sqy, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sqy, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1|p2

#矫正
library(harmony)
sqy<-RunHarmony(sqy,group.by.vars = c("orig.ident"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sqy, 'harmony')
#dim(harmony_embeddings)

#矫正后的降维图和vln
p3 <- DimPlot(object = sqy, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = sqy, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p3|p4

rm(p1,p2,p3,p4,harmony_embeddings,g2m.genes,s.genes)

#UMAP/tsne
sqy <- FindNeighbors(object = sqy, dims = 1:25)       #计算邻接距离
sqy <- FindClusters(object = sqy, resolution = 0.8)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
#sqy <- RunUMAP(object = sqy, dims = 1:15) 
sqy<-RunTSNE(object = sqy, dims = 1:25)

sqy@active.ident<-sqy$seurat_clusters
pdf(file = "tsne_cluster.pdf",width=2.5,height = 2.5)
TSNEPlot(object = sqy, label = TRUE,cols=cors)+
  theme_bw()+
  NoLegend()
dev.off()

pdf(file = "fea_GFRA2_exp.pdf",width=3,height = 3)
FeaturePlot(sqy,features = "GFRA2",label = T,reduction = "tsne")+
  theme_bw()+
  NoLegend()
dev.off()

sqy@active.ident<-sqy$marker_cluster
pdf(file = "Vln_GFRA2_exp.pdf",width = 4,height = 4)
VlnPlot(sqy,features = "GFRA2",cols = cors,pt.size = 0,sort=T)+
  geom_boxplot(outlier.size = 0.1)+
  stat_compare_means(label.x = 1.7,label.y = 2.1)+
  theme_bw()+
  NoLegend()+
  RotatedAxis()+
  ggtitle("GFRA2 High --------> GFRA2 Low")+
  ylab("GFRA2 Expression")+
  xlab("")
dev.off()

##查找每个聚类的差异基因
logFCfilter=0.5
adjPvalFilter=0.05
sqy.markers <- FindAllMarkers(object = sqy,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="cluster_markers.csv")

top <- sqy.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top,file = "cluster_markers_top.csv")

#绘制marker在各个cluster的热图
pdf(file = "doheatmap_cluster_top.pdf",width = 10, height = 10)
DoHeatmap(object = sqy, features = top$gene) + NoLegend()
dev.off()

sqy@active.ident<-sqy$seurat_clusters
#标记细胞群
sqy<-RenameIdents(sqy,
                  "0"="Low",
                  "1"="Low",
                  "2"="Low",
                  "3"="High",
                  "4"="Low",
                  "5"="High",
                  "6"="Low",
                  "7"="High",
                  "8"="Low",
                  "9"="High",
                  "10"="Low")
sqy$GFRA2_group<-sqy@active.ident

sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  "0"="0 CNTN5+",
                  '1'='1 STXBP6+',
                  '2'='2 TENM2+',
                  '3'='3 COBLL1+',
                  '4'='4 CALCR+',
                  '5'='5 ANKFN1+',
                  '6'='6 C8H8orf34+',
                  '7'='7 HS3ST4+',
                  '8'='8 PHACTR1',
                  '9'='9 KIAA1217+',
                  '10'='10 CHN2+')
sqy$marker_cluster<-sqy@active.ident
sqy@active.ident<-sqy$marker_cluster
pdf(file = "tsne_marker_cluster.pdf",width=4,height = 2.5)
TSNEPlot(object = sqy, label = FALSE,cols=cors)+
  theme_bw()
dev.off()

#多巴胺代谢相关基因
Dopa_meta<-c('MAOA','MAOB','COMT','TH','DDC','CACNA1A','CACNA1B','DRD2','SLC6A3','SLC18A1','SLC18A2')
all_genes<-c('GFRA2','DCX','MAP2',Dopa_meta)
sqy@active.ident<-sqy$GFRA2_group
pdf("fea_多巴胺代谢相关基因.pdf",width = 21,height = 7)
f1=FeaturePlot(sqy,features = all_genes[1],ncol = 1,reduction = "tsne",label = T,
            cols = c('lightgrey','#FFD700','orange','red'),
            pt.size = 0.3)+NoLegend()
f2=FeaturePlot(sqy,features = all_genes[2],ncol =1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f3=FeaturePlot(sqy,features = all_genes[3],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f4=FeaturePlot(sqy,features = all_genes[4],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f5=FeaturePlot(sqy,features = all_genes[5],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f6=FeaturePlot(sqy,features = all_genes[6],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f7=FeaturePlot(sqy,features = all_genes[7],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f8=FeaturePlot(sqy,features = all_genes[8],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f9=FeaturePlot(sqy,features = all_genes[9],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f10=FeaturePlot(sqy,features = all_genes[10],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f11=FeaturePlot(sqy,features = all_genes[11],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f12=FeaturePlot(sqy,features = all_genes[12],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f13=FeaturePlot(sqy,features = all_genes[13],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
f14=FeaturePlot(sqy,features = all_genes[14],ncol = 1,reduction = "tsne",label = T,
               cols = c('lightgrey','#FFD700','orange','red'),
               pt.size = 0.3)+NoLegend()
(f1|f2|f3|f4|f5|f6|f7)/(f8|f9|f10|f11|f12|f13|f14)
dev.off()

#计算多巴胺突触通路得分
library("KEGGREST") 
listDatabases()  
gs<-keggGet('hsa04728')#Dopaminergic synapse
#获取通路中gene信息 
gs[[1]]$GENE 
#查找所有基因 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
pathways <- genes[1:length(genes)%%3 ==2] 
pathways <- data.frame(pathways)  
pathways<-as.list(pathways)
sqy<-AddModuleScore(object = sqy,features = pathways,name = "Dopaminergic_synapse")#更改path名字
rm(pathways,i,gene)

#看GFRA2和多巴胺突触的关系
pdf(file = "cor_GFRA2和多巴胺突触的关系.pdf",width = 18,height = 7)
c1<-FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "Dopaminergic_synapse1",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("Dopaminergic_synapse")
  
c2=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "DCX",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("DCX")

c3=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "MAP2",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("MAP2")

c4=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "MAOA",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("MAOA")

c5=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "MAOB",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("MAOB")

c6=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "COMT",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("COMT")

c7=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "TH",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("TH")

c8=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "DDC",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("DDC")

c9=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "CACNA1A",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("CACNA1A")

c10=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "CACNA1B",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("CACNA1B")

c11=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "DRD2",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("DRD2")

c12=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "SLC6A3",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("SLC6A3")
  
c13=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "SLC18A1",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("SLC18A1")
  
c14=FeatureScatter(sqy,feature1 = "GFRA2",feature2 = "SLC18A2",group.by = "seurat_clusters",cols = cors)+
  stat_cor()+
  geom_smooth(method = lm)+
  NoLegend()+
  ggtitle("SLC18A2")

(c1|c4|c5|c6|c7|c8)/(c9|c10|c11|c12|c13|c14)
dev.off()

#细胞死亡
pdf("vln_高低GFRA2_caspase.pdf",width = 3,height = 4)
VlnPlot(sqy,features = "Caspase_cascade1",group.by = "GFRA2_group",cols = cors,pt.size = 0)+
  geom_boxplot()+
  stat_compare_means()+
  NoLegend()
dev.off()
pdf("vln_高低GFRA2_BCL2.pdf",width = 3,height = 4)
VlnPlot(sqy,features = "BCL2",group.by = "GFRA2_group",cols = cors,pt.size = 0)+
  geom_boxplot()+
  stat_compare_means()+
  NoLegend()
dev.off()
pdf("vln_高低GFRA2_自噬.pdf",width = 3,height = 4)
VlnPlot(sqy,features = "Autophagy",group.by = "GFRA2_group",cols = cors,pt.size = 0)+
  geom_boxplot()+
  stat_compare_means()+
  NoLegend()
dev.off()


#看GFRA2群之间的差异基因
sqy@active.ident<-factor(sqy$GFRA2_group)
logFCfilter=0
adjPvalFilter=1
sqy.markers <- FindAllMarkers(object = sqy,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="GFRA2_cluster_markers.csv")

top <- sqy.markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)
write.csv(top,file = "GFRA2_cluster_markers_top.csv")

library(ggplot2)
library(ggrepel)
cut_off_pvalue = 0.05
cut_off_logFC = 0.5
vol_input<-subset(sig.markers,sig.markers$cluster=="High")
vol_input<-vol_input[,c(7,2,5)]
colnames(vol_input)<-c("Gene","logFC","Pvalue")
vol_input$change = ifelse(vol_input$Pvalue < cut_off_pvalue & abs(vol_input$logFC) >= cut_off_logFC, 
                          ifelse(vol_input$logFC> cut_off_logFC ,'Up','Down'),
                          'Not.sig')
vol_input$nega_logP=-log10(as.numeric(vol_input$Pvalue))
vol_input<-subset(vol_input,vol_input$change!="NA")

pdf(file = "vol_input.pdf",height = 4,width = 5)
ggplot(
  #设置数据
  vol_input, 
  aes(x = logFC, 
      y = nega_logP, 
      colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-as.numeric(cut_off_logFC),as.numeric(cut_off_logFC)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2FC",
       y="-log10 (Adj.p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+
  geom_text_repel(
    data = subset(vol_input, vol_input$Gene%in%top$gene),
    aes(label = Gene),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
dev.off()


###### 2 monocle ######
library(monocle)
dir.create("monocle")
setwd("monocle")

#选定要拟时序分析的细胞种类
mmt<-sqy

#把TAM，monocyte等sub细胞分类信息给R
mmt$cell_type_val<-Idents(mmt)
mmt@meta.data[1:5,]

#寻找高变基因，作为拟时序降维的基因
mmt<-FindVariableFeatures(mmt,nfeatures = 2000)
#>这一步，不光是可以使用FindVariableFeatures寻找多变基因
#>也可以用前面的用于pca降维的2000个多变基因作为降维基因
#>二者都要试试，看哪个更符合生物学进程和假说，就用哪个
#>生信需要反复调整参数，而不是去跑一个流程就行
#>?FindVariableFeatures

#提取数据
matrix<-as.matrix(mmt@assays$RNA@counts)  #提取最原始的count数据
dim(matrix) #看多少行多少列
#matrix[1:500,1:10]


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
sc_cds2 <- orderCells(sc_cds2,root_state = 3)  #把cell的顺序排出来、
saveRDS(sc_cds2,"monocle_for_plot.rds")
#beepr::beep(1)


plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T,)#拟时序，颜色表示假时间

#降维图形绘制
#plot_cell_trajectory(sc_cds2, markers = "SIAH2",show_branch_points = T)  #拟时序，颜色表示细胞sub种类
#>看这一步做出的图平不平滑，是不是一段线段只有基本上一种sub细胞
#>如果不是的话记得修改前面的参数
#sc_cds2 <- orderCells(sc_cds2,root_state = 8)
#pData(sc_cds_2)[1:5,]
#p1<-plot_cell_trajectory(sc_cds2, color_by = "RNA_snn_res.0.8",show_branch_points = F)#+facet_wrap(~RNA_snn_res.0.8)

pdf("pca_细胞种类.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "GFRA2_group",show_branch_points = T)+
  scale_color_manual(values = c(cors2[2],cors2[1]))+
  theme_bw()#+facet_wrap(~Site) #拟时序，颜色表示细胞state
dev.off()

pdf("pca_假时间.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = T,)+
  theme_bw()#拟时序，颜色表示假时间
dev.off()

pdf("pca_分群.pdf",width = 5.5,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "marker_cluster",show_branch_points = T,)+
  scale_color_manual(values = cors2)+
  theme_bw()#拟时序，颜色表示假时间
dev.off()

pdf("pca_GFRA2.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("GFRA2"),
                     use_color_gradient = T,
                     show_branch_points = F, 
                     show_tree = T,show_backbone = T)+
  theme_bw()#拟时序，颜色表示marker表达量
dev.off()

pdf("pca_DCX.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("DCX"),
                     use_color_gradient = T,
                     show_branch_points = F, 
                     show_tree = T,
                     show_backbone = T)+
  theme_bw()#拟时序，颜色表示marker表达量
dev.off()

pdf("pca_MAP2.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("MAP2"),
                     use_color_gradient = T,
                     show_branch_points = F, 
                     show_tree = T,
                     show_backbone = T)+
  theme_bw()#拟时序，颜色表示marker表达量
dev.off()

pdf("pca_TH.pdf",width = 5,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("TH"),
                     use_color_gradient = T,
                     show_branch_points = F, 
                     show_tree = T,
                     show_backbone = T)+
  theme_bw()#拟时序，颜色表示marker表达量
dev.off()

save.image("24.3.7.RData")
#240行24.3.7



#画想画的基因的热图
pdf("heatmap_DaNs功能_拟时序热图.pdf",width = 3,height = 3)
to_be_plot <- row.names(subset(fData(sc_cds2), gene_short_name %in% all_genes))
cds_subset1 <- sc_cds2[to_be_plot,]
plot_pseudotime_heatmap(cds_subset1,show_rownames = T,norm_method = c("log", "vstExprs"),num_clusters = 4,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()

#散点图（横轴，假时间；纵轴，自定义marker基因的表达量）
F_map2cor<-
  plot_genes_in_pseudotime(sc_cds2["MAP2",],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor(label.y = 150)+
  geom_smooth(method = "loess")+
  ylim(20,150)+
  theme_bw()+
  NoLegend()

F_dcxcor<-
  plot_genes_in_pseudotime(sc_cds2["DCX",],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor(label.y = 1.1)+
  theme_bw()+
  geom_smooth(method = "loess")+
  NoLegend()+
  ylab("")+
  xlab("")

F_gfra2cor<-
  plot_genes_in_pseudotime(sc_cds2["GFRA2",],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor(label.y = 9)+
  theme_bw()+
  geom_smooth(method = "loess")+
  ylim(0,10)+
  NoLegend()+
  ylab("")

F_calb2cor<-
  plot_genes_in_pseudotime(sc_cds2["CALB2",],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
    stat_cor()+
  theme_bw()+
    geom_smooth(method = "loess")+
  NoLegend()+
  xlab("")#Calretinin
  
pdf("成熟和幼稚标志物cor.pdf",width = 5,height = 5)
(F_calb2cor|F_dcxcor)/(F_map2cor|F_gfra2cor)
dev.off()

#寻找随着假时间变化的基因
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 2)
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))


#记得先看看想研究的基因是不是随时间变化的，再做后续分析
#输出随时间变化的基因
write.csv(sc_cds2[sig_gene_names,],"pseudotime_genes.csv")


#绘制随着假时间变化的基因表达的热图
pdf("pseudoheatmap.pdf",width = 5,height = 5)
pseudoplot<-plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                    num_clusters = 10,###进行调节分群，后面也要改
                                    cores = 10,
                                    show_rownames = F,return_heatmap = T)
dev.off()


#获得每个亚群的基因名称，以便于后续富集分析
clusters <- cutree(pseudoplot$tree_row, k = 10)###进行调节分群，前面也要改
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering,"cor_time_gene_cluster.csv")
#beep(1)
Dopa_meta_clustering<-clustering[Dopa_meta,]
Dopa_meta_clustering<-cbind("Dopa_meta_genes"=Dopa_meta,"Cluster"=Dopa_meta_clustering)
write.csv(Dopa_meta_clustering,"dopa_meta_clustering.csv",row.names = F,quote = F)


#看多巴胺代谢的拟时序
pdf(file = "cor_pse_分解.pdf",width = 3,height = 4)
plot_genes_in_pseudotime(sc_cds2[c("MAOA",'MAOB','COMT'),],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor()+
  theme_bw()+
  geom_smooth(method = "loess")+
  NoLegend()+
  ggtitle("Breakdown")
dev.off()

pdf(file = "cor_pse_合成.pdf",width = 3,height = 4)
plot_genes_in_pseudotime(sc_cds2[c("TH",'DDC'),],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor()+
  theme_bw()+
  geom_smooth(method = "loess")+
  NoLegend()+
  ggtitle("Biosynthesis")+
  ylim(5,100)
dev.off()

pdf(file = "cor_pse_释放.pdf",width = 3,height = 4)
plot_genes_in_pseudotime(sc_cds2[c("CACNA1A",'CACNA1B'),],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor()+
  theme_bw()+
  geom_smooth(method = "loess")+
  NoLegend()+
  ggtitle("Release")+
  ylim(10,100)
dev.off()

pdf(file = "cor_pse_重吸收.pdf",width = 3,height = 4)
plot_genes_in_pseudotime(sc_cds2[c("DRD2",'SLC6A3'),],color_by = "Pseudotime",trend_formula = y~x/100,cell_size = -1)+
  stat_cor()+
  theme_bw()+
  geom_smooth(method = "loess")+
  NoLegend()+
  ggtitle("Reabsorption")+
  ylim(1,40)
dev.off()


##### 3 反卷积 #####

library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library('e1071')
library(future)
library(CIBERSORT)

setwd("~/sqy/24.3.6 pd/GSE178265_MACACA/")
dir.create("cibersort")
setwd("cibersort/")


#### GFRA2+cell丰度
setwd("~/sqy/24.3.6 pd/GSE178265_MACACA/cibersort")
dir.create("GFRA2")
setwd("GFRA2/")

#导入单细胞数据,形成系数表
Idents(sqy) <- sqy$GFRA2_group#细胞类型是哪一列
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
write.csv(results,"result_cibersort.csv",quote = F,row.names = T)

library(reshape2)
library(ggpubr)
library(corrplot)

rt<-results[,1:(ncol(results)-3)]

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

#绘制箱线图
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
data$Type<-gsub("Treat","PD",data$Type)
data$Type<-gsub("Control","Normal",data$Type)

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
pdf(file="immune.diff.cluster.pdf", width=5, height=5)
print(boxplot)
dev.off()



#### cluster丰度
dir.create("cluster")
setwd("cluster/")

rt<-read.csv("result_cibersort.csv",sep = ",",header = T,row.names = 1)
data<-data.frame(row.names = rownames(rt))
data$High_GFRA2<-rt$X3.COBLL1.+rt$X5.ANKFN1.+rt$X7.HS3ST4.+rt$X9.KIAA1217.
data$Low_GFRA2<-rt$X0.CNTN5.+rt$X1.STXBP6.+rt$X2.TENM2.+rt$X4.CALCR.+rt$X6.C8H8orf34.+rt$X8.PHACTR1+rt$X10.CHN2.
write.csv(data,"GFRA2_fraction.csv",row.names = T,quote = F)

library(reshape2)
library(ggpubr)
library(corrplot)

rt=data

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

#绘制箱线图
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
data$Type<-gsub("Treat","PD",data$Type)
data$Type<-gsub("Control","Normal",data$Type)

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
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("p<0.001", "p<0.01", "p<0.05", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.GFRA2.pdf", width=5, height=5)
print(boxplot)
dev.off()

#看GFRA2和GFRA2+DaNs的关系
GFRA2_exp<-read.table("merge.normalize.txt",sep = "\t",row.names = 1,header = T)["GFRA2",]
GFRA2_exp<-t(GFRA2_exp)
GFRA2_exp_GFRA2_DaNs<-cbind(GFRA2_exp,rt)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
GFRA2_exp_GFRA2_DaNs<-cbind("Type"=Type,GFRA2_exp_GFRA2_DaNs)
GFRA2_exp_GFRA2_DaNs_PD<-subset(GFRA2_exp_GFRA2_DaNs,GFRA2_exp_GFRA2_DaNs$Type=="Treat")

pdf(file = "GFRA2_GFRA2DaNs_cor.pdf",width =4,height = 4)
p2<-ggplot(GFRA2_exp_GFRA2_DaNs_PD, aes(x=GFRA2, y=High_GFRA2)) + geom_point()+
  geom_smooth(method = lm)+
  stat_cor(method = "spearman")+
  theme_bw()
print(p2)
dev.off()





######################### 4 GSEA ########################

# 1 提取分析的通路
library(devtools)
#install_github("arc85/singleseqgset")##https://arc85.github.io/singleseqgset/articles/singleseqgset.html
library(singleseqgset)
library(dplyr)
#?msigdbr
msigdbr_collections()%>%print(n=33) #查找能用的基因集的category和subcategory,在下一句话中更改gmt文件

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
expr.mat=sqy@assays$RNA@data
cluster.ids=sqy$marker_cluster
logfc.data <- logFC(cluster.ids=cluster.ids,expr.mat=expr.mat)  #大矩阵可能报错，内存不够
#logfc.data <- logFC(cluster.ids=sqy$Seurat_harmony[1:1000],expr.mat=expr.mat[,c(1:1000)])
#names(logfc.data)
gse.res <- wmw_gsea(expr.mat=expr.mat,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
#rm(logfc.data)

names(gse.res)
res.stats <- gse.res[["GSEA_statistics"]]
#res.stats[1:5,1:5]  #前五行，前五列

ac<-data.frame(cluster=factor(colnames(res.stats))) #变成作图需要的矩阵格式
unique(ac$cluster)
rownames(ac)=colnames(res.stats)

ann_colors = list(
  cluster = c("0 CNTN5+"=cors[1],
              "1 STXBP6+"=cors[2],
              "2 TENM2+"=cors[3],
              "4 CALCR+"=cors[5],
              "3 COBLL1+"=cors[4],
              "6 C8H8orf34+"=cors[7],
              "8 PHACTR1"=cors[9],
              "10 CHN2+"=cors[11],
              "5 ANKFN1+"=cors[6],
              "7 HS3ST4+"=cors[8],
              "9 KIAA1217+"=cors[10]
              )
  #GeneClass = c("CD4+" = "#7570B3", "CD8+" = "#E7298A","NK" = "#66A61E")
) #更改颜色，复制前面的即可

#cors

#dev.off()
library(stringr)
rownames(res.stats)<-str_replace_all(rownames(res.stats),pattern = "KEGG_",replacement ="") #把"HALLMARK_"，"KEGG_"等前缀变成空格
rownames(res.stats)<-str_replace_all(rownames(res.stats),pattern = "_",replacement =" ") #把"_"变成空格
pdf("pheatmap_enrichment.pdf",width = 7,height = 30)
pheatmap::pheatmap(res.stats,
                   fontsize_row = 8,
                   annotation_col = ac,
                   annotation_legend = F,
                   annotation_colors = ann_colors,
                   cluster_rows = T,
                   cluster_cols = F,
                   scale = "row")
dev.off()
#write.csv(res.stats,file="gsea.result.NES.keratinocyte.csv")
#sb的r报错的话直接导出自己画热图，爷nb爷就是这么任性
write.csv(res.stats,"gsea.result.NES.csv")
write.csv(gse.res$GSEA_p_values,"gsea.result.pval.csv")

useful_keggpath<-read.table("有用的keggpath.txt",header = F,sep = "\t")[,1]
useful_res.stats<-res.stats[useful_keggpath,]
pdf("pheatmap_enrichment_useful.pdf",width = 7,height = 7)
pheatmap::pheatmap(useful_res.stats,
                   fontsize_row = 8,
                   annotation_col = ac,
                   annotation_legend = F,
                   annotation_colors = ann_colors,
                   cluster_rows = T,
                   cluster_cols = T,
                   scale = "row")
dev.off()

metabolism_keggpath<-read.table("代谢keggpath.txt",header = F,sep = "\t")[,1]
metabolism_res.stats<-res.stats[metabolism_keggpath,]
col_anno<-data.frame(row.names = levels(sqy$marker_cluster),"Group"=c("Low",'Low',"low",'High','Low',"High","Low","High","Low","High","Low"))
pdf("pheatmap_enrichment_metabolism.pdf",width = 7,height = 7)
pheatmap::pheatmap(metabolism_res.stats,
                   fontsize_row = 8,
                   annotation_col = ac,
                   annotation_legend = F,
                   annotation_colors = ann_colors,
                   cluster_rows = T,
                   cluster_cols = T,
                   scale = "row",
                   cutree_rows = 2,
                   cutree_cols = 2)
dev.off()


