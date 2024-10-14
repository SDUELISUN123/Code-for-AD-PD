setwd("~/sqy/24.3.6 pd/GSE175814_AD/Seurat")


###### 1 导入原来的数据并重新处理 #####
setwd("~/sqy/24.3.6 pd/GSE175814_AD/Seurat")

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

plan("multicore", workers = 8) ###set the compute core
options(future.globals.maxSize = 500000 * 1024^2)#500G

sqy<-readRDS("rawdata_cbind.rds")

library(devtools)
library(harmony)
sqy=CreateSeuratObject(counts = sqy,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^MT-")#鼠源的换成mt

#质控
VlnPlot(object = sqy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")  #ncol是指每行放几张图
sqy=subset(x = sqy, subset = nFeature_RNA > 50 & percent.mt < 10)    #对数据进行过滤

#测序深度的相关性图
plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot1|plot2
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
#pc=30

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
sqy <- FindNeighbors(object = sqy, dims = 1:30)       #计算邻接距离
sqy <- FindClusters(object = sqy, resolution = 0.8)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
#sqy <- RunUMAP(object = sqy, dims = 1:15) 
sqy<-RunTSNE(object = sqy, dims = 1:30)

sqy@active.ident<-sqy$seurat_clusters
pdf(file = "tsne_cluster.pdf",width=4,height = 4)
TSNEPlot(object = sqy, label = TRUE,cols=cors)+
  theme_bw()+
  NoLegend()
dev.off()

pdf(file = "fea_GFRA2_exp.pdf",width=3,height = 3)
FeaturePlot(sqy,features = "GFRA2",label = T,reduction = "tsne")+
  theme_bw()+
  NoLegend()
dev.off()

pdf(file = "Vln_GFRA2_exp.pdf",width = 6,height = 4)
VlnPlot(sqy,features = "GFRA2",cols = cors,pt.size = 0)+
  geom_boxplot(outlier.size = 0.1)+
  stat_compare_means(label.x = 4)+
  NoLegend()
dev.off()

sqy@active.ident<-sqy$orig.ident
sqy<-RenameIdents(sqy,
                  "AD1"="AD",
                  "AD2"="AD",
                  "Control1"="Control",
                  "Control2"="Control")
sqy$group<-sqy@active.ident

sqy@active.ident<-sqy$seurat_clusters
DotPlot(sqy,features = c("SCN2A","MAP2",
                         "TH","SLC6A3",
                         "SLC17A6",
                         'GAD2',"GAD1",
                         "GRIK1",
                         "SLC17A7","CAMK2A",
                         "GFRA2",
                         "GFAP",
                         "SLC6A4",
                         "CHAT",
                         "MNX1",
                         "MOG","MAG",
                         "APBB1IP","CD86","CD14",
                         "CD3D","CD3E","NKG7",
                         "VWF",
                         "DCX"
                         ))+
  RotatedAxis()+
  NoLegend()


#一级注释
sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  "0"="Oligodendrocyte",
                  '1'='Oligodendrocyte',
                  '2'='Neuron',
                  '3'='Microglia',
                  '4'='Others',
                  '5'='Neuron',
                  '6'='Others',
                  '7'='Neuron',
                  '8'='Astrocyte',
                  '9'='Neuron',
                  '10'='Neuron',
                  '11'='Astrocyte',
                  '12'='EC',
                  '13'='Neuron',
                  '14'='Astrocyte',
                  '15'='Microglia',
                  '16'='Neuron',
                  '17'='Neuron',
                  '18'='Neuron',
                  '19'='Neuron',
                  '20'='Neuron',
                  '21'='Neuron',
                  '22'='Neuron',
                  '23'='Neuron',
                  '24'='Neuron',
                  '25'='Neuron',
                  '26'='Neuron',
                  '27'='Microglia',
                  '28'='Neuron',
                  '29'='Microglia'
)
sqy$celltype<-sqy@active.ident

pdf(file = "tsne_cluster_celltype.pdf",width=5,height = 4)
TSNEPlot(object = sqy, label = F,cols=cors)+
  theme_bw()
dev.off()


F1<-(FeaturePlot(sqy,features = 'GFRA2')+
    NoLegend()|
FeaturePlot(sqy,features = 'VWF')+
  NoLegend()|
FeaturePlot(sqy,features = 'APBB1IP')+
  NoLegend())
F2<-(FeaturePlot(sqy,features = 'MOG')+
  NoLegend()|
FeaturePlot(sqy,features = 'MAP2')+
  NoLegend()|
FeaturePlot(sqy,features = 'GFAP')+
  NoLegend())
pdf("fea_marker.pdf",width = 9,height = 6)
F1/F2
dev.off()

pdf("fea_GFRA2.pdf",width = 4,height = 4)
FeaturePlot(sqy,features = 'GFRA2',cols = c("lightgrey","red"))+
  theme_bw()+
  NoLegend()
dev.off()

#二级注释——神经元
sqy@active.ident<-sqy$seurat_clusters
DotPlot(sqy,features = c("SCN2A","MAP2",
                         "TH","SLC6A3",
                         'GAD2',"GAD1",
                         "SLC17A7","CAMK2A",
                         "GFRA2",
                         "SLC6A4",
                         "CHAT",
                         "MNX1",
                         "DCX"
))+
  RotatedAxis()+
  NoLegend()
sqy<-RenameIdents(sqy,
                  "0"="Oligodendrocyte",
                  '1'='Oligodendrocyte',
                  '2'='Neuron_inhibitory',
                  '3'='Microglia',
                  '4'='Others',
                  '5'='Neuron_excitatory',
                  '6'='Others',
                  '7'='Neuron_excitatory',
                  '8'='Astrocyte',
                  '9'='Neuron_excitatory',
                  '10'='Neuron_excitatory',
                  '11'='Astrocyte',
                  '12'='EC',
                  '13'='Neuron_inhibitory',
                  '14'='Astrocyte',
                  '15'='Microglia',
                  '16'='Neuron_inhibitory',
                  '17'='Neuron_excitatory',
                  '18'='Neuron_excitatory',
                  '19'='Neuron_excitatory',
                  '20'='Neuron_inhibitory',
                  '21'='Neuron_inhibitory',
                  '22'='Neuron_excitatory',
                  '23'='Neuron_inhibitory',
                  '24'='Neuron_excitatory',
                  '25'='Neuron_excitatory',
                  '26'='Neuron_excitatory',
                  '27'='Microglia',
                  '28'='Neuron_excitatory',
                  '29'='Microglia'
)
sqy$celltype_minor<-sqy@active.ident
pdf(file = "tsne_cluster_celltypeminor.pdf",width=5,height = 4)
TSNEPlot(object = sqy, label = F,cols=cors)+
  theme_bw()
dev.off()



pdf("fea_marker_neuron.pdf",width = 9,height = 3)
FeaturePlot(sqy[,sqy$celltype%in%"Neuron"],features = 'GFRA2')+
  NoLegend()|
  FeaturePlot(sqy[,sqy$celltype%in%"Neuron"],features = 'GAD1')+
  NoLegend()|
  FeaturePlot(sqy[,sqy$celltype%in%"Neuron"],features = 'CAMK2A')+
  NoLegend()
dev.off()


pdf("vln_GFRA2_exp_celltypeminor.pdf",width = 4,height = 4)
VlnPlot(sqy,features = "GFRA2",group.by = "celltype_minor",cols = cors,pt.size = 0)+
  geom_boxplot(outlier.size = 0)+
  stat_compare_means(label.x = 1.3,label.y = 3.8)+
  theme_bw()+
  NoLegend()+
  RotatedAxis()
dev.off()

pdf("vln_兴奋神经元中GFRA2表达情况.pdf",width = 4,height = 4)
VlnPlot(sqy[,sqy$celltype_minor%in%"Neuron_excitatory"],features = "GFRA2",group.by = "seurat_clusters",cols = cors,pt.size = 0)+
  geom_boxplot(outlier.size = 0)+
  stat_compare_means(label.x = 2,label.y = 3.8)+
  NoLegend()
dev.off()

pdf("vln_兴奋神经元中GFRA2表达情况_AD和con之间.pdf",width = 3,height = 4)
VlnPlot(sqy[,sqy$celltype_minor%in%"Neuron_excitatory"],features = "GFRA2",cols = cors,pt.size = -1,group.by = "group",sort=T)+
  geom_boxplot(outlier.size = -1)+
  stat_compare_means()+
  theme_bw()+
  NoLegend()+
  RotatedAxis()+
  ggtitle("Excitatory Neurons Only")+
  ylab("GFRA2 Expression")
dev.off()

saveRDS(sqy,"GSE175814注释完毕.rds")




#### 2 兴奋型神经元重新聚类 ####

sqy@active.ident<-sqy$celltype_minor
pdf("cor_GFRA2_CAMK2A.pdf",width = 3,height = 3)
FeatureScatter(sqy,feature1 = "GFRA2",feature2 = c("CAMK2A"),cols = cors,pt.size = 0.3,)+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_bw()+
  NoLegend()+
  ggtitle("")
dev.off()

sqy@active.ident<-sqy$celltype_minor
pdf("cor_GFRA2_SLC17A7.pdf",width = 3,height = 3)
FeatureScatter(sqy,feature1 = "GFRA2",feature2 = c("SLC17A7"),cols = cors,pt.size = 0.3,)+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_bw()+
  NoLegend()+
  ggtitle("")
dev.off()

Neuron_excitatory<-sqy[,sqy$celltype_minor%in%"Neuron_excitatory"]

#寻找每个亚群的标志物，并进行重新注释
Neuron_excitatory@active.ident<-Neuron_excitatory$seurat_clusters

logFCfilter=0
adjPvalFilter=0.05
sqy.markers <- FindAllMarkers(object = Neuron_excitatory,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="Neuron_excita_cluster_markers.csv")
top10 <- sqy.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "Neuron_excita_cluster_markers_top10.csv")

#绘制marker在各个cluster的热图
pdf(file = "doheatmap_cluster_top10.pdf",width = 10, height = 10)
DoHeatmap(object = Neuron_excitatory, features = top10$gene) + NoLegend()
dev.off()

#重新命名亚群
Neuron_excitatory@active.ident<-Neuron_excitatory$seurat_clusters
Neuron_excitatory<-RenameIdents(Neuron_excitatory,
                                '5'='5 AL117329.1+',
                                '7'='7 CALM1+',
                                '9'='9 AC109466.1+',
                                '10'='10 TSHZ2+',
                                '17'='17 TOX+',
                                '18'='18 GFRA2_lo XIST+',
                                '19'='19 ADARB2+',
                                '22'='22 ADAMTSL1+',
                                '24'='24 KIAA1217+',
                                '25'='25 TRPS1+',
                                '26'='26 GPC5+',
                                '28'='28 GFRA2_hi TRPM3+')
Neuron_excitatory$marker_cluster<-Neuron_excitatory@active.ident

Neuron_excitatory@active.ident<-Neuron_excitatory$marker_cluster

pdf(file = "tsne_cluster_Neuron_excitatory.pdf",width=5,height = 4)
TSNEPlot(object = Neuron_excitatory, label = F,cols=cors)+
  theme_bw()
dev.off()

pdf("vln_兴奋神经元中GFRA2表达情况_重新聚类.pdf",width = 4.5,height = 3)
VlnPlot(Neuron_excitatory,features = "GFRA2",cols = cors,pt.size = -1,group.by = "marker_cluster",sort=T)+
  geom_boxplot(outlier.size = -1)+
  stat_compare_means(label.x = 2,label.y = 3.5)+
  theme_bw()+
  NoLegend()+
  ggtitle("GFRA2 High --------> GFRA2 Low")+
  ylab("GFRA2 Expression")+
  rotate_x_text(angle = 45)+
  xlab("")
dev.off()

pdf("vln_兴奋神经元中XIST表达情况_重新聚类.pdf",width = 4.5,height = 3)
VlnPlot(Neuron_excitatory,features = "XIST",cols = cors,pt.size = -1,group.by = "marker_cluster",sort=T)+
  geom_boxplot(outlier.size = -1)+
  stat_compare_means(label.x = 2,label.y = 4)+
  theme_bw()+
  NoLegend()+
  ggtitle("XIST High --------> XIST Low")+
  ylab("XIST Expression")+
  rotate_x_text(angle = 45)+
  xlab("")
dev.off()
#18群的GFRA2表达最低，前面证明了他是AD中特有的群
#28群的GFRA2最高，并且后面证明这群在AD中减少

pdf(file = "tsne_cluster_Neuron_excitatory_split_group.pdf",width=8,height = 4)
TSNEPlot(object = Neuron_excitatory, label = F,cols=cors,split.by="group")+
  theme_bw()
dev.off()
#18群是AD相关的兴奋性神经元

pdf(file = "tsne_cluster_Neuron_excitatory_group_group.pdf",width=5,height = 4)
TSNEPlot(object = Neuron_excitatory, label = F,cols=cors,group.by="group")+
  theme_bw()
dev.off()

pdf("fea_GFRA2_XIST_exp_excitatory_neurons.pdf",width = 3,height = 6)
F1=FeaturePlot(Neuron_excitatory,features = c("GFRA2"),reduction = "tsne",cols = c("lightgrey","red"),ncol = 1)+
  theme_bw()+
  NoLegend()
F2=FeaturePlot(Neuron_excitatory,features = c("XIST"),reduction = "tsne",cols = c("lightgrey","red"),ncol = 1)+
  theme_bw()+
  NoLegend()
F1/F2
dev.off()

saveRDS(sqy,'GSE175814注释完毕.rds')
save.image("24.3.14.RData")





####### 3 monocle excitatory neuron #####
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(monocle)
library(ggsci)
library(beepr)
library(ggplot2)
library(ggpubr)

plan("multicore", workers = 4) ###set the compute core
options(future.globals.maxSize = 15000 * 1024^2)
#getwd()

setwd("E:\\科研\\230712 糖足胞葬\\2.6 singlecell2\\1 foot\\融合\\monocle")

dir.create("monocle")
setwd("monocle")

#选定要拟时序分析的细胞种类
mmt<-Neuron_excitatory

#把TAM，monocyte等sub细胞分类信息给R
mmt$seurat_clusters<-Idents(mmt)
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
sc_cds2 <- orderCells(sc_cds2,
                      #root_state = 3,#可以更改其实分支
                      )  #把cell的顺序排出来、
saveRDS(sc_cds2,"monocle_for_plot.rds")
#beepr::beep(1)


pdf("pca_假时间.pdf",width = 4,height = 4)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = T,)+
  theme_bw()+
  theme(legend.position = c(0.2,0.25))#拟时序，颜色表示假时间
dev.off()

pdf("pca_分群.pdf",width = 6,height = 6)
plot_cell_trajectory(sc_cds2, color_by = "group",show_branch_points = F,)+
  scale_color_manual(values = c(cors[2],cors[1]))+
  facet_wrap(~marker_cluster)+
  theme_bw()+
  theme(legend.position = "top",)#拟时序，颜色表示假时间
dev.off()
#18群（AD特有，在最末，GFRA2最低）
#28群（Normal特有，在最初，GFRA2最高）

pdf("pca_GFRA2.pdf",width = 4,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("GFRA2"),use_color_gradient = T,show_branch_points = F, show_tree = T,show_backbone = T)+
  theme_bw()+
  theme(legend.position = c(0.25,0.25))+
  scale_color_gradient2(low="navy",mid="white",high="#B00000")#拟时序，颜色表示marker表达量
dev.off()
pdf("pca_DCX.pdf",width = 4,height = 4)
plot_cell_trajectory(sc_cds2, markers = c("DCX"),use_color_gradient = T,show_branch_points = F, show_tree = T,show_backbone = T)+
  theme_bw()+
  theme(legend.position = c(0.25,0.25))+
  scale_color_gradient2(low="navy",mid="white",high="#B00000")#拟时序，颜色表示marker表达量
dev.off()

#查找和拟时间有关的基因
expressed_genes=row.names(subset(fData(sc_cds2),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(sc_cds2[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]

saveRDS(sc_cds2, file = "monocle_for_plot.rds")
write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

#BEAM节点分析
BEAM_res <- BEAM(sc_cds2, branch_point = 1, progenitor_method = "duplicate",cores =10)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

#节点热图
tmp1=plot_genes_branched_heatmap(sc_cds2[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 10, #这些基因被分成几个group
                                 cores = 10,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 hmcols = NULL, #默认值
                                 #hmcols = colorRampPalette(c("navy","white","firebrick3"))(100),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)

pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()

#获得每个亚群的基因名称，以便于后续富集分析
clusters<-tmp1$annotation_row
write.csv(clusters,"每个基因在哪个cluster.csv",row.names = T,quote = F)
want_cluster<-as.character(clusters["GFRA2",])
want_cluster
cluster_want_genes<-subset(clusters,clusters$Cluster==want_cluster)%>%rownames(.)
write.csv(cluster_want_genes,"和目标基因在一个cluster的基因.csv",row.names = F,quote = F,col.names = F)
#放到微生信里富集分析（clusterprofiler）

pdf("cors_genes_兴奋神经元标志物.pdf",width = 3,height = 4)
plot_genes_branched_pseudotime(sc_cds2[c("GFRA2","SLC17A7","CAMK2A"),],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               branch_labels = c("Cell fate 1", "Cell fate 2"),
                               cell_size = -1,)+
  theme_bw()
dev.off()




####### 4 cibersort #######
setwd("~/sqy/24.3.6 pd/GSE175814_AD/Seurat")
dir.create("cibersort")
setwd("cibersort/")



library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library('e1071')
library(future)
library(CIBERSORT)



#导入单细胞数据,形成系数表
Idents(sqy) <- sqy$celltype_minor#细胞类型是哪一列
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
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
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
data$Type<-gsub("treat","AD",data$Type)
data$Type<-gsub("con","Normal",data$Type)


#绘制箱线图
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=c("#FF5555", "#0088FF"),)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.cluster.pdf", width=5, height=5)
print(boxplot)
dev.off()


#看GFRA2和GFRA2+Exciatatory的关系
GFRA2_exp<-read.table("merge.normalize.txt",sep = "\t",row.names = 1,header = T)["GFRA2",]
GFRA2_exp<-t(GFRA2_exp)
GFRA2_exp_cells<-cbind(GFRA2_exp,rt)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
GFRA2_exp_cells<-cbind("Type"=Type,GFRA2_exp_cells)
GFRA2_exp_cells<-as.data.frame(GFRA2_exp_cells)


GFRA2_exp_cells_AD<-subset(GFRA2_exp_cells,GFRA2_exp_cells$Type%in%"treat")
GFRA2_exp_cells_AD$Neuron_excitatory<-as.numeric(GFRA2_exp_cells_AD$Neuron_excitatory)
GFRA2_exp_cells_AD$GFRA2<-as.numeric(GFRA2_exp_cells_AD$GFRA2)
pdf(file = "GFRA2_Excitatory_cor_AD.pdf",width =4,height = 4)
p2<-ggplot(GFRA2_exp_cells_AD, aes(x=GFRA2, y=Neuron_excitatory)) + geom_point()+
  geom_smooth(method = loess)+
  stat_cor(method = "spearman")+
  theme_bw()+
  ylim(0,0.7)
print(p2)
dev.off()


GFRA2_exp_cells_con<-subset(GFRA2_exp_cells,GFRA2_exp_cells$Type%in%"con")

GFRA2_exp_cells_con$Neuron_excitatory<-as.numeric(GFRA2_exp_cells_con$Neuron_excitatory)
GFRA2_exp_cells_con$GFRA2<-as.numeric(GFRA2_exp_cells_con$GFRA2)
pdf(file = "GFRA2_Excitatory_cor_con.pdf",width =4,height = 4)
p2<-ggplot(GFRA2_exp_cells_con, aes(x=GFRA2, y=Neuron_excitatory)) + geom_point()+
  geom_smooth(method = loess)+
  stat_cor(method = "spearman")+
  theme_bw()+
  ylim(0,0.85)
print(p2)
dev.off()


##看各个亚群
setwd("~/sqy/24.3.6 pd/GSE175814_AD/Seurat/cibersort")
dir.create("每个cluster的丰度")
setwd("每个cluster的丰度/")


#导入单细胞数据,形成系数表
Idents(sqy) <- Neuron_excitatory$marker_cluster#细胞类型是哪一列
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
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
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
data$Type<-gsub("treat","AD",data$Type)
data$Type<-gsub("con","Normal",data$Type)
#data=subset(data,data$Immune%in%c("18 GFRA2_lo XIST+","28 GFRA2_hi TRPM3+"))

#绘制箱线图
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=c("#FF5555", "#0088FF"),outlier.shape =NULL)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.cluster.pdf", width=5, height=5)
print(boxplot)
dev.off()






###### 5 cellchat #####

setwd("~/sqy/24.3.6 pd/GSE175814_AD/Seurat")
sqy<-readRDS("GSE175814注释完毕.rds")

dir.create("cellchat")
setwd("cellchat")


library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(stringr)
library(CellChat)
library(patchwork)

plan("multicore", workers = 10) ###set the compute core

cors<-ggsci::pal_igv(alpha = 0.3)(51)

sqy@active.ident<-factor(sqy$celltype_minor)
table(sqy$celltype_minor)

pdf("dot_ridge_GFRA2_exp_excitatory.pdf",width = 6,height = 4)
DotPlot(sqy[,sqy$celltype_minor%in%"Neuron_excitatory"],features = "GFRA2",group.by = "seurat_clusters")|
RidgePlot(sqy[,sqy$celltype_minor%in%"Neuron_excitatory"],features = "GFRA2",group.by = "seurat_clusters",sort=T,cols = cors)
dev.off()

pdf("tsne_cluster_excitatory.pdf",width = 4,height = 4)
TSNEPlot(object = sqy[,sqy$celltype_minor%in%"Neuron_excitatory"],
         group.by="seurat_clusters",
         cols=cors,label=T)+
  theme_bw()+
  NoLegend()
dev.off()

#准备cellchat用的Seurat文件
sqy <- readRDS("~/sqy/24.3.6 pd/GSE175814_AD/Seurat/GSE175814注释完毕.rds")
Neuron_excitatory <- readRDS("~/sqy/24.3.6 pd/GSE175814_AD/Seurat/Neuron_excitatory.rds")
sqy<-sqy[,!sqy$celltype%in%c("Neuron","EC","Others")]
sqy$cell_type<-sqy$celltype
Neuron_excitatory<-Neuron_excitatory[,Neuron_excitatory$seurat_clusters%in%c("18","28")]
Neuron_excitatory$cell_type<-Neuron_excitatory$marker_cluster
sqy2<-merge(x=sqy,y=Neuron_excitatory)
table(sqy2$cell_type)
sqy<-sqy2
rm(Neuron_excitatory,sqy2)
gc()

#cellchat
table(sqy$cell_type)#看各种细胞类型的名字和数目，以便后边改（要去掉没有的细胞）
#sqy<-sqy[,sqy$cell_type!=c("RBCs")]  #去掉没有的细胞，防止报错，慎重覆盖变量！
normal.input <- GetAssayData(sqy, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(sqy$cell_type)#用哪一列当细胞类型
meta <- data.frame(group = labels, row.names = rownames(sqy@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_") #把斜杠替换成下划线，以免报错
#下句话记得把上面各种细胞类型粘上，然后把"/"改成"_"，并且不能放数目为0的细胞，防报错
meta$group <- factor(meta$group,levels=c("Astrocyte","Microglia","Oligodendrocyte","18 GFRA2_lo XIST+","28 GFRA2_hi TRPM3+"))
normal_cellchat <- createCellChat(object = normal.input, meta = meta, group.by = "group")
table(normal_cellchat@idents)  #看有没有缺失值NA和数目为0的细胞
saveRDS(normal_cellchat,file="cellchat_original.rds")
rm(normaldata,normal.input,labels,meta)


#预处理要跑十分钟左右，去休息一下吧！
normal_cellchat@DB <- CellChatDB.human#物种：人类,可以改为“CellChatDB.mouse”
#View(CellChatDB.human$interaction)  #看有哪些互作
normal_cellchat <- subsetData(normal_cellchat) # subset the expression data of signaling genes for saving computation cost
normal_cellchat <- identifyOverExpressedGenes(normal_cellchat)  #确定每个细胞亚群中的过表达基因
normal_cellchat <- identifyOverExpressedInteractions(normal_cellchat) #寻找过表达的interaction
normal_cellchat <- projectData(normal_cellchat, PPI.human)  #向ppi投射
normal_cellchat <- computeCommunProb(normal_cellchat,raw.use = T) #算interaction的可能性
normal_cellchat <- filterCommunication(normal_cellchat, min.cells = 10)  #去除interaction很少的细胞
normal_cellchat <- computeCommunProbPathway(normal_cellchat)  #计算通路
normal_cellchat <- netAnalysis_computeCentrality(normal_cellchat, slot.name = "netP") #信号网络的一个拓扑学分析

##保存
saveRDS(normal_cellchat,file="normal_cellchat.rds")

####################################### 结果

normal_cellchat <- aggregateNet(normal_cellchat)
groupSize <- as.numeric(table(normal_cellchat@idents))
groupSize
table(normal_cellchat@idents) #看有多少细胞



####看所有细胞的通讯权重和数目（总览通讯数目）
#netVisual_circle(normal_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = F, label.edge= F, title.name = "Number of interactions")  #看细胞通讯的number数目

pdf(file="net_Interaction weights_strength.pdf",width = 4,height = 5)
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength") #看细胞通讯的weight权重
dev.off()


##总览所有细胞的所有通路的发出和接收（多信号，多细胞）
pdf("总览所有细胞的所有通路的发出和接收.pdf",width = 5.5,height = 7)
ht1 <- netAnalysis_signalingRole_heatmap(normal_cellchat, 
                                         pattern = "outgoing",
                                         height = 13,
                                         width = 3.5,
                                         font.size =5,)
ht2 <- netAnalysis_signalingRole_heatmap(normal_cellchat, 
                                         pattern = "incoming",
                                         height = 13,
                                         width = 3.5,
                                         font.size =5,)
ht1 + ht2
dev.off()

####查看特定细胞之间的所有信号发出和接受关系（细胞对，信号对）
levels(normal_cellchat@idents)  #看每类细胞的序号，方便下文进行选择 
pdf(file = "特定细胞之间的所有信号.pdf",width = 8,height = 16)
nb1<-netVisual_bubble(normal_cellchat, sources.use = c(4,5), targets.use = c(1,2), remove.isolate = FALSE)#看第1类细胞对第2-10类细胞的通讯
nb2<-netVisual_bubble(normal_cellchat, sources.use = c(1,2), targets.use = c(4,5), remove.isolate = FALSE) #这个图完全可以放主图
nb1|nb2
dev.off()


####### 6 看和细胞死亡的关系 #####


