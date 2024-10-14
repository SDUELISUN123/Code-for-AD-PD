#### 1 sva ####
#引用包
library(limma)
library(sva)
setwd("E:\\科研\\24.2.29 pan_neuro_MR\\1 pd\\06.sva")      #设置工作目录

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
  if(file=="merge.preNorm.txt"){next}
  if(file=="merge.normalize.txt"){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  if(file=="merge.preNorm.txt"){next}
  if(file=="merge.normalize.txt"){next}
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  
  #数据合并
  if(i==1){
    allTab=rt[interGenes,]
  }else{
    allTab=cbind(allTab, rt[interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#输出合并后的表达数据
outTab=rbind(geneNames=colnames(allTab), allTab)
write.table(outTab, file="merge.preNorm.txt", sep="\t", quote=F, col.names=F)

#对合并后数据进行批次矫正，输出批次矫正后的表达数据
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)





#### 2 PCA ####
#引用包
library(ggplot2)
library(ggpubr)
setwd("E:/科研/24.2.29 pan_neuro_MR/1 pd/1 sva")   #设置工作目录
cors<-ggsci::pal_igv(alpha = 0.4)(51)

#定义PCA分析的函数
bioPCA=function(inputFile=null, outFile=null, titleName=null){
  #读取输入文件,提取数据
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  data=t(rt)
  Project=gsub("(.*?)\\_.*", "\\1", rownames(data))    #获取GEO数据库研究的id
  
  #PCA分析
  data.pca=prcomp(data)
  pcaPredict=predict(data.pca)
  PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)
  
  #绘制图形
  pdf(file=outFile, width=5.5, height=4.25)
  p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape="Type", 
               ellipse=T, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
               size=2, main=titleName, legend="right")+
    theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
    theme_bw()+
    scale_color_manual(values = cors)
  print(p1)
  dev.off()
}

#调用函数, 绘制批次矫正前的图形
bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf", titleName="Before batch correction")
#调用函数, 绘制批次矫正后的图形
bioPCA(inputFile="merge.normalize.txt", outFile="PCA.normalzie.pdf", titleName="After batch correction")





#### 3 boxplot showing batch effect removing ####
library(tidyverse)

setwd("E:/科研/24.2.29 pan_neuro_MR/2 AD/0 sva/")

pre <- read.table("merge.preNorm.txt", header = T, sep = "\t", row.names = 1)
after <- read.table("merge.normalize.txt", header = T, sep = "\t", row.names = 1)

pdf("pre.pdf",width = 20,height = 5)
boxplot(pre,las=2)
dev.off()

pdf("after.pdf",width = 20,height = 5)
boxplot(after,las=2)
dev.off()




#### 4 DEGs filtering, volcano plotting, and heat-map plotting ####
cors<-ggsci::pal_npg()(10)
#引用包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

logFCfilter=0.25          #logFC的过滤条件
adj.P.Val.Filter=0.05      #矫正后p值的过滤条件
inputFile="merge.normalize.txt"      #表达数据文件
setwd("E:\\科研\\24.2.29 pan_neuro_MR\\1 pd\\2 batchdiff")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#获取样品的分组信息
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,order(Type)]       #根据样品的分组信息对样品进行排序
Project=gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))    #获得GEO数据研究的id
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data)=gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

#差异分析
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("Control","Treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(Treat-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#输出所有基因的差异情况
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#输出显著的差异基因
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#输出差异基因表达量
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#绘制差异基因热图
geneNum=50     #定义展示基因的数目
diffUp=diffSig[diffSig$logFC>0,]
diffDown=diffSig[diffSig$logFC<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=data[c(geneUp,geneDown),]
#准备注释文件
names(Type)=colnames(data)
Type=as.data.frame(Type)
Type=cbind(Project, Type)
#输出图形
pdf(file="heatmap.pdf", width=5, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5.5,
         fontsize_col=8)
dev.off()


#定义显著性
rt=read.table("all.txt", header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("blue2", "grey","red2"))+
  labs(title = " ")+
  theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))+
  theme_bw()
#输出图形
pdf(file="vol.pdf", width=4, height=3.5)
print(p)
dev.off()





#### 5 MR ####
#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureFile="exposure.F.csv"         #暴露数据文件
outcomeID="finn-b-PDSTRICT_EXMORE"        #结局数据id(需修改)
outcomeName="Parkinson Disease"       #图形中展示疾病的名称
#setwd("C:\\biowolf\\geneMR\\12.MR")     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据
outcomeData=extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)
write.csv(outcomeData, file="outcome.csv", row.names=F)

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult=mr(dat)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("C:\\biowolf\\geneMR\\13.IVWfilter")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)
#提取IVW方法pvalue<0.05的基因
ivw=data.frame()
for(geneName in unique(rt$exposure)){
	geneData=rt[rt$exposure==geneName,]
	#提取五种方法OR方向一致的基因
	if(nrow(geneData)==5){
		if(geneData[geneData$method=="Inverse variance weighted","pval"]<0.05){
			if(sum(geneData$or>1)==nrow(geneData) | sum(geneData$or<1)==nrow(geneData)){
				ivw=rbind(ivw, geneData)
			}
		}
	}
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除多效性pvalue小于0.05的基因
pleRT=pleRT[pleRT$pval>0.05,]
immuneLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% immuneLists,]
write.csv(outTab, file="IVW.filter.csv", row.names=F)




#### 6 Venn ####
library(VennDiagram)      #引用包
setwd("E:\\科研\\24.2.29 pan_neuro_MR\\1 pd\\5 差异基因和mr基因交集")     #设置工作目录
upList=list()        #上调基因的列表
downList=list()      #下调基因的列表

#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
upRT=rt[rt$logFC>0,]          #提取上调基因的数据
downRT=rt[rt$logFC<0,]        #提取下调基因的数据
upGenes=as.vector(upRT[,1])         #提取上调基因名称
downGenes=as.vector(downRT[,1])     #提取下调基因名称
upList[["DEG_up"]]=upGenes
downList[["DEG_down"]]=downGenes

#读取孟德尔随机化分析结果文件
rt=read.csv("IVW.filter.csv", header=T, sep=",", check.names=F)
rt=rt[rt$method=="Inverse variance weighted",]
upRT=rt[rt$or>1,]          #提取高风险基因的数据
downRT=rt[rt$or<1,]        #提取低风险基因的数据
upGenes=unique(upRT[,"exposure"])         #提取高风险基因名称
downGenes=unique(downRT[,"exposure"])     #提取低风险基因名称
upList[["MR_or>1"]]=upGenes
downList[["MR_or<1"]]=downGenes

#绘制上调venn图
venn.plot=venn.diagram(upList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.1)
pdf(file="up.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#绘制下调venn图
venn.plot=venn.diagram(downList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.1)
pdf(file="down.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集基因的文件
upInterGenes=Reduce(intersect, upList)        #交集的上调基因
upTab=cbind(upInterGenes, "up")
downInterGenes=Reduce(intersect, downList)    #交集的下调基因
downTab=cbind(downInterGenes, "down")
interTab=rbind(upTab, downTab)
colnames(interTab)=c("Gene", "Type")
#输出交集基因的属性文件
write.table(file="interGenes.Type.txt", interTab, sep="\t", quote=F, col.names=T, row.names=F)
#输出交集基因的列表文件
write.table(file="interGenes.List.txt", interTab[,1], sep="\t", quote=F, col.names=F, row.names=F)





#### 7 boxplot and ROC ####
library(limma)

cors<-ggsci::pal_igv()(51)

#读取表达数据文件
rt=read.table("merge.normalize.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集基因的表达量
geneRT=read.table('GFRA2.txt', header=F, sep="\t", check.names=F)
data=t(data[as.vector(geneRT[,1]),])
rownames(data)=geneRT[,1]
data<-as.data.frame(t(data))

#获取样品的分组信息(对照组和实验组)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(data))
rt=cbind(as.data.frame(data), Type)
rt$Type<-gsub("Treat","PD",rt$Type)
rt$Type<-gsub("Control","Normal",rt$Type)


pdf("GFRA2exp.pdf",width = 2,height = 3)
p <- 
  ggplot(data=rt,aes(x=Type,
                     y=GFRA2,
                     color=Type))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = cors)+
  theme_classic()+ 
  theme(text = element_text(size=16)) + 
  #ylim(0.0,1.3)+
  ylab(paste0("GFRA2"," ","Expression"))+
  stat_compare_means(method = "wilcox")+
  theme_bw()+
  NoLegend()
print(p)
dev.off()

library(pROC)
pdf("ROC.pdf",width = 4,height = 4)
roc1<- roc(rt$Type, as.numeric(rt[,1]))
plot(roc1,#可以用smooth函数变得光滑
     col="red",
     #legacy.axes=T,#更改y轴格式
     #print.auc=TRUE,#显示auc面积
     #print.thres=TRUE,#添加节点和95ci
     
)
legend<-paste0("AUC=",round(auc(roc1),3))
legend("bottomright", legend=legend,
       col="red",lty=1)
dev.off()





#### 8 circlize ####
library(circlize)      #引用包
geneFile="interGenes.List.txt"      #基因列表文件
posFile="geneREF.txt"               #基因位置信息文件
setwd("E:\\科研\\24.2.29 pan_neuro_MR\\1 pd\\6 circlize")     #设置工作目录

#读取基因位置信息文件
genepos=read.table(posFile, header=T, sep="\t", check.names=F)
colnames(genepos)=c('genename','chr','start','end')
genepos=genepos[,c('chr','start','end','genename')]
row.names(genepos)=genepos[,'genename']

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
genepos=genepos[as.vector(geneRT[,1]),]
bed0=genepos

#绘制图形
pdf(file="circlize.pdf", width=6, height=6)
#初始化圈图
circos.clear()
circos.initializeWithIdeogram(species="hg38", plotType=NULL)
#展示每条染色体的注释信息
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col=rand_color(24))
  circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height=0.15, bg.border = NA)
#绘制基因组的图形
circos.genomicIdeogram(species = "hg38", track.height=mm_h(6))
#在染色体相应位置上标注基因的名称
circos.genomicLabels(bed0, labels.column=4, side = "inside", cex=0.8)
circos.clear()
dev.off()




#### 9 GSEA ####


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

gene="GFRA2"      #基因的名称
expFile="merge.normalize.txt"          #表达数据文件
gmtFile="c2.cp.kegg.Hs.symbols.gmt"    #基因集文件
setwd("E:/科研/24.2.29 pan_neuro_MR/1 pd/9 gsea")      #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=meanH-meanL
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
#根据pvalue<0.05对数据进行过滤,得到显著富集的结果
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.KEGG.txt",sep="\t",quote=F,row.names = F)

#绘制高表达组富集的图形
termNum=5     #设置展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
  pdf(file="GSEA.highExp.KEGG.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}

#绘制低表达组富集的图形
termNum=5     #设置展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
  pdf(file="GSEA.lowExp.KEGG.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}


gene="GFRA2"      #基因的名称
expFile="merge.normalize.txt"          #表达数据文件
gmtFile="c5.go.symbols.gmt"    #基因集文件
#setwd("E:/科研/24.2.29 pan_neuro_MR/1 pd/9 gsea")      #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=meanH-meanL
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
#根据pvalue<0.05对数据进行过滤,得到显著富集的结果
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.GO.txt",sep="\t",quote=F,row.names = F)

#绘制高表达组富集的图形
termNum=5     #设置展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
  pdf(file="GSEA.highExp.GO.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}

#绘制低表达组富集的图形
termNum=5     #设置展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
  pdf(file="GSEA.lowExp.GO.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}





#### 10 cibersort ####
#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

rm(list=ls())
inputFile="merge.normalize.txt"      #表达数据文件
setwd("E:/科研/24.2.29 pan_neuro_MR/1 pd/10 cibersort")#设置工作目录
library(CIBERSORT)     #引用包
gc()

#免疫细胞浸润分析
results <- cibersort(sig_matrix = "ref.txt", mixture_file = inputFile,perm = 1000,QN = T)
results <- as.data.frame(results)
results2 <-results[,1:(ncol(results)-3)]
results3 <- subset(results, results$`P-value`<0.05)
write.csv(results2,"result2.csv",quote = F,row.names = T)
write.csv(results3,"result3.csv",quote = F,row.names = T)





#### 11 ssGSEA ####
#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="merge.normalize.txt"         #表达输入文件
gmtFile="immune.gmt"            #免疫数据集文件
clusterFile="Cluster.txt"    #m6A分型结果文件
#setwd("C:\\biowolf\\m6A\\18.ssGSEA")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取m6A分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,"Cluster",drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("Cluster"))
colnames(data)=c("m6Acluster", "Immune", "Fraction")

#绘制箱线图
cors<-ggsci::pal_igv(alpha = 1)(2)
bioCol=c(cors[2],cors[1])
bioCol=bioCol[1:length(levels(factor(data[,"m6Acluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="m6Acluster",
            xlab="", 
            ylab="Immune infiltration",
            legend.title="GFRA2",
            palette=bioCol)
p=p+rotate_x_text(50)
#输出图形文件
pdf(file="boxplot.pdf", width=6, height=6)
p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()





#### 12 ceRNA network ####
library(dplyr)


setwd("E:\\科研\\24.2.29 pan_neuro_MR\\1 pd\\13 ceNetwork")

want<-read.table("gene.txt",header = F)
want<-as.vector(want$V1)

#读取miRNA-gene对
miRanda<-read.table("miRanda.tsv",sep = "\t",header = F)
colnames(miRanda)=c("miRNA","gene")
miRDB<-read.table("miRDB.tsv",sep = "\t",header = F)
colnames(miRDB)=c("miRNA","gene")
TargetScan<-read.table("TargetScan.tsv",sep = "\t",header = F)
colnames(TargetScan)=c("miRNA","gene")
miRanda_want_mi<-subset(miRanda,miRanda$gene==want)
miRDB_want_mi<-subset(miRDB,miRDB$gene==want)
TargetScan_want_mi<-subset(TargetScan,TargetScan$gene==want)
all_want_mi<-intersect(miRanda_want_mi$miRNA,miRDB_want_mi$miRNA)
all_want_mi<-intersect(all_want_mi,TargetScan_want_mi$miRNA)
all_want_mi<-all_want_mi[!duplicated(all_want_mi)]
miRNA<-all_want_mi

#读取lncRNA-miRNA对
spongeScan<-read.table("spongeScan.tsv",header=T)
all_want_miRNAs<-as.vector(miRNA)
spongeScan_want_ln<-subset(spongeScan,spongeScan$miRNAs%in%all_want_miRNAs)
lncRNA<-spongeScan_want_ln$lncRNA


#合并
##node种类
want<-as.data.frame(want)
want$Type<-"mRNA"
colnames(want)[1]="Node"
miRNA<-as.data.frame(miRNA)
miRNA$Type<-"miRNA"
colnames(miRNA)[1]="Node"
lncRNA<-as.data.frame(lncRNA)
lncRNA$Type<-"lncRNA"
colnames(lncRNA)[1]="Node"
node_matrix<-rbind(want,miRNA,lncRNA)
write.table(node_matrix,"Node.txt",sep = "\t",quote = F,row.names = F)

##interaction_pair
mRNA_miRNA_pair<-data.frame("Node1"=rep("GFRA2",length(all_want_mi)),"Node2"=all_want_mi)
mRNA_miRNA_pair$Interaction<-"mRNA"
lncRNA_miRNA_pair<-spongeScan_want_ln
colnames(lncRNA_miRNA_pair)=c("Node1","Node2")
lncRNA_miRNA_pair$Interaction<-"lncRNA"
network_matrix<-rbind(mRNA_miRNA_pair,lncRNA_miRNA_pair)
write.table(network_matrix,"Network.txt",sep = "\t",quote = F,row.names = F)

# Further plot the network in cytoscape. 



