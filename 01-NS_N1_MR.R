

rm(list = ls());gc()

# dir.create('./Result/01-NS_N1_MR')
setwd('./Result/01-NS_N1_MR')

color.lib <- c("#E31A1C", "#55c2fc", "#A6761D", "#F1E404", "#33A02C", "#1F78B4",
               "#FB9A99", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#F4B3BE",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
               "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
               "#634795", "#005B1D", "#26418A", "#CB8A93", "#B2DF8A", "#E22826",
               "#A6CEE3", "#F4D31D", "#F4A11D", "#82C800", "#8B5900", "#858ED1",
               "#FF72E1", "#CB50B2", "#007D9B", "#26418A", "#8B495F", "#FF394B")


#----------------------------------------------------------------------------------
#  Step 1: cell type anno
#----------------------------------------------------------------------------------


# 发现关键细胞亚群

library(Seurat)
# 读取数据

dir_name <- c(paste0('CT',1:2),paste0('S',1:2),paste0('NS',1:2))
datalist <- list()
for (i in 1:length(dir_name)){
  # i <- 3
  dir.10x <- paste0("../../Data/",dir_name[i])
  if(dir_name[i] %in% c('CT3','CT4') ){
    my.data <- Read10X(data.dir = dir.10x)[[1]]
  } else{
    my.data <- Read10X(data.dir = dir.10x)
  }
  datalist[[i]] <- CreateSeuratObject(counts = my.data,
                                      project = dir_name[i],
                                      min.cells = 10,
                                      min.features = 500)
}
names(datalist) <- dir_name

# avoid CheckDuplicateCellNames() warning
datalist <- pbapply::pblapply(1:length(datalist),FUN = function(x){
  data <- RenameCells(datalist[[x]],add.cell.id = dir_name[x])
})

scRNA <- merge(datalist[[1]],y=datalist[2:length(datalist)])
scRNA@meta.data$tissue_type=scRNA@meta.data$orig.ident
scRNA@meta.data$tissue_type=stringr::str_remove(scRNA@meta.data$tissue_type,'[0-9]')
as.data.frame(scRNA@assays$RNA@counts[1:10, 1:2])
head(scRNA@meta.data, 10)
table(scRNA@meta.data$orig.ident)
dim(scRNA)


# 质控
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA))
HB.genes <- rownames(scRNA@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
#head(scRNA@meta.data)
library(ggplot2)
col.num <- length(levels(scRNA@active.ident))
# 过滤前
VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
        cols =rainbow(col.num),
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 4) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave('01-vlnplot_before_qc.pdf',width = 15,height = 5)

##设置质控标准，比较随意
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=4000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
        cols =rainbow(col.num),
        pt.size = 0.1,
        ncol = 4) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave('02-vlnplot_after_qc.pdf',width = 15,height = 5)

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#降维聚类#########################
library(Seurat)
library(tidyverse)
library(patchwork)


scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
ggsave('03-dimplot_before_pca.pdf',width = 5.5,height = 5)

### harmony去批次
#BiocManager::
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")
ggsave('04-dimplot_after_pca.pdf',width = 5.5,height = 5)

ElbowPlot(scRNA_harmony,reduction = 'harmony')

# 一定要指定“harmony”！
scRNA <- FindNeighbors(scRNA_harmony, dims = 1:10, reduction = "harmony")
scRNA <- FindClusters(scRNA,resolution = 0.1)
scRNA <- RunUMAP(scRNA, dims = 1:10,reduction = 'harmony')

# 去批次成功
DimPlot(scRNA,split.by = 'tissue_type')
ggsave('05-dimplot_by_tiisue_type.pdf',width = 15,height = 5.5)

rm(scRNA_harmony)
#BiocManager
library(SingleR)
# 人用下面
refdata <- SingleR::HumanPrimaryCellAtlasData()

# 鼠用下面
#refdata <- SingleR::MouseRNAseqData()

library(Seurat)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)


scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(scRNA$celltype)
# Idents(scRNA)=scRNA$celltype


DimPlot(scRNA,split.by = 'tissue_type',group.by = 'celltype',cols = color.lib[1:20])
ggsave('06-dimplot_by_singleR.pdf',width = 15,height = 5.5)


library(ggthemes)
library(ggpubr)
library(cowplot)
# 人工注释--------
# find marker
scRNA <- BuildClusterTree(object = scRNA)
all.markers <- FindAllMarkers(object = scRNA, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1)
all.markers <- all.markers[which(all.markers$p_val_adj < 0.05 & all.markers$avg_log2FC > 0), ]
all.markers <- all.markers[which(all.markers$pct.1 > 0.25), ]
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unique(top10$gene)
p <- DotPlot(scRNA, features = gene.list, dot.scale = 8, cols = c("#DDDDDD", "#003366" ), col.min = -2) + RotatedAxis()
p <- p + theme_few() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
p <- p + theme(axis.text.y = element_text(size = 20))
p <- p + scale_size(range = c(1, 7))
p <- p + gradient_color(c("#EEEEEE","#ffb459","#e8613c","#b70909"))
p1 <- DimPlot(scRNA,group.by = 'RNA_snn_res.0.1',cols = color.lib[1:20],label = T)
p2 <- plot_grid(p1,p,nrow = 1,rel_widths = c(1,3))
p2
ggsave('07-top.markers.pdf', width = 25, height = 7)


gene.list.2 <-c("MS4A1", #B cells
                "CD14","LYZ", #CD14+ monocytes
                "IL7R","CCR7","CD27",#CD4+ T cells
                "CD8A",#CD8+ T cells
                "FCER1A","CST3",'CD123','GZMB',#DCs
                'GYPB','AHSP',#erythroid precursors
                'FCGR3A','MS4A7',#FCGR3A+ monocytes
                'JAML','SERPINB',#neutrophils
                'GNLY','NKG7',#NK cells
                'PPBP'#platelets
)
p <- DotPlot(scRNA, features = gene.list.2, dot.scale = 8, cols = c("#DDDDDD", "#003366" ), col.min = -2) +
  RotatedAxis()+
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))+
theme(axis.text.y = element_text(size = 20))+
scale_size(range = c(1, 7))+
gradient_color(c("#EEEEEE","#ffb459","#e8613c","#b70909"))
p1 <- DimPlot(scRNA,group.by = 'RNA_snn_res.0.1',cols = color.lib[1:20],label = T)
p2 <- plot_grid(p1,p,nrow = 1,rel_widths = c(1,3))
p2
ggsave('08-literature.markers.2.pdf', width = 20, height = 7)


# annote cell
scRNA$celltype_manual <- scRNA$RNA_snn_res.0.1
scRNA$celltype_manual <- recode(scRNA$celltype_manual,
                                '0' = 'T cells',
                                '1' = 'B cells',
                                '2' = 'Monocytes',
                                '3' = 'NK cells',
                                '4' = 'Platelets',
                                '5' = 'NK cells',
                                '6' = 'Monocytes')
Idents(scRNA) <- scRNA$celltype_manual
table(scRNA@active.ident)

# Expression for each cluster
FeaturePlot(object = scRNA, features = c('IL7R','MS4A1','CD14','NKG7','PPBP'),
                 cols = c("#CCCCCC", "red"), pt.size = 0.3, ncol = 3,
                 reduction = "umap")
ggsave('08-featurePlot.pdf', width = 14, height = 12)

VlnPlot(object = scRNA, features = c('IL7R','MS4A1','CD14','NKG7','PPBP'),
             pt.size = 0, cols = color.lib, slot = "data", ncol = 3)
ggsave( "08-VlnPlot.pdf", width = 16, height = 10)

# final umap
DimPlot(scRNA,split.by = 'tissue_type',group.by = 'celltype_manual',cols = color.lib[1:20])
ggsave('09-dimplot_final.pdf',width = 15,height = 5.5)


## 保存数据
saveRDS(scRNA,file = 'scRNA_anno.RDS')


## 细胞比例图
library(reshape2)
library(ggplot2)
prop_df <- table(scRNA@meta.data$tissue_type,scRNA@meta.data$celltype_manual)%>%as.data.frame()
prop_df <- ddply(prop_df,'Var1',transform,percent_Freq=Freq/sum(Freq)*100)
prop_df <- ddply(prop_df,'Var1',transform,lable=cumsum(percent_Freq)-0.5*percent_Freq)
prop_df$ll <- paste0(round(prop_df$percent_Freq/100,2)*100,'%')

sample_color <-  c("#E31A1C", "#55c2fc", "#A6761D", "#F1E404", "#33A02C")


ggplot(prop_df,aes(Var1,percent_Freq,fill=Var2))+
  geom_bar(stat = 'identity',width = 0.76,color = '#f3f4f4')+
  geom_text(aes(label=ll),size=3.8,
            position = position_stack(vjust = 0.5),
            color='white')+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )+coord_flip()
ggsave('10-prop_celltype_tissuetype.pdf',width = 15,height = 5.5)






#----------------------------------------------------------------------------------
#  Step 2: extract cell
#----------------------------------------------------------------------------------



### 根据原文 我们直接提取T亚群 -----------
scRNA_T=subset(scRNA,celltype_manual=='T cells')
library(Seurat)

# 提T细胞亚群,重新降维聚类
scRNAsub<- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)


## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')


scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:15)
scRNAsub <- FindClusters(scRNAsub,resolution = 0.5)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:15)

DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type')

scRNAsub = subset(scRNAsub,ident = c(0:5,7))

scRNAsub<- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)
## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNAsub,reduction = 'harmony')
scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:15)
scRNAsub <- FindClusters(scRNAsub,resolution = 0.5)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:15)
DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type')

scRNAsub = subset(scRNAsub,ident = c(0:6))

scRNAsub<- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)
## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNAsub,reduction = 'harmony')
scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:15)
scRNAsub <- FindClusters(scRNAsub,resolution = 0.5)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:15)
DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type')
# T细胞亚群注释参考


# https://blog.csdn.net/weixin_52505487/article/details/126687526
Idents(scRNAsub)=scRNAsub$seurat_clusters
T_marker <- c("CCR7","LEF1", "TCF7",'SELL','KLF2', #CD4_Naive
                "ANXA1", "CXCR4", "IL2", #CD4_EM
                "BCL6", "CXCR5","ICA1", #CD4_FH
                'IL23R',"CCR6",'CAPG','RORC','IL17A', #TH17
                'FOXP3','IL2RA','IL1R2',#CD4_REG
              'CD8A','GZMK',  # CD8_EM
              'GZMA','CCL5',  #CD8_CM
               'HAVCR2','PDCD1','LAG3', # CD8_exhau
                'EPCAM','CD19','CD3E')

T_marker <- c('CD3D','CD4' ,#CD4 T
              'CD3D','CD8A','CD8B',#CD8 T
              'CD3D','NCR3',#NKT
              'FGFBP2','NKG7',# NK
              'MKI67','STMN1',#Pro-T
              'CD4','CCR7',#CD4+ Tn
              'LTB','IL7R',#CD4+ Tpm
              'IL2RA','FOXP3',#Treg
              'CD4','CD28',#CD4+ CD28+ T
              'CD8A','CCR7',#CD8+ Tn
              'CCL5','GZMB',#CD8+ Te
              'GZMA','GZMK'#CD8+ Tem
)

T_marker <- c('CD3D', 'CD3E' ,'CD4' ,#CD4 T
              'CCR7', 'SELL', 'CD5', 'HLA-B',#CD4+ naive T cell
              'CD44', 'IFNF', 'S100A4' ,'GPR183',#CD4+ memory T cell
              'CD4' ,'FASLG', 'IFNG' ,'CD44', 'FAS',# CD4+ effector T cell
              'IL2' ,'IFNG',#CD4+ Th1
              'IL4', 'IL5' ,'IL13',#CD4+ Th2
              'IL17A' ,'IL22',#CD4+ Th17
              'IL9',	'IL22',#Treg
              'CXCR5' ,'ICOS', 'PDCD1', 'IL21', 'BCL6',#CD4+ Tfh
              'CD3D', 'CD3E' ,'CD8',#CD8+ T cell
              'CCR7' ,'SELL',#CD8+ naive T cell
              'CD8A', 'CD44', 'FASLG' ,'FAS',#CD8 effector T cell
              'TRGV9', 'TRDV2'#γδ T cell
)

DotPlot(scRNAsub,
        features = unique(T_marker),
        group.by = "seurat_clusters") + coord_flip()


## 不要求非常准确
T_celltype=c(
             'CD4_Naive',
             'CD8_T',
             'CD4_EM',
             'CD4_Naive',
             'CD4_Naive',
             'CD8_CM')

Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clust
names(T_celltype) <- levels(scRNAsub)
scRNAsub<- RenameIdents(scRNAsub, T_celltype)

scRNAsub@meta.data$T_celltype <- Idents(scRNAsub)


#设置idents主要识别标签
Idents(scRNAsub)=scRNAsub@meta.data$T_celltype

colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#bc3c29')
DimPlot(scRNAsub, group.by="T_celltype", label=F, label.size=3,cols = color.lib[7:15],
        pt.size = 0.5,split.by = 'tissue_type')
ggsave('11-dimplot_T_sub.pdf',width = 15,height = 5.5)


DimPlot(scRNAsub,group.by = 'T_celltype',cols = color.lib[7:15],label = T)
ggsave('11-dimplot_T_sub_final.pdf',width = 6,height = 5.5)

# Expression for each cluster
FeaturePlot(object = scRNAsub, features = c('CCR7','CD8A','ANXA1','CCL5'),
            cols = c("#CCCCCC", "red"), pt.size = 0.3, ncol = 4,
            reduction = "umap")
ggsave('11-featurePlot.pdf', width = 20, height = 5)

VlnPlot(object = scRNAsub, features = c('CCR7','CD8A','ANXA1','CCL5'),
        pt.size = 0, cols = color.lib[7:15], slot = "data", ncol = 4)
ggsave( "11-VlnPlot.pdf", width = 20, height = 5)




## 细胞比例图再次

library(reshape2)
library(ggplot2)
library(plyr)

prop_df <- table(scRNAsub@meta.data$tissue_type,scRNAsub@meta.data$T_celltype)%>%as.data.frame()
prop_df <- ddply(prop_df,'Var1',transform,percent_Freq=Freq/sum(Freq)*100)
prop_df <- ddply(prop_df,'Var1',transform,lable=cumsum(percent_Freq)-0.5*percent_Freq)
prop_df$ll <- paste0(round(prop_df$percent_Freq/100,2)*100,'%')

sample_color <-  color.lib[7:11]


ggplot(prop_df,aes(Var1,percent_Freq,fill=Var2))+
  geom_bar(stat = 'identity',width = 0.76,color = '#f3f4f4')+
  geom_text(aes(label=ll),size=3.8,
            position = position_stack(vjust = 0.5),
            color='white')+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )+coord_flip()

ggsave('12-prop_T_celltype_tissuetype.pdf',width = 15,height = 5.5)



saveRDS(scRNAsub,file = 'scRNA_T.RDS')


scRNA_T <- readRDS('scRNA_T.RDS')


#----------------------------------------------------------------------------------
#  Step 3: hub cell analysis
#----------------------------------------------------------------------------------

# 关键细胞亚群的深入分析


## 读取数据
scRNA_T=scRNAsub
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

## 细胞轨迹分析----------------------------

## 载入示例数据
table(scRNA_T$T_celltype)

cellinfo <- scRNA_T@meta.data


## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_T)

## run
sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'PCA',
                           start.clus = c(3,5), shrink = 0.2)



lin1 <- getLineages(sce_slingshot,
                   clusterLabels = "seurat_clusters",
                  start.clus ="0",#可指定起始细胞簇，用处不大
                   end.clus="5",#可指定终点细胞簇,用处不大
                   reducedDim = "UMAP")



dev.off()
## 可视化
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)

## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）
igraph::igraph.options(sparsematrices = FALSE)

## 曲线折线仍选一个即可
# lines(SlingshotDataSet(sce_slingshot), lwd=2,col = '#ebc03e')#,type = 'lineages'
lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = '#b42e20')


legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)
ggsave('13-slingshot.pdf',width = 5.5,height = 5.5)




# 细胞通讯----------------------
# 在SLE和衰老中分别做，看看功能是不是一样的

## cell-cell chat----
scRNA=readRDS('scRNA_anno.RDS')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)

# 加一列！！！
scRNA_other$T_celltype =scRNA_other$celltype_manual
scRNA_chat=merge(scRNA_other,scRNA_T)
rm(scRNA_T)
rm(scRNA_other)
gc()

scRNA_chat_S=subset(scRNA_chat,tissue_type=='S')


## 抽2000细胞做
set.seed(123)
a=sample(1:ncol(scRNA_chat_S),2000)
scRNA_chat_S=scRNA_chat_S[,a]

meta =scRNA_chat_S@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_S@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")

CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= T, sources.use = 'CD4_EM',
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = 'CD4_EM',
                           remove.isolate = FALSE)+coord_flip()
p_bubble

# 14-cellchat_s 6 4

## SLE
scRNA_chat_NS=subset(scRNA_chat,tissue_type=='NS')

set.seed(123)
a=sample(1:ncol(scRNA_chat_NS),2000)
scRNA_chat_NS=scRNA_chat_NS[,a]

meta =scRNA_chat_NS@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_NS@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")

CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= T, sources.use = 'CD4_EM',
                 title.name = "Number of interactions")
dev.off()

pp_bubble= netVisual_bubble(cellchat,
                            sources.use = 'CD4_EM',
                            remove.isolate = FALSE)+coord_flip()
pp_bubble

# 14-cellchat_ns 6 4





#----------------------------------------------------------------------------------
#  Step 4: hub gene
#----------------------------------------------------------------------------------

## 关键基因的孟德尔随机化分析

# S --------

library(dplyr)
# 读取数据
scRNA=readRDS('scRNA_anno.RDS')
scRNAsub=readRDS('scRNA_T.RDS')
library(Seurat)
scRNA_other=subset(scRNA,celltype !='T_cells')
scRNA_CM=subset(scRNAsub,T_celltype == 'CD4_EM')

# S --------
scRNA_other=subset(scRNA,tissue_type=='S')
scRNA_CM=subset(scRNAsub,tissue_type=='S')
##！！ 加一列
gc()
scRNA_other$T_celltype=scRNA_other$celltype_manual
scRNA_compare=merge(scRNA_other,scRNA_CM)
table(scRNA_compare$T_celltype)
## 1. CM相对其他T亚群
df_CM=FindMarkers(scRNAsub,ident.1 = 'CD4_EM',only.pos = T,logfc.threshold = 0.25)
## 2. CM相对其他细胞
df_T=FindMarkers(scRNA_compare,ident.1 = 'CD4_EM',only.pos = T,logfc.threshold = 0.25)
ss=intersect(rownames(df_CM),rownames(df_T))
ss

## 保存，后面要用
saveRDS(ss,file ='key_marker_gene_S.rds') # 50
write.csv(ss,file ='key_marker_gene_S.csv',quote = F)


# NS --------
scRNA_other=subset(scRNA,tissue_type=='NS')
scRNA_CM=subset(scRNAsub,tissue_type=='NS')
##！！ 加一列
gc()
scRNA_other$T_celltype=scRNA_other$celltype_manual
scRNA_compare=merge(scRNA_other,scRNA_CM)
table(scRNA_compare$T_celltype)
## 1. CM相对其他T亚群
df_CM=FindMarkers(scRNAsub,ident.1 = 'CD4_EM',only.pos = T,logfc.threshold = 0.25)
## 2. CM相对其他细胞
df_T=FindMarkers(scRNA_compare,ident.1 = 'CD4_EM',only.pos = T,logfc.threshold = 0.25)
ss=intersect(rownames(df_CM),rownames(df_T))
ss

## 保存，后面要用
saveRDS(ss,file ='key_marker_gene_NS.rds') #36
write.csv(ss,file ='key_marker_gene_NS.csv',quote = F)







#----------------------------------------------------------------------------------
#  Step 5: downstream analysis
#----------------------------------------------------------------------------------


### 下游分析（轨迹、细胞通讯、代谢、bulk）

## 回到单细胞看表达-----------------
scRNA=readRDS('./scRNA_anno.RDS')
load('table1.Rdata')

## 关键基因
gene=unique(table1$Gene)



library(Seurat)
library(viridis)
DotPlot(scRNA,features = gene,cols = c('#dadada','#bc3c29'))


scRNA_T=readRDS('./scRNA_T.RDS')

library(viridis)
FeaturePlot(scRNA_T,features = 'APOBEC3G',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))
FeaturePlot(scRNA_T,features = 'YWHAQ',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))


##  轨迹的基因表达分析----------------------------------

library(mixtools)
#devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
#install.packages('fastglm')
library(SingleCellExperiment)



## time and expression in CD8_CM
scRNA_cm=subset(scRNA_T,T_celltype=='CD8_CM')

cellinfo <- scRNA_cm@meta.data

## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_cm)

## run
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP',
                           start.clus = c(3,5), shrink = 0.2)





dev.off()
## 可视化
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)

## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）
igraph::igraph.options(sparsematrices = FALSE)

## 曲线折线仍选一个即可
#lines(SlingshotDataSet(sce_slingshot), lwd=2,col = 'black')#,type = 'lineages'
lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')


legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)

### 开关基因（驱动基因分析）
allexpdata <- as.matrix(scRNA_cm@assays$RNA@data);dim(allexpdata)
allcells<-colData(sce_slingshot);dim(allcells)

allcells$slingPseudotime_1

cells <- allcells[!is.na(allcells$slingPseudotime_1),];dim(cells)
expdata <- allexpdata[,rownames(cells)];dim(expdata)

#filter genes expressed in less than 5 cells
#过滤少于五个细胞表达的基因
expdata <- expdata[apply(expdata > 0,1,sum) >= 5,];dim(expdata)

rd_UMAP <- Embeddings(object = scRNA_cm, reduction = "umap");dim(rd_UMAP)#原 object = seu3obj.integrated
rd_UMAP <- rd_UMAP[rownames(cells), ];dim(rd_UMAP)
all(rownames(rd_UMAP) == colnames(expdata))

library(mixtools)
library(GeneSwitches)
library(SingleCellExperiment)


## create SingleCellExperiment object with log-normalized single cell data
## 使用对数规范化的单细胞数据创建SingleCellExperiment对象
sce <- SingleCellExperiment(assays = List(expdata = expdata))
## 添加伪时间信息
colData(sce)$Pseudotime <- cells$slingPseudotime_1
## 添加降维，例如 PCA、UMAP、tSNE
reducedDims(sce) <- SimpleList(UMAP = rd_UMAP)
sce_p1 <- sce
### 检查二值化阈值
h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}

## 二值化分析
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2)
# sce_p1 <- binarize_exp(sce_p1, ncores = 3)
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE, zero_ratio = 0.65, ds_cutoff = 0.65)

table(rowData(sce_p1)$prd_quality)

# 过滤出开关基因
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92);dim(sg_allgenes)


sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92,
                                genelists = gs_genelists);dim(sg_gtypes)#, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),]);dim(sg_vis)

## 自己关注的基因
gl <- unique(table1$Gene)
intersect(sg_vis$geneID, gl)
sg_my <- rowData(sce_p1)[gl,];head(sg_my)
sg_my$feature_type <- "Mendelian genes"
sg_vis <- rbind(sg_vis, sg_my)
plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)

## R2大于0，上调型开关基因，R2小于0 ，下调型开关基因
a=sce_p1@assays@data$expdata['APOBEC3G',] ## !!!
b=sce_p1$Pseudotime

df=data.frame(gene=a,time=b)

ggstatsplot::ggscatterstats(data=df,x='time',y='gene')



####细胞通讯--------阴阳性群
gc()

scRNA=readRDS('./scRNA_anno.RDS')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)
gc()
scRNA_T=readRDS('./scRNA_T.RDS')


gc()

##!!!
scRNA_CM=subset(scRNA_T,T_celltype=='CD8_CM')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['APOBEC3G',]>0,'APOBEC3G+CD8_CM','APOBEC3G-CD8_CM')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_CM')

# 加列
scRNA_other$gene_group =scRNA_other$celltype
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_chat=merge(scRNA_CM,c(scRNA_other,scRNA_otherT))

rm(scRNA_CM,scRNA_otherT)
rm(scRNA_T)
rm(scRNA_other)
gc()

##
scRNA_chat_SLE=subset(scRNA_chat,tissue_type=='SLE')

set.seed(123)
a=sample(1:ncol(scRNA_chat_SLE),2000)
scRNA_chat_SLE=scRNA_chat_SLE[,a]

meta =scRNA_chat_SLE@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_SLE@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "gene_group")

CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= T, sources.use = c('APOBEC3G+CD8_CM','APOBEC3G-CD8_CM'),
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('APOBEC3G+CD8_CM','APOBEC3G-CD8_CM'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble


### 代谢相关------------------

scRNA_T=readRDS('./scRNA_T.RDS')

gc()
scRNA_CM=subset(scRNA_T,T_celltype=='CD8_CM')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['APOBEC3G',]>0,'APOBEC3G+CD8_CM','APOBEC3G-CD8_CM')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_CM')

# 加列
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_metab=merge(scRNA_CM,c(scRNA_otherT))

gc()

rm(scRNA_chat)
gc()
scRNA_metab_SLE=subset(scRNA_metab,tissue_type=='SLE')

set.seed(123)
a=sample(1:ncol(scRNA_metab_SLE),2000)
scRNA_metab_SLE=scRNA_metab_SLE[,a]




#### scMetabolism评估巨噬细胞代谢活性
#devtools::install_github("YosefLab/VISION")
#devtools::install_github("wu-yc/scMetabolism")

library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_metab_SLE<-sc.metabolism.Seurat(obj = scRNA_metab_SLE, method = 'AUCell', imputation = F, ncores = 2, metabolism.type = "KEGG")

input.pathway <- rownames(scRNA_metab_SLE@assays[["METABOLISM"]][["score"]])[61:90]
DotPlot.metabolism(obj =scRNA_metab_SLE,
                   pathway = input.pathway, phenotype = "gene_group", norm = "y")


gc()
##差异基因------------------------------------
library(Seurat)
## Idents()
Idents(scRNA_CM)=scRNA_CM$gene_group
df=FindAllMarkers(scRNA_CM,only.pos = T,logfc.threshold =0.25)
write.csv(df,'APOBEC_marker.csv',quote = F)

## 富集分析怎么做,负值csv中的基因列到网站，参见下面教程
#https://mp.weixin.qq.com/s/ClHOFvw3GSM9wvmIPip4VA


## bulk---------------------------------------------
gc()

library(data.table)
rt=fread('./GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt',data.table = F)
rownames(rt)=rt$V1
rt$V1=NULL


#IOBR有一系列依赖包
#https://mp.weixin.qq.com/s/nVziQeInS-4QxVNPQCilVQ
#   devtools::install_github("IOBR/IOBR")
library(IOBR)
#
#gc()
#rm(scRNA_chat_SLE)
#rm(scRNA_CM)
#rm(scRNA_T)
gc()

rt=count2tpm(rt,idType = 'Ensembl',org = 'hsa')

rt=as.data.frame(rt)

max(rt)
rt=log2(rt+1)

a1=grep('SLE',colnames(rt))

exp1=rt[,a1]
exp2=rt[,-a1]

rt=cbind(exp2,exp1)


load('table1.Rdata')

data=rt[unique(table1$Gene),]

### 注释文件
anno=data.frame(row.names =colnames(rt),group=c(rep('Healthy',58),
                                                rep('SLE',62)))
pheatmap::pheatmap(data,cluster_cols = F,
                   scale = 'row',show_colnames = F,annotation_col = anno)

df=data.frame(gene=as.numeric(rt['APOBEC3G',]),group=anno$group)
ggpubr::ggboxplot(data = df, x = 'group',y='gene',color = 'group',palette = 'jco',notch = T,size = 1)+
  stat_compare_means()



#----------------------------------------------------------------------------------

setwd('../..')
getwd()
