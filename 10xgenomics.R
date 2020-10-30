library(Seurat)
library(dplyr)
library(magrittr)

# 载入数据，为cellranger count的输出
list.files("./filtered_feature_bc_matrix/")
raw_data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
# 为6498个细胞中21027个基因的count值矩阵
dim(raw_data)
# [1] 21027  6498

# 创建seurat对象
# 基因至少在min.cells个细胞中表达
# 每个细胞中至少表达min.genes个基因
sc.object <- CreateSeuratObject(counts = raw_data, project = "macaca_brain", min.cells = 3, min.features = 200 )
# dim(sc.object)
# [1] 15282  6498

#画出Count和feature数的分布图
FeatureScatter(object = sc.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sc.object[["percent.mt"]] <- PercentageFeatureSet(sc.object, pattern = "^MT-")
VlnPlot(sc.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#根据feature图删除离群值,选择feature大于375小于5000的细胞
sc.object <- subset(sc.object, subset = nFeature_RNA > 375 & nFeature_RNA < 5000 & percent.mt < 5)
#dim(sc.object)
#[1] 15282  6479

#对数据进行归一化处理
sc.object <- NormalizeData(sc.object, normalization.method = "LogNormalize", scale.factor = 10000)

#分析差异表达的基因，即在数据集中找出在一些细胞中高表达一些细胞中低表达的基因,选取了方差最大的2000个feature做下游的分析
sc.object <- FindVariableFeatures(sc.object, selection.method = "vst", nfeatures = 2000)
#方差最大的十个基因
top10 <- head(VariableFeatures(sc.object), 10)
# 画出这些差异表达基因
plot1 <- VariableFeaturePlot(sc.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data线性转化 
#scale就是将数据归一化为均值为0方差为1的数据
all.genes <- rownames(sc.object)
sc.object <- ScaleData(sc.object, features = all.genes)
#pca降维
sc.object <- RunPCA(sc.object, features = VariableFeatures(object = sc.object))
print(sc.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc.object, dims = 1:2, reduction = "pca")
DimPlot(sc.object, reduction = "pca")

#如何决定数据集的维度如何选择PCs
#在选择几个PC进行下游的分析时，我们可以通过热图查看这些方差的来源
DimHeatmap(sc.object, dims = 1:15, cells = 500, balanced = TRUE)

#用JackStraw随机产生了一个数据集，又跑了一次PCA，构建了一个空集对照，发现显著的PC多集中在低P值的位置
#nu.object <- JackStraw(vt.object, num.replicate = 100)
#nu.object <- ScoreJackStraw(nu.object, dims = 1:20)
#JackStrawPlot(nu.object, dims = 1:15)
#这需要很长的时间，可以用ElowPlot代替
ElbowPlot(sc.object)
#用三个方法选择PCs
#1.观察PCs热图找到异质性的来源，这个热图怎么看
#2.用随机的对照数据集统计，耗时而且无法确定合适的cutoff
#3.ElowPlot直接选择PC的拐点
#选择不同的PCs会对下游的结果产生巨大的影响


#细胞聚类分析
#Seurat采用基于图的聚类方法，简而言之这些方法将单元格嵌入图结构中，
#例如K近邻（KNN）图，在具有相似特征表达模式的单元格之间绘制边，
#然后尝试将该图划分为高度互连的“准斜体”或“社区”。

#首先基于PCA空间中的欧式距离构建一个KNN图
#然后细化他们的边缘权重（相似性
#运用了一个模块优化方法Louvain algorithm (default) or SLM
#选了15个PCs
sc.object <- FindNeighbors(sc.object, dims = 1:15)
sc.object <- FindClusters(sc.object, resolution = 0.5)
#可以看一下前5个细胞属于哪一类
head(Idents(sc.object), 5)

#非线性降维
#例如t-SNE和umap这样的非线性降维方法会考虑数据的各个方面
#将他们降维到一个低维空间里
#上面基于图像的聚类方法应该和这些低维图共定位
sc.object <- RunUMAP(sc.object, dims = 1:15)
DimPlot(sc.object, reduction = "umap")

#找到差异表达基因
#对每个cluster定义一个阳性marker和阴性marker
#一个cluster和所有细胞对比找出差异表达基因
#比如cluster1的marker基因，min.pct参数设置两组细胞中被检测的feature所占的比例
cluster1.markers <- FindMarkers(sc.object, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#找到cluster5中能与cluster0到3区分开的基因
cluster5.markers <- FindMarkers(sc.object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# 找到每个cluster中相比其他剩余所有cluster的marker基因，且只保留positive
sc.object.markers <- FindAllMarkers(sc.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sc.object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#用一些参数对差异表达基因进行test，比如ROC test回归每个marker的统计力，从0到1的数值
cluster1.markers <- FindMarkers(sc.object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#用VlnPlot展示这个基因在所有cluster中的表达概率分布
#featurePlot可以在tSNE或这PCA上标记feature
VlnPlot(sc.object, features = c("C1QL2", "PHLDB2", "GRIK1", "UACA", "SLC1A2", "PDE1A"), pt.size = 0)
FeaturePlot(sc.object, features = c("C1QL2", "GRIK1"))
#Doheatmap可以画出指定基因的表达热图
top5 <- sc.object.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(sc.object, features = top5$gene)
top10 <- sc.object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10, file = "./top10marker", sep = "\t")

#指定细胞类型去识别cluster
#需要凭借背景知识判断marker之后属于哪类细胞
#常用的注释类型包有singleR，数据库有CellMarker
#markers根据已知文献查询得到脑组织常见细胞类型的makers
markers <- c("GAD1","MAG","A2M","SLC1A2","GAD2","SATB2","TAC1","PCDH8", "DRD2","ADORA2A","PENK","PCP4","NECAB2","LMO7","CALB1","PDGFRA","CSPG4","GJA1","MOBP","MOG","RELN","AIF1","CX3CR1","PTPRC","HLA−DRA","A2M")
markers_exist <- intersect(x=sc.object.markers$gene, y = markers)
VlnPlot(sc.object, features = markers_exist, pt.size = 0)
#为了方便改名，另用一个final.object记载改名后的结果
final.object <- RenameIdents(object = sc.object, 
                                  "0" = "undefined",
                                  "1" = "excitatory neuron",
                                  "2" = "Purkinje cell_1",
                                  "3" = "undefined",
                                  "4" = "Purkinje cell_2",
                                  "5" = "inhibitory neurons_1",
                                  "6" = "Astrocyte_1",
                                  "7" = "Oligodendrocyte_3",
                                  "8" = "inhibitory neurons_2",
                                  "9" = "Monocyte",
                                  "10" = "Astrocyte_2",
                                  "11" = "endothelial vascular cells",
                                  "12" = "Oligodendrocyte_4",
                                  "13" = "Oligodendrocyte progenitor cell")

DimPlot(final.object, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)
# 改成缩写，或者在注释出现问题时方便修改
final.object <- RenameIdents(object = sc.object, 
                          "0" = "undefined",
                          "1" = "excitatory neuron",
                          "2" = "Purkinje cell_1",
                          "3" = "undefined",
                          "4" = "Purkinje cell_2",
                          "5" = "inhibitory neurons_1",
                          "6" = "Astro_1",
                          "7" = "OLIG_1",
                          "8" = "inhibitory neurons_2",
                          "9" = "Monocyte",
                          "10" = "Astro_2",
                          "11" = "endothelial vascular cells",
                          "12" = "OLIG_2",
                          "13" = "OPC")
DimPlot(final.object, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, pt.size = 1)

