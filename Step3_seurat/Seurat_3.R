#########################
# Load libraries
#########################
library(Seurat)
library(Matrix)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scater)
library(SeuratDisk)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocManager)
library(celldex)
library(SingleR)
library(cowplot)
library(DoubletFinder)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(stringr)
library(enrichplot)
library(topGO)
#library(GSVA)
library(pathview)
library(DOSE)
library(GSEABase)
library(EnhancedVolcano)
library(enrichR)
library(garnett)
library(rlang)


##############################
## doubletFinder 位置再QC之前
##############################
compare = ScaleData(compare,verbose =F)
compare = RunPCA(compare, verbose = F, npcs = 20)
compare = RunUMAP(compare, dims = 1:20, reduction = "pca",verbose = F)
compare = FindNeighbors(compare, dims = 1:20,verbose = FALSE)
compare = FindClusters(compare, resolution = 0.8, verbose = FALSE)

## 寻找最优pK值
sweep.res.list <- paramSweep_v3(compare, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = ncol(compare)*8*1e-6                     
homotypic.prop <- modelHomotypic(compare$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(compare)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
compare<- doubletFinder_v3(compare, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

# Predict doublet
DF.name = colnames(compare@meta.data)[grepl("DF.classification", colnames(compare@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(compare, group.by = "orig.ident") + NoAxes(),DimPlot(compare, group.by = DF.name) + NoAxes())

# Remove doublet
compare = compare[, compare@meta.data[, DF.name] == "Singlet"]
dim(compare)
# Remove
rm(sweep.res.list,sweep.stats,bcmvn)
gc()


#####################################
## Quality control ##################
#####################################
compare[["percent.mt"]] <- PercentageFeatureSet(compare, pattern = "^MT")
compare[["percent.rb"]] <- PercentageFeatureSet(compare, pattern = "^RP[SL]")

# Detection-based filtering
# consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
selected_c <- WhichCells(compare, expression =  nFeature_RNA >200)
selected_f <- rownames(compare)[Matrix::rowSums(compare) > 3]
compare <- subset(compare, features = selected_f, cells = selected_c)
compare <- subset(compare, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)



#######################################
# normalize data and find feature genes
#######################################
compare <- NormalizeData(compare)
compare <- FindVariableFeatures(compare, verbose = F,selection.method = "vst", nfeatures = 2000)
# Run the standard workflow for visualization and clustering
# 所有基因 处理 平均值为0 标准差为1  高表达基因和低表达同等对待  标准化
# scaladata cell cycle 如果有差异，就只能挑选你想要的
compare <- ScaleData(compare,verbose = F, vars.to.regress = c("nFeature_RNA", "percent.mt"))
#保留变化量 处理2000到20feature 
compare <- RunPCA(compare,npcs =20, verbose = FALSE)
# UMAP TSNE 20维度再将维度2维度 UMAP细胞多，快 全局精准  PC1 PC2
compare <- RunUMAP(compare, reduction = "pca", dims = 1:20,verbose = FALSE)


compare <- FindNeighbors(compare, dims = 1:20,verbose = FALSE)
compare <- FindClusters(compare, resolution = 0.8, verbose = FALSE)
#table(compare@meta.data$RNA_snn_res.0.8)
compare <- RunUMAP(compare, dims = 1:20, do.fast = TRUE)
DimPlot(compare,reduction = "umap",label=T)
saveRDS(compare, file="compare_seurat_1.Rds")

