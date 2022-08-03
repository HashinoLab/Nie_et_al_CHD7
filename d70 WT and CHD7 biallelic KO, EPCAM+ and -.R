
### d70 WT (#D5e) and CHD7 KO/KO (#A11a), FACS sorted EpCAM+ and EpCAM- cells


### Load the libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load in Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D5e_pos <- Read10X("D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/mapping/D5e EpCAM+") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
A11a_pos <- Read10X("D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/mapping/A11a EpCAM+") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
D5e_neg <- Read10X("D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/mapping/D5e EpCAM-") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
A11a_neg <- Read10X("D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/mapping/A11a EpCAM-") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)


#Add metadata
D5e_pos$Genotype <- "D5e"
A11a_pos$Genotype <- "A11a"
D5e_neg$Genotype <- "D5e"
A11a_neg$Genotype <- "A11a"

D5e_pos$FACS <- "Ep+"
A11a_pos$FACS <- "Ep+"
D5e_neg$FACS <- "Ep-"
A11a_neg$FACS <- "Ep-"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Filtering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D5e_pos$mt.per <- PercentageFeatureSet(D5e_pos, pattern = "^MT-")
VlnPlot(D5e_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


A11a_pos$mt.per <- PercentageFeatureSet(A11a_pos, pattern = "^MT-")
VlnPlot(A11a_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


D5e_neg$mt.per <- PercentageFeatureSet(D5e_neg, pattern = "^MT-")
VlnPlot(D5e_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


A11a_neg$mt.per <- PercentageFeatureSet(A11a_neg, pattern = "^MT-")
VlnPlot(A11a_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


D5e_pos <- subset(D5e_pos, subset = nFeature_RNA > 1000 & nFeature_RNA < 8500 & nCount_RNA < 48000 & mt.per < 15)

A11a_pos <- subset(A11a_pos, subset = nFeature_RNA > 1000 & nFeature_RNA < 9500 & nCount_RNA < 80000 & mt.per < 15)

D5e_neg <- subset(D5e_neg, subset = nFeature_RNA > 1300 & nFeature_RNA < 7200 & nCount_RNA < 42000 & mt.per < 15)

A11a_neg <- subset(A11a_neg, subset = nFeature_RNA > 900 & nFeature_RNA < 7500 & nCount_RNA < 50000 & mt.per < 12)


VlnPlot(D5e_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(A11a_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(D5e_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(A11a_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%   Merge Datasets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list_d70pos <- list(D5e_pos, A11a_pos)

d70pos <- merge(list_d70pos[[1]], list_d70pos[[2]], add.cell.ids = c("D5eP", "A11aP"))

head(colnames(d70pos))



list_d70neg <- list(D5e_neg, A11a_neg)

d70neg <- merge(list_d70neg[[1]], list_d70neg[[2]], add.cell.ids = c("D5eN", "A11aN"))

head(colnames(d70neg))


list_d70 <- list(d70pos, d70neg)

d70 <- merge(list_d70[[1]], list_d70[[2]])


rm(D5e_pos, A11a_pos, D5e_neg, A11a_neg, d70pos, d70neg, list_d70, list_d70pos, list_d70neg)

saveRDS(d70, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- Merge before UMAP V220515.rds")


d70 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- Merge before UMAP V220515.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%  Seurat workflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d70 <- SCTransform(d70, vars.to.regress = "mt.per")
d70 <- RunPCA(d70)
DimPlot(d70)

saveRDS(d70, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- Merge after PCA V220515.rds")

d70 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- Merge after PCA V220515.rds")


nPC <- 30
res <- 0.35


d70 <- RunUMAP(d70, dims = 1:nPC)

d70 <- FindNeighbors(d70, dims = 1:nPC)
d70 <- FindClusters(d70, resolution = res)
DimPlot(d70, pt.size = 0.8)

DimPlot(d70, group.by ="Genotype", pt.size = 0.8)
DimPlot(d70, group.by ="FACS", pt.size = 0.8)


rm(nPC, res)

FeaturePlot(d70, features = "nFeature_RNA", pt.size = 1)


saveRDS(d70, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- Merge after UMAP res=0.35 V220516.rds")




#%%%%%%%%%%%%%%%%%%%%%%%%%  Feature Plots, etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Hair cell genes
FeaturePlot(d70, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "FBXO2", "ESPN", "GFI1", "SOX2", 
                                 "CD164L2", "CCER2", "GNG8", "EPCAM", "CDH1", "CDH2"), pt.size = 0.2, ncol = 4)


# Supporting cell and otic lineage genes
FeaturePlot(d70, features = c(
  "FBXO2","SOX2", "EPCAM", "CDH1", "CDH2", "S100A1", "OC90", "CLDN6", "SOX10", "JAG1", "SIX1", "SPARCL1", "CNTN1", 
  "LMX1A", "KRT23","SPP1", "GATA3", "COL9A2", "ARHGEF19", "PLEKHB1", "TBX2", "MARVELD3", "BDNF", "SIX1", "FGF10",
  "ZPLD1", "OTOG"
                              ), pt.size = 0.2, ncol = 5)


# neuronal genes
FeaturePlot(d70, features = c("TUBB3", "NEUROD1", "STMN2", "STMN4", "NEFM", "NEFL", "NEFH", 
                              "ELAVL2", "ELAVL3", "ELAVL4"), pt.size = 0.2, ncol = 3)

# fibroblast genes
FeaturePlot(d70, features = c("MGP", "COL3A1", "COL1A1", "COL1A2", "POSTN", "LGALS1", "LUM"), pt.size = 0.2, ncol = 4)


# muscle genes
FeaturePlot(d70, features = c("ACTC1", "MYLPF", "TNNI1", "MYL1", "TNNC2", "TTN", "MYH3"), pt.size = 0.2, ncol = 3)

# cycling genes
FeaturePlot(d70, features = c("CENPF", "TOP2A", "UBE2C", "NUSAP1", "HMGB2", "CDK1"), pt.size = 0.2, ncol = 3)

# CNC genes
FeaturePlot(d70, features = c("PLP1", "CRYAB", "S100A6", "MPZ", "MEF2C", "CRABP1"), pt.size = 0.2, ncol = 3)


# keratinocyte genes
FeaturePlot(d70, features = c("KRT17", "KRT19", "S100A9", "KRT15", "KRT6A", "KRT13"), pt.size = 0.2, ncol = 3)

# melanocyte genes
FeaturePlot(d70, features = c("DCT", "PMEL", "MLANA", "TYRP1", "TYR", "EDNRB", "MITF"), pt.size = 0.2, ncol = 3)




DotPlot(d70, features = c(
  # hair cell - c5
  "ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "ESPN", "TMPRSS3", "BDNF","CD164L2", "CCER2", "GNG8", "POU4F3",
  # otic - c1 WT SCs, c2 KO SC-like, c11 KO otic-like?
  "FBXO2", "SOX2","OC90", "S100A1", "COL9A2", "TBX2", "PLEKHB1", "MARVELD3", 
  "EPCAM", "CDH1", "CLDN6", "SIX1", "DLX5", "SPARCL1", "BRICD5",
  # keratinocytes? - c6
  "S100A9", "KRT5", "KRT15", "KRT17", "KRT18", 
  # neural crest - c14
  "S100B", "PLP1", "SOX10", "CDH19", "MPZ", "FOXD3",
  # cycling - c12
  "UBE2C", "TOP2A", "CENPF",  "NUSAP1", "MKI67", "BIRC5", 
  # muscle 1 - c9
  "MSC", "MYF5", "TNNT1", "DES", "ACTC1", 
  # muscle 2 - c15
  "MYLPF", "TNNI1", "MYL1", "TNNC2", "MYH3",
  # neuron 1, neuron 2 - c7 KO, c8 WT
  "STMN2", "TUBB3", "DCX", "STMN4", "ELAVL3", 
  "NEUROD1", "POU4F1", "ELAVL2",
  # neuron 3 - c10 WT
  "INSM1", "KIF5C", "GPM6A", "NFIB", "NRN1",
  # glial - c4
  "GFAP", "NTRK2", "ID4", "PTPRZ1", "DMD",
  # chondracytes - c13
  "MATN1", "MATN3", "ACAN", "COL9A1", "COL11A1", 
  # fibroblasts
  "COL3A1", "POSTN", "COL1A1",  "COL6A3", "LUM"
), dot.scale = 10) + RotatedAxis()


# find all markers of cluster xxxxx
cluster0.markers <- FindMarkers(d70, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 100)

write.csv(cluster0.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster0.csv")


cluster1.markers <- FindMarkers(d70, ident.1 = 1, min.pct = 0.25)
write.csv(cluster1.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster1.csv")


cluster2.markers <- FindMarkers(d70, ident.1 = 2, min.pct = 0.25)
write.csv(cluster2.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster2.csv")


cluster3.markers <- FindMarkers(d70, ident.1 = 3, min.pct = 0.25)
write.csv(cluster3.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster3.csv")


cluster4.markers <- FindMarkers(d70, ident.1 = 4, min.pct = 0.25)
write.csv(cluster4.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster4.csv")


cluster5.markers <- FindMarkers(d70, ident.1 = 5, min.pct = 0.25)
write.csv(cluster5.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster5.csv")


cluster6.markers <- FindMarkers(d70, ident.1 = 6, min.pct = 0.25)
write.csv(cluster6.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster6.csv")


cluster7.markers <- FindMarkers(d70, ident.1 = 7, min.pct = 0.25)
write.csv(cluster7.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster7.csv")


cluster8.markers <- FindMarkers(d70, ident.1 = 8, min.pct = 0.25)
write.csv(cluster8.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster8.csv")


cluster9.markers <- FindMarkers(d70, ident.1 = 9, min.pct = 0.25)
write.csv(cluster9.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster9.csv")


cluster10.markers <- FindMarkers(d70, ident.1 = 10, min.pct = 0.25)
write.csv(cluster10.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster10.csv")


cluster11.markers <- FindMarkers(d70, ident.1 = 11, min.pct = 0.25)
write.csv(cluster11.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster11.csv")


cluster12.markers <- FindMarkers(d70, ident.1 = 12, min.pct = 0.25)
write.csv(cluster12.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster12.csv")


cluster13.markers <- FindMarkers(d70, ident.1 = 13, min.pct = 0.25)
write.csv(cluster13.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster13.csv")


cluster14.markers <- FindMarkers(d70, ident.1 = 14, min.pct = 0.25)
write.csv(cluster14.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster14.csv")


cluster15.markers <- FindMarkers(d70, ident.1 = 15, min.pct = 0.25)
write.csv(cluster15.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster15.csv")


cluster16.markers <- FindMarkers(d70, ident.1 = 16, min.pct = 0.25)
write.csv(cluster16.markers, "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/cluster16.csv")



#%%%%%%%%%%%%%%%%%%%%%%%%% Subsetting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#-------------- Subsetting HCs, SCs, and Otic-like mutant cells --------------------

d70_subA <- subset(d70, idents = c(1, 2, 5, 11))
rm(d70)

saveRDS(d70_subA, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subA (1 2 5 11) v220526a.rds")

d70_subA <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subA (1 2 5 11) v220526a.rds")



d70_subA@active.ident <- d70_subA$Genotype %>% as.factor()


DimPlot(d70_subA)

d70_subA@active.ident <- d70_subA$seurat_clusters

DimPlot(d70_subA)



FeaturePlot(d70_subA, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "POU4F3", "ESPN", "TMPRSS3", "BDNF",
                             "FBXO2", "SOX2","OC90", "S100A1", "COL9A2", "TBX2", "SPARCL1", "EPCAM", "CDH1", "CDH2", 
                             "CLDN6", "ATP1B1", "JAG1", "SIX1", "SOX10"
), pt.size = 0.2, ncol = 6)



FeaturePlot(d70_subA, features = c(
  "SOX2", "EPCAM", "SPARCL1", "S100B", "ANXA4", "PCP4"), pt.size = 0.5, ncol = 3)


FeaturePlot(d70_subA, features = c(
  "SOX2", "FBXO2", "EPCAM", "CDH1", "S100A1", "OC90", "CLDN6", "SOX10", "JAG1", "SIX1", "SPARCL1", 
  "SPP1", "TBX2", "SIX1", 
  "ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "POU4F3", "ESPN", "TMPRSS3", "BDNF"), pt.size = 0.5, ncol = 6)


#----------- subsetting WT HCs and CHD7 KO ANXA4 high cells ----------------

d70_subHC <- subset(d70, idents = c(5, 11))
rm(d70)


saveRDS(d70_subHC, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subHC (5 11) v220527a.rds")


d70_subHC <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subHC (5 11) v220527a.rds")



d70_subA@active.ident <- d70_subA$Genotype %>% as.factor()


DimPlot(d70_subA)

d70_subA@active.ident <- d70_subA$seurat_clusters

DimPlot(d70_subA)



#----------- subsetting WT SCs and CHD7 KO ANXA4 low cells ----------------

d70_subSC <- subset(d70, idents = c(1, 2))
rm(d70)


saveRDS(d70_subSC, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subSC (1 2) v220527a.rds")


d70_subSC <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subSC (1 2) v220527a.rds")



d70_subSC@active.ident <- d70_subSC$Genotype %>% as.factor()


DimPlot(d70_subSC)

d70_subA@active.ident <- d70_subA$seurat_clusters

DimPlot(d70_subA)





#%%%%%%%%%%%%%%%%%%%%%%%%%% DEA with DESeq2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")


install.packages("remotes")
remotes::install_github("statOmics/zingeR")




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


library(tidyverse, quietly = TRUE)
library(scRNAseq, quietly = TRUE)

library(SingleCellExperiment, quietly = TRUE)

library(matrixStats, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)

d70_subHC_sce <- SummarizedExperiment(as.matrix(d70_subHC@assays$RNA@counts),
                                   colData = d70_subHC@meta.data)

d70_subSC_sce <- SummarizedExperiment(as.matrix(d70_subSC@assays$RNA@counts),
                                      colData = d70_subSC@meta.data)

saveRDS(d70_subHC_sce, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subHC (5 11) sce v220527a.rds")

saveRDS(d70_subSC_sce, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subSC (1 2) sce v220527a.rds")

rm(d70_subHC) 
rm(d70_subSC) 

filter <- rowSums(assay(d70_subHC_sce)) > 5
d70_subHC_sce <- d70_subHC_sce[filter,]
vars <- assay(d70_subHC_sce) %>% log1p() %>% rowVars()
names(vars) <- rownames(d70_subHC_sce)
vars <- sort(vars, decreasing = T)

d70_subHC_sce <- d70_subHC_sce[names(vars[1:15000]),]
assayNames(d70_subHC_sce)[1] <- "counts"


d70_subHC.colData <- data.frame(Genotype = colData(d70_subHC_sce)$Genotype)
iCounts <- assay(d70_subHC_sce, i = "counts")

d70_subHC.colData$Genotype <- as.factor(d70_subHC.colData$Genotype)

d70_subHC.colData$Genotype <- relevel(d70_subHC.colData$Genotype, ref = "D5e")

design <- model.matrix(~d70_subHC.colData$Genotype)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = d70_subHC.colData, design = ~ Genotype)

memory.limit()

rm(d20_OP_sce)

saveRDS(d20_OP.colData, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP d20_OP.colData V210731.rds")
saveRDS(design, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP design V210731.rds")
saveRDS(iCounts, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP iCounts V210731.rds")


d20_OP.colData <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP d20_OP.colData V210731.rds")
design <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP design V210731.rds")
iCounts <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP iCounts V210731.rds")
library(DESeq2, quietly = TRUE)

# -----------------------------------------------------------

weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                 design = design,
                                 maxit = 500, normalization = "DESeq2_poscounts",
                                 colData = d70_subHC.colData,
                                 designFormula = ~ Genotype, 
                                 verbose = TRUE)

saveRDS(weights, file="D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/RStudio files/d70 D5e A11a EpCAM+- subHC (5 11) weights v220527a.rds")

# To save on IU Carbonate:
saveRDS(weights, file="/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/d70 D5eE5b nT+nT- Merge HC weights V210819.rds")


weights <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP weights V210731.rds")

assays(dse)[["weights"]] <- weights

dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = FALSE) # Can be parallelized on OSX/Linux


full.results <- results(dse, parallel = FALSE)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Analysis/d70 D5e A11a EpCAM+- subHC (5 11) DESeq v220527a.csv")



#  HC Volcano Plots from DESeq2 .csv file   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")



library(EnhancedVolcano)
library(data.table)
library(Seurat)

labeled.genes <- c("PCP4", "CCER2", "ATOH1", "BDNF",
                   "SOX2", "TMPRSS3", "ESPN", "bGHpA", "GFI1", "MYO6", "MYO7A", "PLEKHB1", "FBXO2",
                   "S100A6", "S100A14", "S100A11", "S100A10", "S100A2", "S100A16")




logfc.threshold <- 1
p.threshold <- 1E-10

KO <- "#D6604D"
WT <- "#4393C3"
Highlight <- "black"


d20_OP_DE <- fread(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset HC/d70 D5e A11a EpCAM+- subHC (5 11) DESeq v220527a.csv")


d20_OP_DE <- d20_OP_DE[!is.na(d20_OP_DE$padj),] #Have to remove rows where the padj is NA, otherwise a bunch of stuff breaks

d20_OP_DE <- rbind(d20_OP_DE, d20_OP_DE[(d20_OP_DE$gene %in% labeled.genes), ])
d20_OP_DE[-(1:(nrow(d20_OP_DE)-length(labeled.genes))), 2] <- paste0(d20_OP_DE[[2]][-(1:(nrow(d20_OP_DE)-length(labeled.genes)))], ".1")

color.key <- ifelse(d20_OP_DE$log2FoldChange > logfc.threshold & d20_OP_DE$padj < p.threshold, KO,
                    ifelse(d20_OP_DE$log2FoldChange < -logfc.threshold & d20_OP_DE$padj < p.threshold, WT,
                           "lightgrey"))

color.key[d20_OP_DE$gene %in% c(labeled.genes, paste0(labeled.genes, ".1"))] <- Highlight

names(color.key)[color.key == KO] <- "KO"
names(color.key)[color.key == WT] <- "WT"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == Highlight] <- "Highlight"



custom.VP.jpg <- EnhancedVolcano(d20_OP_DE,
                                 lab = d20_OP_DE$gene,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 1e-10,
                                 ### FCcutoff = 7, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                 pointSize = 4,
                                 colCustom = color.key,
                                 selectLab = labeled.genes,
                                 legendPosition = "below",
                                 labSize = 6.0,
                                 drawConnectors = T,
                                 title = expression(paste("Day 70: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-8, 9)) +
  NoLegend()



custom.VP.eps <- EnhancedVolcano(d20_OP_DE,
                                 lab = d20_OP_DE$gene,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 1e-10,
                                 ### FCcutoff = 7, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                 pointSize = 4,
                                 colCustom = color.key,
                                 selectLab = labeled.genes,
                                 legendPosition = "right",
                                 labSize = 6.0,
                                 drawConnectors = T,
                                 colAlpha = 1, #You can't generate .eps files with transparent color values. This ensures no points are transparent.
                                 title = expression(paste("Day 70 hair cells and HC-like cells: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-8,9)) +
  NoLegend()

ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset HC/d70 HC DESeq2 VolcanoPlot V220527a.jpg",
       plot = custom.VP.jpg,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = jpeg())


ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset HC/d70 HC DESeq2 VolcanoPlot V220527a.eps",
       plot = custom.VP.eps,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = "eps")





# SC  Volcano Plots from DESeq2 .csv file   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")



library(EnhancedVolcano)
library(data.table)
library(Seurat)
labeled.genes <- c("OTOL1", "BRICD5", "SPARCL1", "PLP1",
                   "COL9A2", "SOX2", "SLITRK6", "OC90", "HES5", "OTOG", "GAS2", "CLDN14", "GJB6", "TBX1",
                   "KRT7", "KRT6B", "KRT17", "ANXA1", "LAMC2", "LAMB3", "CST6", "KRT13", "KRT14")




logfc.threshold <- 1
p.threshold <- 1E-10

KO <- "#D6604D"
WT <- "#4393C3"
Highlight <- "black"


d20_OP_DE <- fread(file = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset SC/d70 D5e A11a EpCAM+- subSC (1 2) DESeq2 V220527a.csv")


d20_OP_DE <- d20_OP_DE[!is.na(d20_OP_DE$padj),] 

d20_OP_DE <- rbind(d20_OP_DE, d20_OP_DE[(d20_OP_DE$gene %in% labeled.genes), ])
d20_OP_DE[-(1:(nrow(d20_OP_DE)-length(labeled.genes))), 2] <- paste0(d20_OP_DE[[2]][-(1:(nrow(d20_OP_DE)-length(labeled.genes)))], ".1")

color.key <- ifelse(d20_OP_DE$log2FoldChange > logfc.threshold & d20_OP_DE$padj < p.threshold, KO,
                    ifelse(d20_OP_DE$log2FoldChange < -logfc.threshold & d20_OP_DE$padj < p.threshold, WT,
                           "lightgrey"))

color.key[d20_OP_DE$gene %in% c(labeled.genes, paste0(labeled.genes, ".1"))] <- Highlight

names(color.key)[color.key == KO] <- "KO"
names(color.key)[color.key == WT] <- "WT"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == Highlight] <- "Highlight"



custom.VP.jpg <- EnhancedVolcano(d20_OP_DE,
                                 lab = d20_OP_DE$gene,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 1e-10,
                                 ### FCcutoff = 7, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                 pointSize = 4,
                                 colCustom = color.key,
                                 selectLab = labeled.genes,
                                 legendPosition = "below",
                                 labSize = 6.0,
                                 drawConnectors = T,
                                 title = expression(paste("Day 70: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-5, 6)) +
  NoLegend()



custom.VP.eps <- EnhancedVolcano(d20_OP_DE,
                                 lab = d20_OP_DE$gene,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 1e-10,
                                 ### FCcutoff = 7, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                 pointSize = 4,
                                 colCustom = color.key,
                                 selectLab = labeled.genes,
                                 legendPosition = "right",
                                 labSize = 6.0,
                                 drawConnectors = T,
                                 colAlpha = 1, #You can't generate .eps files with transparent color values. This ensures no points are transparent.
                                 title = expression(paste("Day 70 hair cells and HC-like cells: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-5,6)) +
  NoLegend()

# home computer:
ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset SC/d70 SC DESeq2 VolcanoPlot V220530a.jpg",
       plot = custom.VP.jpg,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = jpeg())


ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e A11a EpCAM+-/Subset SC/d70 SC DESeq2 VolcanoPlot V220530a.eps",
       plot = custom.VP.eps,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = "eps")

