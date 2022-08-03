
### d70 WT (#D5e) and CHD7 KO/+ (#E5b) POU4F3-2a-ntdTomato+ and POU4F3-2a-ntdTomato- cells


### Load the libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load in Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D5e_pos <- Read10X("C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/mapping/D5e nT+") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
E5b_pos <- Read10X("C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/mapping/E5b nT+") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
D5e_neg <- Read10X("C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/mapping/D5e nT-") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
E5b_neg <- Read10X("C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/mapping/E5b nT-") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)



D5e_pos$Genotype <- "D5e"
E5b_pos$Genotype <- "E5b"
D5e_neg$Genotype <- "D5e"
E5b_neg$Genotype <- "E5b"

D5e_pos$FACS <- "nT+"
E5b_pos$FACS <- "nT+"
D5e_neg$FACS <- "nT-"
E5b_neg$FACS <- "nT-"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Filtering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D5e_pos$mt.per <- PercentageFeatureSet(D5e_pos, pattern = "^MT-")
VlnPlot(D5e_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


E5b_pos$mt.per <- PercentageFeatureSet(E5b_pos, pattern = "^MT-")
VlnPlot(E5b_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


D5e_neg$mt.per <- PercentageFeatureSet(D5e_neg, pattern = "^MT-")
VlnPlot(D5e_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


E5b_neg$mt.per <- PercentageFeatureSet(E5b_neg, pattern = "^MT-")
VlnPlot(E5b_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)



D5e_pos <- subset(D5e_pos, subset = nFeature_RNA > 1500 & nFeature_RNA < 9500 & nCount_RNA < 70000 & mt.per < 30)

E5b_pos <- subset(E5b_pos, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 40000 & mt.per < 30)

D5e_neg <- subset(D5e_neg, subset = nFeature_RNA > 1500 & nFeature_RNA < 8000 & nCount_RNA < 50000 & mt.per < 12)

E5b_neg <- subset(E5b_neg, subset = nFeature_RNA > 1200 & nFeature_RNA < 8000 & nCount_RNA < 45000 & mt.per < 10)


VlnPlot(D5e_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(E5b_pos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(D5e_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
VlnPlot(E5b_neg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%   Merge Datasets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list_d70pos <- list(D5e_pos, E5b_pos)

d70pos <- merge(list_d70pos[[1]], list_d70pos[[2]], add.cell.ids = c("D5eP", "E5bP"))

head(colnames(d70pos))



list_d70neg <- list(D5e_neg, E5b_neg)

d70neg <- merge(list_d70neg[[1]], list_d70neg[[2]], add.cell.ids = c("D5eN", "E5bN"))

head(colnames(d70neg))



list_d70 <- list(d70pos, d70neg)

d70 <- merge(list_d70[[1]], list_d70[[2]])


rm(D5e_pos, E5b_pos, D5e_neg, E5b_neg, d70pos, d70neg, list_d70, list_d70pos, list_d70neg)

saveRDS(d70, file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge V210819.rds")


d70 <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge V210819.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%  Seurat workflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(Seurat)

# Running on IU Carbonate interactive session:
d70 <- SCTransform(d70, vars.to.regress = "mt.per")
d70 <- RunPCA(d70)


saveRDS(d70, file = "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/d70 D5eE5b nT+nT- Merge SCT V210819.rds")


# download the SC Transformed .rds file and load in on local PC:

d70 <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge SCT V210819.rds")



DimPlot(d70)


nPC <- 30
res <- 0.2

d70 <- RunUMAP(d70, dims = 1:nPC)

d70 <- FindNeighbors(d70, dims = 1:nPC)
d70 <- FindClusters(d70, resolution = res)
DimPlot(d70, pt.size = 0.8)

DimPlot(d70, group.by ="Genotype", pt.size = 0.8)
DimPlot(d70, group.by ="FACS", pt.size = 0.8)


rm(nPC, res)

FeaturePlot(d70, features = "nFeature_RNA", pt.size = 1)


saveRDS(d70, file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge UMAP V210819.rds")




#%%%%%%%%%%%%%%%%%%%%%%%%%  Feature Plots, etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Hair cell genes
FeaturePlot(d70, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "FBXO2", "ESPN", "GFI1", "SOX2", 
                                 "CD164L2", "CCER2", "GNG8", "EPCAM", "CDH1", "CDH2"), pt.size = 0.2, ncol = 4)


# Supporting cell and otic lineage genes
FeaturePlot(d70, features = c(
  "FBXO2","SOX2", "EPCAM", "CDH1", "CDH2", "S100A1", "OC90", "CLDN6", "SOX10", "JAG1", "SIX1", "SPARCL1", "CNTN1", 
  "LMX1A", "KRT23","SPP1", "GATA3", "COL9A2", "ARHGEF19", "PLEKHB1", "TBX2", "MARVELD3", "BDNF", "SIX1", "FGF10"
                              ), pt.size = 0.2, ncol = 5)


# Dot plot

DotPlot(d70, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "POU4F3", "ESPN", "TMPRSS3", "BDNF",
                             "FBXO2", "SOX2","OC90", "S100A1", "COL9A2", "TBX2", "SPARCL1", "EPCAM", "CDH1", "CDH2", "CLDN6", "ATP1B1", "JAG1", "SIX1", "SOX10", 
                             "NEFL", "LY6H", "SCG2", "FABP7", "PCSK1N", "MAP1B", "GPC3", "CDH6", "GAP43", "CTXN1", "CNIH2",
                             "TUBB3", "STMN2", "STMN4", "DCX", "ELAVL3", "ELAVL4", "ELAVL2", "INA", "CHGB", "ATP1A3", "NEUROD1",
                             "MGP", "COL3A1", "COL1A2", "POSTN", "LGALS1", "LUM",
                             "COL1A1", "ACTG2", "PRRX1", "COL6A3", "TAGLN", "FN1",
                             "S100B", "CRYM", "CST3", "SERPINE2", "APOE", "SERPINA3",
                             "TNNT1", "ACTC1", "DES", "TPM2", "MYOG", "MYL4", "MYOD1",
                             "PLP1", "CRYAB", "MPZ", "EDNRB", "CRABP1", "S100A6",
                             "CENPF", "TOP2A", "UBE2C", "NUSAP1", "HMGB2", "BIRC5", 
                             "KRT17", "KRT19", "S100A8", "KRT5", "KRT13", "TACSTD2", 
                             "MYLPF", "TNNI1", "MYL1", "TNNC2", "TTN", "ACTA1", "MYH3"
)) + RotatedAxis()


# find all markers of cluster xxxxx
cluster0.markers <- FindMarkers(d70, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 100)

write.csv(cluster0.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster0.csv")


cluster1.markers <- FindMarkers(d70, ident.1 = 1, min.pct = 0.25)
write.csv(cluster1.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster1.csv")


cluster2.markers <- FindMarkers(d70, ident.1 = 2, min.pct = 0.25)
write.csv(cluster2.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster2.csv")


cluster3.markers <- FindMarkers(d70, ident.1 = 3, min.pct = 0.25)
write.csv(cluster3.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster3.csv")


cluster4.markers <- FindMarkers(d70, ident.1 = 4, min.pct = 0.25)
write.csv(cluster4.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster4.csv")


cluster5.markers <- FindMarkers(d70, ident.1 = 5, min.pct = 0.25)
write.csv(cluster5.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster5.csv")


cluster6.markers <- FindMarkers(d70, ident.1 = 6, min.pct = 0.25)
write.csv(cluster6.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster6.csv")


cluster7.markers <- FindMarkers(d70, ident.1 = 7, min.pct = 0.25)
write.csv(cluster7.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster7.csv")


cluster8.markers <- FindMarkers(d70, ident.1 = 8, min.pct = 0.25)
write.csv(cluster8.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster8.csv")


cluster9.markers <- FindMarkers(d70, ident.1 = 9, min.pct = 0.25)
write.csv(cluster9.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster9.csv")


cluster10.markers <- FindMarkers(d70, ident.1 = 10, min.pct = 0.25)
write.csv(cluster10.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster10.csv")


cluster11.markers <- FindMarkers(d70, ident.1 = 11, min.pct = 0.25)
write.csv(cluster11.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster11.csv")


cluster12.markers <- FindMarkers(d70, ident.1 = 12, min.pct = 0.25)
write.csv(cluster12.markers, "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/cluster12.csv")




#%%%%%%%%%%%%%%%%%%%%%%%%%% Subsetting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#-------------- Subsetting HCs --------------------

d70_HC <- subset(d70, idents = c(1, 6, 8))
rm(d70)

saveRDS(d70_HC, file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset HC V210819.rds")

saveRDS(d70_HC, file = "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/d70 D5eE5b nT+nT- Merge subset HC V210819.rds")


d70_HC <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset HC V210819.rds")



d70_HC@active.ident <- d70_HC$Genotype %>% as.factor()


DimPlot(d70_HC)

d70_HC@active.ident <- d70_HC$seurat_clusters

DimPlot(d70_HC)



FeaturePlot(d70_HC, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "POU4F3", "ESPN", "TMPRSS3", "BDNF",
                             "FBXO2", "SOX2","OC90", "S100A1", "COL9A2", "TBX2", "SPARCL1", "EPCAM", "CDH1", "CDH2", 
                             "CLDN6", "ATP1B1", "JAG1", "SIX1", "SOX10"
), pt.size = 0.2, ncol = 6)



#-------------- Subsetting maturing HCs (mHCs) --------------------

d70_mHC <- subset(d70, idents = c(1, 6))
rm(d70)

saveRDS(d70_mHC, file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset mHC V210820.rds")


d70_mHC <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset mHC V210820.rds")



d70_mHC@active.ident <- d70_mHC$Genotype %>% as.factor()


DimPlot(d70_mHC)

d70_mHC@active.ident <- d70_mHC$seurat_clusters

DimPlot(d70_mHC)

FeaturePlot(d70_mHC, features = "nFeature_RNA", pt.size = 1)


FeaturePlot(d70_mHC, features = c("ATOH1", "MYO7A", "MYO6", "PCP4", "ANXA4", "GFI1", "POU4F3", "ESPN", "TMPRSS3", "BDNF",
                                 "FBXO2", "SOX2","OC90", "S100A1", "COL9A2", "TBX2", "SPARCL1", "EPCAM", "CDH1", "CDH2", 
                                 "CLDN6", "ATP1B1", "JAG1", "SIX1", "SOX10"
), pt.size = 0.2, ncol = 6)




#----------------- Subsetting SCs -------------------
d70_SC <- subset(d70, idents = 2)
rm(d70)

saveRDS(d70_SC, file = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset SC V210819.rds")

saveRDS(d70_SC, file = "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/d70 D5eE5b nT+nT- Merge subset SC V210819.rds")

d70_SC <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset SC V210819.rds")



d70_SC@active.ident <- d70_SC$Genotype %>% as.factor()


DimPlot(d70_SC)

d70_SC@active.ident <- d70_SC$seurat_clusters

DimPlot(d70_SC)



FeaturePlot(d70_SC, features = c("OC90", "WFDC2", "OTOL1", "KRT18", "KRT8", "FBXO2", "COL2A1",
                                 "HES1", "EPCAM", "CLDN6", "COL9A2", "SPP1", "SLITRK6", "SOX2", 
                                 "SOX10", "SPARCL1", "PLEKHB1", "TBX2"
), pt.size = 0.2, ncol = 5)




FeaturePlot(d70_SC, features = c(
  "HES6", "FGF20", "CHD7", "CALB1", "PCP4", "CIB2", "SOX2", "SIX1", "BDP1", "AIFM1", "RDX"
), pt.size = 0.2, ncol = 5)


FeaturePlot(d70_SC, features = c(
  "SLC2A1", "BNIP3", "THBS1", "OTOR", "COL11A2", "COL11A1", "OTOA", "OTOS", "OTOG", "CRYM"
  ), pt.size = 0.2, ncol = 5)



#-------------- Subsetting Neurons --------------------

d70_neurons <- subset(d70, idents = c(11))
DimPlot(d70_neurons)
rm(d70)

saveRDS(d70_neurons, file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset neurons V220509.rds")


d70_neurons <- readRDS(file = "C:/Research/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge subset neurons V220509.rds")


d70_neurons@active.ident <- d70_neurons$Genotype %>% as.factor()


DimPlot(d70_neurons)

d70_neurons@active.ident <- d70_neurons$seurat_clusters

DimPlot(d70_neurons)



FeaturePlot(d70_neurons, features = c("TUBB3", "NEUROD1", "NEUROD2", "NEUROG1", "STMN2", "STMN4", "NEFH", "NEFM", "NEFL",
                                 "ELAVL2", "ELAVL3", "ELAVL4"
                                 
), pt.size = 0.2, ncol = 4)





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

d70_neurons_sce <- SummarizedExperiment(as.matrix(d70_neurons@assays$RNA@counts),
                                   colData = d70_neurons@meta.data)



saveRDS(d70_neurons_sce, file = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge neurons sce V220510.rds")


d70_neurons_sce <- readRDS(file = "/N/u/jingnie/Carbonate/R/d70 D5e E5b HC/d70 D5eE5b nT+nT- Merge neurons sce V220510.rds")
print("Seurat object read in")

rm(d70_neurons) 
print("SCE created, Seurat removed")

filter <- rowSums(assay(d70_neurons_sce)) > 5
d70_neurons_sce <- d70_neurons_sce[filter,]
vars <- assay(d70_neurons_sce) %>% log1p() %>% rowVars()
names(vars) <- rownames(d70_neurons_sce)
vars <- sort(vars, decreasing = T)

d70_neurons_sce <- d70_neurons_sce[names(vars[1:15000]),]
assayNames(d70_neurons_sce)[1] <- "counts"


d70_neurons.colData <- data.frame(Genotype = colData(d70_neurons_sce)$Genotype)
iCounts <- assay(d70_neurons_sce, i = "counts")

d70_neurons.colData$Genotype <- as.factor(d70_neurons.colData$Genotype)

d70_neurons.colData$Genotype <- relevel(d70_neurons.colData$Genotype, ref = "D5e")

design <- model.matrix(~d70_neurons.colData$Genotype)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = d70_neurons.colData, design = ~ Genotype)


memory.limit()

rm(d20_OP_sce)

saveRDS(d20_OP.colData, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP d20_OP.colData V210731.rds")
saveRDS(design, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP design V210731.rds")
saveRDS(iCounts, file="D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP iCounts V210731.rds")


d20_OP.colData <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP d20_OP.colData V210731.rds")
design <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP design V210731.rds")
iCounts <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP iCounts V210731.rds")
library(DESeq2, quietly = TRUE)


weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                 design = design,
                                 maxit = 500, normalization = "DESeq2_poscounts",
                                 colData = d70_neurons.colData,
                                 designFormula = ~ Genotype, 
                                 verbose = TRUE)

saveRDS(weights, file="D:/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge neurons weights V220510.rds")


weights <- readRDS(file = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/RStudio files/d70 D5eE5b nT+nT- Merge neurons weights V220510.rds")

assays(dse)[["weights"]] <- weights

dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = FALSE) 

full.results <- results(dse, parallel = FALSE)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/Subset neurons/DESeq2 d70 D5eE5b nT+nT- neurons v220510.csv")



#   Volcano Plots from DESeq2 .csv file   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")



library(EnhancedVolcano)
library(data.table)
library(Seurat)

labeled.genes <- c(
  "NEFM", "KRT8", "KRT18", "SIX1", "OTOL1", "TRIP11", "BDP1", "NEFL", "SOX2", 
  "MSX1", "SLK", "CXCL14", "COL9A3", "DLX5", "COL9A2","SLC9A3R1", "MARVELD2", "AIFM1", "USH1C",
  "HES5", "HES1", "STRC", "OTOGL", "CDH23", "JAG2"
)

logfc.threshold <- 1
p.threshold <- 1E-10

KO <- "#D6604D"
WT <- "#4393C3"
Highlight <- "black"


d20_OP_DE <- fread(file = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/Subset HC/d70 D5eE5b nT+nT- Merge HC DESeq2 V210819.csv")


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
                                 title = expression(paste("Day 70 hair cells: ",
                                                          italic("CHD7"), ""^"+/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-3, 3)) +
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
                                 title = expression(paste("Day 70 hair cells: ",
                                                          italic("CHD7"), ""^"+/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-3,3)) +
  NoLegend()

ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/Subset HC/d70 HC DESeq2 VolcanoPlot V210821.jpg",
       plot = custom.VP.jpg,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = jpeg())


ggsave(filename = "D:/BioInfo/Seurat Analysis/d70 D5e E5b/Subset HC/d70 HC DESeq2 VolcanoPlot V210821.eps",
       plot = custom.VP.eps,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = "eps")




















