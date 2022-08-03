

### d20 WT (#A2e), FACS sorted into PAX2-2a-nGFP+ and PAX2-2a-nGFP- populations


### Load the libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(EnhancedVolcano)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load in Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d20_Gpos <- Read10X("D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/mapping/d20 A2e P2G pos2") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
d20_Gneg <- Read10X("D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/mapping/d20 A2e P2G neg2") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)

d20_Gpos$Genotype <- "P2G+"
d20_Gneg$Genotype <- "P2G-"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Filtering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d20_Gpos$mt.per <- PercentageFeatureSet(d20_Gpos, pattern = "^MT-")
VlnPlot(d20_Gpos, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


d20_Gneg$mt.per <- PercentageFeatureSet(d20_Gneg, pattern = "^MT-")
VlnPlot(d20_Gneg, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)



d20_Gpos <- subset(d20_Gpos, subset = nFeature_RNA > 200 & nFeature_RNA < 9200 & nCount_RNA < 65000 & mt.per < 20)

d20_Gneg <- subset(d20_Gneg, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & nCount_RNA < 55000 & mt.per < 12)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%   Merge Datasets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_d20 <- list(d20_Gpos, d20_Gneg)



d20 <- merge(list_d20[[1]], list_d20[[2]])
rm(d20_Gpos, d20_Gneg, list_d20)

saveRDS(d20, file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge before UMAP v20220514a.rds")

d20 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge before UMAP v20220514a.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%  Seurat workflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d20 <- SCTransform(d20, vars.to.regress = "mt.per")
d20 <- RunPCA(d20)
DimPlot(d20)

nPC <- 30
res <- 0.17


d20 <- RunUMAP(d20, dims = 1:nPC)

d20 <- FindNeighbors(d20, dims = 1:nPC)
d20 <- FindClusters(d20, resolution = res)
DimPlot(d20)

DimPlot(d20, group.by ="Genotype")

FeaturePlot(d20, features = "nFeature_RNA", pt.size = 1)


saveRDS(d20, file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge UMAP res0.17 v20220517a.rds")

d20 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge UMAP res0.17 v20220517a.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%  Feature Plots, etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Otic progenitor genes

FeaturePlot(d20, features = c("FBXO2", "EPCAM", "CDH1", "CDH2", "PAX2", "bGHpA", "PAX8", "SOX2", "SOX10", "JAG1", "AP1M2", "GATA3", "COL9A2", "ESPN", "OC90", "PLEKHB1"), pt.size = 0.2, 
            ncol = 4)

DotPlot(d20, features = c("FBXO2", "EPCAM", "CDH1", "CDH2", "PAX2", "PAX8", "SOX2", "SOX10", "JAG1", "AP1M2", "GATA3", "COL9A2", "ESPN", "OC90", "PLEKHB1")) + RotatedAxis()



#%%%%%%%%%%%%%%%%%%%%%%%%%  Cluster annotation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# find all markers of cluster xxxxx
cluster0.markers <- FindMarkers(d20, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 100)

write.csv(cluster0.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster0.csv")


cluster1.markers <- FindMarkers(d20, ident.1 = 1, min.pct = 0.25)
write.csv(cluster1.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster1.csv")


cluster2.markers <- FindMarkers(d20, ident.1 = 2, min.pct = 0.25)
write.csv(cluster2.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster2.csv")


cluster3.markers <- FindMarkers(d20, ident.1 = 3, min.pct = 0.25)
write.csv(cluster3.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster3.csv")


cluster4.markers <- FindMarkers(d20, ident.1 = 4, min.pct = 0.25)
write.csv(cluster4.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster4.csv")


cluster5.markers <- FindMarkers(d20, ident.1 = 5, min.pct = 0.25)
write.csv(cluster5.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster5.csv")


cluster6.markers <- FindMarkers(d20, ident.1 = 6, min.pct = 0.25)
write.csv(cluster6.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster6.csv")


cluster7.markers <- FindMarkers(d20, ident.1 = 7, min.pct = 0.25)
write.csv(cluster7.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster7.csv")


cluster8.markers <- FindMarkers(d20, ident.1 = 8, min.pct = 0.25)
write.csv(cluster8.markers, "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/Analysis/cluster8.csv")




#%%%%%%%%%%%%%%%%%%%%%%%%%% Subsetting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



d20_OP <- subset(d20, idents = c(0,2))


saveRDS(d20_OP, file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP V210729.rds")

d20_OP <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge OP V210729.rds")



d20_OP@active.ident <- d20_OP$Genotype %>% as.factor()


DimPlot(d20_OP)

d20_OP@active.ident <- d20_OP$seurat_clusters

DimPlot(d20_OP)



FeaturePlot(d20_OP, features = c("DLX5", "TBX1", "TBX2", "SIX1", "HMX2", "FBXO2", "SOX2", "COL11A1"), pt.size = 0.2, ncol = 4)


VlnPlot(d20_OP, features = c("DLX5", "TBX2", "SIX1", "FBXO2", "SOX2"), 
        group.by = "Genotype", pt.size = 0.1, ncol = 5, y.max = 3.7, same.y.lims = TRUE)

VlnPlot(d20_OP, features = c("PAX8", "EMX2", "PCOLCE", "PAX3", "FGF17"), 
        group.by = "Genotype", pt.size = 0.1, ncol = 5, y.max = 3.7, same.y.lims = TRUE)

VlnPlot(d20_OP, features = c("HOXB9", "HOXA10", "HOXB7", "EMX2", "PCOLCE"), 
        group.by = "Genotype", pt.size = 0.1, ncol = 5, y.max = 3.3, same.y.lims = TRUE)


#------- subsetting neural crest cells----------------------------------


d20_NC <- subset(d20, idents = c(6))


saveRDS(d20_NC, file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge subsetNC v20220607a.rds")

d20_NC <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e P2G pos neg/RStudio files/d20 A2e P2G+- merge subsetNC v20220607a.rds")



d20_NC@active.ident <- d20_NC$Genotype %>% as.factor()


DimPlot(d20_NC)

d20_NC@active.ident <- d20_NC$seurat_clusters

DimPlot(d20_NC)


library(RColorBrewer)

FeaturePlot(d20_NC, features = c("SOX2"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("S100B"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("EPCAM"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("bGHpA"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("HOXB9"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("PAX3"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("SOX10"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("FOXD3"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("ERBB3"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("MPZ"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_NC, features = c("PLP1"), pt.size = 4, col = c("lightgray", brewer.pal(n=9, name="RdPu")))





# Split Vln Plots   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

split_d20_OP <- subset(d20_OP, idents = c(0,2))
split_d20_OP@active.ident <- factor(as.factor(split_d20_OP$Genotype), levels = c("A2e", "A11a"))
split_d20_OP@active.ident <- as.factor(rep("Day 20", length(Cells(split_d20_OP))))
split_d20_OP$Condition <- factor(as.factor(split_d20_OP$Genotype), levels = c("A2e", "A11a"))


VlnPlot(split_d20_OP, features = c("TBX2", "SIX1", "SOX2"), 
        split.by = "Genotype",
        split.plot = T,
        pt.size = 0.1, ncol = 3, y.max = 3.7, same.y.lims = TRUE)


VlnPlot(split_d20_OP, features = c("HOXB9", "EMX2", "PCOLCE"), 
        split.by = "Genotype",
        split.plot = T,
        pt.size = 0.1, ncol = 3, y.max = 3.3, same.y.lims = TRUE)











