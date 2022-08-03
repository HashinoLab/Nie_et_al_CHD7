
#d20 WT (#A2e) and CHD7 KO/KO (#A11a), FACS sorted PAX2-2a-nGFP+ cells


### Load the libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(EnhancedVolcano)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load in Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d20_A2e <- Read10X("D:/BioInfo/Seurat Analysis/d20 A2e A11a/mapping/d20 A2e+pA") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
d20_A11a <- Read10X("D:/BioInfo/Seurat Analysis/d20 A2e A11a/mapping/d20 A11a+pA") %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)

#Add metadata
d20_A2e$Genotype <- "A2e"
d20_A11a$Genotype <- "A11a"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Filtering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d20_A2e$mt.per <- PercentageFeatureSet(d20_A2e, pattern = "^MT-")
VlnPlot(d20_A2e, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


d20_A11a$mt.per <- PercentageFeatureSet(d20_A11a, pattern = "^MT-")
VlnPlot(d20_A11a, pt.size = 0.5, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)


d20_A2e <- subset(d20_A2e, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 50000 & mt.per < 20)

d20_A11a <- subset(d20_A11a, subset = nFeature_RNA > 1300 & nFeature_RNA < 8000 & nCount_RNA < 50000 & mt.per < 20)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%   Merge Datasets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_d20 <- list(d20_A2e, d20_A11a)

d20 <- merge(list_d20[[1]], list_d20[[2]])
rm(d20_A2e, d20_A11a, list_d20)

saveRDS(d20, file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge V210729.rds")

d20 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge V210729.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%  Seurat workflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d20 <- SCTransform(d20, vars.to.regress = "mt.per")
d20 <- RunPCA(d20)
DimPlot(d20)

nPC <- 30
res <- 0.06

d20 <- RunUMAP(d20, dims = 1:nPC)

d20 <- FindNeighbors(d20, dims = 1:nPC)
d20 <- FindClusters(d20, resolution = res)
DimPlot(d20)

DimPlot(d20, group.by ="Genotype")

FeaturePlot(d20, features = "nFeature_RNA")

saveRDS(d20, file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge UMAP V210729.rds")

d20 <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/RStudio files/d20 A2e A11a +pA Merge UMAP V210729.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%  Feature Plots, etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FeaturePlot(d20, features = c("FBXO2", "EPCAM", "CDH1", "CDH2", "PAX2", "PAX8", "SOX2", "SOX10", "JAG1", "AP1M2", "GATA3", "COL9A2", "ESPN", "OC90", "PLEKHB1"), pt.size = 0.2, 
            ncol = 4)

DotPlot(d20, features = c("FBXO2", "EPCAM", "CDH1", "CDH2", "PAX2", "PAX8", "SOX2", "SOX10", "JAG1", "AP1M2", "GATA3", "COL9A2", "ESPN", "OC90", "PLEKHB1")) + RotatedAxis()


FeaturePlot(d20, features = c("DLX5", "TBX1", "TBX2", "SIX1", "HMX2", "FBXO2", "SOX2", "COL11A1"), pt.size = 0.2, ncol = 4)


VlnPlot(d20, features = c("DLX5", "TBX1", "TBX2", "SIX1", "HMX2", "FBXO2", "SOX2", "COL11A1"), pt.size = 0.2, ncol = 4)


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


#%%%%%%%%%%%%%%%%%%%% RColorBrewer %%%%%%%%%%%%%%%%%%%%%%

library(RColorBrewer)

FeaturePlot(d20_OP, features = c("DLX5"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("SIX1"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("SOX2"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("SOX10"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("HOXB9"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("FBXO2"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))
FeaturePlot(d20_OP, features = c("COL9A2"), pt.size = 2, col = c("lightgray", brewer.pal(n=9, name="RdPu")))


#%%%%%%%%%%%%%%%%%%%%%%%% DE and GSEA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#CDH7 Bi-allelic knockout Day20 Prosensory DE and GSEA
library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(zinbwave, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(edgeR, quietly = TRUE)
library(iDEA, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)

register(MulticoreParam(8))
print("Cores = Registered")
d20.ps <- readRDS(file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA Merge OP V210729.rds")
print("Seurat object read in")

d20.sce <- SummarizedExperiment(as.matrix(d20.ps@assays$RNA@counts),
                                colData = d20.ps@meta.data)
rm(d20.ps)
print("SCE created, Seurat removed")

filter <- rowSums(assay(d20.sce)) > 5
d20.sce <- d20.sce[filter,]
assay(d20.sce) %>% log1p %>% rowVars -> vars #This line is cursed
names(vars) <- rownames(d20.sce)
vars <- sort(vars, decreasing = T)

d20.sce <- d20.sce[names(vars[1:15000]),]
assayNames(d20.sce)[1] <- "counts"


d20.colData <- data.frame(Genotype = colData(d20.sce)$Genotype)
iCounts <- assay(d20.sce, i = "counts")
design <- model.matrix(~d20.colData$Genotype)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = d20.colData, design = ~ Genotype)
weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                 design = design,
                                 maxit = 500, normalization = "DESeq2_poscounts",
                                 colData = d20.colData,
                                 designFormula = ~ Genotype, 
                                 verbose = TRUE)

saveRDS(weights, file="/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA Merge OP weights V210807.rds")
print("weights calculation completed, check SFTP for file")

assays(dse)[["weights"]] <- weights

dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = TRUE)

res <- results(dse, cooksCutoff = F, parallel = T)
res$log2FoldChange %>% summary()
res$lfcSE %>% summary()

#Save DESeq results
full.results <- results(dse, parallel = T)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA OP DESeq2 V210807.csv")



# %%%%%%%%%%%%%%%%%%%  Volcano Plots from DESeq2 .csv file   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")


library(EnhancedVolcano)
library(data.table)
library(Seurat)


labeled.genes <- c("DLX5", "TBX2", "SIX1", "FBXO2", "SOX2", "HOXB9", "HOXA10", "HOXB7", "PCOLCE", "EMX2", 
                   "DLX6", "OC90", "PAX3", "FGF8")

logfc.threshold <- 1
p.threshold <- 1E-10

KO <- "#D6604D"
WT <- "#4393C3"
Highlight <- "black"


d20_OP_DE <- fread(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a - Alex/DE/Prosensory_DESeq2.csv")
d20_OP_DE$log2FoldChange <- d20_OP_DE$log2FoldChange * -1 #Specific to the Prosensory Day 20 dataset, becuase of a previous error in defining the mutant as the reference dataset.
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
                                 title = expression(paste("Day 20 Prosensory: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-5, 5)) +
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
                                 title = expression(paste("Day 20 Prosensory: ",
                                                          italic("CHD7"), ""^"-/-", 
                                                          " vs. ",
                                                          italic("CHD7"), ""^"+/+"))) + 
  xlim(c(-5,5)) +
  NoLegend()


ggsave(filename = "D:/BioInfo/Seurat Analysis/d20 A2e A11a/Merged A2e + A11a/D20 OP VolcanoPlot v210501.jpg",
       plot = custom.VP.jpg,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = jpeg())


ggsave(filename = "D:/BioInfo/Seurat Analysis/d20 A2e A11a/Merged A2e + A11a/D20 OP VolcanoPlot v210501.eps",
       plot = custom.VP.eps,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       device = "eps")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%% HPC_GSEA with iDEA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#CDH7 Bi-allelic knockout Day20 Prosensory GSEA
library(tidyverse, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(iDEA, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(data.table)

register(MulticoreParam(8))
print("Cores = Registered")

res <- fread(file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA OP DESeq2 V210807.csv")
CompleteSets <- readRDS(file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/CompleteSets.rds")

de.summary <- data.frame(log2FoldChange = res$log2FoldChange,
                         lfcSE2 = res$logFC_SE^2,
                         row.names = res$gene)

de.summary.up <- de.summary[de.summary$log2FoldChange > 0, ]

d20.idea <- CreateiDEAObject(de.summary.up, CompleteSets, num_core = 8)
d20.idea <- iDEA.fit(d20.idea, modelVariant = F)
d20.idea <- iDEA.louis(d20.idea)
write.csv(d20.idea@gsea, file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA OP Up iDEA V210811b.csv")

rm(d20.idea, de.summary.up)


de.summary.down <- de.summary[de.summary$log2FoldChange < 0, ]

d20.idea <- CreateiDEAObject(de.summary.down, CompleteSets, num_core = 8)
d20.idea <- iDEA.fit(d20.idea, modelVariant = F)
d20.idea <- iDEA.louis(d20.idea)
write.csv(d20.idea@gsea, file = "/N/u/jingnie/Carbonate/R/d20 A2e A11a +pA OP DESeq2/d20 A2e A11a +pA OP Down iDEA V210811b.csv")


#------- bubble plot -------------

library(data.table)
library(iDEA)
library(tidyverse)

customGeneSetsInfo <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/customGeneSetsInfo.rds")



plotdata<- fread("D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP Up iDEA V210811b.csv")
bp.data <- data.table::merge.data.table(plotdata, customGeneSetsInfo, by.x = "annot_id", by.y = "gset") 

bp.data$Log10_Pvalue_Louis <- -1*log10(bp.data$pvalue_louis)

# Save the bp.data dataframe here
# Renamed "Up" GSEA results into "DOWN" GSEA + Gene set info file. 
write.csv(bp.data, "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP DOWN iDEA+GeneSetsInfo.csv")


library(data.table)
library(iDEA)
library(tidyverse)

customGeneSetsInfo <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/customGeneSetsInfo.rds")



plotdata<- fread("D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP Down iDEA V210811b.csv")
bp.data <- data.table::merge.data.table(plotdata, customGeneSetsInfo, by.x = "annot_id", by.y = "gset") 


bp.data$Log10_Pvalue_Louis <- -1*log10(bp.data$pvalue_louis)

# Save the bp.data dataframe here
# Renamed "Down" GSEA results into "UP" GSEA + Gene set info file. 
write.csv(bp.data, "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP UP iDEA+GeneSetsInfo.csv")



library(data.table)
library(iDEA)
library(tidyverse)
library(stringr)
library(RColorBrewer)


customGeneSetsInfo <- readRDS(file = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/customGeneSetsInfo.rds")

### The initial DESeq2 analysis seems to switched the up- and down- regulated genes, which results in switched GSEA results.
# The following script converts the "up" GSEA results into "down" files while adding gene set info, and vice versa.

plotdata<- fread("D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/csv files/d20 A2e A11a +pA OP Up iDEA V210811b.csv")
bp.data <- data.table::merge.data.table(plotdata, customGeneSetsInfo, by.x = "annot_id", by.y = "gset") 

bp.data$Category <- droplevels(bp.data$gsetBioName)

#### Bubble plot 

includedCats <- c("GO BIOLOGICAL PROCESS",
                  "GO MOLECULAR FUNCTION",
                  "GO CELLULAR COMPONENT",
                  "REACTOME",
                  "KEGG",
                  "PID",
                  "POSITIONAL",
                  "TRANSCRIPTION FACTORS")

caseCats <- c("GO biological process",
              "GO molecular function",
              "GO cellular component",
              "Reactome",
              "KEGG",
              "PID",
              "Positional",
              "Transcription factors")

bp.data <- bp.data[toupper(bp.data$Category) %in% includedCats, ] 
bp.data$Category <- droplevels(bp.data$Category) %>% factor(caseCats) 
bp.data <- bp.data[order(match(toupper(bp.data$Category), includedCats)), ]
bp.data$IDNum <- row.names(bp.data) %>% as.integer() 
bp.data$Log10_Pvalue_Louis <- -1*log10(bp.data$pvalue_louis)


#Optional Step to save the bp.data dataframe here
# Renamed "Up" GSEA results into "DOWN" GSEA for bubble plot file. 
write.csv(bp.data, "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/csv files/d20 A2e A11a +pA OP DOWN iDEA for bubble plot.csv")



##-------------------------------------------------------------
## Configurable Variables
##-------------------------------------------------------------


# Colors
display.brewer.all()

display.brewer.pal(12, "Set3")

brewer.pal(12, "Set3")


bp.colors <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5")

labeled.genesets <- c("GO_EAR_DEVELOPMENT",
                      "GO_INNER_EAR_MORPHOGENESIS",
                      "GO_SENSORY_ORGAN_DEVELOPMENT",
                      "GO_REGULATION_OF_WNT_SIGNALING_PATHWAY",
                      "GO_MORPHOGENESIS_OF_AN_EPITHELIUM",
                      "GO_EXTRACELLULAR_MATRIX",
                      "GO_GROWTH_FACTOR_BINDING",
                      "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                      "NABA_CORE_MATRISOME",
                      "PID_FGF_PATHWAY",
                      "WNT_SIGNALING",
                      "REACTOME_CELL_JUNCTION_ORGANIZATION",
                      "SIX1_TARGET_GENES",
                      "SOX10_TARGET_GENES",
                      "chr19q13")

included.categories <- c("GO biological process",
                         "GO molecular function",
                         "GO cellular component",
                         "Reactome",
                         "KEGG",
                         "PID",
                         "Positional",
                         "Transcription factors")

data.to.plot <- bp.data

plot.title <- expression(paste("Downregulated genes in ", italic("CHD7"), ""^"KO/KO", " otic progenitors at day 20"))


Sig <- data.to.plot[which(data.to.plot$annot_id %in% labeled.genesets),]


Sig$Term <- tolower(Sig$annot_id)

library(ggplot2)
library("ggrepel")

bp <- ggplot(data.to.plot, aes(x = IDNum, y = Log10_Pvalue_Louis, color = Category)) +
  geom_point(shape = 19, alpha=1, size = 4) + 
  labs(x = "",
       y = expression(paste(bold(-log[10]), bold("("), bolditalic(p), bold("-value)")))) +
  ggtitle(label = plot.title) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 30, face="bold"),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'grey80'),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 40, face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.title=element_text(size=16,face = 'bold'),
        legend.text=element_text(size=16),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_hline(yintercept = 1.82, col = 'black', linetype = 2, size=2) +
  scale_color_manual(values = bp.colors) +
  scale_fill_manual(values = bp.colors) +
  theme(legend.direction = "vertical") +
  theme(legend.position = c(0.8, 0.95)) +
  theme(legend.box = "horizontal") +
  theme(legend.title.align = 0) +
  geom_label_repel(
    data = Sig,
    aes(label = Term),
    col = 'black',
    size = 6,
    nudge_y = 0.6)

png("D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP DOWN BubblePlot v210812a.png", width = 12500, height = 7500, res = 600)
bp
dev.off()


ggsave(filename = "D:/BioInfo/Seurat Analysis/d20 A2e A11a + bGHpA/GSEA/d20 A2e A11a +pA OP DOWN BubblePlot v210812a.eps",
       plot = bp,
       width = 20,
       height = 20/1.66667,
       units = "in",
       dpi = 300,
       device = "eps")


