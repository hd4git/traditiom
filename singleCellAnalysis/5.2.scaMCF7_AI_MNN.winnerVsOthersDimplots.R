require("Seurat")
require("dplyr")
require("ggplot2")
require("scRNAseq")
require("scuttle")
require("celldex")
require("Matrix")
require("stringr")
require("slingshot")
require("reshape2")
require("SeuratWrappers")
require("enrichR")
require("DoubletFinder")
require("RColorBrewer")
require("Seurat")
require("SeuratWrappers")
require("dplyr")

### Winner Vs Others dimplots ###
# set.seed(17)
setwd("/rds/general/user/hdhiman/home/scRNAseq/cellRanger")
files<-c("T0_1", "T0_2", "AI1m1", "AI2m1", "AI1m2", "AI2m2","AIaw1", "AIaw2", "AIaw3", "AIaw4")
singlets.combined <- readRDS("/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/rdsAllAI/singlets.combined.MNN.rds")
# singlets.combined <- singlets.combined.ori
colors<-c("#7c1158", "#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#1F9E89FF", "blue", "grey")
names(colors)<-0:11


#############################
#### Cell cycling scores ####
#############################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined <- CellCycleScoring(singlets.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#####################################
#### Assign barcodes to cell IDs ####
#####################################
cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
colnames(cellecta.10X)<-c("value", "L1", "sample")
cellecta.10X$sampleCell<-paste(cellecta.10X$sample, cellecta.10X$L1, sep="_")
cellecta.10X <- cellecta.10X[cellecta.10X$sample %in% files,]
singlets.combined <- AddMetaData(
    object = singlets.combined,
    metadata = Cells(singlets.combined),
    col.name = 'sampleCell')

singlets.combined.sub <- subset(singlets.combined, sampleCell %in% cellecta.10X$sampleCell )
singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = cellecta.10X[match(singlets.combined.sub[[]]$sampleCell, cellecta.10X$sampleCell),]$value,
    col.name = 'bc')
singlets.combined.sub <- subset(singlets.combined.sub, bc != "complex")

selBCs <- c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096", "complex")
reducedBCs <- singlets.combined.sub[[]]$bc
reducedBCs[!reducedBCs %in% selBCs] <- "others"

singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = reducedBCs,
    col.name = "reducedBCs")

#####################################################
#### subset winning and other barcodes #########
#####################################################
subBCs <- function(x,y){
          subSample <- subset(singlets.combined.sub, orig.ident == x)
          subset(subSample, bc %in% y)
}
subBCothers <- function(x,y){
          subSample <- subset(singlets.combined.sub, orig.ident == x)
          sel <- unique(subSample[[]]$bc)
          subset(subSample, bc %in% sel[!sel %in% c(y, "complex")])
}

subBCList <- list(subBCs("AIaw1", "bc-7704123"),
              subBCs("AIaw2", "bc-6864523"),
              subBCs("AIaw3", "bc-2127641"),
              subBCs("AIaw4", "bc-5122096"),
              subBCs("AI1m1", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCs("AI2m1", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCs("AI1m2", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCs("AI2m2", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096"))
              )
subBCothersList <- list(subBCothers("AIaw1", "bc-7704123"),
              subBCothers("AIaw2", "bc-6864523"),
              subBCothers("AIaw3", "bc-2127641"),
              subBCothers("AIaw4", "bc-5122096"),
              subBCothers("AI1m1", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCothers("AI2m1", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCothers("AI1m2", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096")),
              subBCothers("AI2m2", c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096"))              
              )
# colors<-colors3

pdf("~/manuscript/plots/sca/MNNccRegressed/subBCaw0.4.pdf", height=4, width=4)
lapply(subBCList, function(x){
  DimPlot(x, group.by="seurat_clusters", repel=T, label=F, reduction = "umap", label.size=4) + 
  ggplot2::labs(title = "") +
  scale_color_manual(values=colors) +
  ggplot2::theme(legend.position = "none") 
})
dev.off()



pdf("~/manuscript/plots/sca/MNNccRegressed/subBCothers0.4.pdf", height=4, width=4)
lapply(subBCothersList, function(x){
  DimPlot(x, group.by="seurat_clusters", repel=T, label=F, reduction = "umap", label.size=4) + 
  ggplot2::labs(title = "") +
  scale_color_manual(values=colors) +
  ggplot2::theme(legend.position = "none") 
})
dev.off()

