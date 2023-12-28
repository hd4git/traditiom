require("Seurat")
require("dplyr")
require("ggplot2")
require("scRNAseq")
require("scuttle")
require("celldex")
require("Matrix")
require("DoubletFinder")
require("stringr")
require("slingshot")
require("reshape2")
require("SeuratWrappers")

setwd("/rds/general/user/hdhiman/home/scRNAseq/cellRanger/")

singlets <- readRDS("rds/singlets.cellectaAll.rds")
names(singlets)<-c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")

samples <- names(singlets)[grepl("^T", names(singlets))]
singlets.TAM <- singlets[samples]

singlets.combined <- merge(singlets.TAM[[samples[1]]], y = singlets.TAM[samples[2:10]],
					add.cell.ids = names(singlets.TAM), project = "TRADITIOM_MCF7_TAM")
singlets.combined <- NormalizeData(singlets.combined)
singlets.combined <- FindVariableFeatures(singlets.combined, selection.method = "vst", nfeatures = 2000)
singlets.combined <- RunFastMNN(object.list = SplitObject(singlets.combined, split.by = "orig.ident"))
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindNeighbors(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.3)

colors<-c("0"="#E6D800", "1"="#50E991", "2"="#B1961A", "3"="#FFA300", "4"="#E60049", "5"="#DC0AB4", "6"="#4DAF4A", "7"="#B3D4FF")
png("~/manuscriptFeb23/plots/singlets.combined.TAM.png", height=500, width=500)
DimPlot(singlets.combined, group.by="seurat_clusters", reduction = "umap") + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "bottom",
				axis.title = element_text(size=26),
				axis.text = element_text(size=24),
				legend.text = element_text(size=22)) 
dev.off()

singlets.combined$orig.ident <- factor(singlets.combined$orig.ident, levels=c(names(singlets.TAM)[3:10], names(singlets.TAM)[1:2]))
saveRDS(singlets.combined, "~/manuscriptFeb23/rds/singlets.combined.TAM.rds")


png("~/manuscriptFeb23/plots/singlets.combined.TAMsplit.png", height=900, width=1200)
DimPlot(singlets.combined, group.by="seurat_clusters", split.by="orig.ident",  reduction = "umap", ncol=4, raster="FALSE") + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "bottom",
				axis.title = element_text(size=26),
				axis.text = element_text(size=24),
				legend.text = element_text(size=22)) 
dev.off()

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined.tam <- CellCycleScoring(singlets.combined.tam, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
phaseColors<-c("G1"="green4", "G2M"="orange", "S"="red")

pdf("~/manuscriptFeb23/plots/tam.phase.pdf", height=10, width=10)
DimPlot(singlets.combined.tam, group.by="Phase", split.by="orig.ident", reduction = "umap", ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=phaseColors)  
dev.off()

phaseCounts <- table(singlets.combined.tam[[]][,c("orig.ident", "Phase")])
pctPhaseCounts <- 100*(phaseCounts/rowSums(phaseCounts))
pctPhaseCounts.m <- melt(pctPhaseCounts)
pctPhaseCounts.m$orig.ident <- factor(pctPhaseCounts.m$orig.ident, levels=c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4"))

pdf("~/manuscriptFeb23/plots/pctPhaseCountsTAM.pdf", height=2, width=18)
ggplot(pctPhaseCounts.m, aes(x="", y=value, fill=Phase)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.001) +
  coord_polar("y", start=0, clip="off") +  
  theme_void() +
  scale_fill_manual(values=phaseColors) +
  labs(fill="Clusters") +
  facet_grid(.~orig.ident) +
  theme(legend.position="bottom",
    strip.text = element_text(size = 16), 
    legend.text=element_text(size= 16)) 
dev.off()


