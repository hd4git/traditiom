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

setwd("/rds/general/user/hdhiman/home/scRNAseq/T47D/")
singlets<-readRDS("rds/singlets.cellectaAll.rds")
cellecta.10X.t47d<-readRDS("~/scRNAseq/T47D/rds/cellecta.10XAll.rds")

singlets.combined <- merge(singlets[[1]], y = c(singlets[2:4]),
					add.cell.ids = names(singlets)[1:4], project = "TRADITIOM_T47D")

singlets.combined <- NormalizeData(singlets.combined)
singlets.combined <- FindVariableFeatures(singlets.combined, selection.method = "vst", nfeatures = 2000)
singlets.combined <- RunFastMNN(object.list = SplitObject(singlets.combined, split.by = "orig.ident"))
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:30, seed.use="88")
singlets.combined <- FindNeighbors(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.2)
# singlets.combined <- RunTSNE(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined$orig.ident <- factor(singlets.combined$orig.ident, levels=c("T0_1", "T0_2", "1m1", "2m1", "aw1"))
saveRDS(singlets.combined, "rds/singlets.combined.T47D.info.rds")

colors2 <- c("0"="#E6D800", "1"="#4DAF4A", "2"="#50e991", "3"="#FFA300", "4"="#377EB8", "5"="#B1961A", "6"="#E60049", "7"="#0bb4ff", "8"="#B3D4FF")

pdf("~/manuscriptFeb23/plots/figure4f.t47d.dimplot.pdf", height=3, width=3)
DimPlot(singlets.combined, reduction = "umap") + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors2) +
	ggplot2::theme(legend.position = "none")
dev.off()	
pdf("~/manuscriptFeb23/plots/figure4f.t47d.dimplotSplit.pdf", height=3, width=12)
DimPlot(singlets.combined, split.by="orig.ident", reduction = "umap", ncol=5) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors2)  
dev.off()

####################################################################
############ Dimplot and pie charts with cell cycle phases ########
####################################################################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined.t47d <- CellCycleScoring(singlets.combined.t47d, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
phaseColors<-c("G1"="green4", "G2M"="orange", "S"="red")

pdf("~/manuscriptFeb23/plots/t47d.phase.pdf", height=3, width=12)
DimPlot(singlets.combined.t47d, group.by="Phase", split.by="orig.ident", reduction = "umap", ncol=5) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=phaseColors)  
dev.off()

phaseCounts <- table(singlets.combined.t47d[[]][,c("orig.ident", "Phase")])
pctPhaseCounts <- 100*(phaseCounts/rowSums(phaseCounts))
pctPhaseCounts.m <- melt(pctPhaseCounts[1:4,])

pdf("~/manuscriptFeb23/plots/pctPhaseCountsT47D.pdf", height=3, width=12)
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

saveRDS(singlets.combined, "~/../ephemeral/scRNAseq/cellranger/singlets.combined.t47d.tsne30dim.rds")


singlets.combined$orig.ident <- factor(singlets.combined$orig.ident, levels=names(singlets))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined <- CellCycleScoring(singlets.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

saveRDS(singlets.combined, "~/scRNAseq/T47D/rds/singlets.combined.T47D.rds")
