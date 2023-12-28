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

singlets.combined <- readRDS("/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.4.rds")
singlets.combined.ori <- singlets.combined
singlets.combined <- singlets.combined.ori
DefaultAssay(singlets.combined) <- "mnn.reconstructed"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined.ccScore <- CellCycleScoring(singlets.combined,
        s.features = s.genes,
        g2m.features = g2m.genes,
        set.ident = TRUE)

singlets.combined.cc <- ScaleData(singlets.combined.ccScore, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlets.combined.ccScore))

saveRDS(singlets.combined.cc, "/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAImnn.ccAll.rds")
singlets.combined.cc <- RunPCA(singlets.combined.cc, 
	assay =  "mnn.reconstructed", 
	features = VariableFeatures(singlets.combined.cc, 
		assay =  "mnn.reconstructed", 
		slot = "data"), 
	nfeatures.print = 10)

# cell cycle effects strongly mitigated in PCA
DefaultAssay(singlets.combined.cc) <- "mnn.reconstructed"
singlets.combined.cc <- JackStraw(singlets.combined.cc, num.replicate = 100)
singlets.combined.cc <- ScoreJackStraw(singlets.combined.cc, dims = 1:20)

JackStrawPlot(singlets.combined.cc, dims = 1:15)
ElbowPlot(singlets.combined.cc, ndims=50)

singlets.combined.cc <- FindNeighbors(singlets.combined.cc, dims = 1:10)
singlets.combined.cc.5 <- FindClusters(singlets.combined.cc, resolution = 0.45)
singlets.combined.cc <- RunUMAP(singlets.combined.cc, dims = 1:10, seed.use=17)
singlets.combined.cc$orig.ident<-factor(singlets.combined.cc$orig.ident, levels=c("AI1m1", "AI2m1", "AI1m2", "AI2m2","AIaw1", "AIaw2", "AIaw3", "AIaw4", "T0_1", "T0_2"))	
singlets.combined.cc.4 <- singlets.combined.cc

saveRDS(singlets.combined.cc.5, "/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combined.MNNccRegressedCluster0.45Info.rds")

singlets.combined.cc<-readRDS("~/manuscript/rds/singlets.combined.MNNccRegressedCluster0.45Info.rds")
colors<-c("#FFA300", "#E6D800", "#B1961A", "#1F9E89", "#9B19F5", "#DC0AB4", "#B3D4FF")
names(colors)<-0:6
colors2 <- c("#B1961A", "#FFA300","#1F9E89", "#E6D800","#9B19F5", "#DC0AB4", "#E60049", "#B3D4FF", "grey")
names(colors2)<-0:8
singlets.combined<-singlets.combined.cc

pdf("~/manuscript/plots/sca/MNNccRegressed/singlets.combined.AI.ccRegressedAIafterMNN0.5.pdf", height=12, width=12)
png("~/manuscript/plots/sca/MNNccRegressed/singlets.combined.AI.ccRegressedAIafterMNN0.5.phase2.png", height=500, width=1000)
DimPlot(singlets.combined.cc, group.by="seurat_clusters", repel=T, label=T, reduction = "umap", label.size=8) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors2) +
	ggplot2::theme(legend.position = "bottom") 
# DimPlot(singlets.combined.cc.5, group.by="Phase", repel=T, label=F, reduction = "umap", label.size=4) + 
# 	ggplot2::labs(title = "") +
# 	# scale_color_manual(values=colors[1:10]) +
# 	ggplot2::theme(legend.position = "bottom") 
DimPlot(singlets.combined.cc, group.by="Phase",split.by="Phase", repel=T, label=F, reduction = "umap", label.size=4, ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=phaseColors) +
	ggplot2::theme(legend.position = "bottom") 

DimPlot(singlets.combined.cc, group.by="Phase",split.by="orig.ident", repel=T, label=F, reduction = "umap", label.size=4, ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=phaseColors) +
	ggplot2::theme(legend.position = "bottom") 
DimPlot(singlets.combined.cc, group.by="seurat_clusters", split.by="orig.ident", repel=T, label=F, reduction = "umap", label.size=4, ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors2) +
	ggplot2::theme(legend.position = "bottom") 
dev.off()


singlets.combined <- singlets.combined.cc.5
### Assign cell IDs ###
cellecta.10X<-readRDS("~/manuscript_Oct/rds/cellecta.10XAll.rds")
cellecta.10X$sampleCell<-paste(cellecta.10X$SampleName, cellecta.10X$CellName, sep="_")
files <- as.vector(unique(singlets.combined$orig.ident))
cellecta.10X <- cellecta.10X[cellecta.10X$SampleName %in% files,]
singlets.combined <- AddMetaData(
	  object = singlets.combined,
	  metadata = Cells(singlets.combined),
	  col.name = 'sampleCell')

singlets.combined.sub <- subset(singlets.combined, sampleCell %in% cellecta.10X$sampleCell )
singlets.combined.sub <- AddMetaData(
	  object = singlets.combined.sub,
	  metadata = cellecta.10X[match(singlets.combined.sub[[]]$sampleCell, cellecta.10X$sampleCell),]$LineageBarcode,
	  col.name = 'bc')

selBCs <- c("bc-7704123", "bc-2127641", "bc-6864523", "bc-5122096", "complex")
reducedBCs <- singlets.combined.sub[[]]$bc
reducedBCs[!reducedBCs %in% selBCs] <- "others"

singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = reducedBCs,
    col.name = "reducedBCs")


##############################################################
#### subset barcodes for winner Vs others analysis #########
##############################################################
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

pdf("~/manuscript/plots/sca/MNNccRegressed/subBCaw0.45.pdf", height=4, width=4)
lapply(subBCList, function(x){
  DimPlot(x, group.by="seurat_clusters", repel=T, label=F, reduction = "umap", label.size=4) + 
  ggplot2::labs(title = "") +
  scale_color_manual(values=colors2) +
  ggplot2::theme(legend.position = "none") 
})
dev.off()

pdf("~/manuscript/plots/sca/MNNccRegressed/subBCothers0.45.pdf", height=4, width=4)
lapply(subBCothersList, function(x){
  DimPlot(x, group.by="seurat_clusters", repel=T, label=F, reduction = "umap", label.size=4) + 
  ggplot2::labs(title = "") +
  scale_color_manual(values=colors2) +
  ggplot2::theme(legend.position = "none") 
})
dev.off()


getPct <- function(sampleName, bc){
	x <- subset(singlets.combined.sub, orig.ident == sampleName)
	reducedBCs <- x[[]]$bc
	reducedBCs[!reducedBCs %in% c(bc, "complex")] <- "others"
	x <- AddMetaData(
	    object = x,
	    metadata = reducedBCs,
	    col.name = "reducedBCs")
	x <- subset(x, reducedBCs %in% c(bc, "others"))
	x.winners <- subset(x, reducedBCs %in% bc)
	x.others <- subset(x, reducedBCs == "others")

	countsWinner <- table(x.winners$seurat_clusters)
	countsOthers <- table(x.others$seurat_clusters)
	counts <- table(x$seurat_clusters)
	data.frame("winner"=100*(countsWinner/sum(countsWinner)),
				"others"=100*(countsOthers/sum(countsOthers)))[,c(1,2,4)]
}
selWinners<-c("bc-7704123","bc-6864523","bc-2127641", "bc-5122096")
pctInfoAIaw1 <- getPct("AIaw1", "bc-7704123")
pctInfoAIaw2 <- getPct("AIaw2", "bc-6864523")
pctInfoAIaw3 <- getPct("AIaw3", "bc-2127641")
pctInfoAIaw4 <- getPct("AIaw4", "bc-5122096")
pctInfoList <- list(
					"AIaw1"=pctInfoAIaw1, "AIaw2"=pctInfoAIaw2, "AIaw3"=pctInfoAIaw3, "AIaw4"=pctInfoAIaw4)
pctInfo.m <- melt(pctInfoList)


pdf("~/manuscript/plots/sca/MNNccRegressed/pctInfoWinnerVsOthers.pdf", height=7, width=7)
ggplot(pctInfo.m, aes(x="", y=value, fill=winner.Var1)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.001) +
  coord_polar("y", start=0, clip="off") +  
  theme_void() +
  scale_fill_manual(values=colors2) +
  labs(fill="Clusters") +
  facet_grid(L1~variable) +
  theme(legend.position="bottom",
    strip.text = element_text(size = 16), 
    legend.text=element_text(size= 16)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()


getPctAll <- function(sampleName){
	x <- subset(singlets.combined.sub, orig.ident == sampleName)
	winnerBC <- c("bc-7704123", "bc-6864523", "bc-2127641", "bc-5122096")
	reducedBCs <- x[[]]$bc
	reducedBCs[!reducedBCs %in% c(winnerBC, "complex")] <- "others"
	x <- AddMetaData(
	    object = x,
	    metadata = reducedBCs,
	    col.name = "reducedBCs")
	x <- subset(x, reducedBCs %in% c(winnerBC, "others"))

	bcInfoList <- split(x[[]], x[[]]$reducedBCs)
	bcInfoPct <- lapply(bcInfoList, function(splitBC){
		counts <- table(splitBC$seurat_clusters)
		pct <- 100*(counts/sum(counts))
		pct
		})
	bcInfoPct.m <- melt(bcInfoPct)
	colnames(bcInfoPct.m) <- c("cluster", "pct", "bc")
	bcInfoPct.m
}

pctAllBc <- lapply(files[1:6], getPctAll)
names(pctAllBc) <- files[1:6]
pctAllBc.m <- melt(pctAllBc, id.vars=c("bc", "cluster"), value="pct")
pctAllBc.m$cluster <- as.factor(pctAllBc.m$cluster)
pctAllBc.m <- pctAllBc.m[pctAllBc.m$value>0,]
pctAllBc.m$L1 <- factor(pctAllBc.m$L1, levels=files)
pctAllBc.m$bc <- factor(pctAllBc.m$bc, levels=c("bc-5122096", "bc-7704123", "bc-2127641", "bc-6864523", "others"))

pdf("~/manuscript/plots/sca/MNNccRegressed/pctInfoWinnerVsOthersAllBCs.pdf", height=12, width=15)
ggplot(pctAllBc.m, aes(x="", y=value, fill=cluster)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.00001) +
  coord_polar("y", start=0, clip="off") +  
  theme_void() +
  scale_fill_manual(values=colors2) +
  labs(fill="Clusters") +
  facet_grid(bc~L1) +
  theme(legend.position="bottom",
    strip.text = element_text(size = 16), 
    legend.text=element_text(size= 16)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
dev.off()

