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

#### Plot cell counts and feature counts ###

setwd("/rds/general/user/hdhiman/home/scRNAseq/cellRanger/")

singlets <- readRDS("rds/singlets.cellectaAll.rds")
names(singlets)<-c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")

scheme<-c("T0"="#BFBFBF", 
          "Dormancy"="#FFDC00", 
          "Awakening"="#B30078")
getState<-function(x){
		x$state <- ""
		if(any(grepl("T0", x$L1))){ x[grepl("T0", x$L1),]$state <- "T0"}
		if(any(grepl("m1|m2", x$L1))){ x[grepl("m1|m2", x$L1),]$state <- "Dormancy"}
		if(any(grepl("aw", x$L1))){ x[grepl("aw", x$L1),]$state <- "Awakening"}
		x$state<-factor(x$state, levels=c("T0","Dormancy", "Awakening"))
		x
}

cellCounts <- melt(lapply(singlets, function(x){dim(x)[2]}))
cellCounts$L1 <- factor(cellCounts$L1, levels= names(singlets))
cellCounts$samples <- "MCF7"
cellCounts <- getState(cellCounts)
pdf("~/manuscript/plots/filteredCellCount.pdf", height=7, width=4)
ggplot(cellCounts, aes(x=samples, y=value, label=L1)) +
	geom_violin() +
	# geom_jitter(aes(color=state)) +
	geom_text_repel(aes(color=state)) +
	scale_color_manual(values=scheme) +
	labs(x="", y="Number of cells", color="State") +
	theme_bw() +
	theme(legend.position="bottom")
dev.off() 

cellFeatures <- melt(lapply(singlets, function(x){
	mat <- as.data.frame(GetAssayData(x, assay="RNA", slot = "data"))
	featureCounts <- colSums(mat>0)
	featureCounts
	}))
cellFeatures <- getState(cellFeatures)
cellFeatures$L1 <- factor(cellFeatures$L1, levels=names(singlets)) 
png("~/manuscript/plots/filteredCellFeatures.png", height=300, width=1200)
ggplot(cellFeatures, aes(x=L1, y=log10(value), fill=state)) +
	geom_violin() +
	geom_jitter(size=0.1, width=0.1) +
	# geom_text_repel(aes(color=state)) +
	scale_fill_manual(values=scheme) +
	labs(x="", y="log10(Number of features)", fill="State") +
	theme_bw() +
	theme(legend.position="bottom", 
		axis.text = element_text(size=12),
		legend.text = element_text(size=10))
dev.off() 


# tLow.bcs <- lapply(singlets[1:2], function(x){
# 	mat <- as.data.frame(GetAssayData(x, assay="Custom", slot = "data"))
# 	bcs <- rownames(mat[rowSums(mat>0)>0,])
# 	gsub("-", "_", bcs)
# 	})
# tHigh.bcs <- lapply(1:3, function(x){
# 	unique(rownames(freq[freq[,x]>0,]))		### TRADITIOM High frequency matrix
# 	})
# names(tHigh.bcs) <- colnames(freq)[1:3]


# vennBCs <- plotVenn(append(tLow.bcs, tHigh.bcs), 
# 	sNames=names(append(tLow.bcs, tHigh.bcs)),
# 	nCycles = 15000, 
# 	outFile="~/manuscript/plots/overlapBCs.svg",
# 	fontScale = 1.5)

# tlow <- unique(melt(tLow.bcs)$value)
# tHigh <- unique(melt(tHigh.bcs)$value)
# tlow[!tlow %in% tHigh]
# # 23 cells with bc_988608 with reads in range of 1-12

### Get cellecta barcode counts for each sample ###
bcCounts <- melt(lapply(singlets, function(x){
	mat <- as.data.frame(GetAssayData(x, assay="Custom", slot = "data"))
	dim(mat[rowSums(mat>0)>0,])[1]
	}))
bcCounts$L1 <- factor(bcCounts$L1, levels= names(singlets))
bcCounts$samples <- "MCF7"
bcCounts <- getState(bcCounts)
pdf("~/manuscript/plots/filteredBcCount.pdf", height=7, width=4)
ggplot(bcCounts, aes(x=samples, y=value, label=L1)) +
	geom_violin() +
	# geom_jitter(aes(color=state)) +
	geom_text_repel(aes(color=state)) +
	scale_color_manual(values=scheme) +
	labs(x="", y="Number of BCs", color="State") +
	theme_bw() +
	theme(legend.position="bottom")
dev.off() 

### Merge MCF& AI samples and remove batch effects with FastMNN normalization ###

singlets.combined <- merge(singlets[[1]], y = c(singlets[2:18]),
					add.cell.ids = names(singlets), project = "TRADITIOM_MCF7")
singlets.combined <- NormalizeData(singlets.combined)
singlets.combined <- FindVariableFeatures(singlets.combined, selection.method = "vst", nfeatures = 2000)
singlets.combined <- RunFastMNN(object.list = SplitObject(singlets.combined, split.by = "orig.ident"))
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindNeighbors(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.4)
colors<-c("#7c1158", "#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#1F9E89FF", "blue", "grey")
names(colors)<-0:11

singlets.combined.info<-readRDS("~/manuscript/rds/singlets.combined.all.info.rds")

# pdf("~/manuscript/plots/sca/singlets.combined.All.pdf", height=10, width=10)
# DimPlot(singlets.combined, group.by="seurat_clusters",  reduction = "umap", repel=T, label=T, label.size=6) + 
# 	ggplot2::labs(title = "") +
# 	scale_color_manual(values=colors) +
# 	ggplot2::theme(legend.position = "bottom") 
# dev.off()

# singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.4)
# colors<-c("#7c1158", "#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#1F9E89FF", "blue", "grey")
names(colors)<-0:11
samples<-c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")
singlets.combined$orig.ident<-factor(singlets.combined$orig.ident, levels=c(samples[3:18], samples[1:2]))
pdf("~/manuscript/plots/sca/singlets.combined.All.pdf", height=10, width=10)
DimPlot(singlets.combined, group.by="seurat_clusters", split.by="orig.ident",  reduction = "umap", ncol=4, raster="FALSE") + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "bottom") 
dev.off()

png("~/manuscript/plots/sca/singlets.combined.All.png", height=1000, width=720)
DimPlot(singlets.combined, group.by="seurat_clusters", split.by="orig.ident",  reduction = "umap", ncol=4, raster="FALSE") + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=c(colors)) +
	ggplot2::theme(legend.position = "bottom") 
dev.off()

saveRDS(singlets.combined, "~/../ephemeral/scRNAseq/cellranger/singlets.combined.all.info.rds")

### Assign cellecta barcodes to merged samples object ###
cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
cellecta.10X$sampleCell<-paste(cellecta.10X$sample, cellecta.10X$L1, sep="_")
cellecta.10X <- cellecta.10X[cellecta.10X$sample %in% samples,]
singlets.combined <- AddMetaData(
    object = singlets.combined,
    metadata = Cells(singlets.combined),
    col.name = 'sampleCell')

singlets.combined.sub <- subset(singlets.combined, sampleCell %in% cellecta.10X$sampleCell )
singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = cellecta.10X[match(singlets.combined.sub[[]]$sampleCell, cellecta.10X$sampleCell),]$value,
    col.name = 'bc')

### Pie charts for barcode percentage in each cluster for all samples ###
ccInfo2<-singlets.combined.sub[[]][,c("orig.ident","bc", "seurat_clusters")]
# ccInfo2<-ccInfo2[ccInfo2$bc!="complex",]
bcPerSample<-split(ccInfo2, ccInfo2$orig.ident)
bcPerSampleInfo.m<-melt(lapply(bcPerSample, function(x){
	total<-sum(table(x[,2:3]))
	(table(x[,2:3])/total)*100
	}))
bcPct<-split(bcPerSampleInfo.m, bcPerSampleInfo.m$bc)
bcPct<-melt(lapply(bcPct, function(x){sum(x$value)}))
bcOrder<-bcPct[order(bcPct$value, decreasing=T),]$L1

bcPerSampleInfo.m$L1<-factor(bcPerSampleInfo.m$L1, levels=samples)
bcPerSampleInfo.m$bc<-factor(bcPerSampleInfo.m$bc, levels=bcOrder)
bcPerSampleInfo.m$seurat_clusters<-factor(bcPerSampleInfo.m$seurat_clusters, levels=0:11)

pdf("~/manuscript/plots/sca/singlets.combined.all.BCpct.pdf", height=10, width=23)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct[bcPct$value>10,]$L1,], aes(x="", y=value, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors) +
  facet_grid(L1~bc) +
  theme(legend.position="bottom")
dev.off()
