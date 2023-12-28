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

### Percentage barcodes in each cluster across samples ###

singlets.combined <- readRDS("~/manuscriptFeb23/rds/singlets.combined.TAM.rds")
cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
samples <- as.vector(unique(singlets.combined$orig.ident))
samples <- samples[grepl("^T", samples)]
cellecta.10X$sampleCell <- paste(cellecta.10X$SampleName, cellecta.10X$CellName, sep="_")
cellecta.10X <- cellecta.10X[cellecta.10X$SampleName %in% samples,]
singlets.combined <- AddMetaData(
    object = singlets.combined,
    metadata = Cells(singlets.combined),
    col.name = 'sampleCell')

singlets.combined.sub <- subset(singlets.combined, sampleCell %in% cellecta.10X$sampleCell )
singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = cellecta.10X[match(singlets.combined.sub[[]]$sampleCell, cellecta.10X$sampleCell),]$LineageBarcode,
    col.name = 'bc')
singlets.combined.sub <- subset(singlets.combined.sub, bc != "complex")

ccInfo2<-singlets.combined.sub[[]][,c("orig.ident","bc", "seurat_clusters")]
bcPerSample<-split(ccInfo2, ccInfo2$orig.ident)
bcPerSampleInfo.m<-melt(lapply(bcPerSample, function(x){
	total<-sum(table(x[,2:3]))
	(table(x[,2:3])/total)*100
	}))
bcPct<-split(bcPerSampleInfo.m, bcPerSampleInfo.m$bc)
bcPct.m<-melt(lapply(bcPct, function(x){sum(x$value)}))
bcOrder<-bcPct.m[order(bcPct.m$value, decreasing=T),]$L1

bcPerSampleInfo.m$L1<-factor(bcPerSampleInfo.m$L1, levels=samples)
bcPerSampleInfo.m$bc<-factor(bcPerSampleInfo.m$bc, levels=bcOrder)
bcPerSampleInfo.m$seurat_clusters<-factor(bcPerSampleInfo.m$seurat_clusters, levels=0:11)

bcPerSampleInfo.m <- bcPerSampleInfo.m[bcPerSampleInfo.m$value>0,]
bcColors <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#FFFAC8", "#9A6324")

pdf("~/manuscriptFeb23/plots/singlets.combined.TAM.BCpct.pdf", height=10, width=6)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct.m[bcPct.m$value>10,]$L1,], aes(x=seurat_clusters, y=value, fill=bc)) +
	geom_bar(stat="identity") +
	theme_bw() +
	scale_fill_manual(values = bcColors) +
	facet_grid(L1~., space="free") +
	labs(x="Clusters", y="Percentage of barcode in each cluster", fill="Barcodes") +
	theme(legend.position="bottom",
		strip.background=element_rect(fill="white"))
dev.off()



