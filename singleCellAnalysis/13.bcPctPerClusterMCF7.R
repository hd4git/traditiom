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

singlets.combined <- readRDS("~/manuscript/rds/singlets.combinedAI.MNN0.4.rds")
cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
samples<-c("T0_1", "T0_2", "AI1m1", "AI2m1","AI1m2", "AI2m2", "AIaw1","AIaw2", "AIaw3", "AIaw4")
colors<-c("#ffa300", "#B1961A", "#e6d800", "#50e991", "#9b19f5", "#dc0ab4", "#1F9E89FF", "#e60049", "#b3d4ff", "blue", "#0bb4ff", "grey", "brown")
names(colors)<-0:12

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
# ccInfo2<-ccInfo2[ccInfo2$bc!="complex",]
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
bcPerSampleInfo.m$seurat_clusters<-factor(bcPerSampleInfo.m$seurat_clusters, levels=0:12)

bcPerSampleInfo.m <- bcPerSampleInfo.m[bcPerSampleInfo.m$value>0,]

# multiColors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# bcColors <- sample(multiColors, 13)
bcColors <- c("slateblue1", "tomato2", "chartreuse4", "chocolate1", "tan4", "seagreen", "gainsboro", "paleturquoise3", "cornsilk2", "palegreen4", "seagreen3", "darkslateblue", "firebrick1") 

pdf("~/manuscriptFeb23/plots/bcPerClusterInfo_MCF7.pdf", height = 9, width = 7)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct.m[bcPct.m$value>10,]$L1,], aes(x=seurat_clusters, y=value, fill=bc)) +
	geom_bar(stat="identity") +
	theme_bw() +
	scale_fill_manual(values = bcColors) +
	facet_grid(L1~., space="free") +
	labs(x="Clusters", y="Percentage of barcode in each cluster", fill="Barcodes") +
	theme(legend.position="bottom",
		strip.background=element_rect(fill="white"))
dev.off()




