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
singlets.ori <- singlets
names(singlets) <- c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")
singlets.combined <- merge(singlets[["T0_1"]], 
					y = c(singlets[c("T0_2", 
						names(singlets)[grepl("AI", names(singlets))])]),
					add.cell.ids = c("T0_1", "T0_2", names(singlets)[grepl("AI", names(singlets))]), 
					project = "TRADITIOM_MCF7_AI")
saveRDS(singlets.combined, "/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.rds")
#########################################################
#### Normalize merged data with fastMNN and run UMAP ####
#########################################################
singlets.combined.ori <- singlets.combined
singlets.combined <- NormalizeData(singlets.combined)
singlets.combined <- FindVariableFeatures(singlets.combined, selection.method = "vst", nfeatures = 2000)
singlets.combined <- RunFastMNN(object.list = SplitObject(singlets.combined, split.by = "orig.ident"))
singlets.combined <- FindNeighbors(singlets.combined, reduction = "mnn", dims = 1:30)
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.4)
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:30)

singlets.combined4 <- FindVariableFeatures(singlets.combined4, selection.method = "vst", nfeatures = 2000)
singlets.combined4VF <- subset(singlets.combined4, features=VariableFeatures(singlets.combined4))
mat<-GetAssayData(singlets.combined, assay="mnn.reconstructed", slot = "data")
mat<-Matrix(mat, sparse = T )
Matrix::writeMM(obj = mat, file="~/manuscript/rds/singlets.combinedAIall.mtx")
matVF<-GetAssayData(singlets.combined4VF, assay="mnn.reconstructed", slot = "data")
matVF<-Matrix(matVF, sparse = T )
Matrix::writeMM(obj = matVF, file="~/manuscript/rds/singlets.combinedAIvf.mtx")

# singlets.combined.new <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.5)

saveRDS(singlets.combined, "/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.3.rds")
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.4)
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:30, seed.use=17)
saveRDS(singlets.combined, "/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.4.rds")

singlets.combined3 <- readRDS("/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.3.rds")
# singlets.combined <- readRDS("/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.4.rds")
singlets.combined <- readRDS("~/manuscript/rds/singlets.combinedAI.MNN0.4.rds")

###################
#### Plot UMAP ####
###################
# colors<-c("#7c1158", "#ffa300", "#e6d800", "#50e991", "#9b19f5", "#e60049", "#1F9E89FF", "#dc0ab4", "#b3d4ff", "blue", "#0bb4ff", "grey", "brown")
colors<-c("#ffa300", "#B1961A", "#e6d800", "#50e991", "#9b19f5", "#dc0ab4", "#1F9E89FF", "#e60049", "#b3d4ff", "blue", "#0bb4ff", "grey", "brown")
names(colors)<-0:12
# colors3<-c("#e6d800", "#7c1158", "#50e991", "#9b19f5", "#1F9E89FF", "#dc0ab4", "#b3d4ff", "blue", "#0bb4ff", "grey")
colors3<-c("#e6d800", "#ffa300", "#50e991", "#9b19f5", "#1F9E89FF", "#e60049", "#b3d4ff", "blue", "#0bb4ff", "grey")
names(colors3)<-0:9
filesOrder<-c("AI1m1", "AI2m1", "AI1m2", "AI2m2","AIaw1", "AIaw2", "AIaw3", "AIaw4", "T0_1", "T0_2")

singlets.combined$orig.ident<-factor(singlets.combined$orig.ident, levels=filesOrder)	
singlets.combined3$orig.ident<-factor(singlets.combined3$orig.ident, levels=filesOrder)	
pdf("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.4.pdf", height=7, width=7)
png("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.4.png", height=1000, width=1000, res=150)
DimPlot(singlets.combined, group.by="seurat_clusters", repel=T, label=T, reduction = "umap", label.size=6) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "none", 
				axis.title = element_text(size=16),
				axis.text = element_text(size=14)) +
	guides(fill=guide_legend(nrow=4,byrow=TRUE))
dev.off()
pdf("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.4split.pdf", height=10, width=10)
png("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.4.split.png", height=1800, width=2100, res=150)
DimPlot(singlets.combined, group.by="seurat_clusters", split.by="orig.ident", repel=T, label=F, reduction = "umap", label.size=8, ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "none",, 
				axis.title = element_text(size=20),
				axis.text = element_text(size=18)) 	+
	guides(fill=guide_legend(nrow=4,byrow=TRUE))
dev.off()
pdf("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.3.pdf", height=10, width=10)
png("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.3.png", height=600, width=600)
DimPlot(singlets.combined3, group.by="seurat_clusters", repel=T, label=T, reduction = "umap", label.size=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors3) +
	ggplot2::theme(legend.position = "bottom") +
	guides(fill=guide_legend(nrow=2,byrow=TRUE))
png("~/manuscript/plots/sca/final/singlets.combined.AI.mnn0.3.separate.png", height=1200, width=1200)
DimPlot(singlets.combined3, group.by="seurat_clusters", split.by="orig.ident", repel=T, label=F, reduction = "umap", label.size=4, ncol=4) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors3) +
	ggplot2::theme(legend.position = "bottom") +
	guides(fill=guide_legend(nrow=2,byrow=TRUE))	
dev.off()

#####################################
#### Assign barcodes to cell IDs ####
#####################################
singlets.combined4<-singlets.combined
singlets.combined<-singlets.combined3
singlets.combined<-singlets.combined4

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

bcPerSampleInfo.m$L1<-factor(bcPerSampleInfo.m$L1, levels=files)
bcPerSampleInfo.m$bc<-factor(bcPerSampleInfo.m$bc, levels=bcOrder)
bcPerSampleInfo.m$seurat_clusters<-factor(bcPerSampleInfo.m$seurat_clusters, levels=0:9)

pdf("~/manuscript/plots/sca/final/singlets.combinedAI.allBCpctGT30.pdf", height=10, width=12)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct[bcPct$value>=30,]$L1,], aes(x="", y=value, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors) +
  facet_grid(L1~bc) +
  theme(legend.position="bottom") +
	guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()

pdf("~/manuscript/plots/sca/final/singlets.combinedAI0.3.allBCpctGT30.pdf", height=10, width=12)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct[bcPct$value>=30,]$L1,], aes(x="", y=value, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors3) +
  facet_grid(L1~bc) +
  theme(legend.position="bottom") +
	guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()



######################
#### Find markers ####
######################
allMarkers <- FindAllMarkers(singlets.combined, assay="RNA", logfc.threshold = 0.50, test.use = "roc", only.pos = TRUE)
allMarkers.l <- split(allMarkers, allMarkers$cluster)
enrichedList <- lapply(allMarkers.l[c(1:7, 9:12)], function(x){
	res <- enrichr(rownames(x), "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020
	res$count <- str_split(res$Overlap, "/", simplify=TRUE)[,1]
	res <- res[as.numeric(res$count) >= 3 & res$Adjusted.P.value<=0.05, ]
	res
})

samplePhase <- paste(singlets.combined.ccScore[[]]$orig.ident, singlets.combined.ccScore[[]]$Phase, sep="_")

singlets.combined.ccScore<-AddMetaData(
						  object = singlets.combined.ccScore,
						  metadata = samplePhase,
						  col.name = 'samplePhase')
singlets.combined.ccScore.sub <- subset(singlets.combined.ccScore, Phase == "G2M")


g1DormancyVsAwakening<-FindMarkers(singlets.combined, ident.1 = 0, ident.2 = c(1,2), p_val_adj=0.001, only.pos = TRUE, logfc.threshold = 0.5, min.diff.pct = 0.25)
enrichr(rownames(g1DormancyVsAwakening), "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020
enrichr(rownames(g1DormancyVsAwakening), "TRANSFAC_and_JASPAR_PWMs")[[1]]

singlets.combined.ori<-singlets.combined

Idents(object = singlets.combined.ccScore.sub) <- "samplePhase"
g2mDormancyVsAwakening<-FindMarkers(singlets.combined.ccScore.sub, ident.1 = "AI1m1_G2M", ident.2 = "AIaw4_G2M", p_val_adj=0.001, only.pos = TRUE, logfc.threshold = 0.5, min.diff.pct = 0.25)


pdf("g2mDormancyVsAwakening.pdf", height=50, width=20)
RidgePlot(singlets.combined.ccScore, features=rownames(g2mDormancyVsAwakening), ncol=3) 
dev.off()


enrichr(rownames(g2mDormancyVsAwakening), "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020


#####################################
#### Assign barcodes to cell IDs ####
#####################################
cellecta.10X<-readRDS("rds/cellecta.10XAll.rds")
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

bcPerSampleInfo.m$L1<-factor(bcPerSampleInfo.m$L1, levels=files)
bcPerSampleInfo.m$bc<-factor(bcPerSampleInfo.m$bc, levels=bcOrder)
bcPerSampleInfo.m$seurat_clusters<-factor(bcPerSampleInfo.m$seurat_clusters, levels=0:9)

pdf("plotsAllAImnn/manuscript/singlets.combined.AI.allBCpctGT20.pdf", height=10, width=12)
ggplot(bcPerSampleInfo.m[bcPerSampleInfo.m$bc %in% bcPct[bcPct$value>=20,]$L1,], aes(x="", y=value, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors) +
  facet_grid(L1~bc) +
  theme(legend.position="bottom")
dev.off()

ccInfo2<-singlets.combined.sub[[]][,c("orig.ident","bc", "seurat_clusters")]
# ccInfo2<-ccInfo2[ccInfo2$bc!="complex",]
bcPerCluster<-split(ccInfo2, ccInfo2$seurat_clusters)
bcPerClusterInfo.m<-melt(lapply(bcPerCluster, function(x){
	total<-sum(table(x[,2]))
	(table(x[,2])/total)*100
	}))

bcPerClusterInfo.m<-bcPerClusterInfo.m[order(bcPerClusterInfo.m$value, decreasing=T),]
library(randomcoloR)
n <- length(as.vector(unique(bcPerClusterInfo.m[bcPerClusterInfo.m$value>1,]$Var1)))

colors<- rep("grey", length(unique(bcPerClusterInfo.m$Var1)))
names(colors)<-unique(bcPerClusterInfo.m$Var1)
colors[c("bc-2127641", "bc-150411", "bc-7704123", "bc-5122096", "bc-6864523", "complex")]<-c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "blue", "black")
colors[names(colors) %in% as.vector(unique(bcPerClusterInfo.m[bcPerClusterInfo.m$value>1,]$Var1)) & !names(colors) %in% c("bc-2127641", "bc-150411", "bc-7704123", "bc-5122096", "bc-6864523","complex")]<-distinctColorPalette(15)

bcOrder<-as.vector(unique(bcPerClusterInfo.m[order(bcPerClusterInfo.m$value, decreasing=T),]$Var1))
bcPerClusterInfo.m$Var1<-factor(bcPerClusterInfo.m$Var1, levels=bcOrder)
# bcPerClusterInfo.m$L1<-factor(bcPerClusterInfo.m$L1, levels=c(1,7,3, 0, 5,6,8,9,4,2))
pdf("plotsAllAImnn/singlets.combined.AI.BCclusterPct.pdf", height=10, width=23)
ggplot(bcPerClusterInfo.m, aes(x="", y=value, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.1) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors) +
  labs(fill="Barcodes") +
  facet_grid(.~L1) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))
dev.off()


