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

#### Create seurat object for T47D and get stats ###
setwd("/rds/general/user/hdhiman/home/scRNAseq/T47D/data")
files<-c("T0_1", "T0_2", "1m1", "2m1", "aw1")
data<-lapply(files, function(x){
	fn<-paste("filtered_feature_bc_matrix_T47D", x, sep="_")
	data10X<-Read10X(data.dir = fn)
	sample<-CreateSeuratObject(counts = data10X$'Gene Expression', project = x, min.cells = 3, min.features = 200)
	sample
	})
names(data)<-files
data10X<-lapply(files, function(x){
	fn<-paste("filtered_feature_bc_matrix_T47D", x, sep="_")
	Read10X(data.dir = fn)
})
names(data10X)<-files

pctMT<-lapply(data, function(x){
	x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
	x
	})
setwd("/rds/general/user/hdhiman/home/scRNAseq/T47D/")

pdf("plots/violin_unfiltered_T47D.pdf")
lapply(names(pctMT), function(x){
        VlnPlot(pctMT[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) &
        ggplot2::labs(subtitle = x)
	})
dev.off()

### Filter cells to remove low quality cells and doublets ###
subsetData<-function(x, nFeature){
	x<-subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature & percent.mt < 20)
}
maxFeature<-as.numeric(c("T0_1"=7000, "T0_2"=7000, 
	"1m1"=5500, "2m1"=6000, 
	"aw1"=7000))


subsetDataInfo<-lapply(1:length(pctMT), function(x){
	pctMT[[x]]<-subsetData(pctMT[[x]], maxFeature[x])
	pctMT[[x]]
	})
names(subsetDataInfo)<-files
dataScaled<-lapply(subsetDataInfo, function(x){
	x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
	x <- ScaleData(x, features = rownames(x))
	x <- RunPCA(x, features = VariableFeatures(object = x))
	})

#############################################
#### Remove doublets with doublet finder ####
#############################################

pKinfo<-lapply(names(dataScaled), function(x){
	sweep.res.list <- paramSweep_v3(dataScaled[[x]], PCs = 1:10, sct = FALSE)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats)

	pK=as.numeric(as.character(bcmvn$pK))
	BCmetric=bcmvn$BCmetric
	pK_choose = pK[which(BCmetric %in% max(BCmetric))]
	dev.off()
	
	list("pK"=pK, "BCmetric"=BCmetric, "pK_choose"=pK_choose)
	})
names(pKinfo)<-files
saveRDS(pKinfo, "/rds/general/user/hdhiman/home/scRNAseq/T47D/rds/pKinfo_T47D.rds")


pdf("plots/pkInfo_t47D.pdf")
lapply(names(pKinfo), function(x){
	pK<-pKinfo[[x]][[1]]
	BCmetric<-pKinfo[[x]][[2]]
	pK_choose<-pKinfo[[x]][[3]]
	par(mar=c(5,4,4,8)+1, cex.main=1.2, font.main=2)
	plot(x = pK, y = BCmetric, pch = 16, type="b", col = "blue", lty=1)
	abline(v=pK_choose, lwd=2, col='red', lty=2)
	title(x)
	text(pK_choose, max(BCmetric), as.character(pK_choose), pos = 4, col = "red")	
})
dev.off()


pKvalues<-c("T0_1"=0.28, "T0_2"=0.23, 
	"1m1"=0.005, "2m1"=0.26, 
	"aw1"=0.14)
getDoubletInfo<-function(x, pKchosen){
	annotations <- x@meta.data$ClusteringResults
	homotypic.prop <- modelHomotypic(annotations)
	nExp_poi <- round(<corrected value>*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = pKchosen, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
}

getDoubletInfoAdj<-function(x, pKchosen){
	annotations <- x@meta.data$ClusteringResults
	homotypic.prop <- modelHomotypic(annotations)
	nExp_poi <- round(<corrected value>*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = pKchosen, nExp = nExp_poi, reuse.pANN = nExp_poi.adj, sct = FALSE)
}

dataDoubletInfo<-lapply(1:length(dataScaled), function(x){
	dataScaled[[x]]<-getDoubletInfo(dataScaled[[x]], pKvalues[x])
	dataScaled[[x]]
	})
### subset singlets ####
singlets<-lapply(dataDoubletInfo, function(x){
	x <- AddMetaData(
	  object = x,
	  metadata = x@meta.data[,6],
	  col.name = 'DF.classifications'
	)
	x<-subset(x, subset = DF.classifications=="Singlet")
	x
})


names(singlets)<-files
saveRDS(singlets, "rds/singlets.rds")

### Assing cellecta barcodes to singlets ####
singlets<-readRDS("rds/singlets.rds")
singlets.ori<-singlets
singlets<-lapply(names(data10X), function(x){
	cellecta<-data10X[[x]]$'Custom'
	cells.10X<-Cells(singlets[[x]])	
	cellecta.selected<-cellecta[,cells.10X]
	singlets[[x]][['Custom']] = CreateAssayObject(counts = cellecta.selected)	
	singlets[[x]]
})
names(singlets)<-files
saveRDS(singlets, "rds/singlets.cellectaAll.rds")
singlets<-readRDS("rds/singlets.cellectaAll.rds")
cellecta.10XAll<-lapply(files, function(x){
	cellectaBCsInfo<-singlets[[x]][['Custom']]
	cellectaBCsInfo.m<-as.matrix(GetAssayData(cellectaBCsInfo))
	melt(cellectaBCsInfo.m)
})

names(cellecta.10XAll)<-files
#### Annotate cells with unique barcodes ####
cellectaBCsDim2<-lapply(files, function(nm){
	cellectaBCsInfo<-singlets[[nm]][['Custom']]
	cellectaBCsInfo.df<-as.data.frame(GetAssayData(cellectaBCsInfo))
	gt0<-cellectaBCsInfo.df[rowSums(cellectaBCsInfo.df>0)>0, colSums(cellectaBCsInfo.df>0)>0]
	gt0.m<-melt(as.matrix(gt0))
	gt0.m<-gt0.m[gt0.m$value>0,]
	gt0.list<-split(gt0.m, gt0.m$Var2)
	filtered<-lapply(gt0.list, function(x){
			ifelse(dim(x)[1]>=1 & min(x$value)<max(x$value), 
				as.vector(x[x$value>=as.numeric(quantile(x$value, 0.9)),]$Var1), 
				ifelse(dim(x)[1]>=1 & min(x$value)==max(x$value) & quantile(x$value)[5]>3, 
					as.vector(x$Var1), 
					"complex"))			
		})
	melt(filtered)
})
names(cellectaBCsDim2)<-files

cellecta.10X<-do.call(rbind, lapply(files, function(x){
	cellectaBCsDim2[[x]]$sample<-x
	cellectaBCsDim2[[x]]
	}))
colnames(cellecta.10X)<-c("LineageBarcode", "CellName", "SampleName")
saveRDS(cellecta.10X, "~/scRNAseq/T47D/rds/cellecta.10XAll.rds")

####################################################################
####### Merge samples and remove batche effects with fastMNN  ######
####################################################################


singlets.combined <- merge(singlets[[1]], y = c(singlets[2:5]),
					add.cell.ids = names(singlets), project = "TRADITIOM_T47D")
singlets.combined <- NormalizeData(singlets.combined)
singlets.combined <- FindVariableFeatures(singlets.combined, selection.method = "vst", nfeatures = 2000)
singlets.combined <- RunFastMNN(object.list = SplitObject(singlets.combined, split.by = "orig.ident"))
singlets.combined <- RunUMAP(singlets.combined, reduction = "mnn", dims = 1:20)
singlets.combined <- FindNeighbors(singlets.combined, reduction = "mnn", dims = 1:20)
singlets.combined <- FindClusters(singlets.combined, reduction = "mnn", resolution = 0.1)

colors<-c("#7c1158", "#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#1F9E89FF", "blue", "grey", "black")
names(colors)<-0:12

singlets.combined$orig.ident <- factor(singlets.combined$orig.ident, levels=files)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined <- CellCycleScoring(singlets.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("plots/singlets.combined.0.1.pdf", height=10, width=10)
DimPlot(singlets.combined, group.by="seurat_clusters",  reduction = "umap", repel=T, label=T, label.size=6) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors[1:6]) +
	ggplot2::theme(legend.position = "bottom") 
DimPlot(singlets.combined, group.by="seurat_clusters", split.by="orig.ident", reduction = "umap", repel=T, label=F, label.size=6, ncol=2) + 
	ggplot2::labs(title = "") +
	scale_color_manual(values=colors[1:6]) +
	ggplot2::theme(legend.position = "bottom") 
DimPlot(singlets.combined, group.by="Phase", split.by="orig.ident", reduction = "umap", repel=T, label=F, label.size=6, ncol=2) + 
	ggplot2::labs(title = "") +
	# scale_color_manual(values=colors) +
	ggplot2::theme(legend.position = "bottom") 
dev.off()




