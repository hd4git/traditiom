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

data<-lapply(files, function(x){
	fn<-paste(x, paste("filtered_feature_bc_matrix", x, sep="_"), sep="/")
	data10X<-Read10X(data.dir = fn)
	sample<-CreateSeuratObject(counts = data10X$'Gene Expression', project = x, min.cells = 3, min.features = 200)
	sample
	})
names(data)<-files

pctMT<-lapply(data, function(x){
	x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
	x
	})
setwd("/rds/general/user/hdhiman/home/scRNAseq/cellRanger/")

VlnPlot(singlets.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


subsetData<-function(x, nFeature){
	x<-subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature & percent.mt < 20)
}
maxFeature<-as.numeric(c("T0_1"=6500, "T0_2"=8500, 
	"TAM1m1"=7500, "TAM2m1"=9000, "AI1m1"=5500, "AI2m1"=6000, 
	"TAM1m2"=7500, "TAM2m2"=8000, "AI1m2"=8000, "AI2m2"=9500, 
	"TAMaw1"=9500, "TAMaw2"=5750, "TAMaw3"=7000, "TAMaw4"=7500, 
	"AIaw1"=9000, "AIaw2"=9500,"AIaw3"=6500, "AIaw4"=7000))


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
getDoubletInfo<-function(x, pKchosen){
	annotations <- x@meta.data$ClusteringResults
	homotypic.prop <- modelHomotypic(annotations)
	nExp_poi <- round(0.075*nrow(x@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = pKchosen, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
}
pKvalues<-c("T0_1"=0.21, "T0_2"=0.29, 
	"TAM1m1"=0.04, "TAM2m1"=0.21, "AI1m1"=0.005, "AI2m1"=0.22, 
	"TAM1m2"=0.2, "TAM2m2"=0.27, "AI1m2"=0.09, "AI2m2"=0.29, 
	"TAMaw1"=0.30, "TAMaw2"=0.25, "TAMaw3"=0.30, "TAMaw4"=0.15, 
	"AIaw1"=0.28, "AIaw2"=0.3, "AIaw3"=0.26, "AIaw4"=0.3)

dataDoubletInfo<-lapply(1:length(dataScaled), function(x){
	dataScaled[[x]]<-getDoubletInfo(dataScaled[[x]], pKvalues[x])
	dataScaled[[x]]
	})
singlets<-lapply(dataDoubletInfo, function(x){
	x <- AddMetaData(
	  object = x,
	  metadata = x@meta.data[,6],
	  col.name = 'DF.classifications'
	)
	x<-subset(x, subset = DF.classifications=="Singlet")
	x
})


names(singlets)<-c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")
saveRDS(singlets, "rds/singletsAll.rds")
singlets<-lapply(names(data), function(x){
	cellecta<-data[[x]]$'Custom'
	cells.10X<-Cells(singlets[[x]])	
	cellecta.selected<-cellecta[,cells.10X]
	singlets[[x]][['Custom']] = CreateAssayObject(counts = cellecta.selected)	
	singlets[[x]]
})
saveRDS(singlets, "rds/singlets.cellectaAll.rds")

