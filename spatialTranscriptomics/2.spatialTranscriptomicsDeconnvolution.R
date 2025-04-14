require("NanoStringNCTools")
require("GeomxTools")
require("GeoMxWorkflows")
require("ggplot2")
require("knitr")
require("dplyr")
require("ggforce")
require("scales")
require("stringr")
require("reshape2")
require("SpatialDecon")
require("GeomxTools")
require("celldex")

setwd("~/spatialTranscriptomics/")

target_patientData <- readRDS("rds/target_patientData.rds")
target_patientData <- readRDS("rds/target_patientData_all.rds")

bg <- SpatialDecon::derive_GeoMx_background(norm = assayDataElement(target_patientData, elt = "log_q"),
                             probepool = fData(target_patientData)$Module,
                             negnames = c("NegProbe-WTX"))

hpca.se <- HumanPrimaryCellAtlasData()
hpca.mat <- assay(hpca.se)
hpca.mat.m<- melt(hpca.mat)

hpca.se.info <- as.data.frame(colData(hpca.se))
hpca.mat.m$label <- hpca.se.info[match(hpca.mat.m$Var2, rownames(hpca.se.info)),]$label.main
hpca.exp <- reshape2::dcast(hpca.mat.m[,-2], Var1~label, mean)
hpca.exp.mat <- hpca.exp[,-1]
rownames(hpca.exp.mat) <- hpca.exp[,1]
hpca.exp.mat <- as.matrix(hpca.exp.mat)

res <- runspatialdecon(object = target_patientData,
                      norm_elt = "log_q",
                      raw_elt = "exprs",
                      X = hpca.exp.mat,
                      align_genes = TRUE)
# res <- runspatialdecon(object = target_patientData,
#                       norm_elt = "log_q",
#                       raw_elt = "exprs",
#                       X = safeTME,
#                       align_genes = TRUE)

cellProportions <- as.data.frame(melt(t(pData(res)$prop_of_all)))
cellProportions <- cellProportions[cellProportions$value>0,]
cellProportions$sampleName <- pData(res)[match(cellProportions$Var2, rownames(pData(res))), ]$sampleName
cellProportions$segment <- pData(res)[match(cellProportions$Var2, rownames(pData(res))), ]$segment
cellProportions$roi <- pData(res)[match(cellProportions$Var2, rownames(pData(res))), ]$roi

# colors <- sample(color, 21)
# names(colors) <- ct
# saveRDS(colors, "rds/colors.rds")
colors <- readRDS("rds/colors.rds")
# colors <- colors[1:18]
# names(colors) <- as.vector(unique(cellProportions$Var1))
newNames <- read.table("rds/rename.txt", header=F)
cellProportions$sampleNameNew <- newNames[match(cellProportions$sampleName, newNames$V1),]$V2
cellProportions$roiNew <- paste(cellProportions$sampleNameNew, cellProportions$roi, sep=".")
cellProportions$roiNew <- gsub("Tumor core ", "core", cellProportions$roiNew)
cellProportions$roiNew <- gsub("Tumor edge ", "edge", cellProportions$roiNew)
cellProportions$roiNew <- gsub("Lympho ", "lympho", cellProportions$roiNew)
cellProportions$roiNew <- gsub("00", "", cellProportions$roiNew)

cellProportions$sampleNameNew <- factor(cellProportions$sampleNameNew, levels=c("D1L", "D2", "D3", "S1L", "S1R", "S2", "S3", "R1R"))
cellProportions$segment <- factor(cellProportions$segment, levels=c("CK+", "CD45+", "Rest"))

getReducedCellTypes <-	function(x){
	x$label<-""
	x[x$Var1 %in% c("Epithelial_cells", "Neuroepithelial_cell"),]$label<-"Epithelial_cells"
	x[x$Var1 %in% c("Embryonic_stem_cells", "iPS_cells", "Tissue_stem_cells"),]$label<-"iPS|Embryonic|Tissue_stem_cells"
	x[x$Var1 == "Chondrocytes",]$label<-"Chondrocytes"
	x[x$Var1 == "Hepatocytes",]$label<-"Hepatocytes"
	x[x$Var1 == "Fibroblasts",]$label<-"Fibroblasts"
	x[x$Var1 == "Keratinocytes",]$label<-"Keratinocytes"
	x[x$Var1 == "Gametocytes",]$label<-"Gametocytes"
	x[x$Var1 == "Neurons",]$label<-"Neurons"
	x[x$Var1 == "CMP",]$label<-"CMP"
	x[x$Var1 %in% c("B_cell", "Macrophage", "Monocyte", "Neutrophils", "NK_cell", "T_cells", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", "DC", "Myelocyte"),]$label<-"Immune_cells"
	# splitInfo<-split(x, x$label)
	# as.matrix(vapply(splitInfo, function(y){
	# 		sum(y$Freq)
	# 	}, 1))
	x
}
cellProportions2 <- getReducedCellTypes(cellProportions)
cellProportions2List <- split(cellProportions2, paste(cellProportions2$roiNew, cellProportions2$segment, sep="_"))
cellProportions2reduced <- lapply(cellProportions2List, function(x){
	y <- split(x, x$label)
	lapply(y, function(info){
		sum(as.numeric(info$value))
	})
})
names(cellProportions2reduced) <- names(cellProportions2List)
cellProportions2reduced.m <- melt(cellProportions2reduced)
cellProportions2reduced.m$roiNew <- str_split(cellProportions2reduced.m$L1, "_", simplify=TRUE)[,1]
cellProportions2reduced.m$segment <- str_split(cellProportions2reduced.m$L1, "_", simplify=TRUE)[,2]
cellProportions2reduced.m$sampleNameNew <- str_split(cellProportions2reduced.m$roiNew, "\\.", simplify=TRUE)[,1]
cellProportions2reduced.m$sampleNameNew <- factor(cellProportions2reduced.m$sampleNameNew, levels=c("D1L", "D2", "D3", "S1L", "S1R", "S2", "S3", "R1R"))
cellProportions2reduced.m$segment <- factor(cellProportions2reduced.m$segment, levels=c("CK+", "CD45+", "Rest"))

pdf("plots/deconvolution_ms_all.pdf", height=7, width=18)
ggplot(cellProportions, aes(x=roiNew, y=value, fill=Var1)) +
	geom_bar(stat="identity") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90),
		legend.pos = "bottom", 
		strip.background = element_rect(fill="white")) +
	labs(x="", y="Proportion of cell type", fill="Cell types") +
	scale_fill_manual(values=colors) +
	# scale_fill_brewer(palette = "Set1") +
	facet_grid(.~segment+sampleNameNew, space="free", scales="free") 
	# facet_wrap(sampleName+roi~segment, ncol=5, nrow=17, scales="free")
	# facet_grid(sampleName+roi~segment, space="free", scales="free") 
dev.off()


pdf("plots/deconvolution_ms_all2.pdf", height=7, width=18)
ggplot(cellProportions2reduced.m, aes(x=roiNew, y=value, fill=L2)) +
	geom_bar(stat="identity") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90),
		legend.pos = "bottom", 
		strip.background = element_rect(fill="white")) +
	labs(x="", y="Proportion of cell type", fill="Cell types") +
	scale_fill_manual(values=c("Epithelial_cells" = "#4DAF4A", "Immune_cells" = "#FFFF33", "Chondrocytes" = "#E41A1C", "Hepatocytes" = "#377EB8", "Gametocytes" = "#984EA3", "Fibroblasts" = "#FF7F00", "Neurons" = "#A65628", "Keratinocytes" = "#F781BF", "CMP" = "#999999", "iPS|Embryonic|Tissue_stem_cells" = "grey")) +
	# scale_fill_brewer(palette = "Set1") +
	facet_grid(.~segment+sampleNameNew, space="free", scales="free") 
	# facet_wrap(sampleName+roi~segment, ncol=5, nrow=17, scales="free")
	# facet_grid(sampleName+roi~segment, space="free", scales="free") 
dev.off()



