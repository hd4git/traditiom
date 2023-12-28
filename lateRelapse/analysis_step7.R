#module anaconda3/personal
#source activate jvCalls
 
require("VariantAnnotation")
require("stringr")
require("ggplot2")
require("reshape2")
require("circlize")
require("viridis")

setwd("/rds/general/user/hdhiman/home/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header")
filteredFiles <- list.files(pattern = "\\_filtered_filtered2.vcf.gz$")

### Get allele frequency and consequence for each variant in all samples ###
dataInfo <- lapply(filteredFiles, function(ip){
	test <- readVcf(ip, "hg38")
	testInfo <- geno(test)$AF
	colnames(testInfo) <- "AF"
	testGene <- info(test)$CSQ
	testGeneInfo <- lapply(seq_len(length(testGene)), function(x){
		# message(x)
		splitInfo <- as.data.frame(str_split(testGene[[x]], "\\|", simplify=T))
		splitInfo <- splitInfo[splitInfo[,4]!="",]
		splitInfo <- splitInfo[splitInfo[,3]=="MODERATE" | splitInfo[,3]=="HIGH",]
		data.frame("Gene"=splitInfo[,4])
		# aa <- str_split(splitInfo[,12], ":", simplify=T)[,2]
		# unique(data.frame("Gene"=splitInfo[,4], "AA"=splitInfo[,16]))
		}) 

	testGeneInfo.m <- do.call(rbind, lapply(testGeneInfo, function(x){
		data.frame("Gene"=paste(unique(x$Gene), collapse=";"), "AA"=paste(unique(x$AA), collapse="; "))
		# data.frame("Gene"=paste(unique(x$Gene), collapse="; "))
		}))
	testInfo <- cbind(testInfo,testGeneInfo.m)
	testInfo <- testInfo[!grepl("chrM|chrX", rownames(testInfo)),]
	# testInfo$Mutation <- paste(testInfo$Gene, testInfo$AA, sep=":")
	testInfo$AF <- as.numeric(testInfo$AF)
	testInfo <- testInfo[order(testInfo$AF, decreasing=TRUE),]
	testInfo
	})
sampleList2 <- lapply(filteredFiles, function(x){
	substr(x, 1, 5)
})

names(dataInfo) <- as.vector(melt(sampleList2)[,1])


dataInfo <- lapply(names(dataInfo), function(x){
	dataInfo[[x]]$sample <- x 
	dataInfo[[x]]
	})
setwd("/rds/general/user/hdhiman/home/lateRelapse_targetedPanel")

dataInfo.m <- do.call(rbind, dataInfo)
dataInfo.m$Variant <- rownames(dataInfo.m)
dataInfo.m$AF <- as.numeric(as.character(dataInfo.m$AF))

geneTable <- as.data.frame(table(dataInfo.m$Gene))
dataInfo.m$count <- geneTable[match(dataInfo.m$Gene, geneTable$Var1),]$Freq
dataInfo.m <- dataInfo.m[dataInfo.m$count >1,]
test <- split(dataInfo.m,dataInfo.m$Gene)
testInfo <- melt(lapply(test, function(x){
	data.frame("sumAF"=sum(x$AF), "maxAF"=max(x$AF), "count"=unique(x$count))
	}), id.vars=c("sumAF", "maxAF", "count"))
dataInfo.m$sumAF <- testInfo[match(dataInfo.m$Gene, testInfo$L1),]$sumAF
dataInfo.m$maxAF <- testInfo[match(dataInfo.m$Gene, testInfo$L1),]$maxAF
dataInfo.m <- dataInfo.m[order(dataInfo.m$sumAF, dataInfo.m$maxAF, dataInfo.m$count), ]
dataInfo.m$Gene <- factor(dataInfo.m$Gene, levels=unique(dataInfo.m$Gene))
test2 <- split(dataInfo.m,dataInfo.m$sample)
testInfo2 <- melt(lapply(test2, function(x){
	data.frame("sumAF"=sum(x$AF), "maxAF"=max(x$AF), "count"=unique(x$count))
	}), id.vars=c("sumAF", "maxAF", "count"))
dataInfo.m$sumAF2 <- testInfo[match(dataInfo.m$Gene, testInfo2$L1),]$sumAF
dataInfo.m$maxAF2 <- testInfo[match(dataInfo.m$Gene, testInfo2$L1),]$maxAF
dataInfo.m$count2 <- testInfo[match(dataInfo.m$Gene, testInfo2$L1),]$count
dataInfo.m$sample <- factor(dataInfo.m$sample, levels=unique(dataInfo.m[order(dataInfo.m$sumAF2, dataInfo.m$maxAF2, -dataInfo.m$count2), ]$sample))

relapseSiteInfo <- read.table("rds/relapseSite.txt", header=TRUE, sep="\t")

dataInfo.m$Location <- relapseSiteInfo[match(dataInfo.m$sample, relapseSiteInfo$Patient),]$Location
dataInfo.m$RelapseSite <- relapseSiteInfo[match(dataInfo.m$sample, relapseSiteInfo$Patient),]$RelapseSite
dataInfo.m$TimeToRelapse <- relapseSiteInfo[match(dataInfo.m$sample, relapseSiteInfo$Patient),]$TimeToRelapse

dataInfoMatrix <- acast(dataInfo.m[order(dataInfo.m$AF, decreasing = TRUE),], Gene~sample, value.var="AF", fun.aggregate = function(x) x[1])
write.table(dataInfoMatrix, "rds/dataInfoMatrix.txt", sep="\t", quote=F, row.names=T, col.names=T)

### Upload druver gene information ###
intogenDrivers <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes.tsv", header=TRUE, sep="\t")
intogenDriversBRCA <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes_BRCA.tsv", header=TRUE, sep="\t")
resDrivers <- c("ESR1", "TP53", "AKT1", "RIC8A", "RB1", "NF1", "FRG1", "KMT2C", "NCOR1")

### Upload dN/dS analysis information for heatmap ###
infoDataDNDS <- readRDS("rds/infoData.rds")
require("ComplexHeatmap")
getPlotInfo <- function(x){
set.seed(1) 
	dataPlot <- x
	dataPlot[is.na(dataPlot)] = 0
dataPlot <- dataPlot[order(100*(rowSums(dataPlot>0)/48), decreasing=TRUE), order(colMaxs(dataPlot), decreasing=TRUE)]
colAnnDownAll <- HeatmapAnnotation("TimeToRelapse" = relapseSiteInfo[relapseSiteInfo$Patient %in% colnames(dataPlot),2], 
								"RelapseSite" = relapseSiteInfo[relapseSiteInfo$Patient %in% colnames(dataPlot),6], 
								which = "column", 
								col=list("TimeToRelapse" = colorRamp2(c(10, 35), c("white", "blue")),
									"RelapseSite" = c("Local" = "grey", "Distal" = "black")))
rowAnnAll <- HeatmapAnnotation("%Samples"=100*(rowSums(dataPlot>0)/48), 
	which = "row", annotation_name_side = "top", annotation_name_rot = 0, 
	col=list("%Samples" = colorRamp2(c(0, 40), c("white", "gold"))))
dataPlot[dataPlot == 0] <- NA
qGlobal_cv <- infoDataDNDS[[1]][[1]]$sel_cv[,c("gene_name", "qglobal_cv")]
qGlobal_cv <- qGlobal_cv[qGlobal_cv$gene_name %in% rownames(dataPlot),]
qGlobal_cv <- qGlobal_cv[match(rownames(dataPlot), qGlobal_cv$gene_name),]
qGlobal_cv[,2][is.na(qGlobal_cv[,2])] = 0 
is_sig = qGlobal_cv[,2] < 0.1
pch = rep("*", length(qGlobal_cv[,2]))
pch[!is_sig] = NA
rowAnnAll2 <- HeatmapAnnotation("dN/dS\nqGlobal_cv"=anno_simple(qGlobal_cv[,2], 
	pch = pch ,
	col = colorRamp2(c(0, 1), c("red", "white"))),
which = "row", annotation_name_side = "top", annotation_name_rot = 0)

# define the legend for pvalue
lgd1 = Legend(title = "dN/dS\nqGlobal_cv", 
	col_fun = colorRamp2(c(0, 1), c("red", "white"))) 
# define the legend for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.1")

ht <- Heatmap(dataPlot, name = "AF",
	col = viridis(100, begin = 0.2, end = 1, direction = -1, option = "G"),
    row_order = rownames(dataPlot), 
    column_order = colnames(dataPlot),
	na_col = "white",
	show_column_names = FALSE,
	show_row_names = TRUE,
	row_names_side ="left",
	show_column_dend = FALSE,
	show_row_dend = FALSE,
	left_annotation=rowAnnAll,
	right_annotation=rowAnnAll2,
	bottom_annotation=colAnnDownAll)
draw(ht, annotation_legend_list = list(lgd1, lgd_sig))

}

###### All ######
dataAll <- dataInfoMatrix

pdf("plots/varHeatmapPASSgermlineAll2.pdf", height=12, width=12)	
getPlotInfo(dataAll)
dev.off()

###### Intogen ######
dataIntogen <- dataInfoMatrix[rownames(dataInfoMatrix) %in% intogenDrivers$Symbol, ]
pdf("plots/varHeatmapPASSgermlineIntogen2.pdf", height=12, width=12)	
getPlotInfo(dataIntogen)
dev.off()

###### IntogenBRCA ######

dataIntogenBRCA <- dataInfoMatrix[rownames(dataInfoMatrix) %in% intogenDriversBRCA$Symbol, ]
pdf("plots/varHeatmapPASSgermlineIntogenBRCA2.pdf", height=10, width=12)	
getPlotInfo(dataIntogenBRCA)
dev.off()

###### Resistance drivers ######
dataRes <- dataInfoMatrix[rownames(dataInfoMatrix) %in% resDrivers, ]

pdf("plots/varHeatmapPASSgermlineRes2.pdf", height=3, width=12)	
getPlotInfo(dataRes)
dev.off()
