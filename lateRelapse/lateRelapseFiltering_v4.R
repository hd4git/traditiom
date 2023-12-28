cd ~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr

while read f; 
do 
grep "#" $f\_tumor-only_mutect2_filtered_chr.vcf | grep -vE "##contig=<ID=chr[A-Za-z0-9]*_|##contig=<ID=HLA|##contig=<ID=chrEBV" > $f\_PASSgermline.vcf
grep -v "#" $f\_tumor-only_mutect2_filtered_chr.vcf | grep -E "PASS|germline" >> $f\_PASSgermline.vcf
bgzip $f\_PASSgermline.vcf 
tabix -p vcf $f\_PASSgermline.vcf.gz
done < <(ls | awk "{split($0,a,"_"); print a[1]"_"a[2]}")

while read f; 
do 
bcftools view $f\.vcf.gz --regions-file ~/wgs/oncomineCoveredByProbesHg38.bed > $f\_targeted.vcf 
bgzip $f\_targeted.vcf 
tabix -p vcf $f\_targeted.vcf.gz 
done < <(ls *vcf.gz | awk "{split($0, a, "."); print a[1]}")



require("VariantAnnotation")
require("stringr")
require("reshape2")
setwd("~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr")
files <- list.files(pattern = "\\PASSgermline_targeted.vcf.gz$")

dataFilteredList <- lapply(files, function(x){
	sample <- substr(x, 1, 5)
	data <- read.table(x, header=F, sep="\t")
	AD <- str_split(data$V10, ":", simplify=T)[,2]
	ADref <- as.numeric(str_split(AD, ",", simplify=T)[,1])
	ADalt <- as.numeric(str_split(AD, ",", simplify=T)[,2])
	ADsum <- ADref + ADalt
	AF <- str_split(data$V10, ":", simplify=T)[,3]
	F1R2 <- str_split(data$V10, ":", simplify=T)[,5]
	F2R1 <- str_split(data$V10, ":", simplify=T)[,6]
	ALT_F1R2 <- as.numeric(str_split(F1R2, ",", simplify=T)[,2])
	ALT_F2R1 <- as.numeric(str_split(F2R1, ",", simplify=T)[,2])
	ALT_F1R2_F2R1 <- ALT_F1R2+ALT_F2R1
	data$ADsum <- ADsum
	data$AF <- AF
	data$ALT_F1R2_F2R1 <- ALT_F1R2_F2R1
	dataFiltered <-  data[data$ADsum >=20 & data$ALT_F1R2_F2R1 >=4 & data$AF >=0.1,]
	dataFiltered
})

sampleList <- lapply(files, function(x){
	substr(x, 1, 5)
})

names(dataFilteredList) <- melt(sampleList)[,1]
setwd("~/lateRelapse_targetedPanel/")
lapply(names(dataFilteredList), function(x){
	write.table(dataFilteredList[[x]], file = paste(paste("rds/PON_PASSgermline_chr_filtered", x, sep="/"), "vcf", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
	})

cd ~/lateRelapse_targetedPanel/rds

while read f; 
do 
zgrep "#" ~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr/$f\_PASSgermline_targeted.vcf.gz  > ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header/$f\.vcf
cat PON_PASSgermline_chr_filtered/$f\.vcf >> ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header/$f\.vcf
done < <(ls PON_PASSgermline_chr_filtered | awk "{split($0,a,"."); print a[1]}")


cd ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header
while read f; 
do 
bgzip $f
tabix -p vcf $f\.gz
done < <(ls)

while read f; 
do 
bcftools view $f\.vcf.gz --regions-file ~/wgs/oncomineCoveredByProbesHg38.bed > $f\_targeted.vcf 
bgzip $f\_targeted.vcf 
tabix -p vcf $f\_targeted.vcf.gz 
done < <(ls *.vcf.gz | awk "{split($0, a, "."); print a[1]}")


source activate genomics
while read f; 
do 
vcf=$f\.vcf.gz 
ofile=$f\_vep.vcf.gz
ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205/Homo_sapiens_assembly38.fasta
vep --everything \
-i $vcf \
-o $ofile \
--vcf \
--offline \
--fasta $genome_fasta \
--compress_output bgzip \
--allele_number \
--dir_cache ~/wgs/GRCh38/vep/
done < <(ls *targeted.vcf.gz | awk "{split($0, a, "."); print a[1]}")

while read f; 
do 
zgrep -vE "germline;" $f\.vcf.gz | zgrep -vE ";germline" > $f\_filtered.vcf.gz
done < <(ls *vep.vcf.gz | awk "{split($0, a, "."); print a[1]}")

while read f; 
do 
zgrep "#" $f\.vcf.gz > $f\_filtered2.vcf
zgrep -E "MODERATE|HIGH\|" $f\.vcf.gz >> $f\_filtered2.vcf
bgzip $f\_filtered2.vcf
tabix -p vcf $f\_filtered2.vcf.gz
done < <(ls *_filtered.vcf.gz | awk "{split($0, a, "."); print a[1]}")




source activate jvCalls
 
require("VariantAnnotation")
require("stringr")
require("ggplot2")
require("reshape2")
require("circlize")
require("viridis")

setwd("/rds/general/user/hdhiman/home/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header")
filteredFiles <- list.files(pattern = "\\_filtered_filtered2.vcf.gz$")
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


intogenDrivers <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes.tsv", header=TRUE, sep="\t")
intogenDriversBRCA <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes_BRCA.tsv", header=TRUE, sep="\t")
resDrivers <- c("ESR1", "TP53", "AKT1", "RIC8A", "RB1", "NF1", "FRG1", "KMT2C", "NCOR1")
dataInfoMatrix <- acast(dataInfo.m[order(dataInfo.m$AF, decreasing = TRUE),], Gene~sample, value.var="AF", fun.aggregate = function(x) x[1])


write.table(dataInfoMatrix, "rds/dataInfoMatrix.txt", sep="\t", quote=F, row.names=T, col.names=T)
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

# see how we define the legend for pvalue
lgd1 = Legend(title = "dN/dS\nqGlobal_cv", 
	col_fun = colorRamp2(c(0, 1), c("red", "white"))) 
	# at = c(0, 1, 2, 3), 
    # labels = c("1", "0.1", "0.01", "0.001"))
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.1")
# these two self-defined legends are added to the plot by `annotation_legend_list`

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

