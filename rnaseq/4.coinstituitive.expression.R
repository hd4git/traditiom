require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")
require("ggrepel")
require("pheatmap")
require("ggplotify")
require("viridis")
require("cowplot")

setwd("~/rnaseq/results/traditiom")
fpkm<-readRDS("rds/fpkm.rds")
pa.up<-read.table("rds/pa.up.new.txt", header=F, stringsAsFactors=F)[,1]
pa.down<-read.table("rds/pa.down.new.txt", header=F, stringsAsFactors=F)[,1]
info.m<-read.table("rds/info.m.txt", header=T, stringsAsFactors=F)
cellCycle <- c("ABL1", "ANAPC1", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "ATM", "ATR", "BUB1", "BUB1B", "BUB3", "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNH", "CDC14A", "CDC14B", "CDC16", "CDC20", "CDC23", "CDC25A", "CDC25B", "CDC25C", "CDC26", "CDC27", "CDC45", "CDC6", "CDC7", "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "CHEK1", "CHEK2", "CREBBP", "CUL1", "DBF4", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "EP300", "ESPL1", "FZR1", "GADD45A", "GADD45B", "GADD45G", "GSK3B", "HDAC1", "HDAC2", "MAD1L1", "MAD2L1", "MAD2L2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MDM2", "MYC", "ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "PCNA", "PKMYT1", "PLK1", "PRKDC", "PTTG1", "PTTG2", "RAD21", "RB1", "RBL1", "RBL2", "RBX1", "SFN", "SKP1", "SKP1P2", "SKP2", "SMAD2", "SMAD3", "SMAD4", "SMC1A", "SMC1B", "SMC3", "STAG1", "STAG2", "TFDP1", "TFDP2", "TGFB1", "TGFB2", "TGFB3", "TP53", "TTK", "WEE1", "WEE2", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ", "ZBTB17")
mat<-readRDS("rds/lnCounts.rds")
deg.info.list<-readRDS("rds/deg.info.list2.rds")
deg.info.dormancyList <- deg.info.list[names(deg.info.list)[grepl("AI30|AI60|AI90|T30|T60", names(deg.info.list))]]
deg.common.dormancyList <- Reduce(intersect,deg.info.dormancyList)

enrichMarkers <- function(x){
	res <- enrichr(x, "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020
	res$count <- str_split(res$Overlap, "/", simplify=TRUE)[,1]
	res[as.numeric(res$count) >= 3 & res$Adjusted.P.value<=0.05, ]
}

matrix  <- mat - rowMeans(mat)
sd<-rowSds(matrix)
matrix<-matrix/sd
samplesDormancy <- colnames(matrix)[grepl("AI30|AI60|AI90|T30|T60", colnames(matrix))]

##### Dormancy Up
matrix.sub1 <- matrix[rowSums(matrix[,samplesDormancy]>0)==10 & rowSums(matrix[,c("POT.1", "POT.2", "POT.3")]<0)==3,]
matrix.sub1 <- matrix.sub1[rownames(matrix.sub1) %in% deg.common.dormancyList,]
sub1 <- matrix.sub1[,samplesDormancy]
sub1 <- sub1[rowSums(sub1>=quantile(sub1, 0.1))>6,]
rowZ1 <- rowSums(sub1)
sub1 <- sub1[order(rowSums(sub1), decreasing = TRUE),]
dormancySign1 <- rownames(sub1)
dormancySignwoCC1 <- dormancySign1[!dormancySign1 %in% cellCycle]
# cov1 <- rowSds(sub1)/rowMeans(sub1)
# sub1.top50 <- names(head(cov1[order(cov1)], n=50L))

sub1.top50 <- head(sub1[order(rowSds(sub1)/rowMeans(sub1), -1*(rowSums(sub1))),], n=50L)


##### Dormancy Down
matrix.sub2 <- matrix[rowSums(matrix[,samplesDormancy]<0)==10 & rowSums(matrix[,c("POT.1", "POT.2", "POT.3")]>0)==3,]
matrix.sub2 <- matrix.sub2[rownames(matrix.sub2) %in% deg.common.dormancyList,]
sub2 <- matrix.sub2[,samplesDormancy]
sub2 <- sub2[rowSums(sub2<=quantile(sub2, 0.9))>6,]
rowZ2 <- rowSums(sub2)
sub2 <- sub2[order(rowSums(sub2), decreasing = TRUE),]
dormancySign2 <- rownames(sub2)
dormancySignwoCC2 <- dormancySign2[!dormancySign2 %in% cellCycle]
# cov2 <- rowSds(sub2)/rowMeans(sub2)
# sub2.top50 <- names(head(cov2[order(cov2)], n=50L))
sub2.top50 <- head(sub2[order(rowSds(sub2)/rowMeans(sub2), -1*(abs(rowSums(sub2)))),], n=50L)


getHeatmapInfo <- function(dormancySignwoCC, sub){
		sel <- dormancySignwoCC
		matrix<-mat[rownames(mat) %in% sel, ]
		matrix  <- matrix - rowMeans(matrix)
		sd<-rowSds(matrix)
		matrix<-matrix/sd
		matrix <- matrix[dormancySignwoCC,]
		matrix.m<-melt(matrix)

		info<-getState(matrix.m)
		info<-info[info$state!="Latency",]
		info$Var1 <- factor(info$Var1, levels=rownames(sub))
		info$state2 <- info$state
		info[info$state=="Dormancy" & grepl("^AI", info$Var2), ]$state2 <- "AI.Dormancy"
		info[info$state=="Dormancy" & grepl("^T", info$Var2), ]$state2 <- "T.Dormancy"
		info[info$state=="Awakening" & grepl("^AI", info$Var2), ]$state2 <- "AI.Awakening"
		info[info$state=="Awakening" & grepl("^T", info$Var2), ]$state2 <- "T.Awakening"
		info[info$state=="TEPs" & grepl("^AI", info$Var2), ]$state2 <- "AI.TEPs"
		info[info$state=="TEPs" & grepl("^T", info$Var2), ]$state2 <- "T.TEPs"
		info$state2 <- factor(info$state2, levels=unique(info$state2))
		info
}
getState<-function(x){
		x$state <- "Awakening"
		if(any(grepl("POT", x$Var2))){ x[grepl("POT", x$Var2),]$state <- "POT"}
		if(any(grepl("^U", x$Var2))){ x[grepl("^U", x$Var2),]$state <- "Untreated"}
		if(any(grepl("^AI7\\.|^AI14\\.|^T7\\.|^T14\\.", x$Var2))){ x[grepl("^AI7\\.|^AI14\\.|^T7\\.|^T14\\.", x$Var2),]$state <- "Latency"}
		if(any(grepl("^AI30\\.|^AI60\\.|^AI90\\.|^T30\\.|^T60\\.", x$Var2))){ x[grepl("^AI30\\.|^AI60\\.|^AI90\\.|^T30\\.|^T60\\.", x$Var2),]$state <- "Dormancy"}
		if(any(grepl("_TEP$", x$Var2))){ x[grepl("_TEP$", x$Var2),]$state <- "TEPs"}
		# x$state<-factor(x$state, levels=c("POT","Untreated", "Latency", "Dormancy", "Awakening", "TEPs"))
		x
}

getMat <- function(dormancySignwoCC, sub){
		sel <- dormancySignwoCC
		matrix<-mat[rownames(mat) %in% sel, ]
		matrix
}

infoUpMat <- getMat(dormancySignwoCC1, sub1)
infoDownMat <- getMat(dormancySignwoCC2, sub2)
write.table(infoUpMat, "rds/infoUpMat.txt", sep="\t", quote=F, row.names=T, col.names=T)
write.table(infoDownMat, "rds/infoDownMat.txt", sep="\t", quote=F, row.names=T, col.names=T)

infoUp <- getHeatmapInfo(dormancySignwoCC1, sub1)
infoUp$enriched <- "DormancyUp"
infoDown <- getHeatmapInfo(dormancySignwoCC2, sub2)
infoDown$enriched <- "DormancyDown"
infoUp50 <- getHeatmapInfo(rownames(sub1.top50), sub1)
infoUp50$enriched <- "DormancyUp"
infoDown50 <- getHeatmapInfo(rownames(sub2.top50), sub2)
infoDown50$enriched <- "DormancyDown"

info <- rbind(infoUp, infoDown)
info$enriched <- factor(info$enriched, levels=c("DormancyUp", "DormancyDown"))

info50 <- rbind(infoUp50, infoDown50)
info50$enriched <- factor(info50$enriched, levels=c("DormancyUp", "DormancyDown"))


pdf("~/manuscriptFeb23/plots/figure3c.dormancy.pdf", height=5, width=5)
ggplot(info[!info$Var2 %in% c("U170.3"),], aes(x=Var2, y=Var1, fill=value)) +
	geom_tile() +
	theme_bw() +
	scale_fill_gradient2(low = 'yellow', mid = 'white', high = '#50157A') +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=8), 
		strip.background=element_rect(fill="white"), 
		legend.position="bottom",
		# axis.text.y=element_text(size=5),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.ticks.x=element_line(size=0.1),
		axis.title=element_blank(),
		strip.text.x = element_text(size = 8),
		legend.key.size = unit(0.3, 'cm'),
		legend.text = element_text(size=6),
		legend.title = element_text(size=8),
		panel.spacing = unit(0.02, "lines")) + 
	labs( x="Samples", y="Genes", fill="Z score") +
	facet_grid(enriched~state2, scales="free", space="free")
dev.off()

pdf("~/manuscriptFeb23/plots/figure3c.dormancy50.pdf", height=10, width=7)
ggplot(info50[!info50$Var2 %in% c("U170.3"),], aes(x=Var2, y=Var1, fill=value)) +
	geom_tile() +
	theme_bw() +
	scale_fill_gradient2(low = 'yellow', mid = 'white', high = '#50157A') +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=8), 
		strip.background=element_rect(fill="white"), 
		legend.position="bottom",
		axis.text.y=element_text(size=5),
		# axis.text.y=element_blank(),
		# axis.ticks.y=element_blank(),
		axis.ticks.x=element_line(size=0.1),
		axis.title=element_blank(),
		strip.text.x = element_text(size = 8),
		legend.key.size = unit(0.3, 'cm'),
		legend.text = element_text(size=6),
		legend.title = element_text(size=8),
		panel.spacing = unit(0.02, "lines")) + 
	labs( x="Samples", y="Genes", fill="Z score") +
	facet_grid(enriched~state2, scales="free", space="free")
dev.off()


#### GSEA 
dormancyUp <- enrichMarkers(dormancySignwoCC1)
dormancyDown <- enrichMarkers(dormancySignwoCC2)

dormancyPAup <- dormancySignwoCC1[dormancySignwoCC1 %in% pa.up] # 52 genes
dormancyPAdown <- dormancySignwoCC2[dormancySignwoCC2 %in% pa.down] # 318 genes

dormancyUp$enriched <- "DormancyUp"
dormancyDown$enriched <- "DormancyDown"
dormancyDown$Adjusted.P.value <- -1*dormancyDown$Adjusted.P.value

enrichedList <- rbind(dormancyUp, dormancyDown)
overlapInfo <- as.data.frame(str_split(enrichedList$Overlap, "/", simplify=TRUE))
enrichedList$ratio <- as.numeric(overlapInfo[,1])/as.numeric(overlapInfo[,2])
enrichedList$Term <- factor(enrichedList$Term, levels = enrichedList[order(enrichedList$Adjusted.P.value),]$Term)
enrichedList$enriched <- factor(enrichedList$enriched, levels = c("DormancyUp", "DormancyDown"))

pdf("~/manuscriptFeb23/plots/figure3c.dormancyEnriched.pdf", height=5, width=6)
ggplot(enrichedList, aes(x=Adjusted.P.value, y=Term, size=ratio, fill=enriched)) +
	geom_point(shape=21, color="black") +
	theme_bw() +
	facet_grid(enriched~., space="free", scale="free") +
	scale_fill_manual(values = c("DormancyUp"="#50157A", "DormancyDown"="yellow")) +
	labs(fill="Signature", size="OverlapRatio") +
	theme(axis.title.y = element_blank(),
		strip.background=element_rect(fill="white"), 
		legend.position = "bottom") +
	guides(fill=guide_legend(nrow=2,byrow=TRUE), 
		size=guide_legend(nrow=2,byrow=TRUE))
dev.off()

enrichMarkersGOST <- function(x){
	gprofiler2.addon::gost_custom_gmt(x,
	correction_method = "false_discovery_rate", 
	custom_gmt = "rds/h.all.v2023.1.Hs.symbols.pa.gmt",
	user_threshold = 0.05)[[1]]
}

dormancyUp <- enrichMarkersGOST(dormancySignwoCC1)
dormancyDown <- enrichMarkersGOST(dormancySignwoCC2)


dormancyUp$enriched <- "DormancyUp"
dormancyDown$enriched <- "DormancyDown"
dormancyDown$p_value <- -1*dormancyDown$p_value

enrichedList <- rbind(dormancyUp, dormancyDown)
enrichedList$ratio <- as.numeric(enrichedList$intersection_size)/as.numeric(enrichedList$term_size)
enrichedList$term_id <- str_replace(enrichedList$term_id, "HALLMARK_", "")
enrichedList$term_id <- factor(enrichedList$term_id, levels = unique(enrichedList[order(enrichedList$p_value),]$term_id))
enrichedList$enriched <- factor(enrichedList$enriched, levels = c("DormancyUp", "DormancyDown"))


pdf("~/manuscriptFeb23/plots/figure3c.dormancyEnrichedGOST2.pdf", height=4.5, width=3)
ggplot(enrichedList, aes(x=p_value, y=term_id, size=ratio, fill=enriched)) +
	geom_point(shape=21, color="black") +
	theme_bw() +
	facet_grid(enriched~., space="free", scale="free") +
	scale_fill_manual(values = c("DormancyUp"="#50157A", "DormancyDown"="yellow")) +
	labs(fill="Signature", size="OverlapRatio") +
	theme(axis.title.y = element_blank(),
		strip.background=element_rect(fill="white"), 
		legend.position = "none", 
		axis.text = element_blank(),
		axis.title = element_blank()) +
	guides(fill=guide_legend(nrow=2,byrow=TRUE), 
		size=guide_legend(nrow=2,byrow=TRUE))
dev.off()

