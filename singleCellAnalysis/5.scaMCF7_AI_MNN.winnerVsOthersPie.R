require("Seurat")
require("dplyr")
require("ggplot2")
require("scRNAseq")
require("scuttle")
require("celldex")
require("Matrix")
require("stringr")
require("reshape2")
require("SeuratWrappers")
require("enrichR")
require("DoubletFinder")
require("RColorBrewer")

### Winner Vs Others pie charts ###

setwd("~/manuscript")
singlets.combined <- readRDS("/rds/general/user/hdhiman/ephemeral/scRNAseq/cellranger/singlets.combinedAI.MNN0.4.rds")
colors <- c("#ffa300", "#B1961A", "#e6d800", "#50e991", "#9b19f5", "#dc0ab4", "#1F9E89FF", "#e60049", "#b3d4ff", "blue", "#0bb4ff", "grey", "brown")
names(colors)<-0:12
colors3 <- c("#e6d800", "#ffa300", "#50e991", "#9b19f5", "#1F9E89FF", "#e60049", "#b3d4ff", "blue", "#0bb4ff", "grey")
names(colors3)<-0:9
files <- c("T0_1", "T0_2","AI1m1", "AI2m1", "AI1m2", "AI2m2","AIaw1", "AIaw2", "AIaw3", "AIaw4")
filesOrder <- c("AI1m1", "AI2m1", "AI1m2", "AI2m2","AIaw1", "AIaw2", "AIaw3", "AIaw4", "T0_1", "T0_2")

singlets.combined4$orig.ident<-factor(singlets.combined4$orig.ident, levels=filesOrder)	

singlets.combined<-singlets.combined4
##########################
#### cell cycle phase ####
##########################
			 	
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined <- CellCycleScoring(singlets.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#####################################
#### Assign barcodes to cell IDs ####
#####################################
# setwd("/rds/general/user/hdhiman/home/scRNAseq/cellRanger/")
# cellecta.10X<-readRDS("rds/cellecta.10XAll.rds")
cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
colnames(cellecta.10X)<-c("value", "L1", "sample")
# files <- filesOrder
cellecta.10X$sampleCell<-paste(cellecta.10X$sample, cellecta.10X$L1, sep="_")
cellecta.10X <- cellecta.10X[cellecta.10X$sample %in% files,]


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
singlets.combined.sub <- subset(singlets.combined.sub, bc != "complex")

singlets.combined.sub$orig.ident <- factor(singlets.combined.sub$orig.ident, levels=files)

getPct <- function(sampleName, bc){
	x <- subset(singlets.combined.sub, orig.ident == sampleName)
	reducedBCs <- x[[]]$bc
	reducedBCs[!reducedBCs %in% c(bc, "complex")] <- "others"
	x <- AddMetaData(
	    object = x,
	    metadata = reducedBCs,
	    col.name = "reducedBCs")
	x <- subset(x, reducedBCs %in% c(bc, "others"))
	x.winners <- subset(x, reducedBCs %in% bc)
	x.others <- subset(x, reducedBCs == "others")

	countsWinner <- table(x.winners$seurat_clusters)
	countsOthers <- table(x.others$seurat_clusters)
	counts <- table(x$seurat_clusters)
	data.frame("winner"=100*(countsWinner/sum(countsWinner)),
				"others"=100*(countsOthers/sum(countsOthers)))[,c(1,2,4)]
}
selWinners<-c("bc-7704123","bc-6864523","bc-2127641", "bc-5122096")
pctInfoAIaw1 <- getPct("AIaw1", "bc-7704123")
pctInfoAIaw2 <- getPct("AIaw2", "bc-6864523")
pctInfoAIaw3 <- getPct("AIaw3", "bc-2127641")
pctInfoAIaw4 <- getPct("AIaw4", "bc-5122096")
pctInfoList <- list(
					"AIaw1"=pctInfoAIaw1, "AIaw2"=pctInfoAIaw2, "AIaw3"=pctInfoAIaw3, "AIaw4"=pctInfoAIaw4)
pctInfo.m <- melt(pctInfoList)
colors<-colors3
pdf("~/manuscript/plots/sca/MNNccRegressed/pctInfoWinnerVsOthers.pdf", height=7, width=7)
ggplot(pctInfo.m, aes(x="", y=value, fill=winner.Var1)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.001) +
  coord_polar("y", start=0, clip="off") +  
  theme_void() +
  scale_fill_manual(values=colors) +
  labs(fill="Clusters") +
  facet_grid(L1~variable) +
  theme(legend.position="bottom",
    strip.text = element_text(size = 16), 
    legend.text=element_text(size= 16)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()


getPctAll <- function(sampleName){
	x <- subset(singlets.combined.sub, orig.ident == sampleName)
	winnerBC <- c("bc-7704123", "bc-6864523", "bc-2127641", "bc-5122096")
	reducedBCs <- x[[]]$bc
	reducedBCs[!reducedBCs %in% c(winnerBC, "complex")] <- "others"
	x <- AddMetaData(
	    object = x,
	    metadata = reducedBCs,
	    col.name = "reducedBCs")
	x <- subset(x, reducedBCs %in% c(winnerBC, "others"))

	bcInfoList <- split(x[[]], x[[]]$reducedBCs)
	bcInfoPct <- lapply(bcInfoList, function(splitBC){
		counts <- table(splitBC$seurat_clusters)
		pct <- 100*(counts/sum(counts))
		pct
		})
	bcInfoPct.m <- melt(bcInfoPct)
	colnames(bcInfoPct.m) <- c("cluster", "pct", "bc")
	bcInfoPct.m
}

pctAllBc <- lapply(files[1:6], getPctAll)
names(pctAllBc) <- files[1:6]
pctAllBc.m <- melt(pctAllBc, id.vars=c("bc", "cluster"), value="pct")
pctAllBc.m$cluster <- as.factor(pctAllBc.m$cluster)
pctAllBc.m <- pctAllBc.m[pctAllBc.m$value>0,]
pctAllBc.m$L1 <- factor(pctAllBc.m$L1, levels=files)
pctAllBc.m$bc <- factor(pctAllBc.m$bc, levels=c("bc-5122096", "bc-7704123", "bc-2127641", "bc-6864523", "others"))

pdf("~/manuscript/plots/sca/MNNccRegressed/pctInfoWinnerVsOthersAllBCs.pdf", height=12, width=15)
ggplot(pctAllBc.m, aes(x="", y=value, fill=cluster)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.00001) +
  coord_polar("y", start=0, clip="off") +  
  theme_void() +
  scale_fill_manual(values=colors) +
  labs(fill="Clusters") +
  facet_grid(bc~L1) +
  theme(legend.position="bottom",
    strip.text = element_text(size = 16), 
    legend.text=element_text(size= 16)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
dev.off()


