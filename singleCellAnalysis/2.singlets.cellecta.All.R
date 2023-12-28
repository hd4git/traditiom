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
#####################################
#### Assign barcodes to cell IDs ####
#####################################
singlets <- readRDS("/rds/general/user/hdhiman/home/scRNAseq/cellRanger/rds/singlets.cellectaAll.rds")
names(singlets) <- c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")
files <- names(singlets)
cellecta.10XAll<-lapply(files, function(x){
	cellectaBCsInfo<-singlets[[x]][['Custom']]
	cellectaBCsInfo.m<-as.matrix(GetAssayData(cellectaBCsInfo))
	melt(cellectaBCsInfo.m)
})
##### Get stats for cells with cellecta barcodes ###########

stats <- lapply(cellecta.10XAll, function(test){
	total <- length(unique(test$Var2))
	withBCs <- length(unique(test[test$value>0,]$Var2))
	freq <- as.data.frame(table(test[test$value>0,]$Var2))$Freq
	list("cellFreq" = withBCs/total, "bcFreq" = freq)
})
names(stats) <- names(singlets)

statsCellsWithBCs <- melt(lapply(stats, function(x){x[[1]]})) 
statsCellsWithBCs$info <- "CellsWithBCs"

pdf("~/manuscript/plots/statsCellsWithBCs.pdf", height=8, width=5)
ggplot(statsCellsWithBCs, aes(x=info, y=value, color=value, label=L1)) +
	geom_violin() +
	geom_jitter(position=position_jitter(0.2)) +
	geom_text_repel(data=statsCellsWithBCs[statsCellsWithBCs$value<0.7,], size=5, color="black") +
	# ylim(0,1) +
	scale_color_viridis(direction=-1) +
	labs(x="", y="Frequency of cells with BCs", color="Freq") +
	theme_bw() +
	theme(legend.position="bottom",
		text = element_text(size = 15))
dev.off()

statsBCfreq <- melt(lapply(stats, function(x){x[[2]]})) 
statsBCfreq$L1 <- factor(statsBCfreq$L1, levels=names(singlets))
pdf("~/manuscript/plots/statsBCfreq.pdf", height=3, width=15)
ggplot(statsBCfreq, aes(x=L1, y=value)) +
	geom_violin() +
	theme_bw() +
	theme(text = element_text(size = 15)) +
	labs(x="Samples", y="Number of BCs per cell")
dev.off()

#### Annotate cells with unique cellecta barcode IDs ##########
names(cellecta.10XAll)<-files

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
saveRDS(cellecta.10X, "~/manuscript/rds/cellecta.10XAll.rds")

