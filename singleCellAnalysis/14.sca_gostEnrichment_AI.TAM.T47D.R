require("monocle3")
require("Seurat")
require("dplyr")
require("org.Hs.eg.db")
require("garnett")
require("Matrix")
require("ggplot2")
require("viridis")
require("RColorBrewer")
require("gprofiler2")

setwd("~/manuscript")
singlets.combined.ai <- readRDS("~/manuscript/rds/singlets.combinedAI.MNN0.4.rds")
singlets.combined.tam <- readRDS("~/manuscriptFeb23/rds/singlets.combined.TAM.rds")
singlets.combined.t47d <- readRDS("~/scRNAseq/T47D/rds/singlets.combined.T47D.info.rds")

singlets.combined.list <- list("AI" = singlets.combined.ai, 
								"TAM" = singlets.combined.tam, 
								"T47D" = singlets.combined.t47d)

markerist <- lapply(singlets.combined.list, function(x){
	FindAllMarkers(x, test.use="roc", min.pct = 0.25, only.pos = TRUE)
})

geneLists <- lapply(markerist, function(x){
	geneLists.l<-split(x, x$cluster)
	geneLists.gs<-lapply(geneLists.l, function(x){x$gene})
	names(geneLists.gs)<-names(geneLists.l)
	geneLists.gs
})


geneLists.l<-split(markerist[[1]], markerist[[1]]$cluster)
geneLists.gs<-lapply(geneLists.l, function(x){x$gene})
names(geneLists.gs)<-names(geneLists.l)
geneLists.gmt<-lapply(names(geneLists.gs), function(x){
	genes<-paste(geneLists.gs[[x]], collapse="	")
	data.frame(paste("cluster", x, sep=""), paste("cluster", x, sep=""), genes)
	})
geneLists.gmt<-do.call(rbind, geneLists.gmt)
write.table(geneLists.gmt,file ="~/manuscriptFeb23/rds/geneListsAI.gmt",sep="\t",quote=F,row.names=F,col.names=F)


enrichMarkersGOST <- function(x){
	gprofiler2.addon::gost_custom_gmt(x,
	correction_method = "false_discovery_rate", 
	custom_gmt = "~/manuscriptFeb23/rds/geneListsAI.gmt",
	user_threshold = 0.05)[[1]]
}

enrichedTAM <- lapply(geneLists[[2]], function(x){
		if(length(x)>0)
		{enrichMarkersGOST(x)}
})
enrichedT47D <- lapply(geneLists[[3]], function(x){
		if(length(x)>0)
		{enrichMarkersGOST(x)}
})
names(enrichedTAM) <- paste("TAM",0:7, sep="_")
names(enrichedT47D) <- paste("T47D",0:8, sep="_")

enrichedLists1 <- do.call(rbind, enrichedTAM)
enrichedLists1$sample <- "TAM"
enrichedLists2 <- do.call(rbind, enrichedT47D)
enrichedLists2$sample <- "T47D"
enriched.m <- rbind(enrichedLists1[enrichedLists1$intersection_size>=3,], enrichedLists2[enrichedLists2$intersection_size>=3,])
enriched.m$sampleCluster <- str_split(rownames(enriched.m), "\\.", simplify=TRUE)[,1]
enriched.m$intersectionRatio <- enriched.m$intersection_size/enriched.m$query_size
enriched.m$term_id <- str_replace(enriched.m$term_id, "cluster", "AI_")
enriched.m$term_id <- factor(enriched.m$term_id, levels=rev(paste("AI", c(3,6,8,2,1,10,0, 9,5,7,11,12, 4), sep="_")))
enriched.m$sampleCluster <- factor(enriched.m$sampleCluster, levels=c(paste("T47D", c(1,2,7,8,4,0,5,3,6), sep="_"),
	paste("TAM", c(1,6,7,0,2,3,4,5), sep="_")))
enriched.m$sample <- factor(enriched.m$sample, levels=c("TAM", "T47D"))
colors2 <- c("T47D_0"="#E6D800", "T47D_1"="#4DAF4A", "T47D_2"="#50e991", "T47D_3"="#FFA300", "T47D_4"="#377EB8", "T47D_5"="#B1961A", "T47D_6"="#E60049", "T47D_7"="#0bb4ff", "T47D_8"="#B3D4FF",
	"TAM_0"="#E6D800", "TAM_1"="#50E991", "TAM_2"="#B1961A", "TAM_3"="#FFA300", "TAM_4"="#E60049", "TAM_5"="#DC0AB4", "TAM_6"="#4DAF4A", "TAM_7"="#B3D4FF")

pdf("~/manuscriptFeb23/plots/gost_AI.TAM.T47D.pdf", height=3, width=5)
ggplot(enriched.m, aes(x=sampleCluster, y=term_id, color=sampleCluster, size=-log10(p_value))) +
	geom_point() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		axis.title = element_blank(),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="none") +
	scale_color_manual(values=colors2) +
	facet_grid(.~sample, space="free", scale="free")
dev.off()

pdf("~/manuscriptFeb23/plots/gost_AI.TAM.T47D_legend.pdf", height=7, width=12)
ggplot(enriched.m, aes(x=sampleCluster, y=term_id, color=sampleCluster, size=-log10(p_value))) +
	geom_point() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		axis.title = element_blank(),
		strip.background = element_rect(fill="#FFFFFF")) +
	scale_color_manual(values=colors2) +
	facet_grid(.~sample, space="free", scale="free")
dev.off()



