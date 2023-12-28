
singlets.combined <- readRDS("~/manuscript/rds/singlets.combined.MNNccRegressedCluster0.45Info.rds")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
singlets.combined <- CellCycleScoring(singlets.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


cellecta.10X<-readRDS("~/manuscript/rds/cellecta.10XAll.rds")
colnames(cellecta.10X)<-c("value", "L1", "sample")
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
selWinners<-c("bc-7704123","bc-6864523","bc-2127641", "bc-5122096")
singlets.combined.sub <- subset(singlets.combined.sub, bc != "complex")  

# info<-singlets.combined.sub[[]]
# info$cyclingWinners<-"G1.others"
# info[info$bc %in% selWinners & info$Phase %in% c("G2M", "S"),]$cyclingWinners<-"G2MS.winners"
# info[info$bc %in% selWinners & info$Phase == "G1",]$cyclingWinners<-"G1.winners"
# info[!info$bc %in% selWinners & info$Phase %in% c("G2M", "S"),]$cyclingWinners<-"G2MS.others"
# singlets.combined.sub <- AddMetaData(
#     object = singlets.combined.sub,
#     metadata = info$cyclingWinners,
#     col.name = 'cyclingWinners')

info2<-singlets.combined.sub[[]]
info2$stage<-""
info2[info2$orig.ident %in% c("T0_1", "T0_2"),]$stage="T0"
info2[info2$orig.ident %in% c("AI1m1", "AI2m1", "AI1m2", "AI2m2"),]$stage="Dormancy"
info2[info2$orig.ident %in% c("AIaw1", "AIaw2", "AIaw3", "AIaw4"),]$stage="Awakening"

singlets.combined.sub <- AddMetaData(
    object = singlets.combined.sub,
    metadata = info2$stage,
    col.name = 'stage')
singlets.combined.sub$stage<-factor(singlets.combined.sub$stage, levels=c("T0", "Dormancy", "Awakening"))

### Subset G2M and S phase cells with winning barcodes in awakening and all cells in dormancy ###

singlets.combined.subG2MSwinners2 <- subset(singlets.combined.sub, (Phase %in% c("G2M", "S") & bc %in% selWinners & stage == "Awakening")
     | (Phase %in% c("G2M", "S") & stage == "Dormancy") )
Idents(object = singlets.combined.subG2MSwinners2) <- "seurat_clusters"
getDEGs2<-function(test, ref){
    markers <-  FindMarkers(singlets.combined.subG2MSwinners2, assay="RNA", ident.1 = test, ident.2 = ref, p_val_adj=0.001, only.pos = TRUE) 
    res <- enrichr(rownames(markers), "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020
    res$count <- str_split(res$Overlap, "/", simplify=TRUE)[,1]
    res[as.numeric(res$count) >= 3 & res$Adjusted.P.value<=0.05, ]
}

dormancyClusters<-c(0,3)
awakeningClusters<-c(4)

dormVsAw <- getDEGs2(dormancyClusters, awakeningClusters)
awVsdorm <- getDEGs2(awakeningClusters, dormancyClusters)

enrichedGeneSets<-list(dormVsAw, awVsdorm)
names(enrichedGeneSets)<-c("DormancyVsAwakening", "AwakeningVsDormancy")
enrichedGeneSets.m<-melt(enrichedGeneSets)
enrichedGeneSets.m<-enrichedGeneSets.m[enrichedGeneSets.m$variable=="Adjusted.P.value",]
# enrichedGeneSets.m[enrichedGeneSets.m$L1 %in% c("T0VsDorm", "awVsdorm","t0VsAw"),]$L1<-(enrichedGeneSets.m[enrichedGeneSets.m$L1 %in% c("T0VsDorm", "awVsdorm","t0VsAw"),]$L1)*-1
enrichedGeneSets.m$count<-as.numeric(enrichedGeneSets.m$count)
# enrichedGeneSets.m$L1<-factor(enrichedGeneSets.m$L1, levels=c("dormVsAw", "dormVsT0",   "awVsT0","awVsdorm","T0VsDorm", "t0VsAw"))
enrichedGeneSets.df<-acast(enrichedGeneSets.m, Term~L1)
enrichedGeneSets.df[is.na(enrichedGeneSets.df)]=0
enrichedGeneSets.df<-as.data.frame(enrichedGeneSets.df)
# enrichedGeneSets.m$Term<-factor(enrichedGeneSets.m$Term, levels= rownames(enrichedGeneSets.df[order(rowSums(enrichedGeneSets.df[,1:3]>0), rownames(enrichedGeneSets.df), decreasing=T),]))
enrichedGeneSets.m$Term<-factor(enrichedGeneSets.m$Term, levels= rownames(enrichedGeneSets.df[order(enrichedGeneSets.df[,1], decreasing=T),]))

enrichedGeneSets.m[enrichedGeneSets.m$L1=="DormancyVsAwakening",]$value <- (enrichedGeneSets.m[enrichedGeneSets.m$L1=="DormancyVsAwakening",]$value)*(-1)
enrichedGeneSets.m$Term <- factor(enrichedGeneSets.m$Term, levels=unique(enrichedGeneSets.m[order(enrichedGeneSets.m$value, decreasing=TRUE),]$Term))
enrichedGeneSets.m[enrichedGeneSets.m$value<0,]$L1 <- "Dormancy"
enrichedGeneSets.m[enrichedGeneSets.m$value>0,]$L1 <- "Awakening"
# enrichedGeneSets.m[enrichedGeneSets.m$Term=="Estrogen Response Early",]$L1 <- "Both"
enrichedGeneSets.m$L1<-factor(enrichedGeneSets.m$L1, levels=c("Dormancy","Awakening"))

pdf("~/manuscriptFeb23/plots/enrichedGeneSets.G2MScells.clWise.pdf", height=3, width=5)
ggplot(enrichedGeneSets.m, aes(x=value, y=Term, size=count, fill=L1)) +
    geom_point(shape=21, color="black") +
    theme_bw() +
    scale_fill_manual(values=c("Dormancy"="yellow", "Awakening"="#50157A")) +
    theme(legend.position="none", 
        legend.key.width = unit(0.8, 'cm'),
        # axis.text.x = element_text(angle = 45),
        legend.text = element_text(size=8)) +
    labs(x="Adjusted P.value", y="",fill="stage", size="Overlapping\nGenes")
    # +
    # facet_grid(L1~., space="free", scale="free") +
    # guides(size = guide_legend(nrow = 2),
    #     fill = guide_legend(nrow = 2),byrow=TRUE)
dev.off()

pdf("~/manuscriptFeb23/plots/enrichedGeneSets.G2MScells.clWise.legend.pdf")
ggplot(enrichedGeneSets.m, aes(x=value, y=Term, size=count, fill=L1)) +
    geom_point(shape=21, color="black") +
    theme_bw() +
    scale_fill_manual(values=c("Dormancy"="yellow", "Awakening"="#50157A")) +
    theme(legend.position="bottom", 
        legend.key.width = unit(0.8, 'cm'),
        # axis.text.x = element_text(angle = 45),
        legend.text = element_text(size=8)) +
    labs(x="Adjusted P.value", y="",fill="stage", size="Overlapping\nGenes")+
    # facet_grid(L1~., space="free", scale="free") +
    guides(size = guide_legend(nrow = 2),
        fill = guide_legend(nrow = 2),byrow=TRUE)
dev.off()



colors2 <- c("#B1961A", "#FFA300","#1F9E89", "#E6D800","#9B19F5", "#DC0AB4", "#E60049", "#B3D4FF", "grey")
names(colors2)<-0:8

samples<-c("T0_1", "T0_2", "TAM1m1", "TAM2m1", "AI1m1", "AI2m1", "TAM1m2", "TAM2m2", "AI1m2", "AI2m2", "TAMaw1","TAMaw2", "TAMaw3", "TAMaw4", "AIaw1","AIaw2", "AIaw3", "AIaw4")
DimPlot(singlets.combined.subG2MSwinners2, group.by="seurat_clusters", split.by="orig.ident",  reduction = "umap", ncol=4, raster="FALSE") + 
    ggplot2::labs(title = "") +
    scale_color_manual(values=colors2) +
    ggplot2::theme(legend.position = "bottom") 
dev.off()




