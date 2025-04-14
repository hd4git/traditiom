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
require("enrichR")
require("viridis")
require("enrichR")
require("ggrepel") 
require("pheatmap") 

setwd("~/spatialTranscriptomics/")

sigDEGs <- readRDS("rds/sigDEGs.RDS") 
target_patientData <- readRDS("rds/target_patientData_all.rds")
# write.table(melt(sigDEGs[[1]], id.vars=colnames(sigDEGs[[1]][[1]])), "rds/sigDEGsPOS.txt", sep="\t",quote=F,row.names=F,col.names=T)
# write.table(melt(sigDEGs[[2]], id.vars=colnames(sigDEGs[[2]][[1]])), "rds/sigDEGsNEG.txt", sep="\t",quote=F,row.names=F,col.names=T)

dbs <- listEnrichrDbs()
selDB.TF<-c("ChEA_2022", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TF_Perturbations_Followed_by_Expression", "TRANSFAC_and_JASPAR_PWMs", "Epigenomics_Roadmap_HM_ChIP-seq", "TF-LOF_Expression_from_GEO", "ENCODE_Histone_Modifications_2015", "Transcription_Factor_PPIs")
selDB.PW <- c("Reactome_2022", "BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "NCI-Nature_2015", "Panther_2015", "Kinase_Perturbations_from_GEO_down", "Kinase_Perturbations_from_GEO_up", "GO_Molecular_Function_2021")

sigDEGsPosEnriched <- lapply(sigDEGs[[1]], function(x){
                            info <- lapply(c(selDB.TF, selDB.PW), function(y){
                                res <- enrichr(x$Gene, y)
                                if(dim(res[[y]])[1]>0)
                                {   res <- res[[y]]
                                    res$count <- str_split(res$Overlap, "/", simplify=TRUE)[,1]
                                    res <- res[as.numeric(res$count) >= 3 & res$Adjusted.P.value<=0.05, ]
                                    if(dim(res)[1]>0){res$libraryName <- y}
                                    res
                                    }})
                            do.call(rbind, info)
                            })
compPos <- c("S1R_S1R:S1L_CK" ,"R1R_R1R:S1R_CK" ,"R1R_R1R:S1L_CK" ,"S2_S2:D2_CK" ,"S3_S3:D3_CK" ,"R1R_R1R:S1L_CD" ,"S1R_S1R:S1L_Rest" ,"R1R_R1R:S1L_Rest")
sigDEGsPosEnriched2 <- lapply(1:length(names(sigDEGsPosEnriched)), function(x){
        if(dim(sigDEGsPosEnriched[[x]])[1] >0)
        {sigDEGsPosEnriched[[x]]$comp <- compPos[x]
            }
            sigDEGsPosEnriched[[x]]})
sigDEGsPosEnriched.df <- do.call(rbind, sigDEGsPosEnriched2)


sigDEGsNegEnriched <- lapply(sigDEGs[[2]], function(x){
                            info <- lapply(c(selDB.TF, selDB.PW), function(y){
                                res <- enrichr(x$Gene, y)
                                if(dim(res[[y]])[1]>0)
                                {   res <- res[[y]]
                                    res$count <- str_split(res$Overlap, "/", simplify=TRUE)[,1]
                                    res <- res[as.numeric(res$count) >= 3 & res$Adjusted.P.value<=0.05, ]
                                    if(dim(res)[1]>0){res$libraryName <- y}
                                    res
                                    }})
                            do.call(rbind, info)
                            })
compNeg <- c("S1L_S1R:S1L_CK" ,"S1R_R1R:S1R_CK" ,"S1L_R1R:S1L_CK" ,"D2_S2:D2_CK" ,"D3_S3:D3_CK" ,"S1L_R1R:S1L_CD" ,"S1L_S1R:S1L_Rest" ,"S1L_R1R:S1L_Rest")
sigDEGsNegEnriched2 <- lapply(1:length(names(sigDEGsNegEnriched)), function(x){
        if(dim(sigDEGsNegEnriched[[x]])[1] >0)
        {sigDEGsNegEnriched[[x]]$comp <- compNeg[x]
            }
            sigDEGsNegEnriched[[x]]})
sigDEGsNegEnriched.df <- do.call(rbind, sigDEGsNegEnriched2)


sigDEGsPosEnriched.df$factor = 1
sigDEGsNegEnriched.df$factor = -1
sigDEGsEnriched.all <- rbind(sigDEGsPosEnriched.df, sigDEGsNegEnriched.df)
sigDEGsEnriched.all$type <- ""
sigDEGsEnriched.all[sigDEGsEnriched.all$libraryName %in% selDB.TF, ]$type <- "TFs"
sigDEGsEnriched.all[sigDEGsEnriched.all$libraryName %in% selDB.PW, ]$type <- "Pathways"
mainCompInfo <- str_split(sigDEGsEnriched.all$comp, "_", simplify=T)[,2:3]
mainCompInfo <- paste(mainCompInfo[,1], mainCompInfo[,2], sep="_")
sigDEGsEnriched.all$mainComp <- mainCompInfo

sigDEGsEnriched.all <- sigDEGsEnriched.all[sigDEGsEnriched.all$type == "Pathways" | (sigDEGsEnriched.all$type == "TFs" & grepl("Human", sigDEGsEnriched.all$Term)),]
overlap <- str_split(sigDEGsEnriched.all$Overlap, "/", simplify=TRUE)

sigDEGsEnriched.all$countRatio <- -log10(as.numeric(overlap[,1])/as.numeric(overlap[,2]))
sigDEGsEnriched.all$Adjusted.P.value2 <- sigDEGsEnriched.all$Adjusted.P.value * sigDEGsEnriched.all$factor
sigDEGsEnriched.all$Term2 <- paste(sigDEGsEnriched.all$libraryName, sigDEGsEnriched.all$Term, sep=":")
sigDEGsEnriched.all.list <- split(sigDEGsEnriched.all, sigDEGsEnriched.all$mainComp)

saveRDS(sigDEGsEnriched.all.list, "rds/sigDEGsEnriched.all.list.rds")
sigDEGsEnriched.all.list<-readRDS("rds/sigDEGsEnriched.all.list.rds")
# pdf("plots/enrichment.pdf", height=10, width=20)
# lapply(sigDEGsEnriched.all.list, function(x){
scheme <- c("-1"="purple", "1"="green")
plotEnrichment <- function(x,y){
    overlap <- str_split(x$Overlap, "/", simplify=TRUE)
    x$countRatio <- as.numeric(overlap[,1])/as.numeric(overlap[,2])
    x$Term2 <- str_split(x$Term2, ":", simplify=TRUE)[,2]
    x$Term2<-factor(x$Term2, levels=unique(x[order(x$Adjusted.P.value2,decreasing=T),]$Term2))
    x$factor <- as.factor(x$factor)
    # labelInfo <- unique(x$mainComp)
    # labelInfo <- str_split(labelInfo, ":|_", simplify=T)
    ggplot(x, aes(x=Adjusted.P.value2, y=Term2, size=countRatio, color=factor)) +
        geom_point() +
        theme_bw() +
        # facet_grid(type~., space="free", scales="free") +
        xlim(-0.05,0.05) +
        # scale_color_viridis(limits=c(-0.05,0.05)) +
        scale_color_manual(values=y, breaks=waiver()) +
        labs(
            # title=x$mainComp, 
            # x=paste(labelInfo[2], "Adjusted.P.value", labelInfo[1], sep="                                "), 
            y="", size="OverlapRatio", color="Comparison") +
        theme(legend.position="bottom", 
            strip.background = element_rect(fill="white"),
            legend.key.width = unit(0.5, 'cm')) +
        guides(size = guide_legend(nrow = 3),
            color = guide_legend(nrow = 2),byrow=TRUE)
    }
# scheme1<-c("0"="black", "diagnostic_A.L"="#0084A5", 
#     "surgery_A.L"="#660044", "surgery_A.R"="#B30078", 
#     "relapse_A.R"="#7D21C0", 
#     "diagnostic_B"="#00C4F5", "surgery_B"="#FF0AAD", 
#     "diagnostic_C"="#85E7FF", "surgery_C"="#FF70CF")
# scheme2<-c("CK+"="red", "CD45+"="green", "Rest"="grey", "0"="black")


pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[5], "_MSigDB.pdf", sep=""), height=3, width=7)
plotEnrichment(sigDEGsEnriched.all.list[[5]][sigDEGsEnriched.all.list[[5]]$libraryName=="MSigDB_Hallmark_2020",], c("-1"="#660044", "1"="#B30078"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[7], "_MSigDB.pdf", sep=""), height=3, width=7)
plotEnrichment(sigDEGsEnriched.all.list[[7]][sigDEGsEnriched.all.list[[7]]$libraryName=="MSigDB_Hallmark_2020",], c("-1"="#660044", "1"="#B30078"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[8], "_MSigDB.pdf", sep=""), height=3, width=7)
plotEnrichment(sigDEGsEnriched.all.list[[8]][sigDEGsEnriched.all.list[[8]]$libraryName=="MSigDB_Hallmark_2020",], c("-1"="#660044", "1"="#B30078"))    
dev.off()    




plotEnrichment(sigDEGsEnriched.all[sigDEGsEnriched.all$libraryName == "MSigDB_Hallmark_2020",], c("-1"="#B30078", "1"="#7D21C0"))    
dev.off()    


pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[1], ".pdf", sep=""), height=7, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[1]], c("-1"="#B30078", "1"="#7D21C0"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[2], ".pdf", sep=""), height=20, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[2]], c("-1"="#B30078", "1"="#7D21C0"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[3], ".pdf", sep=""), height=5, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[3]], c("-1"="#B30078", "1"="#7D21C0"))    
dev.off()    


pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[4], ".pdf", sep=""), height=5, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[4]], c("-1"="#B30078", "1"="#7D21C0"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[5], ".pdf", sep=""), height=7, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[5]], c("-1"="#660044", "1"="#B30078"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[6], ".pdf", sep=""), height=10, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[6]], c("-1"="#660044", "1"="#B30078"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[7], ".pdf", sep=""), height=10, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[7]], c("-1"="#00C4F5", "1"="#FF0AAD"))    
dev.off()    

pdf(paste("plots/DE/",names(sigDEGsEnriched.all.list)[8], ".pdf", sep=""), height=10, width=15)
plotEnrichment(sigDEGsEnriched.all.list[[8]], c("-1"="#85E7FF", "1"="#FF70CF"))    
dev.off()    


