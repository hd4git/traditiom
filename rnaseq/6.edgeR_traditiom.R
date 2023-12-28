require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")
require("tximport")
require("edgeR")
require("ggrepel")
require("pheatmap")
require("viridis")
require("org.Hs.eg.db")
require("clusterProfiler")
require("gprofiler2")

### DE of TEP vs awakening ###
samples<-list.files("pseudoAlignment")
samples<-samples[!samples %in% c("R51", "R52", "R53", "R54", "R55", "R56", "R57", "R58")]

info.m<-readRDS("rds/sample_details.rds")
info.m<-info.m[samples,]

dds<-readRDS("rds/dds.sampleTime.rds")

counts<-assays(dds)[["counts"]]
colnames(counts)<-info.m$Sample
counts.aw<-counts[,grep("AI\\.|^T\\.", colnames(counts))]
sel.order<-c("T.alpha","T.alpha_TEP", "T.beta","T.beta_TEP", "T.gamma","T.gamma_TEP", "T.delta","T.delta_TEP", "T.epsilon", "T.epsilon_TEP", "T.zeta", "T.zeta_TEP", "T.eta",  "T.theta", "AI.alpha", "AI.alpha_TEP", "AI.beta", "AI.beta_TEP",  "AI.gamma", "AI.gamma_TEP", "AI.delta", "AI.epsilon")
counts.aw<-counts.aw[,sel.order]
group<-factor(1:22)
dgList <- DGEList(counts=counts.aw, genes=rownames(counts.aw), group=group)
countsPerMillion <- cpm(dgList)


############################

runEdgeR<-function(x){
	counts.aw<-counts[,grep(x, colnames(counts))]
	group<-factor(1:2)
	dgList <- DGEList(counts=counts.aw, genes=rownames(counts.aw), group=1:2)
	countsPerMillion <- cpm(dgList, normalized.lib.sizes=FALSE)

	keep <- filterByExpr(dgList)
	dgList <- dgList[keep, , keep.lib.sizes=FALSE]
	dgList <- calcNormFactors(dgList, method="TMM")

	design <- model.matrix(~group)
	# dgList <- estimateDisp(dgList, design)
	dgList <- estimateGLMCommonDisp(dgList, method="deviance", robust=TRUE, subset=NULL)
	# et <- exactTest(dgList, dispersion=0.2^2)
	fit <- glmFit(dgList, design)
	lrt <- glmLRT(fit)

	res<-lrt$table
	res$Name<-rownames(res)
	res$fcSign<-sign(res$logFC)
	res$logP<--log10(res$PValue)
	res$metric<-res$logP/res$fcSign
	rankFile<-res[,c("Name", "metric")]
	rankFile<-rankFile[order(rankFile$metric, decreasing=T),]
	outFile<-paste(x, "edger.rnk", sep="_")
	outFile2<-paste(x, "edger_sig.txt", sep="_")
	outFile3<-paste(x, "edger_go_gp.txt", sep="_")
	outFile4<-paste(x, "edger_go_cp.txt", sep="_")
	write.table(rankFile,file=paste("reports_edger", outFile, sep="/"), quote=F, sep="\t", row.names=F)

	sig<-glmTreat(fit, coef=2, lfc = 1)
	sig<-topTags(lrt,  adjust.method = "BH", sort.by = "PValue", p.value = 0.001, n=dim(lrt)[1])
	write.table(sig,file=paste("reports_edger", outFile2, sep="/"), quote=F, sep="\t", row.names=F)

	gostres <- gost(query = rownames(sig), 
                organism = "hsapiens", significant = TRUE,  
                user_threshold = 0.005, correction_method = "g_SCS", sources = c("GO:BP"))
	go<-as.data.frame(gostres$result)
	write.table(as.matrix(go), file=paste("reports_edger", outFile3, sep="/"), quote=F, sep="\t", row.names=F)

	ggo <- groupGO(gene     = rownames(sig),
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

}

runEdgeR("AI.alpha")
runEdgeR("AI.beta")
runEdgeR("AI.gamma")

runEdgeR("T.alpha")
runEdgeR("T.beta")
runEdgeR("T.gamma")
runEdgeR("T.delta")
runEdgeR("T.epsilon")
runEdgeR("T.zeta")



