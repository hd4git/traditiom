require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")
require("tximport")
require("edgeR")
require("rhdf5")
require("ensembldb")
require("EnsDb.Hsapiens.v86")
require("ggrepel")
require("pheatmap")

samples<-list.files("pseudoAlignment")
samples<-samples[!samples %in% c("R51", "R52", "R53", "R54", "R55", "R56", "R57", "R58")]

info.m<-readRDS("rds/sample_details.rds")
info.m<-info.m[samples,]

tx2gene<-readRDS("rds/tx2gene.rds")
txi.kallisto.tsv<-readRDS("rds/txi.rds")


info.m$sampleTime <- factor(info.m$sampleTime, levels = c("POT","UT30", "UT120", "UT170", "T7", "T14", "T30", "T60", "T.alpha", "T.beta", "T.gamma", "T.delta", "T.epsilon", "T.zeta", "T.theta", "T.eta", "T.alpha_TEP", "T.beta_TEP", "T.gamma_TEP", "T.delta_TEP", "T.epsilon_TEP", "T.zeta_TEP", "AI7", "AI14", "AI30", "AI60", "AI90", "AI.alpha", "AI.beta", "AI.gamma", "AI.delta", "AI.epsilon", "AI.alpha_TEP", "AI.beta_TEP", "AI.gamma_TEP", "AI.delta_TEP", "AI.epsilon_TEP"))
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, info.m, ~sampleTime)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
#Filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds.filter<- dds[keep,]
dds<-dds.filter

saveRDS(dds, "rds/dds.sampleTime.rds")
saveRDS(dds, "rds/dds.condition.rds")
######################################################

dds<-readRDS("rds/dds.sampleTime.rds")
##############################
######## DE host vs exp ###
##############################
# dds$type <- factor(dds$type, levels = c("host", "exp"))
design(dds) <- formula(~ sampleTime)

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
getDE<-function(contrast, test, control)
{
	de_results <- lfcShrink(DDS_batch, coef=contrast, type="apeglm")
	summary(de_results)

	#Export sig DE list
	# sig_de_results <- subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
	sig_de_results <-subset(de_results,  abs(log2FoldChange)>1 & padj < 0.01)
	sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
	summary(sig_de_results)

	fileName<-paste("reports2", contrast, sep="/")
	sig_de_results<-read.table(paste(fileName, "txt", sep="."), sep="\t", header=T, stringsAsFactors=F)

	write.table(as.data.frame(sig_de_results),file =paste(fileName, "txt", sep="."),sep="\t",quote=F,row.names=T,col.names=T)

	### Using contrast for getting ranks for GSEA 
	de_results<-results(DDS_batch, 
	                    contrast = c("sampleTime",test,control),
	                    lfcThreshold=0, 
	                    independentFiltering =T, 
	                    pAdjustMethod="BH", 
	                    alpha = 0.01)
	#Export sig DE list
	summary(de_results)

	sig_de_results <-subset(de_results,  abs(log2FoldChange)> 1 & padj < 0.01)
	sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
	summary(sig_de_results)

	res<-as.data.frame(de_results)
	res<-res[order(res$stat,decreasing=TRUE),]
	res2<-subset(res, select="stat")
	rownames(res2)<-toupper(rownames(res2))
	write.table(res2,file =paste(fileName, "rnk", sep="."),sep="\t",quote=F,row.names=T,col.names=F)
}

getDE("sampleTime_UT30_vs_POT", "UT30", "POT")
getDE("sampleTime_UT120_vs_POT", "UT120", "POT")
getDE("sampleTime_UT170_vs_POT", "UT170", "POT")

getDE("sampleTime_T30_vs_POT", "T30", "POT")
getDE("sampleTime_T60_vs_POT", "T60", "POT")

getDE("sampleTime_AI30_vs_POT", "AI30", "POT")
getDE("sampleTime_AI60_vs_POT", "AI60", "POT")
getDE("sampleTime_AI90_vs_POT", "AI90", "POT")









