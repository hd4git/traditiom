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
require("sleuth")
require("rgl")
## Import samples
setwd("/rds/general/user/hdhiman/rnaseq/results/traditiom/pseudoAlignment")
samples<-list.files("pseudoAlignment")
samples<-samples[!samples %in% c("R51", "R52", "R53", "R54", "R55", "R56", "R57", "R58")]
info<-read.table("../src/sample_info.txt", 
                header=T, 
                stringsAsFactors=F,
                sep="\t")
info<-info[c(1:50, 59),]
str_info<-as.data.frame(str_split(info$Sample, "\\.", simplify=T), stringsAsFactors=F)
info$condition<-""
info[grep("T", info$Sample),]$condition<-"T"
info[grep("AI", info$Sample),]$condition<-"AI"
info[grep("^U", info$Sample),]$condition<-"U"
info[grep("POT", info$Sample),]$condition<-"POT"
info$time<-""
info[grep("POT", info$Sample),]$time<-"0"
info[grep("7", info$Sample),]$time<-"7"
info[grep("14", info$Sample),]$time<-"14"
info[grep("30", info$Sample),]$time<-"30"
info[grep("60", info$Sample),]$time<-"60"
info[grep("90", info$Sample),]$time<-"90"
info[grep("120", info$Sample),]$time<-"120"
# info[grep("150", info$Sample),]$time<-"150"
info[grep("170", info$Sample),]$time<-"170"
info[info$Sample %in% c("AI.alpha", "AI.delta"),]$time<-"120"
info[info$Sample %in% c("AI.beta", "AI.gamma", "AI.epsilon"),]$time<-"150"
info[info$Sample %in% c("T.theta", "T.eta"),]$time<-"90"
info[info$Sample %in% c("T.alpha", "T.beta", "T.gamma", "T.delta", "T.epsilon"),]$time<-"120"
info[info$Sample %in% c("T.theta", "T.eta"),]$time<-"90"
info[info$Sample %in% c("T.alpha", "T.beta", "T.gamma", "T.delta", "T.epsilon", "T.zeta"),]$time<-"120"
info[grep("TEP", info$Sample),]$time<-"180"
info$state<-""
info[grep('POT', info$Sample),]$state<-"POT"
info[grep('U', info$Sample),]$state<-"Untreated"
info[grep('7\\.|14\\.', info$Sample),]$state<-"Latency"
info[grep('T30\\.|AI30\\.|T60\\.|AI60\\.|AI90', info$Sample),]$state<-"Dormancy"
info[grep('AI\\.|T\\.', info$Sample),]$state<-"Awakening"

info$sampleTime<-str_info[,1]
info$rep<-str_info[,2]
info.m<-info
saveRDS(info, "rds/sample_details.rds")

###########################################################
############# Import counts with tximport #################
###########################################################
files<-as.data.frame(
	do.call(rbind, lapply(
		samples, function(x){paste(x, "abundance.tsv", sep="/")
	})), stringsAsFactors=F)$V1


## Creating tx2gene 
edb <- EnsDb.Hsapiens.v86
Tx <- transcripts(edb)
tx2gene<-genes(edb,
          columns = c("tx_id", "gene_name"),
          return.type = "DataFrame")
saveRDS(tx2gene, "rds/tx2gene.rds")

## Loading info
setwd("pseudoAlignment")
txi.kallisto.tsv <- tximport(files, type = "kallisto", 
                                    tx2gene = tx2gene, 
                                    ignoreAfterBar = TRUE, 
                                    ignoreTxVersion=T)
colnames(txi.kallisto.tsv$counts)<-samples
saveRDS(txi.kallisto.tsv, "../rds/txi.rds")

txi.kallisto.tsv<-readRDS("rds/txi.rds")
## Importing info to deseq dataset
rownames(info.m) <- info.m$ID
info.m<-info.m[samples,]
info.m$condition<-factor(info.m$condition, levels=c("POT", "UT", "T", "AI"))
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, info.m,  ~condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
#Filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds.filter<- dds[keep,]
dds<-dds.filter


dds<-readRDS("~/rnaseq/results/traditiom/rds/dds.condition.rds")
vsd <- vst(dds)


#### Plot PCA for top 500 variant genes
pcaData <- plotPCA(vsd, intgroup = c("ID", "condition", "Sample", "time", "rep", "state"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

scheme<-c("POT"="#BFBFBF", 
          "Untreated"="#000000", 
          "Latency"="#0084A5", 
          "Dormancy"="#FFDC00", 
          "Awakening"="#B30078", 
          "TEP"="#50157A")
scheme.shape<-c("POT"=15, "UT"=16, "AI"=17, "T"=18)
pcaData[grep("UT30", pcaData$Sample),]$state<-"Untreated"
pcaData[grep("TEP", pcaData$Sample),]$state<-"TEP"
pcaData[grep("POT", pcaData$Sample),]$state<-"POT"

greek<-pcaData[grep("AI\\.|T\\.", pcaData$Sample),]$Sample
greek<-greek[-grep('POT', greek)]
relabel<-read.table("../src/relabel.txt", sep="\t",header=F, stringsAsFactors=F)

pcaData$greek<-pcaData$Sample
pcaData[pcaData$Sample %in% greek,]$greek<-relabel$V2


pcaData$state<-factor(pcaData$state, levels=c("POT", 
                                              "Untreated", 
                                              "Latency", 
                                              "Dormancy", 
                                              "Awakening", 
                                              "TEP"))

pcaData2<-na.omit(pcaData)
pcaData2<-pcaData2[-grep("UT|POT", pcaData2$Sample),]

pdf("~/rnaseq/results/traditiom/plots/fgure5a.pdf", height=3, width=2)
ggplot(na.omit(pcaData), 
      aes(y = PC1, 
          x = PC2, 
          shape = condition,  
          label=greek)) +
      geom_point(size = 2, aes(color = state)) +  
      # geom_text_repel(data=na.omit(pcaData), 
      #               aes(label=greek), parse=T,  
      #               color="#000000", 
      #               size=1.5,
      #               max.overlaps=100) +
      labs(y=paste0("PC1: ", percentVar[1], "% variance"), 
          x=paste0("PC2: ", percentVar[2], "% variance"),
          color="State", shape="Condition") +
      theme_bw() +
      scale_shape_manual(values=scheme.shape) +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      theme(legend.position="bottom", 
        text = element_text(size=10),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size=3),
        legend.title = element_text(size=3)) +
      guides(shape=guide_legend(nrow=2,byrow=TRUE))
 dev.off() 

order.samples<-c("POT.1", "POT.2", "POT.3", "U30.1", "U30.2", "U120.1", "U120.2", "U170.1", "U170.2", "U170.3",
  "T7.1", "T7.2", "T14.1", "T14.2", "AI7.1", "AI7.2", "AI14.1", "AI14.2",  "T30.1", "T30.2", "T60.1", "T60.2", "AI30.1", "AI30.2", "AI60.1", "AI60.2", "AI90.1", "AI90.2",  
  "T.alpha",  "T.beta", "T.gamma", "T.delta",  "T.epsilon", "T.zeta", "T.theta", "T.eta", "AI.alpha", "AI.beta",  "AI.gamma", "AI.delta", "AI.epsilon", 
  "T.alpha_TEP", "T.beta_TEP","T.gamma_TEP", "T.delta_TEP","T.epsilon_TEP", "T.zeta_TEP", "AI.alpha_TEP", "AI.beta_TEP","AI.gamma_TEP", "AI.epsilon_TEP")
vsd.counts<-assay(vsd)
colnames(vsd.counts)<-info.m$Sample
vsd.counts<-vsd.counts[,order.samples]
write.table(vsd.counts, "rds/vsd.counts.txt", sep="\t", quote=F, row.names=T, col.names=T)

dds.counts<-assay(dds)
colnames(dds.counts)<-info.m$Sample
dds.counts<-dds.counts[,order.samples]
write.table(dds.counts, "rds/dds.counts.txt", sep="\t", quote=F, row.names=T, col.names=T)






     
ggplot(na.omit(pcaData2), 
      aes(x = PC1, 
          y = PC2, 
          shape = condition,  
          label=greek)) +
      geom_point(size=3, aes(color = state)) +   
      geom_text_repel(data=na.omit(pcaData2), 
                      aes(label=greek), parse=T,  
                      color="#000000", 
                      size=4, 
                      alpha=1, 
                      max.overlaps=500) +
      labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
          y=paste0("PC2: ", percentVar[2], "% variance"),
          color="State") +
      theme_bw() +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      scale_shape_manual(values=scheme.shape) +
      theme(legend.position="bottom") 
ggplot(na.omit(pcaData), 
      aes(x = PC1, 
          y = PC2, 
          shape = condition,  
          label=greek)) +
      geom_point(size = 3, aes(color = state)) +  
      geom_text_repel(data=na.omit(pcaData), 
                    aes(label=greek), parse=T,  
                    color="#000000", 
                    size=4, 
                    alpha=1, 
                    max.overlaps=100) +
      labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
          y=paste0("PC2: ", percentVar[2], "% variance"),
          color="State") +
      theme_bw() +
      scale_shape_manual(values=scheme.shape) +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      theme(legend.position="bottom",
          axis.title=element_blank(),
          axis.text=element_blank()) 
ggplot(na.omit(pcaData2), 
      aes(x = PC1, 
          y = PC2, 
          shape = condition,  
          label=greek)) +
      geom_point(size=3, aes(color = state)) +   
      geom_text_repel(data=na.omit(pcaData2), 
                      aes(label=greek), parse=T,  
                      color="#000000", 
                      size=4, 
                      alpha=1, 
                      max.overlaps=500) +
      labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
          y=paste0("PC2: ", percentVar[2], "% variance"),
          color="State") +
      theme_bw() +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      scale_shape_manual(values=scheme.shape) +
      theme(legend.position="bottom",
          axis.title=element_blank(),
          axis.text=element_blank()) 
dev.off()

############################################################


pdf("../plots/pca-vsd.traditiom.pdf")
ggplot(na.omit(pcaData), aes(x = PC1, 
					y = PC2, 
					# fill = state,
          color=state,
          shape = condition,  
					label=Sample)) +
  geom_point(size = 3) +  
  # scale_colour_manual(values=cond_colours) + 
  geom_text_repel(data=na.omit(pcaData), aes(label=Sample), parse=T,  color="#000000", size=2.5, alpha=1, max.overlaps=100) +
  stat_ellipse(geom="polygon", aes(fill = state),
                # color="NA", 
                 alpha = 0.2, 
                 show.legend = FALSE, 
                 level = 0.95) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
  		y=paste0("PC2: ", percentVar[2], "% variance"),
  		color="State") +
      theme_bw() +
      scale_shape_manual(values=scheme.shape) +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
  	theme(legend.position="bottom") 

ggplot(na.omit(pcaData), aes(x = PC1, 
          y = PC2, 
          # fill = state,
          color=state,
          # shape = condition,  
          label=Sample)) +
  geom_point(size = 3) +  
  # scale_colour_manual(values=cond_colours) + 
  geom_text_repel(data=na.omit(pcaData), aes(label=Sample), parse=T,  color="#000000", size=2.5, alpha=1, max.overlaps=100) +
  stat_ellipse(geom="polygon", aes(fill = state),
                # color="NA", 
                 alpha = 0.2, 
                 show.legend = FALSE, 
                 level = 0.95) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
      y=paste0("PC2: ", percentVar[2], "% variance"),
      color="State") +
      theme_bw() +
      scale_shape_manual(values=scheme.shape) +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
    theme(legend.position="bottom") 



dev.off()


pdf("../plots/pca-vsd.traditiom.sub.pdf")
ggplot(na.omit(pcaData2), aes(x = PC1, 
          y = PC2, 
          color = state, 
          shape = condition,  
          label=Sample)) +
  geom_point(size=3) +  
  stat_ellipse(geom="polygon", aes(fill = state),
                # color="NA", 
                 alpha = 0.2, 
                 show.legend = FALSE, 
                 level = 0.95) +  
  geom_text_repel(data=na.omit(pcaData2), aes(label=Sample), parse=T,  color="#000000", size=2.5, alpha=1, max.overlaps=500) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
      y=paste0("PC2: ", percentVar[2], "% variance"),
      color="State") +
      theme_bw() +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      scale_shape_manual(values=scheme.shape) +
    theme(legend.position="bottom") 

ggplot(na.omit(pcaData2), aes(x = PC1, 
          y = PC2, 
          color = state, 
          # shape = condition,  
          label=Sample)) +
  geom_point(size=3) +  
  stat_ellipse(geom="polygon", aes(fill = state),
                # color="NA", 
                 alpha = 0.2, 
                 show.legend = FALSE, 
                 level = 0.95) +  
  geom_text_repel(data=na.omit(pcaData2), aes(label=Sample), parse=T,  color="#000000", size=2.5, alpha=1, max.overlaps=500) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
      y=paste0("PC2: ", percentVar[2], "% variance"),
      color="State") +
      theme_bw() +
      scale_color_manual(values=scheme) +
      scale_fill_manual(values=scheme) +
      scale_shape_manual(values=scheme.shape) +
    theme(legend.position="bottom") 
dev.off()

###############################################################################################
## Get estimated counts per transcript
counts.list<-lapply(samples, function(x){read.table(paste(x, "abundance.tsv", sep="/"), header=T, stringsAsFactors=F)})
names(counts.list)<-samples

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ target_id, value.var="est_counts"))
counts[is.na(counts)] <- 0
colnames(counts)[26]<-"R378B" 			### Since sample R378B was named R378II earlier
colnames(counts)<-info.m[match(colnames(counts), info.m$value),]$sample
counts<-counts[,c(paste("pre", 1:23, sep="."), paste("post", 1:23, sep="."))]
saveRDS(counts, "rds/est_counts.rds")
counts<-readRDS("rds/est_counts.rds")

