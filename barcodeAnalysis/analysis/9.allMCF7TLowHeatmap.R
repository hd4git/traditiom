require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("ggpubr")
require("viridis")
require("pheatmap")
require("matrixStats")

####################################
######## Genomic barcodes ##########
####################################
## Import samples
setwd("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered")
info<-read.table("rds/sampleInfo_allMCF7TLow.txt", header=T)
rownames(info)<-info[,2]
samples<-list.files("allMCF7TLow")
counts.list<-lapply(samples, function(x){read.table(paste("allMCF7TLow",x, sep="/"), header=F)})
ids<-gsub("bc_mat_|\\.txt|\\.counts", "", samples)
names(counts.list)<-ids

# for(i in names(counts.list)) {counts.list[[i]]$id<-i}
# counts.list2<-do.call(rbind, counts.list)
counts.list2<-melt(counts.list)
counts<-t(acast(counts.list2, L1 ~ V1, value.var="value"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts)
counts<-counts[,info$sampleID]
colnames(counts)<-info$sampleName

countsG<-counts[,grepl(".g", colnames(counts))]
countsSC<-counts[,grepl(".sc", colnames(counts))]

countsG<-countsG[rowSums(countsG>10)>0,]
countsSC<-countsSC[rowSums(countsG>3)>0,]

countsNew<-merge(countsG, countsSC, by=0, all=TRUE)
countsNew[is.na(countsNew)] = as.numeric(0) 
counts<-countsNew[,-1]
rownames(counts)<-countsNew[,1]
counts<-counts[rowSums(counts>0)>0,]
counts<-counts[,rownames(info)]

total<-colSums(counts)
freq<-t(t(counts)/total)
saveRDS(freq, "rds/freqMCF7TLow.rds")
freq <- readRDS("rds/freqMCF7TLow.rds")
#############################
##### Survival dynamics #####
#############################
surviving <- data.frame("sample"=names(colSums(counts>0)),
						"counts"=as.data.frame(colSums(counts>0))[,1])
surviving$type <- str_split(surviving$sample, "\\.", simplify=T)[,2]
surviving$sample2 <- str_split(surviving$sample, "\\.", simplify=T)[,1]
surviving$condition <- ""
surviving[grepl("POT", surviving$sample),]$condition <- "POT"
surviving[grepl("T0", surviving$sample),]$condition <- "T0"
surviving[grepl("m1|m2", surviving$sample),]$condition <- "Dormancy"
surviving[grepl("aw", surviving$sample),]$condition <- "Awakening"
surviving[grepl("TEP", surviving$sample),]$condition <- "Relapse"

surviving$treatment <- ""
surviving[grepl("POT", surviving$sample),]$treatment <- "POT"
surviving[grepl("TAM", surviving$sample),]$treatment <- "TAM"
surviving[grepl("AI", surviving$sample),]$treatment <- "-E2"
surviving[grepl("T0", surviving$sample),]$treatment <- "POT"

surviving$treatment <- factor(surviving$treatment, levels=c("POT", "TAM", "-E2"))
surviving$sample2 <- factor(surviving$sample2, levels=c(
						paste("POT", 1:3, sep="_"), 
						paste("T0", 1:2, sep="_"), 
						paste("T0", paste("TAM", 1:3, sep=""), sep="_"), 
						paste("T0", paste("AI", 1:3, sep=""), sep="_"),
						"TAM1m1", "TAM2m1", "TAM1m2", "TAM2m2", 
						"AI1m1", "AI2m1", "AI1m2", "AI2m2",
						paste("TAMaw", 1:10, sep=""), 
						paste(paste("TAMaw", 5:10, sep=""), "TEP", sep="_"),
						paste("AIaw", 1:10, sep=""), 
						paste(paste("AIaw", 5:10, sep=""), "TEP", sep="_")))

extras<-surviving[surviving$condition %in% c("Awakening", "Relapse") & grepl("5|6|7|8|9|10", surviving$sample),]$sample
surviving4AWs<-surviving[!surviving$sample %in% extras,]

scheme <- c("POT"="grey","T0"="grey", "Untreated"="#000000", "Dormancy"="#FFDC00", "Awakening"="#B30078", "Relapse"="#50157A")
pdf("~/manuscript/plots/survivingBCsMCF7low4AWs.pdf", height=10, width =10)
ggplot(surviving4AWs[surviving4AWs$type=="sc",], aes(x=sample2, y=counts, fill=condition)) + 
geom_bar(stat="identity", position="dodge2", color="black") + 
geom_text(stat="identity", aes(label=counts), position=position_stack(vjust = 0.94), vjust = -1) +
theme_bw() +  
labs(x=" ", y="Number of barcodes", fill="Condition") + 
theme(text = element_text(size=10),
	axis.text = element_text(size=10),
	axis.text.x = element_text(size=10, angle = 90),  
	legend.position="none", strip.background=element_rect(fill="white")) +
scale_fill_manual(values=scheme) +
facet_wrap(treatment ~ ., ncol=1, scales="free_x") 
ggplot(surviving4AWs[surviving4AWs$type=="g",], aes(x=sample2, y=counts, fill=condition)) + 
geom_bar(stat="identity", position="dodge2", color="black") + 
geom_text(stat="identity", aes(label=counts), position=position_stack(vjust = 0.94), vjust = -1) +
theme_bw() +  
labs(x=" ", y="Number of barcodes", fill="Condition") + 
theme(text = element_text(size=10),
	axis.text = element_text(size=10),
	axis.text.x = element_text(size=10, angle = 90),  
	legend.position="none", strip.background=element_rect(fill="white")) +
scale_fill_manual(values=scheme) +
facet_wrap(treatment ~ ., ncol=1, scales="free_x") 
dev.off()

###################
##### Heatmap #####
###################
freqGT5pct<-freq[rowSums(freq>0.05)>0,]
freqGT5pct[freqGT5pct==0]<-NA
freqGT5pct<-freqGT5pct[order(rowMaxs(freqGT5pct[,1:3]), decreasing=T),]
winners<-c("bc_8749691", "bc_9683153", "bc_2065240", "bc_1758787", "bc_150411", "bc_5189530", "bc_7099037", "bc_5112179", "bc_7305702", "bc_7704123","bc_6864523","bc_2127641","bc_5122096", "bc_9625985", "bc_159071")
selWinners <- c("TAMaw1"="bc_150411", "TAMaw2"="bc_5189530", "TAMaw3"="bc_7099037", "TAMaw4"="bc_5112179",
	"AIaw1"="bc_7704123","AIaw2"="bc_6864523","AIaw3"="bc_2127641","AIaw4"="bc_5122096")

freqGT5pct <- freqGT5pct[winners,]

info$condition<-factor(info$condition, levels=c("POT", "T0", "TAM", "AI"))
info$flask<-factor(info$flask, levels=1:10)
info$sequencing<-factor(info$sequencing, levels=c("g", "sc"))
info$phase<-factor(info$phase, levels=c("T0", "Dormancy", "Awakening"))

annoColor<-list(
	# condition=c("POT"="#6E6F72", "T0"="#ACADAF", "TAM"="#008FCC", "AI"="#6CD4FF"), 
	phase=c("T0"="#ACADAF","Dormancy"="#FFDC00", "Awakening"="#B30078"), 
	sequencing=c("g"="#144ECB", "sc"="#F49FBC"))
pdf("~/manuscript/plots/heatmap_allMCF7TLow5pct.pdf", height=5, width=12)	
pheatmap(freqGT5pct, 
	color = viridis(100, begin = 0.2, end = 1, direction = -1, option = "B"),
	na_col = "white",
	annotation_colors = annoColor,
	border_color = NA,
	angle_col = 90,
	gaps_col = c(3, 5,8,12,14,16,18,20,22,24,26,28,30, 32, 35,39,41,43,44,46,48,50,51,53,54),
	annotation_col = info[,c(5,3,6)], 
	cluster_rows=F, 
	cluster_cols=F)
dev.off()

selWinners2 <- c("TAMaw1"="bc_150411", 
				"TAMaw2"="bc_5189530", 
				# "TAMaw3"="bc_7099037", 
				# "TAMaw4"="bc_5112179",
				"AIaw1"="bc_7704123",
				"AIaw2"="bc_6864523",
				"AIaw3"="bc_2127641",
				"AIaw4"="bc_5122096")
selSamples <- colnames(freqGT5pct)[c(1:5, 36:46)]
pdf("~/manuscript/plots/heatmap_MCF7TLow5pct4AWs2.pdf", height=3, width=6)	
# pheatmap(freqGT5pct[selWinners,!grepl("5|6|7|8|9|10", colnames(freqGT5pct))], 
pheatmap(freqGT5pct[selWinners2,selSamples], 
	color = viridis(100, begin = 0.2, end = 1, direction = -1, option = "B"),
	na_col = "white",
	annotation_colors = annoColor,
	border_color = NA,
	angle_col = 90,
	# gaps_col = c(3, 5,8,12,14,16,18,20,23, 25,27,29,31,32),
	gaps_col = c(3, 5,7,9, 11, 13, 14, 16),
	# annotation_col = info[,c(5,3,6)], 
	annotation_col = info[,c(5,6)], 
	cluster_rows=F, 
	cluster_cols=F)
dev.off()






freqGT1pct<-freq[rowSums(freq>0.01)>0,]
freqGT1pct[freqGT1pct==0]<-NA
freqGT1pct<-freqGT1pct[order(rowMaxs(freqGT1pct[,1:3]), decreasing=T),]

pdf("plots/heatmap_allMCF7TLow.pdf", height=7, width=12)	
pheatmap(freqGT1pct, 
	color = viridis(100, begin = 0.2, end = 1, direction = -1, option = "B"),
	na_col = "white",
	annotation_colors = annoColor,
	border_color = NA,
	angle_col = 90,
	gaps_col = c(3, 5,8,12,14,16,18,20,22,24,26,28,30, 32, 35,39,41,43,44,46,48,50,51,53,54),
	annotation_col = info[,c(5,3,6)], 
	cluster_rows=F, 
	cluster_cols=F)
dev.off()


