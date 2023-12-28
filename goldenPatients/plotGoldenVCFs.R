require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("stringr")
require("viridis")
require("ggpubr")
require("GenomicRanges")
require("rtracklayer")
require("ggrepel")
require("stringr")

setwd("~/wgs/golden")
###################################################################################################################
drivers<-read.table("rds/intogenDriverGenes.txt", header=T, stringsAsFactors=F)
drivers.BRCA<-read.table("rds/intogenDriverGenesBRCA.txt", header=T, stringsAsFactors=F)
damaging<-read.table("rds/damaging.txt", header=F, stringsAsFactors=F)
resistanceDrivers <- c("ESR1", "TP53", "AKT1", "RIC8A", "RB1", "NF1", "FRG1", "KMT2C", "NCOR1")

###################################################################################################################
a13<-readRDS("rds/files/Golden/Golden_CS13_1763_A1_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
a15<-readRDS("rds/files/Golden/Golden_CS15_3459_A5_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
b13<-readRDS("rds/files/Golden/Golden_CS13_1763_B1_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
b15<-readRDS("rds/files/Golden/Golden_CS15_3459_B3_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
a18<-readRDS("rds/files/Golden/Golden_CS18_T_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
###################################################################################################################
listA<-list("A13"=a13,
			"A15"=a15,
			"B13"=b13,
			"B15"=b15,
			"A18"=a18)
listA.mutations<-lapply(listA, function(x){
	data<-x[[1]]$mutations	
	str_info<-as.data.frame(str_split(data$INFO, "\\|", simplify=T), stringsAsFactors=F)
	data$type<-str_info[,2]
	data$impact<-str_info[,3]
	data$gene<-str_info[,4]
	data
})
names(listA.mutations)<-names(listA)


listA.af<-lapply(listA.mutations, function(x){
	coord<-paste(x$from, x$to, sep="_")
	pos<-paste(x$chr, coord, sep=":")
	mut<-paste(x$ref, x$alt, sep=">")
	var<-paste(pos, mut, sep=":")
	data.frame("Variant"=var, "AF"=x$VAF, "gene"=x$gene, "impact"=x$impact, "type"=gsub("_", ".", x$type))
})
names(listA.af)<-names(listA)

df.A<-as.data.frame(acast(melt(listA.af), Variant+gene+impact+type ~ L1, value.var="value", fill=0), stringsAsFactors=F)

df.a13.a15<-df.A[,c("A13", "A15")]
df.a13.a15<-df.a13.a15[rowSums(df.a13.a15>0)>0,]
df.b13.b15<-df.A[,c("B13", "B15")]
df.b13.b15<-df.b13.b15[rowSums(df.b13.b15>0)>0,]
df.a15.a18<-df.A[,c("A15", "A18")]
df.a15.a18<-df.a15.a18[rowSums(df.a15.a18>0)>0,]

listA.sel<-list("A13.A15"=df.a13.a15,
				"B13.B15"=df.b13.b15,
				"A15.A18"=df.a15.a18)

listA.sel.gene<-lapply(listA.sel, function(x){
	x$gene<-as.data.frame(str_split(rownames(x), "_", simplify=T), stringsAsFactors=F)[,3]
	x$impact<-as.data.frame(str_split(rownames(x), "_", simplify=T), stringsAsFactors=F)[,4]
	x$type<-as.data.frame(str_split(rownames(x), "_", simplify=T), stringsAsFactors=F)[,5]
	x$driver<-""
	x$driver.brca<-""
	x$resistance<-""
	check1<-x$gene %in% drivers$Symbol
	check2<-x$gene %in% drivers.BRCA$Symbol 
	check3<-x$gene %in% resistanceDrivers
	
	x[check1,]$driver<-x[check1,]$gene
	x[check2,]$driver.brca<-x[check2,]$gene
	x[check3,]$resistance<-x[check3,]$gene
	x$sample<-""
	x[x[,1]>0 & x[,2]==0,]$sample<-colnames(x)[1]
	x[x[,1]==0 & x[,2]>0,]$sample<-colnames(x)[2]
	x[x[,1]>0 & x[,2]>0,]$sample<-"Both"
	colnames(x)[1:2]<-c("sample1", "sample2")
	x
})

saveRDS(listA.sel.gene, "rds/patientGoldenNew.rds")
listA.sel.gene <- readRDS("rds/patientGoldenNew.rds")
listA.sel.gene.new <- listA.sel.gene
listA.sel.gene<-readRDS("rds/patientGolden.rds")
sel.color1<-c("PreA.R"="#0084A5", "PostA.R"="#B30078", "Both"="grey40")
sel.color2<-c("PreA.L"="#0084A5", "PostA.L"="#B30078", "Both"="grey40")
sel.color3<-c("PostA.R"="#B30078", "RelapseA.R"="#7D21C0", "Both"="grey40")

set1<-listA.sel.gene[[1]]
set2<-listA.sel.gene[[2]]
set3<-listA.sel.gene[[3]]

set1[set1$sample=="A13",]$sample<-"PreA.R"
set1[set1$sample=="A15",]$sample<-"PostA.R"

set2[set2$sample=="B13",]$sample<-"PreA.L"
set2[set2$sample=="B15",]$sample<-"PostA.L"

set3[set3$sample=="A15",]$sample<-"PostA.R"
set3[set3$sample=="A18",]$sample<-"RelapseA.R"

plotFreq<-function(x, a, b, sel.color){
	mutInfo <- rownames(x)
	splitMutInfo <- str_split(mutInfo, "_", simplify=TRUE)
	splitMutInfo[,5] <- str_replace_all(splitMutInfo[,5], "\\.", "_")
	mutInfoNew <- apply(splitMutInfo[,c(1,2,3,5)], 1 , paste , collapse = "_" )
	p <- ggplot(x, aes(x=sample1, y=sample2)) + 
		geom_point(data=x, 
					shape = 21, 
					stroke=0, 
					fill="gray50", 
					color="gray50", 
					alpha=0.3, 
					size=2) +
		# geom_point(data=x[rownames(x)%in%deleterious$V1 & x$driver.brca!="" & !grepl("synonymous|intron|intergenic|upstream|downstream", x$type) & rowSums(x[1:2]>=0.1)>0 ,], 
		geom_point(data=x[mutInfoNew %in% damaging$V1 & x$driver.brca!="" & x$impact %in% c("MODERATE", "HIGH") & rowSums(x[1:2]>=0.1)>0 ,], 					
					shape = 21, 
					stroke=0.5, 
					aes(fill=sample), 
					color="black", 
					size=3) +
		# geom_label_repel(data=x[rownames(x)%in%deleterious$V1 & x$gene %in% resistanceDrivers & rownames(x)%in%deleterious$V1 & x$driver.brca!="" & !grepl("synonymous|intron|intergenic|upstream|downstream", x$type) & rowSums(x[1:2]>=0.1)>0,], 
		geom_label_repel(data=x[mutInfoNew %in% damaging$V1 & x$gene %in% resistanceDrivers & grepl("MODERATE|HIGH", x$impact) & rowSums(x[1:2]>=0.1)>0,], 
						aes(label = driver.brca, fill=sample, color=sample), 
						nudge_x = 0.2, 
						nudge_y = 0.2,
						# color="black", 
						alpha=0.5,
						size=8) +
		theme_bw() + 
		# geom_text_repel(data=x[rownames(x)%in%deleterious$V1 & !x$gene %in% resistanceDrivers & rownames(x)%in%deleterious$V1 & x$driver.brca!="" & !grepl("synonymous|intron|intergenic|upstream|downstream", x$type) & rowSums(x[1:2]>=0.1)>0,], 
		geom_text_repel(data=x[mutInfoNew %in% damaging$V1 & x$driver.brca!="" & !x$gene %in% resistanceDrivers & grepl("MODERATE|HIGH", x$impact) & rowSums(x[1:2]>=0.1)>0,], 
						aes(label = driver.brca, color=sample), 
						max.overlaps = 100,
						nudge_x = 0.2, 
						nudge_y = 0.2, 
						size=7) +
		scale_x_continuous(limits=c(0,1)) +
		scale_y_continuous(limits=c(0,1)) +
		labs(x=a, y=b, fill="BreastCancer", color="Resistance") +
		scale_fill_manual(values=sel.color) +
		scale_color_manual(values=sel.color) +
		theme(legend.position="none",
		axis.title = element_text(size=24), 
        axis.text = element_text(size=22))+
			# axis.title=element_blank(),
			# axis.text=element_blank()) +
		guides(color=guide_legend(nrow=3,byrow=TRUE),
			fill=guide_legend(nrow=3,byrow=TRUE))
	ggExtra::ggMarginal(p, type = "histogram", fill="grey")
}

png("plots/patientGolden.delSNV.new4.png", height=400, width=1400)
# png("Rplots.png", height=400, width=1400)
	plotA<-plotFreq(set1,"Diagnostic Biopsy A (Right)", "Surgical Biopsy A (Right)", sel.color1)
	plotB<-plotFreq(set2,"Diagnostic Biopsy A (Left)", "Surgical Biopsy A (Left)", sel.color2)
	plotC<-plotFreq(set3,"Surgical Biopsy A (Right)", "Local Relapse A (Right)", sel.color3)
	ggarrange(plotA, plotB, plotC, nrow=1, ncol=3)
dev.off()

pdf("plots/patientGolden.delSNV.new4.pdf", height=4, width=4)
	plotA<-plotFreq(set1,"Diagnostic Biopsy A (Right)", "Surgical Biopsy A (Right)", sel.color1)
	plotB<-plotFreq(set2,"Diagnostic Biopsy A (Left)", "Surgical Biopsy A (Left)", sel.color2)
	plotC<-plotFreq(set3,"Surgical Biopsy A (Right)", "Local Relapse A (Right)", sel.color3)
	# ggarrange(plotA, plotB, plotC, nrow=1, ncol=3)
	plotA
	plotB
	plotC
dev.off()

#######################################################

plotFreqOD<-function(x, a, b, sel.color){
	mutInfo <- rownames(x)
	splitMutInfo <- str_split(mutInfo, "_", simplify=TRUE)
	splitMutInfo[,5] <- str_replace_all(splitMutInfo[,5], "\\.", "_")
	mutInfoNew <- apply(splitMutInfo[,c(1,2,3,5)], 1 , paste , collapse = "_" )
	p <- ggplot(x, aes(x=sample1, y=sample2)) + 
		geom_point(data=x, 
					shape = 21, 
					stroke=0, 
					fill="gray50", 
					color="gray50", 
					alpha=0.3, 
					size=2) +
		geom_point(data=x[mutInfoNew %in% damaging$V1  & x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) & rowSums(x[1:2]>=0.1)>0 ,], 
					shape = 21, 
					stroke=0.5, 
					aes(fill=sample), 
					color="black", 
					size=3) +
		geom_text_repel(data=x[mutInfoNew %in% damaging$V1  & x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) & rowSums(x[1:2]>=0.1)>0,], 
						aes(label = driver, color=sample), 
						nudge_x = 0.2, 
						nudge_y = 0.2,
						size=7) +
		geom_label_repel(data=x[x$gene %in% resistanceDrivers & mutInfoNew %in% damaging$V1  & x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) & rowSums(x[1:2]>=0.1)>0,], 
						aes(label = driver, fill=sample, color=sample), 
						nudge_x = 0.2, 
						nudge_y = 0.2,
						# color="black", 
						alpha=0.5,
						size=10) +
		theme_bw() + 
		scale_x_continuous(limits=c(0,1)) +
		scale_y_continuous(limits=c(0,1)) +
		labs(x=a, y=b, fill="Samples") +
		scale_fill_manual(values=sel.color) +
		scale_color_manual(values=sel.color) +
		theme(legend.position="none",
		axis.title = element_text(size=24), 
        axis.text = element_text(size=22))
			# axis.title=element_blank(),
			# axis.text=element_blank())  
	ggExtra::ggMarginal(p, type = "histogram", fill="grey")
}


		# axis.title = element_text(size=24), 
  #       axis.text = element_text(size=22),
		# 	axis.text.y=element_blank(),
		# axis.ticks.y=element_blank())

png("plots/patientGolden.OtherDrivers.delSNVnew4.png", height=500, width=1500)
	plotA<-plotFreqOD(set1,"Diagnostic Biopsy A (Right)", "Surgical Biopsy A (Right)", sel.color1)
	plotB<-plotFreqOD(set2,"Diagnostic Biopsy A (Left)", "Surgical Biopsy A (Left)", sel.color2)
	plotC<-plotFreqOD(set3,"Surgical Biopsy A (Right)", "Local Relapse A (Right)", sel.color3)
	ggarrange(plotA, plotB, plotC, nrow=1, ncol=3)
dev.off()

pdf("plots/patientGolden.OtherDrivers.delSNV4.pdf", height=5, width=5)
	plotA<-plotFreqOD(set1,"Diagnostic Biopsy A (Right)", "Surgical Biopsy A (Right)", sel.color1)
	plotB<-plotFreqOD(set2,"Diagnostic Biopsy A (Left)", "Surgical Biopsy A (Left)", sel.color2)
	plotC<-plotFreqOD(set3,"Surgical Biopsy A (Right)", "Local Relapse A (Right)", sel.color3)
	# ggarrange(plotA, plotB, plotC, nrow=1, ncol=3)
	plotA
	plotB
	plotC
dev.off()
