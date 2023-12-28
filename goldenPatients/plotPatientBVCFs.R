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

setwd("~/wgs/golden")
###################################################################################################################
drivers<-read.table("rds/intogenDriverGenes.txt", header=T, stringsAsFactors=F)
drivers.BRCA<-read.table("rds/intogenDriverGenesBRCA.txt", header=T, stringsAsFactors=F)
# deleterious<-read.table("210211-185553_export_variant.tsv", )
damaging<-read.table("rds/damaging.txt", header=F, stringsAsFactors=F)

##################################################################################################################
preB<- readRDS("rds/files/MD/MD_DB_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
postB<-readRDS("rds/files/MD/MD_T_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
preC<- readRDS("rds/files/C/C_DB_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
postC<-readRDS("rds/files/C/C_T_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
preD<- readRDS("rds/files/ML/ML_DB_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
postD<-readRDS("rds/files/ML/ML_T_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
preE<- readRDS("rds/files/E/E_DB_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")
postE<-readRDS("rds/files/E/E_T_AD20_F1R2_6_mutect_calls_filtered_annotated.rds")


###################################################################################################################
listAll<-list("PreB"=preB,
			"PreC"=preC,
			"PreD"=preD,
			"PreE"=preE,
			"PostB"=postB,
			"PostC"=postC,
			"PostD"=postD,		
			"PostE"=postE)
listAll.mutations<-lapply(listAll, function(x){
	data<-x[[2]]$mutations	
	str_info<-as.data.frame(str_split(data$INFO, "\\|", simplify=T), stringsAsFactors=F)
	data$type<-str_info[,2]
	data$impact<-str_info[,3]
	data$gene<-str_info[,4]
	data
})
names(listAll.mutations)<-names(listAll)


listAll.af<-lapply(listAll.mutations, function(x){
	coord<-paste(x$from, x$to, sep="_")
	pos<-paste(x$chr, coord, sep=":")
	mut<-paste(x$ref, x$alt, sep=">")
	var<-paste(pos, mut, sep=":")
	data.frame("Variant"=var, "AF"=x$VAF, "gene"=x$gene, "impact"=x$impact, "type"=gsub("_", ".", x$type))
})
names(listAll.af)<-names(listAll)

df.All<-as.data.frame(acast(melt(listAll.af), Variant+gene+impact+type ~ L1, value.var="value", fill=0), stringsAsFactors=F)

df.preB.postB<-df.All[,c("PreB", "PostB")]
df.preB.postB<-df.preB.postB[rowSums(df.preB.postB>0)>0,]
df.preC.postC<-df.All[,c("PreC", "PostC")]
df.preC.postC<-df.preC.postC[rowSums(df.preC.postC>0)>0,]
df.preD.postD<-df.All[,c("PreD", "PostD")]
df.preD.postD<-df.preD.postD[rowSums(df.preD.postD>0)>0,]
df.preE.postE<-df.All[,c("PreE", "PostE")]
df.preE.postE<-df.preE.postE[rowSums(df.preE.postE>0)>0,]


listAll.sel<-list("preB.postB"=df.preB.postB,
				"preC.postC"=df.preC.postC,
				"preD.postD"=df.preD.postD,
				"preE.postE"=df.preE.postE)
listAll.sel.gene<-lapply(listAll.sel, function(x){
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
	pre<-x[,1]>0 & x[,2]==0
	post<-x[,1]==0 & x[,2]>0
	both<-x[,1]>0 & x[,2]>0
	if(sum(pre)>0){x[pre,]$sample<-colnames(x)[1]}
	if(sum(post)>0){x[post,]$sample<-colnames(x)[2]}
	if(sum(both)>0){x[both,]$sample<-"Both"}
	colnames(x)[1:2]<-c("sample1", "sample2")
	x
})

saveRDS(listAll.sel.gene, "rds/patientBCDE.New.rds")

listAll.sel.gene<-readRDS("rds/patientBCDE.New.rds")
sel.color2<-c("PreB"="#0084A5", "PostB"="#B30078", "Both"="grey40")
sel.color3<-c("PreC"="#0084A5", "PostC"="#B30078", "Both"="grey40")
sel.color4<-c("PreD"="#0084A5", "PostD"="#B30078", "Both"="grey40")
sel.color5<-c("PreE"="#0084A5", "PostE"="#B30078", "Both"="grey40")


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
						nudge_y = 0.1, 
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


png("plots/patientBCDE.delSNVnew4.png", height=350, width=1600)
	plotA<-plotFreq(listAll.sel.gene[[1]],"Diagnostic Biopsy B", "Surgical Biopsy B", sel.color2)
	plotB<-plotFreq(listAll.sel.gene[[2]],"Diagnostic Biopsy C", "Surgical Biopsy C", sel.color3)
	plotC<-plotFreq(listAll.sel.gene[[3]],"Diagnostic Biopsy D", "Surgical Biopsy D", sel.color4)
	plotD<-plotFreq(listAll.sel.gene[[4]],"Diagnostic Biopsy E", "Surgical Biopsy E", sel.color5)
	ggarrange(plotA, plotB, plotC,plotD, nrow=1, ncol=4)
dev.off()

pdf("patientBCDE.delSNVnew4.pdf", height=4, width=4)
	plotFreq(listAll.sel.gene[[1]],"Diagnostic Biopsy B", "Surgical Biopsy B", sel.color2)
	plotFreq(listAll.sel.gene[[2]],"Diagnostic Biopsy C", "Surgical Biopsy C", sel.color3)
	plotFreq(listAll.sel.gene[[3]],"Diagnostic Biopsy D", "Surgical Biopsy D", sel.color4)
	plotFreq(listAll.sel.gene[[4]],"Diagnostic Biopsy E", "Surgical Biopsy E", sel.color5)
	# ggarrange(plotA, plotB, plotC,plotD, nrow=1, ncol=4)
dev.off()

###############################################
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
		geom_point(data=x[mutInfoNew %in% damaging$V1 &  x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) &  rowSums(x[1:2]>=0.1)>0 ,], 
					shape = 21, 
					stroke=0.5, 
					aes(fill=sample), 
					color="black", 
					size=3) +
		geom_text_repel(data=x[mutInfoNew %in% damaging$V1 &  x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) &  rowSums(x[1:2]>=0.1)>0 ,], 
						aes(label = driver, color=sample), 
						# nudge_x = 0.1, 
						# nudge_y = 0.1,
						size=7) +
		geom_label_repel(data=x[x$gene %in% resistanceDrivers & mutInfoNew %in% damaging$V1 &  x$driver!="" & x$driver.brca=="" & grepl("MODERATE|HIGH", x$impact) &  rowSums(x[1:2]>=0.1)>0,], 
						aes(label = driver, fill=sample, color=sample), 
						nudge_x = 0.2, 
						nudge_y = 0.2,
						# color="black", 
						alpha=0.5,
						size=10) +
		theme_bw() + 
		# theme(text = element_text(size=15)) +	
		scale_x_continuous(limits=c(0,1)) +
		scale_y_continuous(limits=c(0,1)) +
		labs(x=a, y=b, fill="Samples") +
		scale_fill_manual(values=sel.color) +
		scale_color_manual(values=sel.color) +
		theme(legend.position="none",
		axis.title = element_text(size=24), 
        axis.text = element_text(size=22))
   #      axis.title=element_blank(),
			# axis.text=element_blank()) 
	ggExtra::ggMarginal(p, type = "histogram", fill="grey")
}



png("plots/patientBCDE.OtherDrivers.delSNVnew4.png", height=350, width=1600)
	plotA<-plotFreqOD(listAll.sel.gene[[1]],"Diagnostic Biopsy B", "Surgical Biopsy B", sel.color2)
	plotB<-plotFreqOD(listAll.sel.gene[[2]],"Diagnostic Biopsy C", "Surgical Biopsy C", sel.color3)
	plotC<-plotFreqOD(listAll.sel.gene[[3]],"Diagnostic Biopsy D", "Surgical Biopsy D", sel.color4)
	plotD<-plotFreqOD(listAll.sel.gene[[4]],"Diagnostic Biopsy E", "Surgical Biopsy E", sel.color5)
	ggarrange(plotA, plotB, plotC,plotD, nrow=1, ncol=4)
dev.off()

pdf("plots/patientBCDE.OtherDrivers.delSNVnew4.pdf", height=5, width=5)
	plotA<-plotFreqOD(listAll.sel.gene[[1]],"Diagnostic Biopsy B", "Surgical Biopsy B", sel.color2)
	plotB<-plotFreqOD(listAll.sel.gene[[2]],"Diagnostic Biopsy C", "Surgical Biopsy C", sel.color3)
	plotC<-plotFreqOD(listAll.sel.gene[[3]],"Diagnostic Biopsy D", "Surgical Biopsy D", sel.color4)
	plotD<-plotFreqOD(listAll.sel.gene[[4]],"Diagnostic Biopsy E", "Surgical Biopsy E", sel.color5)
	# ggarrange(plotA, plotB, plotC,plotD, nrow=1, ncol=4)
	plotA
plotB
plotC
plotD
dev.off()

