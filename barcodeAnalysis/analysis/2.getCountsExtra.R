require("reshape2")
require("ggplot2")
require("stringr")
require("fishplot")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")

## Get frequencies from batch 1 and batch 2
setwd("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered")
data<-readRDS("rds/counts.all.rds")
data<-data[rowSums(data>0)>0,]
bcs<-rownames(data)
total<-colSums(data)
freq<-t(t(data)/total)

## Get frequencies from batch 3
counts.batch3<-readRDS("rds/counts_batch3.rds")
counts.batch3<-counts.batch3[rowSums(counts.batch3>0)>0, ]
total2<-colSums(counts.batch3)
freq2<-t(t(counts.batch3)/total2)

## Get barcodes common to batch1+2 and batch 3
common<-rownames(counts.batch3)[rownames(counts.batch3) %in% rownames(data)]
common.samples<-colnames(data)[colnames(data) %in% colnames(counts.batch3)]
## Get missing samples in batch3
missing.batch3.sample<-colnames(data)[!colnames(data) %in% colnames(counts.batch3)]
## Get missing barcodes in batch3
missing.batch3.bcs<-rownames(data)[!rownames(data) %in% rownames(counts.batch3)]
## Get extra barcodes in batch3
extra.batch3.bcs<-rownames(counts.batch3)[!rownames(counts.batch3) %in% rownames(data)]


############################################
## Get BC info for missing samples in batch 3
data.prev<-data[,missing.batch3.sample]
## Get BC info missing in batch 3 for common samples in batch 1+2
data.prev2<-data[!rownames(data) %in% common, colnames(data) %in% common.samples]
## Get common BC info for common samples in batch 1+2
data.extra<-data[rownames(data) %in% common, colnames(data) %in% common.samples]
## Get common BC info for common samples in batch 3
data.added<-counts.batch3[common,]+data.extra[common,]
## Merge previous and new counts and rearrange
data.extra2<-rbind(data.prev2, data.added)
data2<-merge(data.prev, data.extra2, by="row.names")
data.recal<-data2[,2:117]
rownames(data.recal)<-data2[,1]
data.recal<-data.recal[,colnames(data)]

## Get frequency after adding resequenced data
total.recal<-colSums(data.recal)
freq.recal<-t(t(data.recal)/total.recal)

freq.m<-melt(as.matrix(freq[,c(1:21,24:50)]))
freq.recal.m<-melt(as.matrix(freq.recal[,c(1:21,24:50)]))

## Compare previous frequency distribution with new distribution in each sample
freq.m$version<-"old"
freq.recal.m$version<-"revised.freq"

info<-rbind(freq.m, freq.recal.m)
rename<-read.table("rds/sample.info.all.t0.txt", header=T, stringsAsFactors=F)

info$Var2 <-as.vector(info$Var2)
info$sample<-""
info[info$Var2 %in% rename$sample,]$sample<-rename[match(info[info$Var2 %in% rename$sample,]$Var2, rename$sample),]$sample.renamed
str_info<-as.data.frame(str_split(info$sample, "\\.", simplify=T), stringsAsFactors=F)
info$condition<-str_info$V1
info$flask<-str_info$V2
info$condition<-factor(info$condition, levels=c("POT", "U7", "U14", "U30", "U60", "U90", "U120", "U150", "AI7", "AI14", "AI30", "AI60", "AI90", "AI", "T7", "T14", "T30", "T60", "T"))

pdf("plots/densityPlotCounts.pdf", height=21, width=21)
ggplot(info, aes(x=log10(value), fill=version)) +
	geom_density(alpha=0.4) +
	facet_wrap(condition ~ flask, scales="free") +
	theme_bw() +
	theme(text = element_text(size=20))
dev.off()

## Get counts matrix for additionally identfied BCs
extra<-counts.batch3[extra.batch3.bcs, ]
extra.0<-as.data.frame(matrix(0, nrow = length(extra.batch3.bcs), ncol=length(missing.batch3.sample)))
colnames(extra.0)<-colnames(data)[!colnames(data) %in% colnames(counts.batch3)]
rownames(extra.0)<-extra.batch3.bcs
extra.new<-cbind(extra, extra.0)
extra.new<-extra.new[,colnames(data)]

data.recal<-rbind(data.recal, extra.new)

saveRDS(data.recal, "rds/data.recal.rds")
write.table(data.recal,file ="rds/data.recal.txt",sep="\t",quote=F,row.names=T,col.names=T)
########## Include T0s ########################

counts.recal <- readRDS("rds/data.recal.rds")
counts.batch4<-readRDS("rds/counts_batch4.rds")

counts.all.T0<-merge(counts.recal, counts.batch4, by="row.names", all=TRUE)
counts.all.T0[is.na(counts.all.T0)] <- 0

saveRDS(counts.all.T0, "rds/counts.all.T0.rds")
write.table(counts.all.T0,file ="rds/counts.all.T0.txt",sep="\t",quote=F,row.names=T,col.names=T)


