require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("stringr")
require("viridis")
require("ggrepel")
require("ggpubr")

names.samples<-read.table("rds/rename.all.txt", header=TRUE, stringsAsFactors=F)
winners<-read.table("rds/winners.txt", header=T, stringsAsFactors=F)

counts.all<-readRDS("rds/counts.all.T0.rds")

t0<-counts.all[,c(1:3, 112:131)]
info<-data.frame("reads"=colSums(t0)/1e6, "barcodes"=colSums(t0>0)/1e3, "samples"=colnames(t0))

data<-counts.all[,c(1:62, 112:131)]
data<-data[,-c(22,23)]

total<-colSums(data) 
freq<-as.data.frame(t(t(data)/total), stringsAsFactors=F)
freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")

pot<-90*1e6*(freq[,1:3])
pot$time.1<-(24*13*log10(2))/log10(2*pot$POT.1)
pot$time.2<-(24*13*log10(2))/log10(2*pot$POT.2)
pot$time.3<-(24*13*log10(2))/log10(2*pot$POT.3)

pot1<-pot[, grep("1", colnames(pot))]
pot2<-pot[, grep("2", colnames(pot))]
pot3<-pot[, grep("3", colnames(pot))]

time<-pot[,4:6]
time<-time[rowSums(time>0)>1,]
time.m<-melt(as.matrix(time))

### Doubling time density plot ###
pdf("plots/suppFigure11.pdf", height=3, width=8)
ggplot(time.m[time.m$value<=330,], aes(x=value, color=Var2)) +
	geom_density() +
	theme_bw() +
	scale_x_continuous(breaks=seq(0, 400, by = 10), 
		labels=seq(0, 400, by = 10)) +
	scale_color_discrete(labels=c("POT1", "POT2", "POT3")) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
	labs(x="Time (hours)", y="Density", color="Sample POT")
dev.off()

