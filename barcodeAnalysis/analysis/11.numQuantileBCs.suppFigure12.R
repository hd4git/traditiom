require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("stringr")
require("viridis")
require("ggrepel")
require("ggpubr")
require("matrixStats")
setwd("~/barcode_analysis/results/filtered")
freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")
freq<-freq[,1:68]
freq<-freq[rowSums(freq>0)>0,]

pot<-rowMedians(as.matrix(freq[,1:3]))
names(pot)<-rownames(freq)
pot<-pot[pot>0]
pot.m<-as.data.frame(melt(pot))
pot.m$sample<-"POT"
pot.m$quantile<-""
pot.m[pot.m$value>=quantile(pot)[1] & pot.m$value<=quantile(pot)[2],]$quantile=1
pot.m[pot.m$value>quantile(pot)[2] & pot.m$value<=quantile(pot)[3],]$quantile=2
pot.m[pot.m$value>quantile(pot)[3] & pot.m$value<=quantile(pot)[4],]$quantile=3
pot.m[pot.m$value>quantile(pot)[4] & pot.m$value<=quantile(pot)[5],]$quantile=4

table(pot.m$quantile)
#     1     2     3     4 
# 24101 24054 24058 24061 
length(pot)
# [1] 96274


bc.1<-rownames(pot.m[pot.m$quantile==1,])
bc.2<-rownames(pot.m[pot.m$quantile==2,])
bc.3<-rownames(pot.m[pot.m$quantile==3,])
bc.4<-rownames(pot.m[pot.m$quantile==4,])

colors<-c(rep("#AEB4A9", each=length(bc.1)), 
		rep("#E0C1B3", each=length(bc.2)), 
		rep("#D89A9E", each=length(bc.3)), 
		rep("#C37D92", each=length(bc.4)))
names(colors)<-c(bc.1, bc.2, bc.3, bc.4)


plotFreq<-function(ip){
	ip<-ip[rowSums(ip>0)>0,	]
	ip.m<-melt(as.matrix(ip))
	ip.m<-ip.m[ip.m$value>0,]

	ggplot(ip.m, aes(x=Var2, y=log10(value))) +
		geom_jitter(aes(color=Var1), size=0.5, width=0.3) +
		geom_violin(alpha=0) +
		scale_colour_manual(values=colors) +
		theme_bw() +
		theme(legend.position="none",
			axis.text = element_text(size=16),
			axis.title = element_text(size=18)) +
		labs(x="", y="log10(Frequency)")
}
m1<-freq[,grepl("AI30|T30", colnames(freq))]
m2<-freq[,grepl("AI60|T60", colnames(freq))]

png("plots_MS/m1freq.png", height=300, width=600)
	plotFreq(m1)
dev.off()
png("plots_MS/m2freq.png", height=300, width=600)
	plotFreq(m2)
dev.off()

pdf("plots/m1freq.pdf")
	plotFreq(m1)
dev.off()
pdf("plots/m2freq.pdf")
	plotFreq(m2)
dev.off()


freq.m<-melt(as.matrix(freq))
freq.m<-freq.m[freq.m$value>0,]
selected<-c("AI30.1", "AI30.2", "T30.1", "T30.2", "AI60.1", "AI60.2", "T60.1", "T60.2")
pctQuantileBCs<-lapply(selected, function(x){
	# c("bc.1"=length(unique(freq.m[freq.m$Var1 %in% bc.1 & freq.m$Var2==x,]$Var1))/length(bc.1)*100,
	#   "bc.2"=length(unique(freq.m[freq.m$Var1 %in% bc.2 & freq.m$Var2==x,]$Var1))/length(bc.2)*100,
	#   "bc.3"=length(unique(freq.m[freq.m$Var1 %in% bc.3 & freq.m$Var2==x,]$Var1))/length(bc.3)*100,
	#   "bc.4"=length(unique(freq.m[freq.m$Var1 %in% bc.4 & freq.m$Var2==x,]$Var1))/length(bc.4)*100)
	c("bc.1"=length(unique(freq.m[freq.m$Var1 %in% bc.1 & freq.m$Var2==x,]$Var1))/length(unique(freq.m[freq.m$Var2==x,]$Var1))*100,
	  "bc.2"=length(unique(freq.m[freq.m$Var1 %in% bc.2 & freq.m$Var2==x,]$Var1))/length(unique(freq.m[freq.m$Var2==x,]$Var1))*100,
	  "bc.3"=length(unique(freq.m[freq.m$Var1 %in% bc.3 & freq.m$Var2==x,]$Var1))/length(unique(freq.m[freq.m$Var2==x,]$Var1))*100,
	  "bc.4"=length(unique(freq.m[freq.m$Var1 %in% bc.4 & freq.m$Var2==x,]$Var1))/length(unique(freq.m[freq.m$Var2==x,]$Var1))*100)
})
names(pctQuantileBCs)<-selected
pctQuantileBCs.df<-do.call(rbind, pctQuantileBCs)

pctQuantileBCs.df<-t(pctQuantileBCs.df)
new<-(100-colSums(pctQuantileBCs.df))
pctQuantileBCs.df<-rbind(pctQuantileBCs.df, new)

png("plots/m1freq.png", height=600, width=900)
ggarrange(p, plotFreq(m1), widths = c(0.2, 1))
dev.off()
png("plots/m2freq.png", height=600, width=900)
ggarrange(p, plotFreq(m2), widths = c(0.2, 1))
dev.off()

colors2<-c("white", "#AEB4A9","#E0C1B3","#D89A9E","#C37D92")
names(colors2)<-c("new", "bc.1", "bc.2", "bc.3", "bc.4")

pctQuantileBCs.m<-melt(as.matrix(pctQuantileBCs.df))
pctQuantileBCs.m$Var1<-factor(pctQuantileBCs.m$Var1, levels=names(colors2))

pdf("plots_MS/pieFigure4a.pdf", height=0.7, width=0.7, bg="transparent")
lapply(colnames(pctQuantileBCs.df), function(x){
ggplot(pctQuantileBCs.m[pctQuantileBCs.m$Var2==x,], aes(x="", y=value, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  	scale_fill_manual(values=colors2) +
	theme_void()  +
		theme(legend.position = "none")
})
dev.off()


numQuantileBCs<-lapply(selected, function(x){
	c("bc.1"=length(unique(freq.m[freq.m$Var1 %in% bc.1 & freq.m$Var2==x,]$Var1)),
	  "bc.2"=length(unique(freq.m[freq.m$Var1 %in% bc.2 & freq.m$Var2==x,]$Var1)),
	  "bc.3"=length(unique(freq.m[freq.m$Var1 %in% bc.3 & freq.m$Var2==x,]$Var1)),
	  "bc.4"=length(unique(freq.m[freq.m$Var1 %in% bc.4 & freq.m$Var2==x,]$Var1)))
})
names(numQuantileBCs)<-selected

numQuantileBCs.df<-do.call(rbind, numQuantileBCs)
numQuantileBCs.m<-melt(as.matrix(numQuantileBCs.df))
numQuantileBCs.m$condition<-""
numQuantileBCs.m[grepl("AI", numQuantileBCs.m$Var1),]$condition<-"-E2"
numQuantileBCs.m[grepl("T", numQuantileBCs.m$Var1),]$condition<-"T"
numQuantileBCs.m$quantile<-""
numQuantileBCs.m[grepl("bc.1", numQuantileBCs.m$Var2),]$quantile<-"low"
numQuantileBCs.m[grepl("bc.2", numQuantileBCs.m$Var2),]$quantile<-"lowMid"
numQuantileBCs.m[grepl("bc.3", numQuantileBCs.m$Var2),]$quantile<-"midHigh"
numQuantileBCs.m[grepl("bc.4", numQuantileBCs.m$Var2),]$quantile<-"high"
numQuantileBCs.m$time<-""
numQuantileBCs.m[grepl("30", numQuantileBCs.m$Var1),]$time=30
numQuantileBCs.m[grepl("60", numQuantileBCs.m$Var1),]$time=60

colors3<-c("#AEB4A9","#E0C1B3","#D89A9E","#C37D92")
names(colors3)<-c("low", "lowMid", "midHigh", "high")

pdf("plots_MS/fig4d1.pdf", height=2.5, width=2.5)
ggplot(numQuantileBCs.m, aes(x=time, y=value, shape=condition, color=quantile)) +
	geom_point(position=position_dodge2(width=0.3), size=2) +
	theme_bw() +
	scale_color_manual(values=colors3) +
	labs(x="Time(days)", y="Surviving barcodes", shape="Treatment", color="Quantile") +
	theme(axis.title=element_text(size=10),
		axis.text=element_text(size=8),
		# legend.position="none",
		legend.key.size = unit(0.08, "cm"),
		legend.title = element_text(size=9), 
        legend.text = element_text(size=8))+
	guides(color=guide_legend(nrow=4,byrow=TRUE),
		shape=guide_legend(nrow=2,byrow=TRUE))
dev.off()
pdf("plots_MS/fig4d2.pdf", height=3.5, width=3.5)
ggplot(numQuantileBCs.m, aes(x=time, y=value, shape=condition, color=quantile)) +
	geom_point(position=position_dodge2(width=0.3), size=2) +
	theme_bw() +
	scale_color_manual(values=colors3) +
	labs(x="Time(days)", y="Surviving barcodes", shape="Treatment", color="Quantile") +
	theme(axis.title=element_text(size=10),
		axis.text=element_text(size=8),
		legend.position="bottom",
		legend.key.size = unit(0.08, "cm"),
		legend.title = element_text(size=9), 
        legend.text = element_text(size=8))+
	guides(color=guide_legend(nrow=2,byrow=TRUE),
		shape=guide_legend(nrow=2,byrow=TRUE))
dev.off()


