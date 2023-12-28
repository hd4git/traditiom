require("reshape2")
require("ggplot2")
require("stringr")
require("fishplot")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")
require("viridis")

## Import samples
setwd("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered")
samples<-list.files("allT47D")
ids<-do.call(rbind, str_split(samples, pattern="\\.", n=3))[,1]
sampleInfo <- read.table("rds/allT47DsampleInfo.txt", header=T)

## Get raw counts and frequency
counts.list<-lapply(samples, function(x){read.table(paste("allT47D",x, sep="/"), header=F)})
names(counts.list)<-ids

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts)
counts<-counts[,sampleInfo$sampleID]
colnames(counts)<-sampleInfo$sampleName
counts[is.na(counts)] <- 0
counts<-counts[,1:38]
counts.gt10<-counts[rowSums(counts>10)>0,]
counts.gt20<-counts[rowSums(counts>20)>0,]
counts.gt30<-counts[rowSums(counts>30)>0,]
counts.gt40<-counts[rowSums(counts>40)>0,]
counts.gt50<-counts[rowSums(counts>50)>0,]

counts<-counts.gt10
counts<-counts[rowSums(counts[,1:7]>0)>0,]
total<-colSums(counts)
freq<-t(t(counts)/total)
colnames(counts)[17:18]<-c("AI1m1", "AI2m1")
colnames(freq)[17:18]<-c("AI1m1", "AI2m1")
saveRDS(counts, "rds/counts_TLow_T47D.rds")
counts <- readRDS("rds/counts_TLow_T47D.rds")
counts<-counts[,1:30]
freq<-freq[,1:30]
freq.m<-melt(freq)

saveRDS(freq.m, "rds/freq.m_TLow_T47D.rds")

############################################
####### Surviving BCs barcplots ############
############################################

bcCount <- as.data.frame(colSums(counts>0))
colnames(bcCount) <- "counts"
bcCount$sample <- rownames(bcCount)
bcCount$condition <- ""
bcCount[grepl("^POT", bcCount$sample),]$condition <- "POT"
bcCount[grepl("^T0", bcCount$sample),]$condition <- "T0"
bcCount[grepl("m", bcCount$sample),]$condition <- "Dormancy"
bcCount[grepl("aw", bcCount$sample),]$condition <- "Awakening"
bcCount[grepl("TEP", bcCount$sample),]$condition <- "Relapse"
bcCount[grepl("^UT", bcCount$sample),]$condition <- "Untreated"

bcCount$treatment <- ""
bcCount[grepl("^POT", bcCount$sample),]$treatment <- "POT"
bcCount[grepl("^T0", bcCount$sample),]$treatment <- "POT"
bcCount[grepl("^AI|aw|TEP", bcCount$sample),]$treatment <- "-E2"
bcCount[grepl("^UT", bcCount$sample),]$treatment <- "UT"
bcCount$treatment <- factor(bcCount$treatment, levels=c("POT", "UT", "-E2"))

bcCount$condition <- factor(bcCount$condition, levels=c("POT", "T0", "Untreated", "Dormancy", "Awakening", "Relapse"))
bcCount$replicate <- ""
bcCount[grepl("_", bcCount$sample),]$replicate <- as.data.frame(str_split(bcCount[grepl("_", bcCount$sample),]$sample, "_", simplify=T), stringsAsFactors=F)[,1]
bcCount[grepl("\\.", bcCount$sample),]$replicate <- as.data.frame(str_split(bcCount[grepl("\\.", bcCount$sample),]$sample, "\\.|_", simplify=T), stringsAsFactors=F)[,2]

bcCount$sample <- factor(bcCount$sample, levels=c(paste("POT", 1:3, sep="."), paste("T0", 1:4, sep="."), 
						paste("UT", 1:3, sep="."), paste(paste("UT", 1:3, sep="."), "aw", sep="_"), paste(paste("UT", 1:3, sep="."), "TEP", sep="_"), 
						"AI1m1", "AI2m1", paste(LETTERS[1:10], "aw", sep="_"), paste(LETTERS[1:10], "TEP", sep="_")))
scheme <- c("POT"="grey","T0"="grey", "Untreated"="#000000", "Dormancy"="#FFDC00", "Awakening"="#B30078", "Relapse"="#50157A")
pdf("~/manuscript/plots/survivingBCsT47D.pdf", height=10, width =10)
ggplot(bcCount, aes(x=sample, y=counts, fill=condition)) + 
geom_bar(stat="identity", position="dodge2", color="black") + 
geom_text(stat="identity", aes(label=counts), position=position_stack(vjust = 0.96), vjust = -1) +
theme_bw() +  
labs(x=" ", y="Number of barcodes", fill="Condition") + 
theme(text = element_text(size=9),
	axis.text = element_text(size=9),
	axis.text.x = element_text(size=9, angle = 90),  
	legend.position="none", 
	strip.background=element_rect(fill="white")) +
scale_fill_manual(values=scheme) +
facet_wrap(treatment ~ ., ncol=1, scales="free_x") 
dev.off()

############################
####### Heatmap ############
############################
freq.m <- readRDS("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered/rds/freq.m_TLow_T47D.rds")
freqAll.m <- freq.m
freqGT10pct <- freq[rowSums(freq>0.1)>0,]
freq.m <- melt(freqGT10pct)
freq.m$sample <- freq.m$Var2
freq.m$condition <- ""
freq.m[grepl("^POT", freq.m$sample),]$condition <- "POT"
freq.m[grepl("^T0", freq.m$sample),]$condition <- "T0"
freq.m[grepl("m", freq.m$sample),]$condition <- "Dormancy"
freq.m[grepl("aw|TEP", freq.m$sample),]$condition <- "Awakening"
freq.m[grepl("^UT", freq.m$sample),]$condition <- "Untreated"
sampleInfo$sampleName[17:18]<-c("AI1m1", "AI2m1")
freq.m$sample <- factor(freq.m$sample, levels=sampleInfo$sampleName)
freq.m$condition <- factor(freq.m$condition, levels=c("POT", "T0", "Untreated", "Dormancy", "Awakening", "Relapse"))
scheme <- c("POT"="grey","T0"="grey", "Untreated"="#000000", "Dormancy"="#FFDC00", "Awakening"="#B30078", "Relapse"="#50157A")
freq.m$replicate <- ""
freq.m[grepl("_", freq.m$sample),]$replicate <- as.data.frame(str_split(freq.m[grepl("_", freq.m$sample),]$sample, "_", simplify=T), stringsAsFactors=F)[,1]
freq.m[grepl("\\.", freq.m$sample),]$replicate <- as.data.frame(str_split(freq.m[grepl("\\.", freq.m$sample),]$sample, "\\.|_", simplify=T), stringsAsFactors=F)[,2]

freq.m[freq.m$value==0,]$value<-NA

freq.m$Var1 <- factor(freq.m$Var1, levels=rownames(freqGT10pct[order(rowMax(freqGT10pct[,19:30]), decreasing=F),]))

freq.m$groups<-""
freq.m[grepl("POT", freq.m$sample),]$groups <- "POT"
freq.m[grepl("T0", freq.m$sample),]$groups <- "T0"
freq.m[grepl("_aw|_TEP", freq.m$sample),]$groups <- paste("Aw", str_split(freq.m[grepl("_aw|_TEP", freq.m$sample),]$sample, "_", simplify=T)[,1], sep="_")
freq.m[grepl("UT.1", freq.m$sample),]$groups <- "UT.1"
freq.m[grepl("UT.2", freq.m$sample),]$groups <- "UT.2"
freq.m[grepl("UT.3", freq.m$sample),]$groups <- "UT.3"
freq.m[grepl("m", freq.m$sample),]$groups <- "Dormancy"
freq.m$groups <- factor(freq.m$groups, levels=c("POT", "T0", "UT.1", "UT.2", "UT.3", "Dormancy", "Aw_A", "Aw_B", "Aw_C", "Aw_D", "Aw_E", "Aw_F", "Aw_G", "Aw_H", "Aw_I", "Aw_J"))
freq.m$Var1 <- factor(freq.m$Var1, levels=rev(c("bc_7030761", "bc_8515639", "bc_5366038", "bc_4516270", "bc_3304724", "bc_3656514", "bc_4047334", "bc_911374", "bc_138465", "bc_5764295", "bc_1732792", "bc_517798", "bc_9499388", "bc_4204733", "bc_5071444", "bc_1146399", "bc_3566657", "bc_5078363", "bc_3685591")))
# pdf("~/manuscript/plots/winnersHeatmapT47D.pdf", height=7, width=10)
pdf("~/manuscriptFeb23/plots/winnersHeatmapT47D.pdf", height=7, width=8)
ggplot(freq.m[!freq.m$Var2 %in% unique(freq.m$Var2)[c(1:8, 11, 14:16)],], aes(x=sample, y=Var1, fill=value)) + 
	geom_tile() + 
	# geom_tile(data=freq.min1.m[freq.min1.m$Var1 %in% bcCommon,], aes(x=Var2, y=Var1), color="red") + 	
	scale_fill_viridis(option="A", begin = 0.25, end = 1, direction=-1, na.value ="white", limits = c(0,1)) + 
	theme_bw() + 
	theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=0.5)) +
	labs(x="Samples", y="Barcodes", fill="Frequency > 0.1") +
	theme(text = element_text(size=12), 
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom") +
	# scale_x_discrete(labels=relabel) 
	# facet_wrap(. ~ condition+replicate, nrow=1)
	facet_grid(. ~ groups, scale="free_x", space="free_x")
dev.off()



#################################
####### Violin plots ############
#################################
freq.m<-melt(freq)
freq.m$sample <- freq.m$Var2
freq.m$condition <- ""	
freq.m[grepl("^POT", freq.m$sample),]$condition <- "POT"
freq.m[grepl("^T0", freq.m$sample),]$condition <- "T0"
freq.m[grepl("m", freq.m$sample),]$condition <- "Dormancy"
freq.m[grepl("aw|TEP", freq.m$sample),]$condition <- "Awakening"
# freq.m[grepl("^UT", freq.m$sample),]$condition <- "Untreated"

freq.m$sample <- factor(freq.m$sample, levels=sampleInfo$sampleName)
freq.m$condition <- factor(freq.m$condition, levels=c("POT", "T0", "Untreated", "Dormancy", "Awakening", "Relapse"))
freq.m$replicate <- ""
freq.m[grepl("_", freq.m$sample),]$replicate <- as.data.frame(str_split(freq.m[grepl("_", freq.m$sample),]$sample, "_", simplify=T), stringsAsFactors=F)[,1]
freq.m[grepl("\\.", freq.m$sample),]$replicate <- as.data.frame(str_split(freq.m[grepl("\\.", freq.m$sample),]$sample, "\\.|_", simplify=T), stringsAsFactors=F)[,2]

freq.m$value[is.na(freq.m$value)]=0
freq.m<-freq.m[freq.m$value>0,]
freq.m$groups<-""
freq.m[grepl("POT", freq.m$sample),]$groups <- "POT"
freq.m[grepl("T0", freq.m$sample),]$groups <- "T0"
freq.m[grepl("_aw|_TEP", freq.m$sample),]$groups <- paste("Aw", str_split(freq.m[grepl("_aw|_TEP", freq.m$sample),]$sample, "_", simplify=T)[,1], sep="_")
# freq.m[grepl("UT.1", freq.m$sample),]$groups <- "UT.1"
# freq.m[grepl("UT.2", freq.m$sample),]$groups <- "UT.2"
# freq.m[grepl("UT.3", freq.m$sample),]$groups <- "UT.3"
freq.m[grepl("m", freq.m$sample),]$groups <- "Dormancy"
freq.m$groups <- factor(freq.m$groups, levels=c("POT", "T0", "UT.1", "UT.2", "UT.3", "Dormancy", "Aw_A", "Aw_B", "Aw_C", "Aw_D", "Aw_E", "Aw_F", "Aw_G", "Aw_H", "Aw_I", "Aw_J"))
# freq.m$groups <- factor(freq.m$groups, levels=c("POT", "T0", "UT.1", "UT.2", "UT.3", "Dormancy", "Aw_A", "Aw_B", "Aw_C", "Aw_D", "Aw_E", "Aw_F"))
freq.m <- freq.m[!grepl("UT", freq.m$Var2),]
winners<-c("bc_4047334", "bc_911374", "bc_138465", "bc_5764295", "bc_1732792", "bc_517798", "bc_9499388", "bc_4204733", "bc_5071444", "bc_1146399", "bc_3566657", "bc_5078363", "bc_3685591")
freqSel.m <- freq.m[freq.m$Var1 %in% winners,]
scheme <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#0868ac", "blue")
names(scheme) <- winners


pdf("~/manuscript/plots/dataAll.T47Dhigh.violin.pdf", height=5, width=10)
ggplot(freq.m, aes(x=Var2, y=log10(value))) + 
	geom_violin(fill="grey90", color="grey90") + 
	geom_point(data=freqSel.m, aes(x=Var2, y=log10(value), color=Var1), position = position_dodge2(width=0.3), size=1) +
	scale_color_manual(values=scheme) +
	theme_bw() + 
	labs(x="", y="log10(Frequency)", color="Barcodes") +
	theme(text = element_text(size=10), 
		axis.text = element_text(size=10), 
		axis.text.x = element_text(angle = 90, margin=margin(t = 0.1, r = 0.1, b = 0, l = 0, unit = "cm")),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		# legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text.x = element_text(size = 10)) +
	facet_grid(.~condition, scales="free_x", space="free")	
dev.off()

test<-rownames(freq[rowSums(log10(freq[,1:7])>-4)>0,])
freqTest.m <- freq.m[freq.m$Var1 %in% test,]

pdf("~/manuscript/plots/dataAll.T47DhighTest.violin.pdf", height=5, width=10)
ggplot(freq.m, aes(x=Var2, y=log10(value))) + 
	geom_violin(fill="grey", color="grey") + 
	geom_point(data=freqTest.m, aes(x=Var2, y=log10(value), color=Var1), position = position_dodge2(width=0.3), size=1) +
	# scale_color_manual(values=scheme) +
	theme_bw() + 
	labs(x="", y="log10(Frequency)", color="Barcodes") +
	theme(text = element_text(size=10), 
		axis.text = element_text(size=10), 
		axis.text.x = element_text(angle = 90, margin=margin(t = 0.1, r = 0.1, b = 0, l = 0, unit = "cm")),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		# legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        strip.text.x = element_text(size = 10),) +
	facet_grid(.~condition, scales="free_x", space="free")	
dev.off()


##############################
####### quantiles ############
##############################

pot<-rowMedians(as.matrix(freq[,1:7]))
names(pot)<-rownames(freq)
pot<-pot[pot>0]
pot.m<-as.data.frame(melt(pot))
pot.m$sample<-"POT"
pot.m$quantile<-""
pot.m[pot.m$value>=quantile(pot)[1] & pot.m$value<=quantile(pot)[2],]$quantile=1
pot.m[pot.m$value>quantile(pot)[2] & pot.m$value<=quantile(pot)[3],]$quantile=2
pot.m[pot.m$value>quantile(pot)[3] & pot.m$value<=quantile(pot)[4],]$quantile=3
pot.m[pot.m$value>quantile(pot)[4] & pot.m$value<=quantile(pot)[5],]$quantile=4

bc.1<-rownames(pot.m[pot.m$quantile==1,])
bc.2<-rownames(pot.m[pot.m$quantile==2,])
bc.3<-rownames(pot.m[pot.m$quantile==3,])
bc.4<-rownames(pot.m[pot.m$quantile==4,])

colors<-c(rep("#AEB4A9", each=length(bc.1)), 
		rep("#E0C1B3", each=length(bc.2)), 
		rep("#D89A9E", each=length(bc.3)), 
		rep("#C37D92", each=length(bc.4)))
names(colors)<-c(bc.1, bc.2, bc.3, bc.4)


df <- data.frame("data" = log10(pot.m$value))
dens <- density(df$data)
qi<-quantile(log10(pot.m$value))
secA<-dens$x >= qi[1] & dens$x <= qi[2]
secB<-dens$x > qi[2] & dens$x <= qi[3]
secC<-dens$x > qi[3] & dens$x <= qi[4]
secD<-dens$x > qi[4] & dens$x <= qi[5]


new_df1 <- data.frame(y = c(dens$x[secA], rev(dens$x[secA])),
                      x = c(-dens$y[secA], rev(dens$y[secA])),
                      z = '#AEB4A9')
new_df2 <- data.frame(y = c(dens$x[secB], rev(dens$x[secB])),
                      x = c(-dens$y[secB], rev(dens$y[secB])),
                      z = '#E0C1B3')
new_df3 <- data.frame(y = c(dens$x[secC], rev(dens$x[secC])),
                      x = c(-dens$y[secC], rev(dens$y[secC])),
                      z = '#D89A9E')
new_df4 <- data.frame(y = c(dens$x[secD], rev(dens$x[secD])),
                      x = c(-dens$y[secD], rev(dens$y[secD])),
                      z = '#C37D92')

pdf("~/manuscript/plots/medPotQuantileT47D.pdf", height=2, width=2)
ggplot(rbind(new_df1, new_df2, new_df3, new_df4), aes(x, y, fill = z)) + 
  geom_polygon() +
  scale_fill_identity() +
  scale_x_continuous(breaks = 0, expand = c(1, 1), labels = 'median(POTs,T0s)', name = '') +
  theme_bw() +
  labs(x="", y="log10(Frequency)")
dev.off()


selected <- c("AI1m1", "AI2m1")
selected <- c("UT.1", "UT.2", "UT.3")
numQuantileBCs<-lapply(selected, function(x){
	c("bc.1"=length(unique(freq.m[freq.m$Var1 %in% bc.1 & freq.m$Var2==x,]$Var1)),
	  "bc.2"=length(unique(freq.m[freq.m$Var1 %in% bc.2 & freq.m$Var2==x,]$Var1)),
	  "bc.3"=length(unique(freq.m[freq.m$Var1 %in% bc.3 & freq.m$Var2==x,]$Var1)),
	  "bc.4"=length(unique(freq.m[freq.m$Var1 %in% bc.4 & freq.m$Var2==x,]$Var1)))
})
names(numQuantileBCs)<-selected

numQuantileBCs.df<-do.call(rbind, numQuantileBCs)
numQuantileBCs.m<-melt(as.matrix(numQuantileBCs.df))
numQuantileBCs.m$quantile<-""
numQuantileBCs.m[grepl("bc.1", numQuantileBCs.m$Var2),]$quantile<-"low"
numQuantileBCs.m[grepl("bc.2", numQuantileBCs.m$Var2),]$quantile<-"lowMid"
numQuantileBCs.m[grepl("bc.3", numQuantileBCs.m$Var2),]$quantile<-"midHigh"
numQuantileBCs.m[grepl("bc.4", numQuantileBCs.m$Var2),]$quantile<-"high"
numQuantileBCs.m$time<-""
numQuantileBCs.m[grepl("1m1", numQuantileBCs.m$Var1),]$time=30.1
numQuantileBCs.m[grepl("2m1", numQuantileBCs.m$Var1),]$time=30.2

colors3<-c("#AEB4A9","#E0C1B3","#D89A9E","#C37D92")
names(colors3)<-c("low", "lowMid", "midHigh", "high")

pdf("~/manuscript/plots/medPotQuantileCountsT47D.pdf", height=2.5, width=2)
ggplot(numQuantileBCs.m, aes(x=time, y=value/1000, color=quantile)) +
pdf("~/manuscriptFeb23/plots/medPotQuantileCountsT47D.UT.pdf", height=2.5, width=2.5)
ggplot(numQuantileBCs.m, aes(x=Var1, y=value/1000, color=quantile)) +
	geom_point(position=position_dodge2(width=0.3), size=2) +
	theme_bw() +
	scale_color_manual(values=colors3) +
	labs(x="Time(days)", y="Surviving barcodes (*1000)", color="Quantile") +
	theme(axis.title=element_text(size=10),
		axis.text=element_text(size=8),
		# legend.position="bottom",
		legend.key.size = unit(0.08, "cm"),
		legend.title = element_text(size=6), 
        legend.text = element_text(size=6))+
	guides(color=guide_legend(nrow=4,byrow=TRUE))
dev.off()


################################################
####### quantiles Awakening violins ############
################################################


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


png("~/manuscript/plots/medPotQuantileViolinsT47Dm1.png", height=360, width=360)
plotFreq(freq[,grepl("m", colnames(freq))])
dev.off()
png("~/manuscript/plots/medPotQuantileViolinsT47Dteps.png", height=360, width=1000)
plotFreq(freq[,grepl("TEP", colnames(freq))][,4:9])
dev.off()
png("~/manuscript/plots/medPotQuantileViolinsT47Daws.png", height=360, width=1000)
plotFreq(freq[,grepl("aw", colnames(freq))][,4:9])
dev.off()


plotPies<-function(nm){
	ip<-freq[,nm]
	ip<-ip[ip>0]
	total<-length(ip)
	ip1<-length(ip[names(ip) %in% bc.1])
	ip2<-length(ip[names(ip) %in% bc.2])
	ip3<-length(ip[names(ip) %in% bc.3])
	ip4<-length(ip[names(ip) %in% bc.4])

	info<-c("low"=ip1/total, "midLow"=ip2/total,"midHigh"=ip3/total, "high"=ip4/total)	
	info.m<-melt(as.matrix(info))
	info.m$value<-info.m$value*100
	# ip.m<-melt(as.matrix(ip))
	# ip.m<-ip.m[ip.m$value>0,]
ggplot(info.m, aes(x="", y=value, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  ylim(0,100.01) +
  coord_polar("y", start=0) +  
  theme_void() +
  scale_fill_manual(values=colors3) +
  labs(fill="Barcodes", title=nm) +
  # facet_grid(.~L1) +
  theme(legend.position="none") +
  	guides(fill=guide_legend(nrow=2,byrow=TRUE))
}

pdf("~/manuscript/plots/medPotQuantilePiesT47D.pdf", height=2, width=2)
lapply(colnames(freq), plotPies)
dev.off()
