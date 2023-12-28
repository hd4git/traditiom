require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("viridis")
require("ggpubr")
require("matrixStats")
require("PerformanceAnalytics") 
require("ggrepel")

data <- readRDS("~/barcode_analysis/results/filtered/rds/data.recal.rds")
names.samples<-read.table("~/barcode_analysis/results/filtered/rds/sample.info.all.t0.txt", header=TRUE, stringsAsFactors=F)

data<-data[,-c(22,23)]
data<-data[rowSums(data>10)>0,]

colnames(data)[1:60]<-names.samples$sample.renamed[1:60]
samples<-colnames(data)

total<-colSums(data)
freq<-as.data.frame(t(t(data)/total))
freq<-freq[rowSums(freq>0)>0,]

saveRDS(freq, "~/barcode_analysis/results/filtered/rds/freq.all.rds")
write.table(freq,file ="rds/freq.all.txt", sep="\t", quote=F, row.names=T, col.names=T)

freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")
#############################
##### Upto 2 months##########
#############################
setwd("~/barcode_analysis/results/filtered/")
freq<-freq[,c(1:29, 32,33)]
freq<-freq[rowSums(freq>0)>0,]
freq.min1<-freq[rowMaxs(as.matrix(freq))>=0.01,]
freq.min1[freq.min1 == 0] <- NA

freq.min1.m<-melt(as.matrix(freq.min1))
freq.min1.m$Var2<-factor(freq.min1.m$Var2, levels=colnames(freq))
freq.min1.m$condition<-"POT"
freq.min1.m[grep("U", freq.min1.m$Var2), ]$condition<-"UT"
freq.min1.m[grep("AI", freq.min1.m$Var2), ]$condition<-"-E2"
freq.min1.m[grep("^T", freq.min1.m$Var2), ]$condition<-"T"
freq.min1.m$condition<-factor(freq.min1.m$condition, levels=c("POT", "UT", "T", "-E2"))

bcOrder1<-rev(c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_7309655", "bc_5430295", "bc_5162157", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635"))
bcOrder2<-rev(c("bc_5122096", "bc_159071", "bc_7099037", "bc_7704123", "bc_5189530", "bc_150411", "bc_9683153", "bc_8749691", "bc_2065240", "bc_1758787", "bc_7305702", "bc_2127641", "bc_9625985", "bc_3974028", "bc_2012775", "bc_5460077", "bc_3699478", "bc_6864523"))
bcOrder3<-rev(c("bc_9683153", "bc_9652614", "bc_8749691", "bc_7704123", "bc_7488963", "bc_7099037", "bc_6864523", "bc_5498319", "bc_5189530", "bc_5122096", "bc_3974028", "bc_3645673", "bc_364027", "bc_2184209", "bc_2127641", "bc_2012775", "bc_159071", "bc_150411"))

winnersAll<-as.vector(unique(freq.min1.m$Var1))
bcCommon<-unique(c(winnersAll[winnersAll %in% unique(c(bcOrder1, bcOrder2, bcOrder3))], bcOrder1))

bcOrder<-unique(c(winnersAll[!winnersAll %in% unique(c(bcOrder1, bcOrder2, bcOrder3))], bcOrder1))
freq.min1.m$Var1<-factor(freq.min1.m$Var1, levels=bcOrder)


pdf("plots/winnersHeatmap.untilDormancy.pdf", height=6, width=8)
ggplot(freq.min1.m, aes(x=Var2, y=Var1, fill=value)) + 
	geom_tile() + 
	geom_tile(data=freq.min1.m[freq.min1.m$Var1 %in% bcCommon,], aes(x=Var2, y=Var1), color="red") + 	
	scale_fill_viridis(option="A", begin = 0.25, end = 1, direction=-1, na.value ="white") + 
	theme_bw() + 
	theme(axis.text.x = element_text(size=10, angle = 90)) +
	labs(x="Samples", y="Barcodes", fill="Frequency > 0.01") +
	theme(text = element_text(size=12), 
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom") +
	# scale_x_discrete(labels=relabel) 
	facet_grid(. ~ condition, scale="free_x", space="free_x")
dev.off()



#############################
freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")
samples<-colnames(freq)
samples.sel<-c(samples[grep("POT", samples)], 
				samples[grep('U150|U170', samples)], 
				samples[grep('AI\\.|^T\\.', samples)])
samples.sel<-samples.sel[-grep("_", samples.sel)]
freq<-freq[,samples.sel]
# data.sel<-data[,samples.sel]

# total<-colSums(data.sel)
# freq<-as.data.frame(t(t(data.sel)/total))
freq<-freq[rowSums(freq>0)>0,]

####################################################################
#### Heatmap #### 
freq.min10<-freq[rowMaxs(as.matrix(freq[,-8]))>=0.095 | rownames(freq)=="bc_2810486",]
saveRDS(freq.min10, "~/barcode_analysis/results/filtered/rds/winners.freq.min10.rds")
freq.min10<-readRDS("~/barcode_analysis/results/filtered/rds/winners.freq.min10.rds")
otter.dendro <- as.dendrogram(hclust(d = dist(x = freq.min10), method="ward.D2"))
otter.order <- order.dendrogram(otter.dendro)

freq.min10[freq.min10 == 0] <- NA


freq.min10.m<-melt(as.matrix(freq.min10))
info.split<-as.data.frame(str_split(freq.min10.m$Var2, "\\.|\\_", simplify=TRUE), stringsAsFactors=F)

# freq.min10.m$new<-names.samples[match(freq.min10.m$Var2, names.samples$sample),]$sample.renamed
freq.min10.m$Var2<-factor(freq.min10.m$Var2, levels=c("POT.1", "POT.2", "POT.3", "U150.1", "U150.2","U170.1", "U170.2", "U170.3", 
													"AI.Alpha","AI.Alpha+30", "AI.Beta","AI.Beta+30",  "AI.Gamma", "AI.Gamma+30", "AI.Delta", "AI.Epsilon",  
													"T.Alpha","T.Alpha+30", "T.Beta","T.Beta+30", "T.Gamma","T.Gamma+30", "T.Delta","T.Delta+30", "T.Epsilon", "T.Epsilon+30", "T.Zeta", "T.Zeta+30", "T.Eta",  "T.Theta"))
freq.min10.m$condition<-info.split[,1]
freq.min10.m$flask<-tolower(as.data.frame(str_split(freq.min10.m$Var2, "\\.", simplify=T), stringsAsFactors=F)[,2])
# freq.min10.m[freq.min10.m$condition=="T0",]$condition<-"POT"
freq.min10.m[freq.min10.m$condition=="AI",]$condition<-"-E2"
freq.min10.m[grepl("U", freq.min10.m$condition),]$condition<-"UT"

freq.min10.m$condition<-factor(freq.min10.m$condition, levels=c("POT", "UT", "T", "-E2"))

freq.min10.m[grepl("U170.1", freq.min10.m$Var2),]$flask<-"1.TEP"
freq.min10.m[grepl("U170.2", freq.min10.m$Var2),]$flask<-"2.TEP"
freq.min10.m[grepl("U170.3", freq.min10.m$Var2),]$flask<-"3.TEP"

freq.min10.m$Var1 <- factor(x = freq.min10.m$Var1,
                               levels = rev(c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635","bc_2810486")))
freq.min10.m$flask<-factor(freq.min10.m$flask, levels=c("1","1.TEP", "2", "2.TEP", "3", "3.TEP", "alpha", "alpha+30", "beta", "beta+30", "gamma", "gamma+30", "delta", "delta+30", "epsilon", "epsilon+30", "zeta", "zeta+30", "eta", "theta"))
relabel<-c("1"="1", "2"="2", "3"="3", "1.TEP"="1+20", "2.TEP"="2+20", "3.TEP"="3+20", 
			"alpha"=expression(alpha),
			"beta"=expression(beta),
			"gamma"=expression(gamma),
			"delta"=expression(delta),
			"epsilon"=expression(epsilon),
			"zeta"=expression(zeta),
			"eta"=expression(eta),
			"theta"=expression(theta),
			"alpha+30"=expression(alpha+30),
			"beta+30"= expression(beta+30),
			"gamma+30"=expression(gamma+30),
			"delta+30"=expression(delta+30),
		  "epsilon+30"=expression(epsilon+30),
			 "zeta+30"=expression(zeta+30),
			  "eta+30"=expression(eta+30),
			"theta+30"=expression(theta+30))
singles<-c("AI.Delta", "T.Eta", "T.Theta", "U170.3")


tepbisAIeInfo<-readRDS("~/barcode_analysis/results/filtered/rds/tepbisAIeInfo.rds")
tepbisAIeInfo<-tepbisAIeInfo[tepbisAIeInfo$Var1 %in% c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635","bc_2810486"),]
freq.min10.m<-rbind(freq.min10.m, tepbisAIeInfo)



# pdf("plots/Figure4a.winnersHeatmap.pdf", height=2, width=5)
pdf("plots_MS/Figure3d.winnersHeatmap.pdf", height=2.5, width=5)
ggplot(freq.min10.m[!freq.min10.m$Var2 %in% singles,], aes(x=flask, y=Var1, fill=value)) + 
	geom_tile() + 
	scale_fill_viridis(option="A", begin = 0.25, end = 1, direction=-1, na.value ="white") + 
	theme_bw() + 
	theme(axis.text.x = element_text(size=10, angle = 90)) +
	labs(x="", y="", fill="Frequency > 9.5%") +
	theme(text = element_text(size=12), 
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="none",
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank()) +
	scale_x_discrete(labels=relabel) +
	facet_grid(. ~ condition, scale="free_x", space="free_x")
dev.off()


pdf("plots_MS/Figure3a.winnersHeatmap.pdf", height=2, width=5)
ggplot(freq.min10.m[-grep('\\+|U170', freq.min10.m$Var2),], aes(x=flask, y=Var1, fill=value)) + 
	geom_tile() + 
	scale_fill_viridis(option="A", begin = 0.25, end = 1, direction=-1, na.value ="white") + 
	theme_bw() + 
	theme(axis.text.x = element_text(size=10)) +
	labs(x="", y="", fill="Frequency > 9.5%") +
	theme(axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		text = element_text(size=12), 
		strip.background = element_rect(fill="#FFFFFF"),
		legend.position="none") +
	scale_x_discrete(labels=relabel) +
	facet_grid(. ~ condition, scale="free_x", space="free_x")
dev.off()






