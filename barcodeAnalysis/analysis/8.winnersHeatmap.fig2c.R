require("RColorBrewer")
require("reshape2")
require("ggplot2")
require("stringr")
require("stringr")
require("viridis")
require("ggpubr")

#############################
##### Upto 2 months #########
#############################
setwd("~/barcode_analysis/results/filtered/")

freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")
samples<-colnames(freq)
samples.sel<-c(samples[grep('U150|U170', samples)], 
				samples[grep('AI30|AI60|AI90|AI\\.|T30|T60|^T\\.', samples)])
samples.sel<-samples.sel[-grep("_", samples.sel)]
freq<-freq[,samples.sel]

# data.sel<-data[,samples.sel]
# total<-colSums(data.sel)
# freq<-as.data.frame(t(t(data.sel)/total))

freq <- freq[rowSums(freq>0)>0,]

####################################################################
#### Heatmap #### 
freq.min10<-freq[rowMaxs(as.matrix(freq[,-8]))>=0.095 | rownames(freq)=="bc_2810486",]
# saveRDS(freq.min10, "~/barcode_analysis/results/filtered/rds/winners.freq.min10.rds")
# freq.min10<-readRDS("~/barcode_analysis/results/filtered/rds/winners.freq.min10.rds")

freq.min10[freq.min10 == 0] <- NA


freq.min10.m<-melt(as.matrix(freq.min10))
info.split<-as.data.frame(str_split(freq.min10.m$Var2, "\\.|\\_", simplify=TRUE), stringsAsFactors=F)

# freq.min10.m$new<-names.samples[match(freq.min10.m$Var2, names.samples$sample),]$sample.renamed
freq.min10.m$Var2<-factor(freq.min10.m$Var2, levels=c("POT.1", "POT.2", "POT.3","U30.1", "U30.2","U60.1", "U60.2","U90.1", "U90.2", "U150.1", "U150.2","U170.1", "U170.2", "U170.3", 
													"AI30.1", "AI30.2","AI60.1", "AI60.2","AI90.1", "AI90.2","AI.Alpha","AI.Alpha+30", "AI.Beta","AI.Beta+30",  "AI.Gamma", "AI.Gamma+30", "AI.Delta", "AI.Epsilon",  
													"T30.1", "T30.2","T60.1", "T60.2","T.Alpha","T.Alpha+30", "T.Beta","T.Beta+30", "T.Gamma","T.Gamma+30", "T.Delta","T.Delta+30", "T.Epsilon", "T.Epsilon+30", "T.Zeta", "T.Zeta+30", "T.Eta",  "T.Theta"))
freq.min10.m$condition<-info.split[,1]
freq.min10.m$flask<-tolower(as.data.frame(str_split(freq.min10.m$Var2, "\\.", simplify=T), stringsAsFactors=F)[,2])
# freq.min10.m[freq.min10.m$condition=="T0",]$condition<-"POT"
freq.min10.m$phase <- ""
freq.min10.m[grepl("T30|T60", freq.min10.m$Var2),]$phase<-"T.Dormancy"
freq.min10.m[grepl("AI30|AI60|AI90", freq.min10.m$Var2),]$phase<-"-E2.Dormancy"
freq.min10.m[grepl("^AI\\.", freq.min10.m$Var2),]$phase<-"-E2"
freq.min10.m[grepl("^T\\.", freq.min10.m$Var2),]$phase<-"T"
freq.min10.m[grepl("^U", freq.min10.m$Var2),]$phase<-"UT"



freq.min10.m[grepl("^T", freq.min10.m$condition),]$condition<-"T"
freq.min10.m[grepl("AI", freq.min10.m$condition),]$condition<-"-E2"
freq.min10.m[grepl("^U", freq.min10.m$condition),]$condition<-"UT"

freq.min10.m$condition<-factor(freq.min10.m$condition, levels=c("UT", "T", "-E2"))
dormSamples <- freq.min10.m[grepl("AI30|AI60|AI90|T30|T60", freq.min10.m$Var2),]
freq.min10.m[grepl("AI30|AI60|AI90|T30|T60", freq.min10.m$Var2),]$flask<-substr(dormSamples$Var2, (nchar(as.vector(dormSamples$Var2))+1)-4, nchar(as.vector(dormSamples$Var2)))


freq.min10.m[grepl("U170.1", freq.min10.m$Var2),]$flask<-"1.TEP"
freq.min10.m[grepl("U170.2", freq.min10.m$Var2),]$flask<-"2.TEP"
freq.min10.m[grepl("U170.3", freq.min10.m$Var2),]$flask<-"3.TEP"
selBCs <- c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635","bc_2810486")
freq.min10.m <- freq.min10.m[freq.min10.m$Var1 %in% selBCs,]
freq.min10.m$Var1 <- factor(x = freq.min10.m$Var1,
                               levels = rev(selBCs))
freq.min10.m$flask<-factor(freq.min10.m$flask, levels=c("1","1.TEP", "2", "2.TEP", "3", "3.TEP","30.1", "30.2", "60.1", "60.2", "90.1", "90.2", "alpha", "alpha+30", "beta", "beta+30", "gamma", "gamma+30", "delta", "delta+30", "epsilon", "epsilon+30", "zeta", "zeta+30", "eta", "theta"))
relabel<-c("1"="1", "2"="2", "3"="3", "1.TEP"="1+20", "2.TEP"="2+20", "3.TEP"="3+20", "30.1"="30.1", "30.2"="30.2", "60.1"="60.1", "60.2"="60.2", "90.1"="90.1", "90.2"="90.2", 
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
# singles<-c("AI.Delta", "T.Eta", "T.Theta", "U170.3")
singles<-c("T.Theta", "U170.3")


tepbisAIeInfo<-readRDS("~/barcode_analysis/results/filtered/rds/tepbisAIeInfo.rds")
tepbisAIeInfo<-tepbisAIeInfo[tepbisAIeInfo$Var1 %in% c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635","bc_2810486"),]
tepbisAIeInfo$phase <- "-E2"
freq.min10.m<-rbind(freq.min10.m, tepbisAIeInfo)

freq.min10.m$phase <- factor(freq.min10.m$phase, levels=c("UT", "T.Dormancy", "T", "-E2.Dormancy", "-E2"))

# pdf("plots/Figure4a.winnersHeatmap.pdf", height=2, width=5)
# pdf("plots_MS/Figure3d.winnersHeatmap.pdf", height=2.5, width=5)
pdf("~/manuscriptFeb23/plots/figure2c.winnersHeatmap.pdf", height=2.5, width=5)
ggplot(freq.min10.m[!freq.min10.m$Var2 %in% singles ,], aes(x=flask, y=Var1, fill=value)) + 
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
	facet_grid(. ~ condition+phase, scale="free_x", space="free_x")
dev.off()


