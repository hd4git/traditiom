require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("stringr")
require("viridis")
require("ggpubr")

setwd("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered")
counts.all <- readRDS("rds/data.recal.rds")
final.color<-c("common"="gray90", "common.1"="gray80",  "common.2"="gray70",  "common.3"="gray60",  "common.4"="gray50", "common.5"="gray40",  "common.6"="gray30",  "common.7"="gray20", "common.8"="gray10", 
				"POT"="#0084A5", "Untreated"="#0084A5", "Relapse.UT"="#004152", "Latency"="#0084A5", "Dormancy"="#FFDC00", "Awakening"="#B30078", "Relapse"="#50157A", 
				"UT"="#0084A5", "-E2"="#B30078", "TAM"="#B30078")
rename<-read.table("rds/sample.info.all.t0.txt", header=T, stringsAsFactors=F)


#################################################
#############  Get counts #######################
#################################################
getInfo<-function(x, min){
	counts.bc<-as.data.frame(colSums(x>=min))
	colnames(counts.bc)<-"counts"
	info.total<-counts.bc

	info.total$sample<-rownames(info.total)
	info.total<-info.total[-grep("Dead", info.total$sample), ]

	info.sep<-as.data.frame(str_split(rownames(info.total), "_", simplify=T), stringsAsFactors=F)
	info.sep2<-as.data.frame(str_split(info.sep[,1], "\\.", simplify=T), stringsAsFactors=F)
	info.sep3<-as.data.frame(str_split(info.sep[,2], "\\.", simplify=T), stringsAsFactors=F)

	info.total$condition<-info.sep2[,1]
	info.total[info.total$condition=="T0",]$condition<-"POT"
	info.total$rep<-info.sep2[,2]

	info.total$duration<-info.sep[,2]
	# info.total$duration[1:3]<-"T0"
	info.total$duration2<-paste(info.total$duration, info.total$rep, sep=".")

	info.total$condition<-factor(info.total$condition, levels=c("POT","UT", "AI", "TAM", "BARCD", "BARCE"))
	info.total$duration<-factor(info.total$duration, levels=c("0D","7D", "14D", "30D", "M2", "M3", "M4", "LA", "R", "TEP", "A.TEP", "DH.R", "DH.TEP", "BIS.R", "Bulk.7D",  "Low.7D", "LowMid.7D", "MidHi.7D",  "Hi.7D" ,  "Bulk.1M",  "Low.1M" ,   "LowMid.1M", "MidHi.1M" , "Hi.1M"))
	info.total$duration2<-factor(info.total$duration2, levels=c("0D.1", "0D.2", "0D.3", "7D.1", "7D.2", "14D.1", "14D.2", "30D.1", "30D.2", "M2.1", "M2.2", "M3.1", "M3.2", "M4.1", "M4.2", "LA.1", "R.1", "R.2", "R.3", "R.4", "TEP.1", "TEP.2", "TEP.3", "TEP.4", "A.TEP.1", "A.TEP.2", "DH.R.1", "DH.R.2", "DH.TEP.1", "DH.TEP.2", "BIS.R.1","BIS.R.2", "BIS.R.3", "Bulk.7D.1", "Bulk.7D.2", "Low.7D.1", "Low.7D.2", "LowMid.7D.1", "LowMid.7D.2", "MidHi.7D.1", "MidHi.7D.2", "Hi.7D.1", "Hi.7D.2", "Bulk.1M.1", "Bulk.1M.2", "Low.1M.1", "Low.1M.2", "LowMid.1M.1", "LowMid.1M.2", "MidHi.1M.1", "MidHi.1M.2", "Hi.1M.1", "Hi.1M.2"))
	info.total$min=min
	
	info.sel<-info.total[1:60,]

	info.sel$sample2<-""
	info.sel[info.sel$sample %in% rename$sample,]$sample2<-rename[match(info.sel[info.sel$sample %in% rename$sample,]$sample, rename$sample),]$sample.renamed
	info.sel$condition<-as.vector(info.sel$condition)
	info.sel[info.sel$condition=="AI",]$condition<-"-E2"
	info.sel$condition<-factor(info.sel$condition, levels=c("POT","UT", "TAM", "-E2", "BARCD", "BARCE"))

	info.sel$state<-""
	info.sel[1:3,]$state<-"POT"
	info.sel[4:15,]$state<-"Latency"
	info.sel[c(16:29, 32:35),]$state<-"Dormancy"
	info.sel[c(30:31, 36:39, 42:48),]$state<-"Awakening"
	info.sel[49:60,]$state<-"Relapse"
	info.sel[grep("U", info.sel$sample2), ]$state<-"Untreated"

	str.info<-as.data.frame(str_split(info.sel$sample2, "\\.", simplify=T), stringsAsFactors=F)
	info.sel$sample3<-str.info[,1]
	# info.sel[info.sel$sample3=="U150",]$sample3<-"UT"
	info.sel$rep2<-str.info[,2]
	info.sel$rep2<-tolower(info.sel$rep2)
	info.sel$rep2<-factor(info.sel$rep2, levels=tolower(c("1", "2", "3", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Alpha+30", "Beta+30", "Gamma+30", "Delta+30", "Epsilon+30", "Zeta+30")))
	info.sel$rep3<-info.sel$rep2
	tep.str<-as.data.frame(str_split(info.sel$rep3, "\\+", simplify=T), stringsAsFactors=F)
	info.sel$rep3<-tep.str[,1]
	info.sel[grepl("U170", info.sel$sample2), ]$state<-"Relapse.UT"
	info.sel$sample2<-factor(info.sel$sample2, levels=c("POT.1", "POT.2", "POT.3", "U7.1", "U7.2", "T7.1", "T7.2", "AI7.1", "AI7.2", "U14.1", "U14.2", "T14.1", "T14.2", "AI14.1", "AI14.2", "U30.1", "U30.2", "T30.1", "T30.2", "AI30.1", "AI30.2", "U60.1", "U60.2", "T60.1", "T60.2", "AI60.1", "AI60.2", "U90.1", "U90.2", "AI90.1", "AI90.2", "U120.1", "U120.2", "U150.1", "U150.2", "AI.Alpha", "AI.Beta", "AI.Gamma", "AI.Delta", "AI.Epsilon", "T.Alpha", "T.Beta", "T.Gamma", "T.Delta", "T.Epsilon", "T.Zeta", "T.Eta", "T.Theta", "U170.1", "U170.2", "U170.3", "AI.Alpha+30", "AI.Beta+30", "AI.Gamma+30", "T.Alpha+30", "T.Beta+30", "T.Gamma+30", "T.Delta+30", "T.Epsilon+30", "T.Zeta+30"))
	info.sel$state<-factor(info.sel$state, levels=c("POT", "Untreated", "Relapse.UT", "Latency", "Dormancy", "Awakening", "Relapse"))
	info.sel
}

info1<-getInfo(counts.all, 1)
info3<-getInfo(counts.all, 3)
info5<-getInfo(counts.all, 5)
info10<-getInfo(counts.all, 10)
info15<-getInfo(counts.all, 15)
info<-rbind(info1, info3, info5, info10, info15)

plotBCsBar<-function(mat, scheme){
	ggplot(mat, aes(x=sample2, y=counts, fill=state)) + 
	geom_bar(stat="identity", position="dodge2", color="black") + 
	geom_text(stat="identity", aes(label=counts), position=position_stack(vjust = 0.98), vjust = -1) +
	theme_bw() +  
	labs(x=" ", y="Number of barcodes", fill="Condition") + 
	theme(text = element_text(size=16),
		axis.text = element_text(size=16),
		axis.text.x = element_text(size=16, angle = 45, margin=margin(t = 1, r = 0.75, b = 0, l = 0, unit = "cm")),  
		legend.position="none", strip.background=element_rect(fill="white")) +
	facet_grid(min ~ condition, scale="free_x", space="free_x") +
	scale_fill_manual(values=scheme)
}


plotBCsTEP<-function(mat, scheme){
	ggplot(mat, aes(x=condition, y=counts, color=state, label=rep3)) + 
	geom_text(aes(condition,counts,label=rep3, color=state), parse=T, position=position_dodge2(width =0.7)) +
	theme_bw() +  
	labs(x="Condition", y="log10(Number of barcodes)", color="Condition", shape="Flasks") + 
	theme(
		# text = element_blank(),
			text = element_text(size=12), 
			legend.position="none", 
			strip.background=element_rect(fill="white")) +
	# facet_grid(. ~ min, scale="free_x", space="free_x") +
	scale_color_manual(values=scheme) 
}



plotBC.bar<-plotBCsBar(info, final.color)
pdf("plots/plotBC.awake.diffMin.bar.pdf", height=18, width=30)
plotBC.bar
dev.off()

dataToPlot<-info[info$min==10 & (info$state %in% c("POT", "Awakening") | grepl("U150", info$sample2)),]
dataToPlot$counts<-log10(dataToPlot$counts)
plotBC.awake<-plotBCsTEP(dataToPlot, final.color)
pdf("plots/fig3a.survivingBCs.awake.pdf", height=3, width=5)
plotBC.awake
dev.off()


dataToPlot<-info[info$min==10 & (info$state %in% c("POT", "Awakening", "Relapse") | grepl("U150|U170", info$sample2)),]
dataToPlot$counts<-log10(dataToPlot$counts)
plotBC.awake.TEP<-plotBCsTEP(dataToPlot, final.color)
pdf("plots/fig3aSupp.survivingBCs.awake+TEP.pdf", height=3, width=7)
plotBC.awake.TEP
dev.off()





