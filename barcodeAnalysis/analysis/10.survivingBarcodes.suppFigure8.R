freq<-readRDS("~/barcode_analysis/results/filtered/rds/freq.all.rds")
freq<-freq[,1:60]
freq<-freq[rowSums(freq>0)>0,]

freqAIeTep<-readRDS("/rds/general/project/traditiom/live/work/ele_TEPbis_T47D/results/rds//freq.rds")
samples<-read.table("/rds/general/project/traditiom/live/work/ele_TEPbis_T47D/results/rds//sampleInfo.txt", header=T, stringsAsFactors=F)
freqAIeTep<-freqAIeTep[,18]
freqAIeTep<-as.data.frame(freqAIeTep[freqAIeTep>0])
colnames(freqAIeTep)<-"AI.Epsilon+30"

data<-merge(freq, freqAIeTep, by=0, all=TRUE)
data2<-data
data2[is.na(data2)]=0
dataAll<-data2[,2:62]
rownames(dataAll)<-data2[,1]

getInfo<-function(x){
	counts.bc<-as.data.frame(colSums(x>0))
	colnames(counts.bc)<-"counts"
	info.total<-counts.bc

	info.total$sample<-rownames(info.total)
	info.total$condition<-""
	info.total[grepl("POT", info.total$sample),]$condition<-"POT"
	info.total[grepl("U", info.total$sample),]$condition<-"Untreated"	
	info.total[grepl("AI", info.total$sample),]$condition<-"-E2"
	info.total[grepl("^T", info.total$sample),]$condition<-"T"

	info.total$state<-""
	info.total[grepl("POT", info.total$sample),]$state<-"POT"
	info.total[grepl("U", info.total$sample),]$state<-"Untreated"	
	info.total[grepl("AI7|AI14|T7|T14", info.total$sample),]$state<-"Latency"
	info.total[grepl("AI30|AI60|AI90|T30|T60", info.total$sample),]$state<-"Dormancy"
	info.total[grepl("AI\\.|^T\\.", info.total$sample),]$state<-"Awakening"
	info.total[grepl("\\+30", info.total$sample),]$state<-"TEP"
	info.total
}

info<-getInfo(dataAll)

# info.sel<-info[4:48,]
# rename<-read.table("rds/rename.all.txt", stringsAsFactors=F)
# info.sel$sample2<-""
# info.sel[info.sel$sample %in% rename$V1,]$sample2<-rename[match(info.sel[info.sel$sample %in% rename$V1,]$sample, rename$V1),]$V2

greek<-as.vector(info[grep("AI\\.|^T\\.", info$sample),]$sample)
relabel<-read.table("rds/relabel.txt", sep="\t",header=F, stringsAsFactors=F)
# relabel<-relabel[1:13,]

info$greek<-info$sample
info[match(tolower(relabel$V1), tolower(info$sample)),]$greek<-relabel$V2
info$sample<-factor(info$sample, levels=c("POT.1", "POT.2", "POT.3", "U7.1", "U7.2", "T7.1", "T7.2", "AI7.1", "AI7.2", "U14.1", "U14.2", "T14.1", "T14.2", "AI14.1", "AI14.2", "U30.1", "U30.2", "T30.1", "T30.2", "AI30.1", "AI30.2", "U60.1", "U60.2", "T60.1", "T60.2", "AI60.1", "AI60.2", "U90.1", "U90.2", "AI90.1", "AI90.2", "U120.1", "U120.2", "U150.1", "U150.2", "U170.1", "U170.2", "U170.3", "AI.Alpha", "AI.Beta", "AI.Gamma", "AI.Delta", "AI.Epsilon", "T.Alpha", "T.Beta", "T.Gamma", "T.Delta", "T.Epsilon", "T.Zeta", "T.Theta", "T.Eta", "AI.Alpha+30", "AI.Beta+30", "AI.Gamma+30",  "AI.Epsilon+30","T.Alpha+30", "T.Beta+30", "T.Gamma+30", "T.Delta+30", "T.Epsilon+30", "T.Zeta+30"))

scheme<-c("POT"="#BFBFBF", 
          "Untreated"="#000000", 
          "Latency"="#0084A5", 
          "Dormancy"="#FFDC00", 
          "Awakening"="#B30078", 
          "TEP"="#50157A")
info$state<-factor(info$state, levels=names(scheme))
info$condition<-factor(info$condition, level=c("POT", "Untreated", "T", "-E2"))
plotBCs<-function(scheme){
	ggplot(info, aes(x=sample, y=counts, fill=state)) + 
	geom_bar(stat="identity", position="dodge2", color="black") + 
	# geom_text(stat="identity", aes(label=counts), position=position_stack(vjust = 0.9), vjust=-1) +
	geom_text(stat="identity",size=5, aes(label=counts), position=position_stack(vjust = 0.8), vjust = -1) +
	theme_bw() +  
	labs(x=" ", y="Number of barcodes", fill="State") + 
	theme(axis.text.x = element_text(size=16, angle = 90),
		axis.title = element_text(size=16),
		axis.text.y = element_text(size=14),
		text = element_text(size=12),
		legend.position="none", 
		strip.background=element_rect(fill="white"),
		strip.text.x = element_text(size = 12)) +
	facet_wrap(condition~., scale="free", ncol=1, nrow = 6) +
	scale_fill_manual(values=scheme)
}

# pdf("plots/survivingBCs.recal.pdf", height=10, width=7)
pdf("plots/suppFigure8.pdf", height=11, width=12)
plotBCs(scheme)
dev.off()
