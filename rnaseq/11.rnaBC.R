require("reshape2")
require("ggplot2")
require("stringr")
require("fishplot")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")

## Import samples
setwd("/rds/general/user/hdhiman/home/rnaseq/results/traditiom/bcAlnCounts")
files<-list.files("counts_filtered")
samples<-do.call(rbind, str_split(files, pattern="\\.", n=3))[,1]

samples<-samples[!samples %in% c("R51", "R52", "R53", "R54", "R55", "R56", "R57", "R58")]

info.m<-readRDS("../rds/sample_details.rds")
info.m<-info.m[info.m$ID %in% samples,]

## Get raw counts and number of barcodes
counts.list<-lapply(samples, function(x){read.table(paste("counts_filtered",paste(x, "counts.txt", sep="."), sep="/"), header=F)})
names(counts.list)<-samples

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts)
order.samples<-c("POT.1", "POT.2", "POT.3", "UT30.1", "UT30.2", "UT120.1", "UT120.2", "UT170.1", "UT170.2", "UT170.3",
  "T7.1", "T7.2", "T14.1", "T14.2", "AI7.1", "AI7.2", "AI14.1", "AI14.2",  "T30.1", "T30.2", "T60.1", "T60.2", "AI30.1", "AI30.2", "AI60.1", "AI60.2", "AI90.1", "AI90.2",  
  "T.alpha",  "T.beta", "T.gamma", "T.delta",  "T.epsilon", "T.zeta", "T.theta", "T.eta", "AI.alpha", "AI.beta",  "AI.gamma", "AI.delta", "AI.epsilon", 
  "T.alpha_TEP", "T.beta_TEP","T.gamma_TEP", "T.delta_TEP","T.epsilon_TEP", "T.zeta_TEP", "AI.alpha_TEP", "AI.beta_TEP","AI.gamma_TEP", "AI.epsilon_TEP")
colnames(counts)<-info.m$Sample
counts<-counts[,order.samples]
write.table(counts, "../rds/rnaBCcounts.txt", sep="\t", quote=F, row.names=T, col.names=T)
saveRDS(counts, "../rds/rnaBCcounts.rds")


##################################
counts<-readRDS("../rds/rnaBCcounts.rds")
samples.sel<-c(info.m$Sample[grep("^U|T\\.|AI\\.", info.m$Sample)])
samples.sel[grepl("U", samples.sel)] <- c("UT30.1",  "UT30.2",  "UT120.1", "UT120.2", "UT170.1", "UT170.2", "UT170.3")
data.sel<-counts[,colnames(counts) %in% samples.sel]

total<-colSums(data.sel)
freq<-as.data.frame(t(t(data.sel)/total))
freq<-freq[rowSums(freq>0)>0,]

#### Heatmap #### 
bcs<-rownames(freq)
sel.bcs<-c("bc_8955485", "bc_2148708", "bc_3551399", "bc_8538010", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_2845296", "bc_9623528", "bc_5762012", "bc_5423601", "bc_3699478", "bc_8533083", "bc_2864628", "bc_9461258", "bc_5786343", "bc_1152635")
presentGenome<-bcs[bcs %in% sel.bcs]
highRNA<-rownames(freq[rowMaxs(as.matrix(freq))>=0.1,])
missing<-sel.bcs[!(sel.bcs %in% bcs)]
# [1] "bc_3551399" "bc_8533083"
missing.df<-t(data.frame("bc_3551399"=rep(0, 33), "bc_8533083"=rep(0, 33)))
missing.df <- as.data.frame(missing.df)
colnames(missing.df)<-colnames(freq)
sel.all<-na.omit(unique(c(presentGenome, highRNA)))

extras<-highRNA[!(highRNA %in% sel.bcs)]
# [1] "bc_1788019" "bc_2891982" "bc_5369085" "bc_5372795" "bc_7309655"
# [6] "bc_9565975" "bc_9702545"
order<-unique(rev(c("bc_9565975", "bc_2891982", "bc_5369085", "bc_7309655", "bc_8955485", "bc_2148708", "bc_8538010", "bc_1788019", "bc_2823514", "bc_8588957", "bc_9601771", "bc_3637132", "bc_5197405", "bc_2881711", "bc_9623528", "bc_9702545",  "bc_5423601", "bc_3699478", "bc_2864628", "bc_9461258", "bc_5786343","bc_5372795", "bc_1152635", "bc_2845296", "bc_5762012", "bc_3551399", "bc_8533083")))

freq.sel<-freq[sel.all,]
freq.sel<-rbind(freq.sel, missing.df)


freq.sel.m<-melt(as.matrix(freq.sel))
info.split<-as.data.frame(str_split(freq.sel.m$Var2, "\\.|\\_", simplify=TRUE), stringsAsFactors=F)

freq.sel.m$Var2<-factor(freq.sel.m$Var2, levels=c("POT.1", "POT.2", "POT.3", "UT30.1", "UT30.2", "UT120.1", "UT120.2", "UT170.1", "UT170.2", "UT170.3", "AI.alpha","AI.alpha_TEP", "AI.beta","AI.beta_TEP",  "AI.gamma", "AI.gamma_TEP", "AI.delta", "AI.epsilon",  "AI.epsilon_TEP", 
													"T.alpha","T.alpha_TEP", "T.beta","T.beta_TEP", "T.gamma","T.gamma_TEP", "T.delta","T.delta_TEP", "T.epsilon", "T.epsilon_TEP", "T.zeta", "T.zeta_TEP", "T.eta",  "T.theta"))
freq.sel.m$condition<-info.split[,1]
freq.sel.m[grep("UT", freq.sel.m$Var2),]$condition<-"UT"
freq.sel.m[freq.sel.m$condition=="AI",]$condition<-"-E2"
freq.sel.m$condition<-factor(freq.sel.m$condition, levels=c("POT", "UT", "T", "-E2"))

freq.sel.m$Var1<-factor(freq.sel.m$Var1, levels=order)
freq.sel.m$value[freq.sel.m$value==0]<-NA

pdf("../plots/rnaBC.heatmap2.legend.pdf", height=6, width=8)
ggplot(freq.sel.m, aes(x=Var2, y=Var1, fill=value)) + 
	geom_tile() + 
	geom_tile(data=freq.sel.m[freq.sel.m$Var1 %in% sel.bcs & freq.sel.m$value>0 & !is.na(freq.sel.m$value),], aes(x=Var2, y=Var1), color="grey") + 
	# geom_point(data=freq.sel.m[freq.sel.m$Var1 %in% present.not.sel,], aes(x=Var2, y=Var1), color="grey") +
	# geom_point(data=freq.sel.m[freq.sel.m$Var1 %in% sel.not.present,], aes(x=Var2, y=Var1), color="black") + 
	scale_fill_viridis(option="A", begin = 0.25, end = 1, direction=-1, na.value ="white") + 
	theme_bw() + 
	theme(axis.text.x = element_text(size=10, angle = 90, vjust=0.5, hjust=1)) +
	labs(x="Samples", y="Barcodes", fill="Frequency > 9.5%") +
	theme(text = element_text(size=12), strip.background = element_rect(fill="#FFFFFF"), legend.position="bottom") +
	# scale_x_discrete(labels=relabel) +
	facet_grid(. ~ condition, scale="free_x", space="free_x")
dev.off()

