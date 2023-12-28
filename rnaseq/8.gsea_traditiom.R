require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("ggpubr")
require("tximport")
require("ggrepel")
require("pheatmap")
require("viridis")

files<-list.files("gsea")
str_files<-as.data.frame(str_split(files, "_|\\.", simplify=T), stringsAsFactors=F)

samples<-paste(str_files$V5, str_files$V6, sep=".")

data<-lapply(files, function(x){
	info<-read.table(paste("gsea", x, sep="/"), header=T, sep="\t", stringsAsFactors=F)
	strx<-as.data.frame(str_split(x, "_|\\.", simplify=T), stringsAsFactors=F)
	# info.sig<-info[info$FDR.q.val<=0.25, c(1, 4, 8)]
	info.sig<-info[info$FDR.q.val<=0.10, c(1, 4, 8)]
	info.sig$sample<-paste(strx$V5, strx$V6, sep=".")
	info.sig$NAME<-str_replace_all(info.sig$NAME, "HALLMARK_", "")
	info.sig
})

data.all<-do.call(rbind, data)
data.all$sample<-factor(data.all$sample, levels=c("neg.ai30", "neg.ai60", "neg.ai90", "neg.t30", "neg.t60", "neg.ut30", "neg.ut120", "neg.ut170", "pos.ai30", "pos.ai60", "pos.ai90", "pos.t30", "pos.t60", "pos.ut30", "pos.ut120", "pos.ut170"))

data.mat<-acast(data.all, NAME ~ sample, value.var="FDR.q.val")
data.mat[is.na(data.mat)]=0

neg.mat<-data.mat[,1:8]
otter.dendro.neg <- as.dendrogram(hclust(d = dist(x = neg.mat), method="complete"))
otter.order.neg <- order.dendrogram(otter.dendro.neg)
order.neg<-rownames(neg.mat[otter.order.neg,])

pos.mat<-data.mat[,8:16]
otter.dendro.pos <- as.dendrogram(hclust(d = dist(x = pos.mat), method="complete"))
otter.order.pos <- order.dendrogram(otter.dendro.pos)
order.pos<-rownames(pos.mat[otter.order.pos,])

data.pos<-data.all[grep("pos", data.all$sample),]
data.neg<-data.all[grep("neg", data.all$sample),]

data.pos$NAME <- factor(x = data.pos$NAME,
                            levels = order.pos)
data.neg$NAME <- factor(x = data.neg$NAME,
                            levels = order.neg)


pdf("plots/gsea_pot.ordered.pdf", height=9, width=15)
pos<-ggplot(data.pos, aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="magma",  begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
neg<-ggplot(data.neg, aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="viridis", begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
ggarrange(neg, pos, nrow=1, ncol=2, common.legend = F)
dev.off()


pdf("plots/gsea_pot.pdf", height=9, width=15)
pos<-ggplot(data.all[grep("pos", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="magma",  begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
neg<-ggplot(data.all[grep("neg", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="viridis", begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
ggarrange(neg, pos, nrow=1, ncol=2, common.legend = F)
dev.off()



##############################

files<-list.files("gsea_edger")
str_files<-as.data.frame(str_split(files, "_|\\.", simplify=T), stringsAsFactors=F)

samples<-paste(paste(str_files$V5, str_files$V6, sep="_"), str_files$V7, sep=".")

data<-lapply(files, function(x){
	info<-read.table(paste("gsea_edger", x, sep="/"), header=T, sep="\t", stringsAsFactors=F)
	strx<-as.data.frame(str_split(x, "_|\\.", simplify=T), stringsAsFactors=F)
	info.sig<-info[info$FDR.q.val<=0.25, c(1, 4, 8)]
		# info.sig<-info[info$FDR.q.val<=0.10, c(1, 4, 8)]
		info.sig$sample<-paste(paste(strx$V5, strx$V6, sep="_"), strx$V7, sep=".")
		info.sig$NAME<-str_replace_all(info.sig$NAME, "HALLMARK_", "")
	info.sig
})

data.all<-do.call(rbind, data)
data.all$sample<-factor(data.all$sample, levels=c("neg_ai.alpha", "neg_ai.beta", "neg_ai.gamma", "neg_t.alpha", "neg_t.beta", "neg_t.gamma", "neg_t.delta", "neg_t.epsilon", "neg_t.zeta", "pos_ai.alpha", "pos_ai.beta", "pos_ai.gamma", "pos_t.alpha", "pos_t.beta", "pos_t.gamma", "pos_t.delta", "pos_t.epsilon", "pos_t.zeta"))

data.mat<-acast(data.all, NAME ~ sample, value.var="FDR.q.val")
data.mat[is.na(data.mat)]=0

neg.mat<-data.mat[,1:9]
otter.dendro.neg <- as.dendrogram(hclust(d = dist(x = neg.mat), method="complete"))
otter.order.neg <- order.dendrogram(otter.dendro.neg)
order.neg<-rownames(neg.mat[otter.order.neg,])

pos.mat<-data.mat[,10:18]
otter.dendro.pos <- as.dendrogram(hclust(d = dist(x = pos.mat), method="complete"))
otter.order.pos <- order.dendrogram(otter.dendro.pos)
order.pos<-rownames(pos.mat[otter.order.pos,])

data.pos<-data.all[grep("pos", data.all$sample),]
data.neg<-data.all[grep("neg", data.all$sample),]

data.pos$NAME <- factor(x = data.pos$NAME,
                            levels = order.pos)
data.neg$NAME <- factor(x = data.neg$NAME,
                            levels = order.neg)


pdf("plots/gsea_pot.ordered.pdf", height=9, width=15)
pos<-ggplot(data.pos, aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="magma",  begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
neg<-ggplot(data.neg, aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="viridis", begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
ggarrange(neg, pos, nrow=1, ncol=2, common.legend = F)
dev.off()


pdf("plots/gsea_awakening_tep.pdf", height=9, width=18)
pos<-ggplot(data.all[grep("pos", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="magma",  begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
neg<-ggplot(data.all[grep("neg", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="viridis", begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
ggarrange(neg, pos, nrow=1, ncol=2, common.legend = F)
dev.off()



###############################

files<-list.files("gsea_dormant")
str_files<-as.data.frame(str_split(files, "_|\\.", simplify=T), stringsAsFactors=F)

samples<-paste(paste(str_files$V5, str_files$V6, sep="_"), str_files$V7, sep=".")

data<-lapply(files, function(x){
	info<-read.table(paste("gsea_dormant", x, sep="/"), header=T, sep="\t", stringsAsFactors=F)
	strx<-as.data.frame(str_split(x, "_|\\.", simplify=T), stringsAsFactors=F)
	info.sig<-info[info$FDR.q.val<=0.25, c(1, 4, 8)]
		# info.sig<-info[info$FDR.q.val<=0.10, c(1, 4, 8)]
		info.sig$sample<-paste(paste(strx$V5, strx$V6, sep="_"), strx$V7, sep=".")
		info.sig$NAME<-str_replace_all(info.sig$NAME, "HALLMARK_", "")
	info.sig
})

data.all<-do.call(rbind, data)
data.all$sample<-factor(data.all$sample, levels=c("neg_ai.alpha", "neg_ai.beta", "neg_ai.gamma", "neg_t.alpha", "neg_t.beta", "neg_t.gamma", "neg_t.delta", "neg_t.epsilon", "neg_t.zeta", "pos_ai.alpha", "pos_ai.beta", "pos_ai.gamma", "pos_t.alpha", "pos_t.beta", "pos_t.gamma", "pos_t.delta", "pos_t.epsilon", "pos_t.zeta"))

data.mat<-acast(data.all, NAME ~ sample, value.var="FDR.q.val")
data.mat[is.na(data.mat)]=0

neg.mat<-data.mat[,1:8]
otter.dendro.neg <- as.dendrogram(hclust(d = dist(x = neg.mat), method="complete"))
otter.order.neg <- order.dendrogram(otter.dendro.neg)
order.neg<-rownames(neg.mat[otter.order.neg,])

pos.mat<-data.mat[,9:16]
otter.dendro.pos <- as.dendrogram(hclust(d = dist(x = pos.mat), method="complete"))
otter.order.pos <- order.dendrogram(otter.dendro.pos)
order.pos<-rownames(pos.mat[otter.order.pos,])

data.pos<-data.all[grep("pos", data.all$sample),]
data.neg<-data.all[grep("neg", data.all$sample),]

data.pos$NAME <- factor(x = data.pos$NAME,
                            levels = order.pos)
data.neg$NAME <- factor(x = data.neg$NAME,
                            levels = order.neg)

pdf("plots/gsea_awakening_dormant.pdf", height=9, width=18)
pos<-ggplot(data.all[grep("pos", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="magma",  begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
neg<-ggplot(data.all[grep("neg", data.all$sample),], aes(x=sample, y=NAME, color=FDR.q.val,size=SIZE, aes=FDR.q.val)) +
	geom_point() +
	theme_bw() +
	scale_color_viridis(option="viridis", begin=0.2, end=1) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
		axis.text = element_text(size = 10),
		legend.position="bottom")
ggarrange(neg, pos, nrow=1, ncol=2, common.legend = F)
dev.off()


