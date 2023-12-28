require("reshape2")
require("ggplot2")
require("stringr")
require("RColorBrewer")
require("stringr")
require("viridis")
require("ggrepel")


########## Include T0s ########################
counts.recal <- readRDS("rds/data.recal.rds")
counts.batch4<-readRDS("rds/counts_batch4.rds")

counts.all.T0<-merge(counts.recal, counts.batch4, by="row.names", all=TRUE)
counts.all.T0[is.na(counts.all.T0)] <- 0

saveRDS(counts.all.T0, "rds/counts.all.T0.rds")
write.table(counts.all.T0,file ="rds/counts.all.T0.txt",sep="\t",quote=F,row.names=T,col.names=T)

########## No. of barcodes and reads ########################
counts.all.T0<-readRDS("rds/counts.all.T0.rds")

# counts.all.T0<-counts.all
# counts.all<-counts.all.T0[,2:132]
# rownames(counts.all)<-counts.all.T0[,1]

counts.all<-counts.all.T0

t0<-counts.all[,c(1:3, 112:131)]
info<-data.frame("reads"=colSums(t0)/1e6, "barcodes"=colSums(t0>0)/1e3, "samples"=colnames(t0))

data<-counts.all[,c(1:62, 112:131)]
data<-data[,-c(22,23)]

total<-colSums(data) 
freq<-as.data.frame(t(t(data)/total), stringsAsFactors=F)
freq.t0<-freq[, c(1:3, 61:80)]

info<-data.frame("reads"=colSums(data)/1e6, "barcodes"=colSums(data>0)/1e3, "samples"=colnames(data))
names.samples<-read.table("rds/rename.all.txt", header=TRUE, stringsAsFactors=F)
info$new<-""
info[1:60, ]$new<-names.samples[match(info[1:60, ]$samples, names.samples$old),]$new
info[61:80, ]$new<-as.vector(info[61:80, ]$samples)

info$state<-""
info[1:3,]$state<-"POT"
info[4:15,]$state<-"Latency"
info[c(16:29, 32:35),]$state<-"Dormancy"
info[c(30:31, 36:39, 42:48),]$state<-"Awakening"
info[49:60,]$state<-"TEP"
info[grep("U", info$new), ]$state<-"Untreated"
info[61:80,]$state<-"T0"
info$state<-factor(info$state, levels=c("POT", "T0", "Untreated", "Latency", "Dormancy", "Awakening", "TEP"))
scheme<-c("POT"="#0084A5", "T0"="blue", "Untreated"="#0084A5", "Latency"="#0084A5", "Dormancy"="#FFDC00", "Awakening"="#B30078", "TEP"="#50157A")

############## Correlation in frequency ######################
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
png("corr.t0.png", height=2000, width=2000)
pairs(freq.t0, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      font.labels = 10)
dev.off()

