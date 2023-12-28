require("reshape2")
require("ggplot2")
require("stringr")
require("fishplot")
require("RColorBrewer")
require("DESeq2")
require("ggpubr")

## Import samples
setwd("/rds/general/user/hdhiman/home/barcode_analysis/results/filtered")
samples<-list.files("counts_filtered_batch1")
ids<-do.call(rbind, str_split(samples, pattern="\\.", n=3))[,1]
info<-data.frame(sample_id=paste("S", 1:23, sep=""), desc=c("T0.1_0D", "T0.2_0D", "T0.3_0D", "UT.1_7D", "UT.2_7D", "TAM.1_7D", "TAM.2_7D", "AI.1_7D", "AI.2_7D", "UT.1_14D", "UT.2_14D", "TAM.1_14D", "TAM.2_14D", "AI.1_14D", "AI.2_14D", "UT.1_30D", "UT.2_30D", "TAM.1_30D", "TAM.2_30D", "AI.1_30D", "AI.2_30D", "TAM.1.Dead_14D", "TAM.1.Dead_35D"))

## Get raw counts and number of barcodes
counts.list<-lapply(samples, function(x){read.table(paste("counts_filtered_batch1",x, sep="/"), header=F)})
names(counts.list)<-ids

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts)
counts<-counts[,info$sample_id]
colnames(counts)<-info$desc
counts<-counts[,-13]
saveRDS(counts, "rds/counts_batch1.rds")
write.table(counts,file ="rds/counts_filtered_batch1.txt",sep="\t",quote=F,row.names=T,col.names=T)


########################################################################################################
## Import samples
samples<-list.files("counts_filtered_batch2")
info<-read.table("rds/sample_info.txt", header=T)

## Get raw counts and number of barcodes
counts.list<-lapply(samples, function(x){read.table(paste("counts_filtered_batch2",x, sep="/"), header=F)})
ids<-do.call(rbind, str_split(samples, pattern="\\.", n=3))[,1]
names(counts.list)<-ids

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts, stringsAsFactors=F)
nm<-paste("S", 1:94, sep="")
counts<-counts[,na.omit(colnames(counts)[match(nm, colnames(counts))])]
colnames(counts)<-info[info$id %in%  colnames(counts), ]$sample3

saveRDS(counts, "rds/counts_batch2.rds")
write.table(counts,file ="rds/counts_filtered_batch2.txt",sep="\t",quote=F,row.names=T,col.names=T)

########################################################################################################
## Import samples
samples<-list.files("counts_filtered_batch3")
info<-read.table("rds/sample_info.txt", header=T)

## Get raw counts and number of barcodes
counts.list<-lapply(samples, function(x){read.table(paste("counts_filtered_batch3",x, sep="/"), header=F)})
ids<-do.call(rbind, str_split(samples, pattern="\\.", n=3))[,1]
names(counts.list)<-ids

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts, stringsAsFactors=F)
nm<-paste("S", 1:94, sep="")
counts<-counts[,na.omit(colnames(counts)[match(nm, colnames(counts))])]
colnames(counts)<-info[info$id %in%  colnames(counts), ]$sample3

saveRDS(counts, "rds/counts_batch3.rds")
write.table(counts,file ="rds/counts_filtered_batch3.txt",sep="\t",quote=F,row.names=T,col.names=T)

########################################################################################################
## Import samples for T0s
samples<-list.files("counts_filtered_batch4")
info<-read.table("rds/sample_info.txt", header=T)

## Get raw counts and number of barcodes
counts.list<-lapply(samples, function(x){read.table(paste("counts_filtered_batch4",x, sep="/"), header=F)})
ids<-do.call(rbind, str_split(samples, pattern="\\.", n=3))[,1]
names(counts.list)<-ids

for(i in names(counts.list)) {counts.list[[i]]$id<-i}
counts.list2<-do.call(rbind, counts.list)
counts<-t(acast(counts.list2, id ~ V1, value.var="V2"))
counts[is.na(counts)] <- 0

counts<-as.data.frame(counts, stringsAsFactors=F)
nm<-paste("s", 1:20, sep="")
counts<-counts[,na.omit(colnames(counts)[match(nm, colnames(counts))])]
info<-data.frame("id"=nm, "sample"=c(paste("AI", 1:10, sep=""), paste("TAM", 1:10, sep="")))

colnames(counts)<-info[info$id %in%  colnames(counts), ]$sample

saveRDS(counts, "rds/counts_batch4.rds")
write.table(counts,file ="rds/counts_filtered_batch4.txt",sep="\t",quote=F,row.names=T,col.names=T)


########################################################################################################
counts.batch1<-readRDS("rds/counts_batch1.rds")
counts.batch2<-readRDS("rds/counts_batch2.rds")
counts.batch3<-readRDS("rds/counts_batch3.rds")

counts.all<-merge(counts.batch1, counts.batch2, by="row.names", all=TRUE)
counts.all[is.na(counts.all)] <- 0
colnames(counts.all)[30]<-"TAM.2_14D"

data<-counts.all[, c(2:13, 30, 14:29, 31:117)]
# data<-counts.all[, 2:110]
rownames(data)<-counts.all[,1]
old.sample.names<-c("T0.1_0D", "T0.2_0D", "T0.3_0D", "UT.1_7D", "UT.2_7D", "TAM.1_7D", "TAM.2_7D", "AI.1_7D", "AI.2_7D", "UT.1_14D", "UT.2_14D", "TAM.1_14D", "TAM.2_14D", "AI.1_14D", "AI.2_14D", "UT.1_30D", "UT.2_30D", "TAM.1_30D", "TAM.2_30D", "AI.1_30D", "AI.2_30D", "TAM.1.Dead_14D", "TAM.1.Dead_35D")
colnames(data)[1:23]<-old.sample.names
saveRDS(data, "rds/counts.all.rds")
write.table(data,file ="rds/counts.all.txt",sep="\t",quote=F,row.names=T,col.names=T)


########################################################################################################
## Get counts and number of barcodes with at least 2 reads 
counts.gt2<-data
counts.gt2[counts.gt2 == 1] <- 0
counts.gt2<-counts.gt2[rowSums(counts.gt2>0)>0,]
saveRDS(counts.gt2, "rds/counts.gt2.all.rds")
write.table(counts.gt2,file ="rds/counts.gt2.all.txt",sep="\t",quote=F,row.names=T,col.names=T)

########################################################################################################
## Get counts and number of barcodes with at least 3 reads 
counts.gt3<-data
counts.gt3[counts.gt3 == 1] <- 0
counts.gt3[counts.gt3 == 2] <- 0
counts.gt3<-counts.gt3[rowSums(counts.gt3>0)>0,]
saveRDS(counts.gt3, "rds/counts.gt3.all.rds")
write.table(counts.gt3,file ="rds/counts.gt3.all.txt",sep="\t",quote=F,row.names=T,col.names=T)


