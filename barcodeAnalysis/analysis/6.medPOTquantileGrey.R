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

    xd <- data.frame(density(log10(pot.m$value))[c("x", "y")])
    xd2 <- xd
    xd2$y <- -xd2$y
    xd3 <- rbind(xd, xd2)
    p <- ggplot(xd3, aes(x, y)) + 
      geom_line(data=xd, aes(x, y)) +
      geom_area(data = xd[xd$x <= log10(quantile(pot)[2]),], fill = "#AEB4A9") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[2]) & xd$x <= log10(quantile(pot)[3]),], fill = "#E0C1B3") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[3]) & xd$x <= log10(quantile(pot)[4]),], fill = "#D89A9E") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[4]),], fill = "#C37D92") +
      geom_line(data=xd2, aes(x, y)) +
      geom_area(data = xd2[xd2$x <= log10(quantile(pot)[2]),], fill = "#AEB4A9") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[2]) & xd2$x <= log10(quantile(pot)[3]),], fill = "#E0C1B3") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[3]) & xd2$x <= log10(quantile(pot)[4]),], fill = "#D89A9E") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[4]),], fill = "#C37D92") +
      theme_bw() +
      labs(x="log10(Median frequency)") +
      coord_flip() +
      theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              text=element_text(size=10))     
pdf("plots_MS/medPotQuantile.pdf", height=2, width=0.9)
    p 
dev.off()

    p <- ggplot(xd3, aes(x, y)) + 
      geom_line(data=xd, aes(x, y)) +
      geom_area(data = xd[xd$x <= log10(quantile(pot)[2]),], fill = "#f7f7f7") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[2]) & xd$x <= log10(quantile(pot)[3]),], fill = "#cccccc") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[3]) & xd$x <= log10(quantile(pot)[4]),], fill = "#969696") +
      geom_area(data = xd[xd$x > log10(quantile(pot)[4]),], fill = "#525252") +
      geom_line(data=xd2, aes(x, y)) +
      geom_area(data = xd2[xd2$x <= log10(quantile(pot)[2]),], fill = "#f7f7f7") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[2]) & xd2$x <= log10(quantile(pot)[3]),], fill = "#cccccc") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[3]) & xd2$x <= log10(quantile(pot)[4]),], fill = "#969696") +
      geom_area(data = xd2[xd2$x > log10(quantile(pot)[4]),], fill = "#525252") +
      theme_bw() +
      labs(x="log10(Median frequency)") +
      coord_flip() +
      theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              text=element_text(size=10))     
pdf("plots_MS/medPotQuantileGrey.pdf", height=2, width=0.9)
    p 
dev.off()

