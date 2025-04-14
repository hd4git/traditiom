setwd("/rds/general/user/hdhiman/home/wgs/timon")

require("ggplot2")
require("stringr")
require("viridis")
require("RColorBrewer")

resDrivers <- c("ESR1", "TP53", "AKT1", "RIC8A", "RB1", "NF1", "FRG1", "KMT2C", "NCOR1")

sampleNamesMCF7 <- read.table("~/manuscriptFeb23/rds/wgs_newNomenclature.txt", header=T)
T47D_intogen <- read.table("T47D_plot_driver_data_intogen.csv", header=T, sep=",")
T47D_BRCAdriver <- read.table("T47D_plot_driver_data_brca_driver.csv", header=T, sep=",")
T47D_driver <- read.table("T47D_plot_driver_data.csv", header=T, sep=",")

t47d <- list("T47D_intogen" <- T47D_intogen,
			"T47D_BRCAdriver" <- T47D_BRCAdriver,
			"T47D_driver" <- T47D_driver)

getInfo <- function(ip){
	ip$sample <- substr(ip$sample, 6,nchar(ip$sample))
	ip$phase <- "AI"
	ip[grepl("POT", ip$sample),]$phase <-"POT"
	ip[grepl("^UT", ip$sample),]$phase <-"UT"
	ip$phase <- factor(ip$phase, levels=c("POT", "UT", "AI"))
	ip
}


t47d <- lapply(t47d, getInfo)

plotT47Dmako <- function(data, len) {
	data <- data[data$sample!="UT3_TEP",]
	ggplot(data, aes(x=sample, y=variant, fill=VAF)) +
	geom_tile() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		axis.title=element_blank(),
		strip.text.y = element_blank()
		# strip.text.y.right = element_text(angle = 0)
		) +
	scale_fill_viridis(option = "G", begin = 0.25, end = 1, direction=-1, na.value ="white") +
	scale_y_discrete(labels = function(x) str_wrap(x, width = len)) +
	facet_grid(gene~phase, space = "free", scales="free")
}
plotT47DBuPu <- function(data, len) {
	data <- data[data$sample!="UT3_TEP",]
	ggplot(data, aes(x=sample, y=variant, fill=VAF)) +
	geom_tile() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		axis.title=element_blank(), 
		strip.text.y = element_blank()
		# strip.text.y.right = element_text(angle = 0)
		) +
	scale_fill_distiller(type = "seq", palette = "BuPu", direction=1, na.value ="white") +
	scale_y_discrete(labels = function(x) str_wrap(x, width = len)) +
	facet_grid(gene~phase, space = "free", scales="free")
}

pdf("~/manuscriptFeb23/plots/wgs/T47D_intogen.pdf", height = 7, width=10)
plotT47Dmako(t47d[[1]], 40)
plotT47DBuPu(t47d[[1]], 40)
dev.off()
pdf("~/manuscriptFeb23/plots/wgs/T47D_BRCAdriver.pdf", height=4, width=7)
plotT47Dmako(t47d[[2]], 30)
plotT47DBuPu(t47d[[2]], 30)
dev.off()
pdf("~/manuscriptFeb23/plots/wgs/T47D_driver.pdf", height=4, width=5)
plotT47Dmako(t47d[[3]][t47d[[3]]$gene %in% resDrivers,], 30)
plotT47DBuPu(t47d[[3]][t47d[[3]]$gene %in% resDrivers,], 30)
dev.off()

mcf7 <- readRDS("MCF7_heatmap.rds")


MCF7_intogen <- read.table("MCF7_plot_driver_data_all_intogen_driver.csv", header=T, sep=",")
MCF7_BRCAdriver <- read.table("MCF7_plot_driver_data_brca_intogen.csv", header=T, sep=",")
MCF7_driver <- read.table("MCF7_plot_driver_data_brca_driver.csv", header=T, sep=",")

mcf7 <- list("MCF7_intogen" <- MCF7_intogen,
			"MCF7_BRCAdriver" <- MCF7_BRCAdriver,
			"MCF7_driver" <- MCF7_driver)

getInfoMCF7 <- function(ip){
	ip <- ip[!grepl("new", ip$sample),]
	ip$newName <- ""
	ip$newName <- sampleNamesMCF7[match(ip$sample, sampleNamesMCF7$Old),]$New 
	ip$phase <- "Awakening"
	ip[grepl("POT", ip$newName),]$phase <-"POT"
	ip[grepl("TEP", ip$newName),]$phase <-"TEP"
	ip[grepl("^U", ip$newName),]$phase <-"UT"
	ip$condition <- "POT"
	ip[grepl("^U", ip$newName),]$condition <-"UT"
	ip[grepl("^T", ip$newName),]$condition <-"T"
	ip[grepl("^AI", ip$newName),]$condition <-"AI"
	ip$phase <- factor(ip$phase, levels=c("POT", "UT", "Awakening", "TEP"))
	ip$condition <- factor(ip$condition, levels=c("POT", "UT", "T", "AI"))
	ip
}


mcf7 <- lapply(mcf7, getInfoMCF7)

plotMCF7mako <- function(data, len) {
	data <- data[data$newName != "U170.3",]
	data$newName <- factor(data$newName, levels=c("POT.1", "POT.2", "POT.3", "U170.1", "U170.2", "U170.3", "Ta", "TaTEP", "Tb", "TbTEP", "Tg", "TgTEP", "Td", "TdTEP", "Te", "TeTEP", "Th", "Tq", "Tz", "TzTEP", "AI90.1", "AI90.2", "AIa", "AIaTEP", "AIb", "AIbTEP", "AIg", "AIgTEP", "AId", "AIe"))
	ggplot(data, aes(x=newName, y=variant, fill=VAF)) +
	geom_tile() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		axis.title=element_blank(),
		strip.text.y = element_blank()
		# strip.text.y.right = element_text(angle = 0)
		) +
	scale_fill_viridis(option = "G", begin = 0.25, end = 1, direction=-1, na.value ="white") +
	scale_y_discrete(labels = function(x) str_wrap(x, width = len)) +
	facet_grid(gene~condition, space = "free", scales="free")
}
plotMCF7BuPu <- function(data, len) {
	data <- data[data$newName != "U170.3",]
	data$newName <- factor(data$newName, levels=c("POT.1", "POT.2", "POT.3", "U170.1", "U170.2", "U170.3", "Ta", "TaTEP", "Tb", "TbTEP", "Tg", "TgTEP", "Td", "TdTEP", "Te", "TeTEP", "Th", "Tq", "Tz", "TzTEP", "AI90.1", "AI90.2", "AIa", "AIaTEP", "AIb", "AIbTEP", "AIg", "AIgTEP", "AId", "AIe"))
	ggplot(data, aes(x=newName, y=variant, fill=VAF)) +
	geom_tile() +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle = 90),
		strip.background = element_rect(fill="#FFFFFF"), 
		legend.position="bottom",
		axis.title=element_blank(), 
		strip.text.y = element_blank()
		# strip.text.y.right = element_text(angle = 0)
		) +
	scale_fill_distiller(type = "seq", palette = "BuPu", direction=1, na.value ="white") +
	scale_y_discrete(labels = function(x) str_wrap(x, width = len)) +
	facet_grid(gene~condition, space = "free", scales="free")
}

pdf("~/manuscriptFeb23/plots/wgs/MCF7_intogen.pdf", height = 18, width=10)
plotMCF7mako(mcf7[[1]], 40)
plotMCF7BuPu(mcf7[[1]], 40)
dev.off()
pdf("~/manuscriptFeb23/plots/wgs/MCF7_BRCAdriver.pdf", height=9, width=7)
plotMCF7mako(mcf7[[2]], 30)
plotMCF7BuPu(mcf7[[2]], 30)
dev.off()
pdf("~/manuscriptFeb23/plots/wgs/MCF7_driver.pdf", height=4, width=7)
plotMCF7mako(mcf7[[3]][mcf7[[3]]$gene %in% resDrivers,], 30)
plotMCF7BuPu(mcf7[[3]][mcf7[[3]]$gene %in% resDrivers,], 30)
dev.off()




