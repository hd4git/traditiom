require("VariantAnnotation")
require("stringr")
require("ggplot2")
require("reshape2")
require("circlize")
require("viridis")
require("dndscv")

setwd("~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header")
filteredFiles <- list.files(pattern = "\\_filtered.vcf.gz$")


mutations <- lapply(filteredFiles, function(ip){
	test <- readVcf(ip, "hg38")
	names(test)
})
sampleList <- lapply(filteredFiles, function(x){
	substr(x, 1, 5)
})

names(mutations) <- as.vector(melt(sampleList)[,1])
mutations.m <- melt(mutations)
mutSplit <- str_split(mutations.m$value, ":|_|/", simplify=TRUE)

mutations.df <- data.frame("sampleID" = mutations.m$L1, 
							"chr" = mutSplit[,1],
							"pos" = mutSplit[,2],
							"ref" = mutSplit[,3],
							"mut" = mutSplit[,4])
mutations.df$chr <- str_replace(mutations.df$chr, "chr", "")
setwd("~/lateRelapse_targetedPanel/rds")
load("RefCDS.rda")
targetList <- read.table("targetList.txt", header=FALSE)[,1]
targetList <- targetList[!targetList %in% c("APOBEC3B-AS1", "CASC11", "CBR3-AS1", "CDKN2A", "CDKN2B-AS1", "ensembl", "EP300-AS1", "FABP5P3", "GNAS-AS1", "havana", "havana_tagene", "IRS4-AS1", "MIR4673", "MIR4713HG", "MIR4728", "MIR4733HG", "MIR6755", "MTOR-AS1", "NDUFA6-DT", "NTRK3-AS1", "RPL21P4", "SPEN-AS1", "SYNM-AS1", "TBX3-AS1", "TTC36-AS1", "XPC-AS1")]

getInfo <- function(x){
	dndsout <- dndscv(mutations.df, refdb=RefCDS, outmats=T, gene_list=x, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
	ci <- geneci(dndsout)
	sigCI <- ci[(ci[,4]>1 & ci[,6]>1) | (ci[,5]>1 & ci[,7]>1),]

	info <- dndsout$sel_cv
	sel <- data.frame(info$gene_name, info$n_mis/info$n_syn, info$qglobal_cv)
	sigQGlobal <- info[info$qglobal_cv<=0.1, ]
	list("main"=dndsout, "ci"=ci, "sigCI"=sigCI, "sigQGlobal"=sigQGlobal)
}
RefCDSgenes <- melt(lapply(1:length(RefCDS), function(x){RefCDS[x][[1]]$gene_name}))$value
intogenDrivers <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes.tsv", header=TRUE, sep="\t")[,1]
intogenDriversBRCA <- read.table("~/T47D_targetedPanel/IntOGen-DriverGenes_BRCA.tsv", header=TRUE, sep="\t")[,1]
resDrivers <- c("ESR1", "TP53", "AKT1", "RIC8A", "RB1", "NF1", "FRG1", "KMT2C", "NCOR1")

infoAll <- getInfo(targetList)
infoIntogen <- getInfo(intogenDrivers[intogenDrivers %in% RefCDSgenes])
infoIntogenBRCA <- getInfo(intogenDriversBRCA[intogenDriversBRCA %in% RefCDSgenes])
infoResistance <- getInfo(resDrivers[resDrivers %in% RefCDSgenes])

infoData <- list("infoAll"=infoAll,
				"infoIntogen"=infoIntogen, 
				"infoIntogenBRCA"=infoIntogenBRCA,
				"infoResistance"=infoResistance)
saveRDS(infoData, "infoData.rds")

ggplot(data, aes(x, y)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))








# ci.m$type <- str_split(ci.m$variable, "_", simplify=TRUE)[,2]
ci.m$type <- str_split(ci.m$variable, "_", simplify=TRUE)[,1]
ggplot(ci.m, aes(x=variable, y=value)) +
	geom_violin() +
	geom_point(size=0.1, position = position_jitter(seed = 1, width = 0.2)) +
	theme_bw() +
	facet_wrap(type~., scales="free")
dev.off()	


##   sampleID chr      pos ref mut
