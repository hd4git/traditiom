require("VariantAnnotation")
require("stringr")
require("reshape2")
setwd("~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr")
files <- list.files(pattern = "\\PASSgermline_targeted.vcf.gz$")

### Get allele depth (AD) and F1R2 and F2R1 counts ###
### Filter for variants with total AD >= 20; 
### F1R2 and F2R1 of alternate allele >= 4;
### allele frequency >= 0.1

dataFilteredList <- lapply(files, function(x){
	sample <- substr(x, 1, 5)
	data <- read.table(x, header=F, sep="\t")
	AD <- str_split(data$V10, ":", simplify=T)[,2]
	ADref <- as.numeric(str_split(AD, ",", simplify=T)[,1])
	ADalt <- as.numeric(str_split(AD, ",", simplify=T)[,2])
	ADsum <- ADref + ADalt
	AF <- str_split(data$V10, ":", simplify=T)[,3]
	F1R2 <- str_split(data$V10, ":", simplify=T)[,5]
	F2R1 <- str_split(data$V10, ":", simplify=T)[,6]
	ALT_F1R2 <- as.numeric(str_split(F1R2, ",", simplify=T)[,2])
	ALT_F2R1 <- as.numeric(str_split(F2R1, ",", simplify=T)[,2])
	ALT_F1R2_F2R1 <- ALT_F1R2+ALT_F2R1
	data$ADsum <- ADsum
	data$AF <- AF
	data$ALT_F1R2_F2R1 <- ALT_F1R2_F2R1
	dataFiltered <-  data[data$ADsum >=20 & data$ALT_F1R2_F2R1 >=4 & data$AF >=0.1,]
	dataFiltered
})

sampleList <- lapply(files, function(x){
	substr(x, 1, 5)
})

names(dataFilteredList) <- melt(sampleList)[,1]
setwd("~/lateRelapse_targetedPanel/")
lapply(names(dataFilteredList), function(x){
	write.table(dataFilteredList[[x]], file = paste(paste("rds/PON_PASSgermline_chr_filtered", x, sep="/"), "vcf", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
	})
