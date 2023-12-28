library(pbapply)
library(THmisc)
library(dplyr)
library(VariantAnnotation)
library(ggplot2)
library(cowplot)
theme_set(cowplot::theme_cowplot())
source("~/tools/T-Heide/functions/mp_trees.R")
pboptions(type="txt") 

idx = commandArgs(trailingOnly = TRUE)

## Options
data_dir = "/rds/general/user/hdhiman/ephemeral/reWGS_MCF7/results/all/vcf_data"
out_dat_dir = "/rds/general/user/hdhiman/ephemeral/reWGS_MCF7/results/all/filtered_vcf_data2"
out_dat_dir_all = "/rds/general/user/hdhiman/ephemeral/reWGS_MCF7/results/all/filtered_vcf_data_raw2"

# create dirs
dir.create(out_dat_dir, FALSE, TRUE)
dir.create(out_dat_dir_all, FALSE, TRUE)

# functions
do_test = function(x, ...) {
  nv = matrix(unlist(c(geno(x)$NV)), NROW(x), NCOL(x), dimnames = dimnames(x))
  nr = matrix(unlist(c(geno(x)$NR)) - unlist(c(geno(x)$NV)), NROW(x), NCOL(x), dimnames = dimnames(x))
  
  nv = split(nv,  seq_len(NROW(nv)))
  nr = split(nr,  seq_len(NROW(nr)))
  
  # drop NA entries
  keep = mapply(function(x, y) !(is.na(x) | is.na(y) | (x+y) == 0), nv, nr, SIMPLIFY = FALSE)
  nv = mapply("[", nv, keep)
  nr = mapply("[", nr, keep)
  tabs = mapply(cbind, nr, nv)
  
  test = function(d, large_b=1e5) {
    tryCatch({
      chisq.test(d)$p.value
    }, warning=function(e) {
      p = chisq.test(d, simulate.p.value = TRUE, B=1e3)$p.value
      if (p < 0.01) p = chisq.test(d, simulate.p.value = TRUE, B=large_b)$p.value
      return(p)
    }, error = function(e) {
      print(e)
      return(NA)
    })
  }
  
  p = rep(1, length(tabs))
  do_test = sapply(tabs, function(x) NROW(x) > 0 & !all(x[,1] == 0) & !all(x[,2] == 0))
  p[do_test] = unlist(pbapply::pblapply(tabs[do_test], test, ...))
  
  return(p)
}


# detect and sort vcfs:
vcf_files = data_dir %>% 
  list.files("[.]rds$", full.names = TRUE) %>% 
  magrittr::set_names(basename(.))


# load vcfs file by file
if (length(idx)) vcf_files = vcf_files[as.numeric(idx)]
print(vcf_files)

for (i in seq_along(vcf_files)) {
  data = readRDS(vcf_files[i])
  out_file = file.path(out_dat_dir, paste0(names(vcf_files)[i]))
  out_file_raw = file.path(out_dat_dir_all, paste0(names(vcf_files)[i]))
  if (file.exists(out_file)) next()
  info(data)$p = do_test(data)
  saveRDS(data, out_file_raw)
  saveRDS(data[info(data)$p <= 0.01,], out_file)
}


data_dir = "/rds/general/user/hdhiman/home/wgs/timon/MCF7_platypus/filtered_vcf_data2"
vcf_files = data_dir %>% 
  list.files("[.]rds$", full.names = TRUE) %>% 
  magrittr::set_names(basename(.))

for (i in seq_along(vcf_files)) {
  nm <- str_split(vcf_files[i], "\\.", simplify=TRUE)[,2]
  data <- readRDS(vcf_files[i])
  mat <- as.data.frame(info(data))
  # loc <- as.data.frame(str_split(rownames(mat), ":|_|/", simplify=TRUE))
  # colnames(loc) <- c("Chr","Pos", "Ref", "Alt")
  # mat <- cbind(loc, mat)
  write.table(rownames(mat), file = paste("filtered_vcf_txt", paste(paste("MCF7_platypus", nm, sep="_"), "txt", sep="."), sep="/"), sep="\t", quote=F, row.names=F, col.names=T)
}
setwd("~/wgs/timon/MCF7_platypus/filtered_vcf_files/")
for (i in seq_along(vcf_files)) {
  nm <- str_split(vcf_files[i], "\\.", simplify=TRUE)[,2]
  data <- readRDS(vcf_files[i])
  writeVcf(data, file = paste(paste("MCF7_platypus", nm, sep="_"),"vcf", sep="."))
}

cd ~/wgs/timon/MCF7_platypus/filtered_vcf_files/
cat MCF7_platypus_chr1_vep.vcf > ../MCF7_platypus.vcf
while read f;
do 
grep -v "#" $f >> ../MCF7_platypus.vcf
done < <(ls | sort --version-sort | head -23 | tail -22)

cut -f 1-9,13-37,39- MCF7_platypus.vcf > MCF7_platypus_woNew.vcf





grep -v "#" T47D_platypus_joined_mutation_calls_PASS_sorted.vcf | awk '{print $1":"$2"\t"$0}' > T47D_platypus_joined_mutation_calls_PASS_sortedAll.vcf
awk 'NR==FNR{a[$1]=$0; next} ($1 in a) {print $0}' <(awk '{split($0, a, "_"); print a[1]}' T47D_patypus_filtered.txt) T47D_platypus_joined_mutation_calls_PASS_sortedAll.vcf | cut -f 2-> T47D_platypus_filtered_oriAll.txt
