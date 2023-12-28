library(pbmcapply)
library(THmisc)
library(dplyr)
library(VariantAnnotation)
library(ggplot2)
library(cowplot)
theme_set(cowplot::theme_cowplot())
source("~/tools/T-Heide/functions/mp_trees.R")

idx = commandArgs(trailingOnly = TRUE)

## Options
data_dir = "/rds/general/user/hdhiman/ephemeral/T47D_targetedPanel/data/results/final/vep"
out_dat_dir = "/rds/general/user/hdhiman/ephemeral/T47D_targetedPanel/data/results/final/rds"
reference_sample = c("T47D_POT1", "T47D_POT2", "T47D_POT3")
vcf_suffix = "[.]vcf([.]gz)?$"
dir.create(out_dat_dir, FALSE, TRUE)

# detect and sort vcfs:
vcf_files = data_dir %>% 
  list.files(vcf_suffix, full.names = TRUE) %>% 
  magrittr::set_names(gsub(vcf_suffix, "", basename(.)))

# load vcfs file by file
if (length(idx)) vcf_files = vcf_files[as.numeric(idx)]
print(vcf_files)

for (i in seq_along(vcf_files)) {
  data = THmisc::load_vcf_file(f=vcf_files[i], annot = FALSE)
  out_file = file.path(out_dat_dir, paste0(names(vcf_files)[i], ".rds"))
  saveRDS(data, out_file)
}

