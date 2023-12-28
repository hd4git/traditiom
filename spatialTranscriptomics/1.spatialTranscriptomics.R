require("NanoStringNCTools")
require("GeomxTools")
require("GeoMxWorkflows")
require("ggplot2")
require("knitr")
require("dplyr")
require("ggforce")
require("scales")
require("stringr")
require("reshape2")

setwd("~/spatialTranscriptomics/")

datadir <- file.path("~/spatialTranscriptomics/data")

DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- file.path("~/spatialTranscriptomics/data/pkcs/Hs_R_NGS_WTA_v1.0.pkc")
SampleAnnotationFile <- file.path("~/spatialTranscriptomics/data/annotation/LabWorksheet.xlsx")
# load data
demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                           pkcFiles = PKCFiles,
                           phenoDataFile = SampleAnnotationFile,
                           phenoDataSheet = "Template",
                           phenoDataDccColName = "Sample_ID",
                           # protocolDataColNames = c("aoi", "roi"),
                           experimentDataColNames = c("panel"))
setwd("~/spatialTranscriptomics/")
saveRDS(demoData, "rds/data.rds")

demoData<-readRDS("rds/data.rds")
data <- demoData
dim(data)
# Features  Samples
#    18815      173
demoData <-data

## Modules Used
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

rtsInfo <- melt(exprs(data))
rtsInfoList <- split(rtsInfo, rtsInfo$Var2)
rtsCounts <- melt(lapply(rtsInfoList, function(x){sum(x[,3])}))

# pData(data)[["NTC"]] = rtsCounts[match(rownames(pData(data)), rtsCounts$L1),]$value
# pData(demoData)[["NTC"]] = rtsCounts[match(rownames(pData(demoData)), rtsCounts$L1),]$value


ntc.fData <- fData(data)

rtsNtc084Info <- rtsInfo[rtsInfo$Var2 == "DSP-1001660009084-H-A01.dcc",]
rtsNtc086Info <- rtsInfo[rtsInfo$Var2 == "DSP-1001660009086-G-A01.dcc",]

ntc084.fData <- ntc.fData
ntc086.fData <- ntc.fData


ntc084.fData$value <- rtsNtc084Info[match(ntc.fData$RTS_ID, rtsNtc084Info[,1]),]$value
ntc086.fData$value <- rtsNtc086Info[match(ntc.fData$RTS_ID, rtsNtc086Info[,1]),]$value

ntc084.fData <- ntc084.fData[order(ntc084.fData$value, decreasing=T),]
ntc086.fData <- ntc086.fData[order(ntc086.fData$value, decreasing=T),]

ntc084.fData <- ntc084.fData[ntc084.fData$value>0,]
ntc086.fData <- ntc086.fData[ntc086.fData$value>0,]

highRTS0.84 <- rownames(ntc084.fData[ntc084.fData$value>=10,])
highRTS0.86 <- rownames(ntc086.fData[ntc086.fData$value>=10,])

highRTS <- unique(c(highRTS0.84, highRTS0.86))

## Sample Overview
count_mat <- count(pData(demoData), slideName, sampleName, segment, aoi)
test_gr <- gather_set_data(count_mat, 2:3)
test_gr <- na.omit(test_gr)
plot Sankey
pdf("plots/sample_overview.pdf", height=15, width=7)
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = segment), alpha = 0.5, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(color = "white", size = 2) +
    theme_classic(base_size = 12) + 
    theme(legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) +
    labs(x = "", y = "") 
dev.off()



## QC
#Before we begin, we will shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
# Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <- list(minSegmentReads = 1000, # Minimum number of reads (1000)
                  percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
                  percentStitched = 80,   # Minimum % of reads stitched (80%)
                  percentAligned = 80,    # Minimum % of reads aligned (80%)
                  percentSaturation = 50, # Minimum sequencing saturation (50%)
                  minNegativeCount = 10,  # Minimum negative control counts (10)
                  maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
                  # minNuclei = 20,         # Minimum # of nuclei estimated (100)
                  minArea = 5000)         # Minimum segment area (5000)
demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)        
pData(demoData)[is.na(pData(demoData))]=0

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
QCResults <- QCResults[!is.na(QCResults$LowArea),]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))

col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() + guides(fill = "none") +
        facet_wrap(as.formula(paste("~", fill_by)), nrow = 5) +
        labs(x = annotation, y = "Segments, #", title = annotation)
    if(!is.null(scale_trans)) {
        plt <- plt +
            scale_x_continuous(trans = scale_trans)
    }
    plt
}

# calculate the negative geometric means for each module
negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]

# sData(demoData)$segment<-factor(sData(demoData)$segment, levels=c("CK+","CD45+", "Rest", "Control"))
pdf("plots/QC_histogram.pdf")
QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
QC_histogram(sData(demoData), "Aligned (%)", col_by, 80)
QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
QC_histogram(sData(demoData), "area", col_by, 5000, scale_trans = "log10")
# QC_histogram(sData(demoData), "nuclei", col_by, 20)

for(ann in negCols) {
    plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
    print(plt)
}
dev.off()


# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
ntcData <- demoData[, pData(demoData)$slideName == "No Template Control"]

## Finally we plot all of the QC Summary information in a table.
kable(QC_Summary, caption = "QC Summary Table for each Segment")
## As the final step in Segment QC, we remove flagged segments that do not meet our QC cutoffs.
demoData <- demoData[, QCResults$QCStatus == "PASS"]
# Subsetting our dataset has removed samples which did not pass QC
dim(demoData)
# Table: QC Summary Table for each Segment

# |              | Pass| Warning|
# |:-------------|----:|-------:|
# |LowReads      |  171|       0|
# |LowTrimmed    |  165|       6|
# |LowStitched   |  162|       9|
# |LowAligned    |  153|      18|
# |LowSaturation |  167|       4|
# |LowNegatives  |  116|      55|
# |LowArea       |  141|      30|
# |TOTAL FLAGS   |  105|      66|
# Features  Samples
#    18815      107
qcPassData <- demoData
filteredData <- demoData[!fData(demoData)$RTS_ID %in% highRTS,]
demoData <- filteredData



# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
    subset(demoData, 
           fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
               fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
# Features  Samples
#    18764      107
demoData <- ProbeQCPassed 

# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))
#> [1] 18635

# demoData <- demoData[rownames(pData(demoData)) %in% c("DSP-1001660009084-H-A01.dcc", "DSP-1001660009086-G-A01.dcc"), ]

# collapse to targets
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
#> Features  Samples 
#>    18635      107
exprs(target_demoData)[1:5, 1:2]



# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
        LOQ[, module] <-
            pmax(minLOQ,
                 pData(target_demoData)[, vars[1]] * 
                     pData(target_demoData)[, vars[2]] ^ cutoff)
    }
}
pData(target_demoData)$LOQ <- LOQ

## filtering based on LOQ
LOQ_Mat <- c()
for(module in modules) {
    ind <- fData(target_demoData)$Module == module
    Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                       FUN = function(x) {
                           x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]


## Segment gene detection
# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
    pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
    cut(pData(target_demoData)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
pdf("plots/geneDetectionRate0.pdf")
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = segment)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
dev.off()

## We can also create a table to review what kidney tissue type (DKD vs normal) is going to be impacted by each threshold:

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$sampleName))
# |       | diagnostic_A.L| diagnostic_B| diagnostic_C| NTC_084| relapse_A.R| surgery_A.L| surgery_A.R| surgery_B| surgery_C|
# |<1%    |              0|            0|            0|       0|           0|           0|           0|         0|         0|
# |1-5%   |              0|            1|            0|       1|           5|           1|           2|         8|         4|
# |5-10%  |              0|            3|            0|       0|           9|           6|           3|        19|        14|
# |10-15% |              1|            2|            1|       0|           0|           7|           1|         0|         3|
# |>15%   |              0|            0|            4|       0|           0|           8|           2|         0|         2|
target_demoData <-
    target_demoData[, pData(target_demoData)$GeneDetectionRate >= .05]

dim(target_demoData)
# Features  Samples
#    18635       85

## Re-plot the Sankey diagram showing our current working dataset. This is now a dataset that no longer contains segments flagged by Segment QC or that have low gene detection rates.

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat2 <- count(pData(target_demoData), slideName, sampleName, segment,  aoi)
test_gr2 <- gather_set_data(count_mat2, 2:3)
test_gr2 <- na.omit(test_gr2)
# plot Sankey
pdf("plots/sample_overview_woLowGeneDetectionRateSegments.pdf", height=15, width=7)
# pdf("plots/Rplots.pdf", height=25, width=7)
ggplot(test_gr2, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = segment), alpha = 0.5, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(color = "white", size = 2) +
    theme_classic(base_size = 12) + 
    theme(legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) +
    labs(x = "", y = "") 
dev.off()

## Gene detection rate
 # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Gene of interest detection table
# goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
#          "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
# goi_df <- data.frame(
#     Gene = goi,
#     Number = fData(target_demoData)[goi, "DetectedSegments"],
#     DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

pdf("plots/geneDetectionRate.pdf")
ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")
dev.off()



# Subset to target genes detected in at least 5-10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- 
    target_demoData[fData(target_demoData)$DetectionRate >= 0.05 |
                        fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)
#> Features  Samples 
#>    6065      85

# retain only detected genes of interest
# goi <- goi[goi %in% rownames(target_demoData)]




library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "sampleName"
Stat_data <- 
    data.frame(row.names = colnames(exprs(target_demoData)),
               Segment = colnames(exprs(target_demoData)),
               Annotation = pData(target_demoData)[, ann_of_interest],
               Q3 = unlist(apply(exprs(target_demoData), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(target_demoData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))

pdf("plots/q3VsNegProbe.pdf", height=8, width=12)
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
dev.off()


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")



# visualize the first 10 segments with each normalization method
pdf("plots/normalizationEffect10segments.pdf")
boxplot(exprs(target_demoData)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")
boxplot(assayDataElement(target_demoData[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")
boxplot(assayDataElement(target_demoData[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")
dev.off()


### Dimensionality reduction

library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
    umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),  
         config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
set.seed(42) # set the seed for tSNE as well
tsne_out <-
    Rtsne(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),
          perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]


saveRDS(target_demoData, "rds/target_demoData.rds")
target_demoData <- readRDS("rds/target_demoData.rds")

scheme1<-c("0"="black", "diagnostic_A.L"="#0084A5", "surgery_A.L"="#660044", "surgery_A.R"="#B30078", "relapse_A.R"="#7D21C0", 
    "diagnostic_B"="#00C4F5", "surgery_B"="#FF0AAD", 
    "diagnostic_C"="#85E7FF", "surgery_C"="#FF70CF")
scheme2<-c("CK+"="red", "CD45+"="green", "Rest"="grey", "0"="black")
scheme3 <- c("0"="black", "D1L"="#0084A5", "S1L"="#660044", "S1R"="#B30078", "R1R"="#7D21C0", 
    "D2"="#00C4F5", "S2"="#FF0AAD", 
    "D3"="#85E7FF", "S3"="#FF70CF")

pdf("plots/ckTreatment.pdf", height=3, width=3)
ggplot(pData(target_demoData)[pData(target_demoData)$segment == "CK+",],
       aes(x = UMAP1, y = UMAP2, color = sampleNameNew)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_color_manual(values=scheme3) +
    theme(legend.position="bottom")
ggplot(pData(target_demoData)[pData(target_demoData)$segment == "CK+",],
       aes(x = UMAP1, y = UMAP2, color = sampleNameNew)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_color_manual(values=scheme3) +
    theme(legend.position="none")
dev.off()

pdf("plots/dimensionalityReductionClusters.pdf", height=8, width=10)
ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = sampleName, shape = segment)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_color_manual(values=scheme1) +
    theme(legend.position="bottom")
ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = sampleName, shape = segment)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_color_manual(values=scheme1) +
    theme(legend.position="bottom")
dev.off()
pdf("plots/dimensionalityReductionClusters_separate_labelled.pdf", height=12, width=25)
ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = sampleName, shape = segment, label = roiNew)) +
    geom_point(size = 3) +
    geom_text(hjust = 0, nudge_x = 0.05, size=3) +
    theme_bw() +
    facet_wrap(.~segment) +
    scale_color_manual(values=scheme1) +
    theme(legend.position="bottom", 
        axis.text = element_text(size=15),
        legend.text = element_text(size=15),
        strip.text = element_text(size = 15))
dev.off()

ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = sampleName, shape = segment)) +
    geom_point(size = 3) +
    theme_bw() +
    facet_wrap(.~segment) +
    scale_color_manual(values=scheme1) +
    theme(legend.position="bottom")
dev.off()

library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_demoData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#      CTSD  ANKRD30A   SCGB2A2     NPY1R       LYZ
# 0.3154273 0.3019893 0.2738959 0.2577853 0.2551979

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]

allGenes <- fData(target_demoData)$TargetName
sanityChecks <- c("KRT7", "KRT8", "KRT18", "KRT19", "PTPRC", "MIB1", "ESR1", "FOXA1", "PGR", "GREB1", "MYD88", "ANXA2", "BASP1", "S100A10")
pa.sub <- c("ANXA2", "BASP1", "S100A10")
paUp <- read.table("../src/pa.up.new.txt", header=F)[,1]
paDown <- read.table("../src/pa.down.new.txt", header=F)[,1]
g2m <- c("ABL1", "AMD1", "ARID4A", "ATF5", "ATRX", "AURKA", "AURKB", "BARD1", "BCL3", "BIRC5", "BRCA2", "BUB1", "BUB3", "CASP8AP2", "CBX1", "CCNA2", "CCNB2", "CCND1", "CCNF", "CCNT1", "CDC20", "CDC25A", "CDC25B", "CDC27", "CDC45", "CDC6", "CDC7", "CDK1", "CDK4", "CDKN1B", "CDKN2C", "CDKN3", "CENPA", "CENPE", "CENPF", "CHAF1A", "CHEK1", "CHMP1A", "CKS1B", "CKS2", "CTCF", "CUL1", "CUL3", "CUL4A", "CUL5", "DBF4", "DDX39A", "DKC1", "DMD", "DR1", "DTYMK", "E2F1", "E2F2", "E2F3", "E2F4", "EFNA5", "EGF", "ESPL1", "EWSR1", "EXO1", "EZH2", "FANCC", "FBXO5", "FOXN3", "G3BP1", "GINS2", "GSPT1", "H2AX", "H2AZ1", "H2AZ2", "H2BC12", "HIF1A", "HIRA", "HMGA1", "HMGB3", "HMGN2", "HMMR", "HNRNPD", "HNRNPU", "HOXC10", "HSPA8", "HUS1", "ILF3", "INCENP", "JPT1", "KATNA1", "KIF11", "KIF15", "KIF20B", "KIF22", "KIF23", "KIF2C", "KIF4A", "KIF5B", "KMT5A", "KNL1", "KPNA2", "KPNB1", "LBR", "LIG3", "LMNB1", "MAD2L1", "MAP3K20", "MAPK14", "MARCKS", "MCM2", "MCM3", "MCM5", "MCM6", "MEIS1", "MEIS2", "MKI67", "MNAT1", "MT2A", "MTF2", "MYBL2", "MYC", "NASP", "NCL", "NDC80", "NEK2", "NOLC1", "NOTCH2", "NSD2", "NUMA1", "NUP50", "NUP98", "NUSAP1", "ODC1", "ODF2", "ORC5", "ORC6", "PAFAH1B1", "PBK", "PDS5B", "PLK1", "PLK4", "PML", "POLA2", "POLE", "POLQ", "PRC1", "PRIM2", "PRMT5", "PRPF4B", "PTTG1", "PTTG3P", "PURA", "RACGAP1", "RAD21", "RAD23B", "RAD54L", "RASAL2", "RBL1", "RBM14", "RPA2", "RPS6KA5", "SAP30", "SFPQ", "SLC12A2", "SLC38A1", "SLC7A1", "SLC7A5", "SMAD3", "SMARCC1", "SMC1A", "SMC2", "SMC4", "SNRPD1", "SQLE", "SRSF1", "SRSF10", "SRSF2", "SS18", "STAG1", "STIL", "STMN1", "SUV39H1", "SYNCRIP", "TACC3", "TENT4A", "TFDP1", "TGFB1", "TLE3", "TMPO", "TNPO2", "TOP1", "TOP2A", "TPX2", "TRA2B", "TRAIP", "TROAP", "TTK", "UBE2C", "UBE2S", "UCK2", "UPF1", "WRN", "XPO1", "YTHDC1")

pdf("plots/heatmapSanityChecks.pdf", height=7, width=12)
pheatmap(assayDataElement(target_demoData[sanityChecks, ], elt = "log_q"),
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = FALSE,
         cutree_cols = 2,
         cellwidth = 7,
         cellheight = 15,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = viridis(5, direction = -1),
         annotation_col = pData(target_demoData)[, c("segment", "sampleName")], 
         annotation_colors = list("sampleName"=scheme1, "segment"=scheme2),
         border_color = "grey60")
dev.off()

target_demoData <- readRDS("rds/target_demoData.rds")
target_demoData <- readRDS("rds/target_demoData_all.rds")
sanityCheckExp<-as.data.frame(assayDataElement(target_demoData[sanityChecks, ], elt = "log_q"))
pData(target_demoData) <- merge(pData(target_demoData), t(sanityCheckExp), by=0, all=TRUE)
rownames(pData(target_demoData)) <- pData(target_demoData)[,1]

require("MatrixGenerics")
require("ggpubr")
require("viridis")
allExp <- as.data.frame(assayDataElement(target_demoData, elt = "log_q"))

getZscores <- function(mat){
    matrix<-as.matrix(mat)
    matrix  <- matrix - rowMeans(matrix)
    sd<-rowSds(matrix)
    matrix<-matrix/sd
    matrix
}
allExp.Z <- getZscores(allExp)

pa.sub.exp <- allExp.Z[pa.sub,]
data$PAsub <- colMeans(pa.sub.exp)

paUp.exp <- allExp.Z[rownames(allExp.Z) %in% paUp,]
data$PAup <- colMeans(paUp.exp)

paDown.exp <- allExp.Z[rownames(allExp.Z) %in% paDown,]
data$PAdown <- colMeans(paDown.exp)

enrichSignature <- function(data, sign, name){
    exp <- allExp.Z[rownames(allExp.Z) %in% sign,]
    data[,name] <- colMeans(exp)
    data
}
data <- as.data.frame(pData(target_demoData))


dataSign <- enrichSignature(data, paUp, "DormancyUp")
dataSign <- enrichSignature(dataSign, paDown, "DormancyDown")
dataSign <- enrichSignature(dataSign, g2m, "G2M")

roiRename <- read.table("rds/roiRename.txt")
dataSign$roiRename <- roiRename[match(dataSign$roiNew, roiRename$V1),]$V2

pdf("plots/ckSignatures.pdf", height=4, width=4)
signatures <- lapply(c("DormancyUp", "DormancyDown", "G2M"), function(x){
    ggplot(dataSign[dataSign$segment == "CK+",],
    # aes_string(x = "UMAP1", y = "UMAP2", color = x, label = "roiRename")) +
    aes_string(x = "UMAP1", y = "UMAP2", color = x, shape = "patient")) +
    geom_point(size = 2) +
    # geom_text_repel(size = 2) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    # labs(title = x) +
    # facet_wrap(.~segment) +
    theme(legend.position="none")
    })
dev.off()
treatment <- ggplot(pData(target_demoData)[pData(target_demoData)$segment == "CK+",],
       aes(x = UMAP1, y = UMAP2, color = sampleNameNew, , shape = patient)) +
    geom_point(size = 2) +
    theme_bw() +
    scale_color_manual(values=scheme3) +
    theme(legend.position="none")


pdf("plots/figure1c.unlabelled.pdf", height=6, width=6)
# ggarrange(treatment, signatures[[1]], signatures[[2]], signatures[[3]], ncol=2, nrow=2)
cowplot::plot_grid(treatment + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                    signatures[[1]]+ theme(axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank(),
                              plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                    signatures[[2]] + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                    signatures[[3]] + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                    ncol=2, nrow=2)
dev.off()

pdf("plots/ckSignaturesLegend.pdf", height=4, width=4)
lapply(c("DormancyUp", "DormancyDown", "G2M"), function(x){
    ggplot(dataSign[dataSign$segment == "CK+",],
    aes_string(x = "UMAP1", y = "UMAP2", color = x, label = "roiRename")) +
    geom_point(size = 2) +
    geom_text_repel(size = 2) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    labs(title = x) +
    # facet_wrap(.~segment) +
    theme(legend.position="bottom")
    })
dev.off()




sanityChecksPlots <- lapply(sanityChecks, function(x){
    data[,x] <- as.numeric(data[,x])
    ggplot(data,
       aes_string(x = "UMAP1", y = "UMAP2", color = x, shape = "segment", label = "roi")) +
    geom_point(size = 3) +
    geom_text(hjust = 0, nudge_x = 0.05, size=3) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    # scale_color_gradientn(colors=viridis(100, begin = 0.2, end = 1, direction = -1, option = "B"))
    labs(title = x) +
    theme(legend.position="bottom")
    })
sanityChecksPlotsSplit <- lapply(sanityChecks, function(x){
    data[,x] <- as.numeric(data[,x])
    ggplot(data,
       aes_string(x = "UMAP1", y = "UMAP2", color = x, shape = "segment")) +
    geom_point(size = 3) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    labs(title = x) +
    facet_wrap(.~segment) +
    theme(legend.position="bottom")
    })
pdf("plots/sanityChecksPlots.pdf", height=20, width=22)
ggarrange(plotlist = sanityChecksPlots, ncol = 4, nrow=4)
dev.off()
pdf("plots/sanityChecksPlotsSplit.pdf", height=16, width=22)
ggarrange(plotlist = sanityChecksPlotsSplit, ncol=4, nrow=4)
dev.off()

pdf("plots/sanityChecksPlotsSplitPAsub.pdf", height=10, width=18)
    ggplot(data,
       aes_string(x = "UMAP1", y = "UMAP2", color = "PAup", shape = "segment", label = "roiNew")) +
    geom_point(size = 3) +
    geom_text_repel(size = 2) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    labs(title = "PAup") +
    # facet_wrap(.~segment) +
    theme(legend.position="bottom")

    ggplot(data,
       aes_string(x = "UMAP1", y = "UMAP2", color = "PAdown", shape = "segment", label = "roiNew")) +
    geom_point(size = 3) +
    geom_text_repel(size = 2) +
    theme_bw() +
    scale_color_viridis(direction=-1) +
    labs(title = "PAdown") +
    # facet_wrap(.~segment) +
    theme(legend.position="bottom")
dev.off()

expMat <- as.data.frame(assayDataElement(target_demoData, elt = "log_q"))
write.table(expMat, "rds/expMat.txt", sep="\t",quote=F,row.names=T,col.names=T)

