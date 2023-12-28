target_demoData_all <- readRDS("rds/target_demoData_all.rds")
getDE <- function(selSegment, comp, group, formula, selSamples){
    target_demoData <- target_demoData_all
    # target_demoData <- target_demoData[,pData(target_demoData)$segment == selSegment & pData(target_demoData)$sampleNameNew %in% comp]
    target_demoData <- target_demoData[,pData(target_demoData)$segment == selSegment & pData(target_demoData)$region %in% comp & pData(target_demoData)$sampleNameNew %in% selSamples]
    pData(target_demoData)[["sample"]] <- factor(pData(target_demoData)[["sampleNameNew"]])
    pData(target_demoData)[["roi2"]] <- factor(pData(target_demoData)[["roi"]])

    pData(target_demoData)$region <- factor(pData(target_demoData)$region, c("core", "edge"))
    pData(target_demoData)$side <- factor(pData(target_demoData)$side, c("right", "left"))
    pData(target_demoData)$type <- factor(pData(target_demoData)$type, c("relapse", "surgery", "diagnostic"))
    assayDataElement(object = target_demoData, elt = "log_q") <-
        assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")
    mixedOutmc <- mixedModelDE(target_demoData,
                     elt = "log_q",
                     # modelFormula = ~ testRegion + (1 + testRegion | sample),
                     modelFormula = formula,
                     groupVar = group,
                     # nCores = parallel::detectCores(),
                     nCores = 6,
                     multiCore = FALSE)

        r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
        tests <- rownames(r_test)
        r_test <- as.data.frame(r_test)
        r_test$Contrast <- tests
        
        # use lapply in case you have multiple levels of your test factor to
        # correctly associate gene name with it's row in the results table
        r_test$Gene <- 
            unlist(lapply(colnames(mixedOutmc),
                          rep, nrow(mixedOutmc["lsmeans", ][[1]])))
        r_test$Subset <- selSegment 
        r_test$FDR <- p.adjust(r_test$'Pr(>|t|)', method = "fdr")
        r_test <- r_test[, c("Gene","Contrast", "Estimate", 
                             "Pr(>|t|)", "FDR")]
        results <- r_test


    colnames(results)[4] <- "Pr"

    # Categorize Results based on P-value & FDR for plotting
    results$Color <- "NS or FC < 0.5"
    results$Color[results$Pr < 0.05] <- "P < 0.05"
    results$Color[results$FDR < 0.05] <- "FDR < 0.05"
    results$Color[results$FDR < 0.001] <- "FDR < 0.001"
    results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
    results$Color <- factor(results$Color,
                            levels = c("NS or FC < 0.5", "P < 0.05",
                                       "FDR < 0.05", "FDR < 0.001"))

    # pick top genes for either side of volcano to label
    # order genes for convenience:

    results$invert_P <- (-log10(results$Pr)) * sign(results$Estimate)
    results <- results[, -1*ncol(results)] # remove invert_P from matrix
results
}

plotDE <- function(results, selSegment, comp, group, formula){
    # Graph results
    target_demoData <- target_demoData_all
    target_demoData <- target_demoData[,pData(target_demoData)$segment == selSegment & pData(target_demoData)$sampleNameNew %in% comp]
    pData(target_demoData)[["sample"]] <- factor(pData(target_demoData)[["sampleNameNew"]])

    pData(target_demoData)$side <- factor(pData(target_demoData)$side, c("right", "left"))
    pData(target_demoData)$type <- factor(pData(target_demoData)$type, c("relapse", "surgery", "diagnostic"))
    filename <- paste("plots/DE/",comp[1],"Vs",comp[2],".", selSegment,"2.pdf", sep="")
    pdf(filename, height=20, width=10)
        # facet_wrap(~Subset, scales = "free_y")

    GOI <- unique(subset(results, (abs(results$Estimate)>0.5 & FDR < 0.05) | (abs(results$Estimate)>0.5 & Pr < 0.0005))$Gene)
    pheatmap(log2(assayDataElement(target_demoData[GOI, ], elt = "q_norm")),
             scale = "row", 
             show_rownames = TRUE, show_colnames = FALSE,
             border_color = NA,
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             cutree_cols = 2, cutree_rows = 2,
             main = paste("Segment=",selSegment, ", Comparison=",comp[1], "Vs", comp[2], sep=""),
             # breaks = seq(-3, 3, 0.05),
             # color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
             annotation_col = pData(target_demoData)[, c("type", "region", "sampleNameNew")])
   ggplot(results,
           aes(x = Estimate, y = -log10(Pr),
               color = Color, label = Gene)) +
        geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
        geom_hline(yintercept = -log10(0.05), lty = "dashed") +
        geom_point() +
        labs(title = paste("Segment=",selSegment, ", Comparison=",comp[1], "Vs", comp[2], sep=""),
            x = paste("Enriched in", comp[2], "<- log2(FC) -> Enriched in", comp[1]),
             y = "Significance, -log10(P)",
             color = "Significance") +
        scale_color_manual(values = c('FDR < 0.001' = "dodgerblue",
                                      'FDR < 0.05' = "lightblue",
                                      'P < 0.05' = "orange2",
                                      'NS or FC < 0.5' = "gray"),
                           guide = guide_legend(override.aes = list(size = 4))) +
        scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
        geom_text_repel(data = subset(results, (abs(results$Estimate)>0.5 & FDR < 0.05) | (abs(results$Estimate)>0.5 & Pr < 0.0005) ),
                        size = 3, point.padding = 0.15, 
                        min.segment.length = 0.1, box.padding = 0.2, lwd = 2,
                        max.overlaps = 50
                        ) +
        theme_bw(base_size = 16) +
        theme(legend.position = "bottom") #+
    # ggarrange(plotlist = list(a, b), nrow=2)
}


deS1RvsS1L.ck <- getDE("CK+", c("S1R", "S1L"), "side", ~ side + (1 | sample))
deS1RvsR1R.ck <- getDE("CK+", c("R1R", "S1R"), "type", ~ type + (1 | sample))
deS1LvsR1R.ck <- getDE("CK+", c("R1R", "S1L"), "type", ~ type + (1 | sample))
deD2vsS2.ck <- getDE("CK+", c("S2", "D2"), "type", ~ type + (1 | sample))
deD3vsS3.ck <- getDE("CK+", c("S3", "D3"), "type", ~ type + (1 | sample))
deS1LvsR1R.cd <- getDE("CD45+", c("R1R", "S1L"), "type", ~ type + (1 | sample))
deS1RvsS1L.rest <- getDE("Rest", c("S1R", "S1L"), "side", ~ side + (1 | sample))
deS1LvsR1R.rest <- getDE("Rest", c("R1R", "S1L"), "type", ~ type + (1 | sample))

deCoreVsEdge.ck <- getDE("CK+", c("core", "edge"), "region", ~ region + (1 | sample), c("S2", "S3" ))


degList <- list("S1RvsS1L.ck" = deS1RvsS1L.ck,
                "S1RvsR1R.ck" = deS1RvsR1R.ck,
                "S1LvsR1R.ck" = deS1LvsR1R.ck,
                "D2vsS2.ck" = deD2vsS2.ck,
                "D3vsS3.ck" = deD3vsS3.ck,
                "S1LvsR1R.cd" = deS1LvsR1R.cd,
                "S1RvsS1L.rest" = deS1RvsS1L.rest,
                "S1LvsR1R.rest" = deS1LvsR1R.rest)

deD3vsS3 <- rbind(sigDEGs[[1]][[5]], sigDEGs[[2]][[5]])
plotDE(deD3vsS3, "CK+", c("S3", "D3"), "type", ~ type + (1 | sample))
dev.off()

plotDE(deS1RvsS1L, "CK+", c("S1R", "S1L"), "side", ~ side + (1 | sample))
dev.off()
plotDE(deS1RvsR1R, "CK+", c("R1R", "S1R"), "type", ~ type + (1 | sample))
dev.off()
plotDE(deS1LvsR1R, "CK+", c("R1R", "S1L"), "type", ~ type + (1 | sample))
dev.off()
plotDE(deD2vsS2, "CK+", c("S2", "D2"), "type", ~ type + (1 | sample))
dev.off()
plotDE(deD3vsS3, "CK+", c("S3", "D3"), "type", ~ type + (1 | sample))
dev.off()
plotDE(deS1LvsR1R.cd, "CD45+", c("R1R", "S1L"), "type", ~ type + (1 | sample))
dev.off()
plotDE(deS1RvsS1L.rest, "Rest", c("S1R", "S1L"), "side", ~ side + (1 | sample))
dev.off()
plotDE(deS1LvsR1R.rest, "Rest", c("R1R", "S1L"), "type", ~ type + (1 | sample))
dev.off()

# resultsList <- list("deS1RvsS1L" = deS1RvsS1L, "deS1RvsR1R"=deS1RvsR1R, "deS1LvsR1R"=deS1LvsR1R)
# sigDEGs <- lapply(resultsList, function(results){
#     subset(results, (abs(results$Estimate)>0.5 & FDR < 0.05) | (abs(results$Estimate)>0.5 & Pr < 0.0005) )$Gene    
#     })
sigDEGsPos <- lapply(degList, function(results){
    subset(results, results$Estimate>0.5 & Pr < 0.005)   
    })
sigDEGsNeg <- lapply(degList, function(results){
    subset(results, results$Estimate< -0.5 & Pr < 0.005)  
    })

saveRDS(list("sigDEGsPos"= sigDEGsPos, "sigDEGsNeg"=sigDEGsNeg), "rds/sigDEGs.RDS")

