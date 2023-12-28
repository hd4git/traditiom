pa.up<-read.table("rds/pa.up.new.txt", header=F, stringsAsFactors=F)[,1]
pa.down<-read.table("rds/pa.down.new.txt", header=F, stringsAsFactors=F)[,1]


ordered.sub1.high.info <- read.table("rds/infoDormancyUp.txt", sep="\t")
ordered.sub2.high.info <- read.table("rds/infoDormancyDown.txt", sep="\t")

exp <- rownames(fData(target_demoData[sData(target_demoData)$segment=="CK+"]))
write.table(exp, "rds/geomx.ck.genes.txt", sep="\t", quote=F, row.names=F, col.names=F)
expUp <- exp[exp %in% paUp]
expDown <- exp[exp %in% paDown]

listUp1 <- list("PA.Up" = pa.up,
				"Dormancy.Up" = rownames(ordered.sub1.high.info))
listDown1 <- list("PA.Down" = pa.down,
				"Dormancy.Down" = rownames(ordered.sub2.high.info))
a <- ggvenn(listUp, show_percentage= FALSE, auto_scale = TRUE, fill_color = c("#C085EA",  "#50157A"), set_name_size = 4)
b <- ggvenn(listDown, show_percentage= FALSE, auto_scale = TRUE, fill_color = c("#FFFFAD", "yellow"), set_name_size = 4)
ggarrange(a,b, ncol=1, nrow=2)
dev.off()


listUp2 <- list("Geomx.CK+" = exp,
				"PA.Up" = pa.up,
				"Dormancy.Up" = rownames(ordered.sub1.high.info))
listDown2 <- list("Geomx.CK+" = exp,
				"PA.Down" = pa.down,
				"Dormancy.Down" = rownames(ordered.sub2.high.info))
a <- ggvenn(listUp2, show_percentage= FALSE, fill_color = c("grey", "#C085EA",  "#50157A"), set_name_size = 4)
b <- ggvenn(listDown2, show_percentage= FALSE, fill_color = c("grey", "#FFFFAD", "yellow"), set_name_size = 4)
pdf("geomx.ck.venn.pdf", height = 3, width = 6)
ggarrange(a,b, ncol=2, nrow=1)
dev.off()
