library(DiffBind)
library(ggplot2)

dname       <- "H4K20me3_abcam"
numMaxSites <- 2000
threshold   <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)

# reading the data table
treatments <- dba(sampleSheet = "samples.diffBind.csv")

# remove the blacklist
treatments <- dba.blacklist(treatments, blacklist = DBA_BLACKLIST_HG38, greylist = FALSE)

# counting reads from bam files
treatments <- dba.count(treatments, summits = FALSE, bUseSummarizeOverlaps = TRUE)

# establishing a contrast
treatments <- dba.contrast(treatments, categories = DBA_CONDITION, minMembers = 2)

# performing the differential enrichment analysis
# default method: DESeq2
# alternative method: edgeR
treatments <- dba.analyze(treatments, method = DBA_ALL_METHODS)

system("mkdir -p diffBind.dupMark.summitFALSE")

# plot a correlation heatmap among samples
pdf(paste0("diffBind.dupMark.summitFALSE/diffBind.summitFALSE.", dname, ".pdf"), height = 6.4, width = 6.4)
plot(treatments)
dev.off()

# PCA plot
pdf(paste0("diffBind.dupMark.summitFALSE/PCA.summitFALSE.", dname, ".pdf"), height = 4.8, width = 6.4)
dba.plotPCA(treatments, DBA_CONDITION, label = DBA_ID)
dev.off()

# comparing DESeq2 and edgeR results
pdf("diffBind.dupMark.summitFALSE/comparison.diffBind.summitFALSE.methods.FDR.pdf", width = 6.4, height = 4.8)
for (i in 1:length(threshold)) {
	dba.plotVenn(
		treatments,
		contrast = 1,
		method   = DBA_ALL_METHODS,
		th       = threshold[i],
		main     = paste0("FDR < ", threshold[i])
	)
}
dev.off()

pdf("diffBind.dupMark.summitFALSE/comparison.diffBind.summitFALSE.methods.pvalue.pdf", width = 6.4, height = 4.8)
for (i in 1:length(threshold)) {
	dba.plotVenn(
		treatments,
		contrast = 1,
		method   = DBA_ALL_METHODS,
		bUsePval = TRUE,
		th       = threshold[i],
		main     = paste0("pvalue < ", threshold[i])
    )
}
dev.off()

pdf("diffBind.dupMark.summitFALSE/comparison.diffBind.summitFALSE.DESeq2.pdf", width = 6.4, height = 4.8)
dba.plotMA(treatments,  method = DBA_DESEQ2)
dba.plotMA(treatments,  method = DBA_DESEQ2, bXY = TRUE)
dba.plotBox(treatments, method = DBA_DESEQ2)
dev.off()

pdf("diffBind.dupMark.summitFALSE/comparison.diffBind.summitFALSE.edgeR.pdf", width = 6.4, height = 4.8)
dba.plotMA(treatments,  method = DBA_EDGER)
dba.plotMA(treatments,  method = DBA_EDGER, bXY = TRUE)
dba.plotBox(treatments, method = DBA_EDGER)
dev.off()

# output diffBind tables
treatments.DB.DESeq2 <- data.frame(
							dba.report(
								treatments,
                                method   = DBA_DESEQ2,
                                contrast = 1,
                                th       = 1,
                                bUsePval = TRUE
							)
                        )

treatments.DB.edgeR  <- data.frame(
                            dba.report(
                                treatments,
                                method   = DBA_EDGER,
                                contrast = 1,
                                th       = 1,
                                bUsePval = TRUE
							)
                    	)

write.table(
	treatments.DB.DESeq2,
	paste0("diffBind.dupMark.summitFALSE/DB.summitFALSE.", dname, ".DESeq2.csv"),
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
)

write.table(
	treatments.DB.edgeR,
	paste0("diffBind.dupMark.summitFALSE/DB.summitFALSE.", dname, ".edgeR.csv"),
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
)

## annotation in the output file:
# Conc: mean read concentration over all the samples 
#       (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted)
# Fold: shows the difference in mean concentrations between the two groups, with a positive value indicating
#       increased binding affinity in the treatment group and a negative value indicating increased binding
#       affinity in the control group.

diffPeakCounts.DESeq2           <- data.frame(matrix(0, length(threshold), 2))
colnames(diffPeakCounts.DESeq2) <- c("pvalue", "FDR")
diffPeakCounts.DESeq2           <- cbind(threshold, diffPeakCounts.DESeq2)

for (i in 1:length(threshold)) {
	diffPeakCounts.DESeq2[i, 2] <- length(which(treatments.DB.DESeq2$p.value < threshold[i]))
	diffPeakCounts.DESeq2[i, 3] <- length(which(treatments.DB.DESeq2$FDR     < threshold[i]))
}
diffPeakCounts.DESeq2 <- t(diffPeakCounts.DESeq2)
diffPeakCounts.DESeq2 <- transform(diffPeakCounts.DESeq2, method = "DESeq2")

diffPeakCounts.edgeR           <- data.frame(matrix(0, length(threshold), 2))
colnames(diffPeakCounts.edgeR) <- c("pvalue", "FDR")
diffPeakCounts.edgeR           <- cbind(threshold, diffPeakCounts.edgeR)

for (i in 1:length(threshold)) {
	diffPeakCounts.edgeR[i, 2] <- length(which(treatments.DB.edgeR$p.value < threshold[i]))
	diffPeakCounts.edgeR[i, 3] <- length(which(treatments.DB.edgeR$FDR     < threshold[i]))
}
diffPeakCounts.edgeR <- t(diffPeakCounts.edgeR)
diffPeakCounts.edgeR <- transform(diffPeakCounts.edgeR, method = "edgeR")

diffPeakCounts           <- rbind(diffPeakCounts.DESeq2[2:3, ], diffPeakCounts.edgeR[2:3, ])
colnames(diffPeakCounts) <- c(t(threshold), "methods")

write.table(
	diffPeakCounts,
	paste0("diffBind.dupMark.summitFALSE/diffPeakCounts.summitFALSE.", dname, ".csv"),
	col.names = TRUE,
	row.names = TRUE,
	quote     = FALSE,
	sep       = "\t"
)

system("mkdir -p diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/pdf/")
system("mkdir -p diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/png/")

volcano.x <- ceiling(max(abs(treatments.DB.DESeq2$Fold)))
for (i in 1:length(threshold)) {
	p <- ggplot(treatments.DB.DESeq2, aes(x = Fold, y = -1 * log10(FDR))) +
			geom_point(aes(color = ifelse(FDR < threshold[i], "red", "black"))) +
			scale_colour_manual(labels = c(paste0("> ", threshold[i]), paste0("< ", threshold[i])), values = c("black", "red")) +
			labs(color = "FDR") +
			xlim(-volcano.x, volcano.x) +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
    			panel.grid.major = element_blank(),
    			panel.grid.minor = element_blank(),
				axis.line        = element_line(colour = "black")
			)

	pdf(paste0("diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/pdf/volcano.summitFALSE.", dname, "DESeq2.fdr_", threshold[i] , ".pdf"), height = 6, width = 6)
	print(p)
	dev.off()

	png(paste0("diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/png/volcano.summitFALSE.", dname, "DESeq2.fdr_", threshold[i] , ".png"), height = 600, width = 600)
	print(p)
	dev.off()

	p <- ggplot(treatments.DB.DESeq2, aes(x = Fold, y = -1 * log10(p.value))) +
			geom_point(aes(color = ifelse(p.value < threshold[i], "red", "black"))) +
			scale_colour_manual(labels = c(paste0("> ", threshold[i]), paste0("< ", threshold[i])), values = c("black", "red")) + 
			labs(color = "pvalues") +
			xlim(-volcano.x, volcano.x) +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				panel.grid.major = element_blank(),
    			panel.grid.minor = element_blank(),
    			axis.line        = element_line(colour = "black")
			)

	pdf(paste0("diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/pdf/volcano.summitFALSE.", dname, "DESeq2.pvalue_", threshold[i] , ".pdf"), height = 6, width = 6)
	print(p)
	dev.off()

	png(paste0("diffBind.dupMark.summitFALSE/volcano_plot/DESeq2/png/volcano.summitFALSE.", dname, "DESeq2.pvalue_", threshold[i] , ".png"), height = 600, width = 600)
	print(p)
	dev.off()
}

system("mkdir -p diffBind.dupMark.summitFALSE/volcano_plot/edgeR/pdf/")
system("mkdir -p diffBind.dupMark.summitFALSE/volcano_plot/edgeR/png/")

volcano.x <- ceiling(max(abs(treatments.DB.edgeR$Fold)))
for (i in 1:length(threshold)) {
	p <- ggplot(treatments.DB.edgeR, aes(x = Fold, y = -1 * log10(FDR))) +
			geom_point(aes(color = ifelse(FDR < threshold[i], "red", "black"))) +
			scale_colour_manual(labels = c(paste0("> ", threshold[i]), paste0("< ", threshold[i])), values = c("black", "red")) + 
			labs(color = "FDR") +
			xlim(-volcano.x, volcano.x) +
			theme_bw() +
			theme(
    			panel.border     = element_blank(),
    			panel.grid.major = element_blank(),
    			panel.grid.minor = element_blank(),
    			axis.line        = element_line(colour = "black")
    		)

	pdf(paste0("diffBind.dupMark.summitFALSE/volcano_plot/edgeR/pdf/volcano.summitFALSE.", dname, "edgeR.fdr_", threshold[i] , ".pdf"), height = 6, width = 6)
	print(p)
	dev.off()

	png(paste0("diffBind.dupMark.summitFALSE/volcano_plot/edgeR/png/volcano.summitFALSE.", dname, "edgeR.fdr_", threshold[i] , ".png"), height = 600, width = 600)
	print(p)
	dev.off()

	p <- ggplot(treatments.DB.edgeR, aes(x = Fold, y = -1 * log10(p.value))) +
			geom_point(aes(color = ifelse(p.value < threshold[i], "red", "black"))) +
			scale_colour_manual(labels = c(paste0("> ", threshold[i]), paste0("< ", threshold[i])), values = c("black", "red")) +
			labs(color = "pvalues") +
			xlim(-volcano.x, volcano.x) +
			theme_bw() +
			theme(
    			panel.border     = element_blank(),
    			panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
    			axis.line        = element_line(colour = "black")
    		)

	pdf(paste0("diffBind.dupMark.summitFALSE/volcano_plot/edgeR/pdf/volcano.summitFALSE.", dname, "edgeR.pvalue_", threshold[i] , ".pdf"), height = 6, width = 6)
	print(p)
	dev.off()

	png(paste0("diffBind.dupMark.summitFALSE/volcano_plot/edgeR/png/volcano.summitFALSE.", dname, "edgeR.pvalue_", threshold[i] , ".png"), height = 600, width = 600)
	print(p)
	dev.off()
}

###################################################
pdf(paste0("diffBind.dupMark.summitFALSE/heatmap.summitFALSE.", dname, ".pdf"), height = 9.6, width = 6.4)
corvals <- dba.plotHeatmap(treatments, method = DBA_EDGER,  correlations = FALSE)
corvals <- dba.plotHeatmap(treatments, method = DBA_DESEQ2, correlations = FALSE, maxSites = numMaxSites)
dev.off()

corvals <- data.frame(corvals)
write.table(
	corvals,
	paste0("diffBind.dupMark.summitFALSE/top", numMaxSites, ".bed"),
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
)

peakset <- data.frame(dba.peakset(treatments, bRetrieve = TRUE))
write.table(
	peakset,
	paste0("diffBind.dupMark.summitFALSE/peakset.", dname, ".csv"),
	col.names = FALSE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
)

samples   <- treatments$samples[, 1]
lib.size  <- treatments$norm$DESeq2$lib.sizes
norm.facs <- treatments$norm$DESeq2$norm.facs
FRIP      <- treatments$SN
out.table <- cbind(samples, lib.size, norm.facs, FRIP)

write.table(
	out.table,
	paste0("diffBind.dupMark.summitFALSE/stats.summitFALSE.", dname, ".csv"),
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
)

save.image(paste0("diffBind.dupMark.summitFALSE/diffBind.", dname, ".RData"))

sessionInfo()
