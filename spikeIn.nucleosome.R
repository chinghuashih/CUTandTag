library("gplots")

# Notes:
# EpiCypher considers an antibody with <20% binding to all off-target PTMs specific 
# and suitable for downstream data analysis. For IgG, data is normalized to the sum
# of total barcode reads.

plot_heatmap <- function (input_matrix, out_name) {
	pdf(paste0(out_name, ".pdf"), height = 6, width = 20)
	heatmap.2(
		input_matrix,
		cellnote     = input_matrix,
		notecol      = "black",
		na.rm        = TRUE,
		dendrogram   = 'none',
		Rowv         = FALSE,
		Colv         = FALSE,
		col          = rev(heat.colors(100)),
		trace        = 'none',
		main         = "specificity",
		xlab         = "barcode identify of spike-in",
		key          = TRUE,
		keysize      = 1,
		key.title    = NA,
		key.ylab     = NA,
		key.xlab     = "",
		offsetRow    = 0,
		offsetCol    = 0,
		margins      = c(12, 40),
		density.info = "none",
		symkey       = FALSE,
		symbreaks    = FALSE,
		cexCol       = 2,
		cexRow       = 1.5,
		notecex      = 1
	)
	dev.off()
}

# Note: the order of barcode identity is important for matching counts
barcodeIdentify <- c(
			"Unmodified",
			"H3K4me1",
			"H3K4me2",
			"H3K4me3",
			"H3K9me1",
			"H3K9me2",
			"H3K9me3",
			"H3K27me1",
			"H3K27me2",
			"H3K27me3",
			"H3K36me1",
			"H3K36me2",
			"H3K36me3",
			"H4K20me1",
			"H4K20me2",
			"H4K20me3"
			)

spikeIn <- read.csv("spikeIn/spikeIn.counts.txt", header = TRUE,  sep = "\t")
row_odd <- seq_len(nrow(spikeIn)) %% 2
tmp1    <- spikeIn[which(row_odd == 1), ]
tmp2    <- spikeIn[which(row_odd == 0), ]
spikeIn <- cbind(tmp1[, c(1, 2)], tmp1[, 3], tmp2[, 3], tmp1[, c(4, 5)], tmp2[, c(4, 5)])

colnames(spikeIn) <- c(
	"sample",
	"raw_counts",
	"barcode_A",
	"barcode_B",
	"R1.count_A",
	"R2.count_A",
	"R1.count_B",
	"R2.count_B"
)

barcode <- rep(barcodeIdentify, length(unique(spikeIn$sample)))

spikeIn <- cbind(spikeIn[, 1:2], barcode, spikeIn[, c(3:8)]) 
spikeIn <- transform(spikeIn, barcode_counts = R1.count_A + R1.count_B + R2.count_A + R2.count_B)
spikeIn <- transform(spikeIn, on_target_percent = 0)
spikeIn <- transform(spikeIn, sum_percent       = 0)

specificity           <- data.frame(matrix(NA, length(unique(spikeIn$sample)), 17))
colnames(specificity) <- c("library", barcodeIdentify)

for (i in 1:length(unique(spikeIn$sample))) {
	sample.tmp <- unique(spikeIn$sample)[i]
	x.counts   <- which(spikeIn$sample == sample.tmp)
	sum.tmp    <- sum(spikeIn$barcode_counts[x.counts])
	
	specificity$library[i] <- unique(spikeIn$sample)[i]
	
	for (j in 1:16) {
		k <- x.counts[j]
		spikeIn$sum_percent[k] <- spikeIn$barcode_counts[k] / sum.tmp
		specificity[i, j + 1]  <- spikeIn$sum_percent[k] * 100
	}
}

write.table(
	spikeIn,
	"spikeIn/spikeInSummaryCounts.txt",
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
	)

write.table(
	specificity,
	"spikeIn/specificity.txt",
	col.names = TRUE,
	row.names = FALSE,
	quote     = FALSE,
	sep       = "\t"
	)

m           <- specificity[, 2:17]
rownames(m) <- specificity[, 1]
m           <- as.matrix(m)
m           <- round(m, digits = 3)

plot_heatmap(m, "spikeIn/specificity")


