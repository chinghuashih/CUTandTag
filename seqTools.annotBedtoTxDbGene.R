require("BSgenome");
require("data.table")
require("GenomicFeatures")
require("rtracklayer")

annotBedtoTxDbGene = function(bed, tx, org = c(NULL), prefix = c(NULL), promUp = 2500, promDown = 0, org.keytype = c("ENTREZID"), ...) {

	#Debug
	#bed = tabr;        tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0; org.keytype = "ENTREZID"
	#bed = peaksBed;    tx = tx; org = org;  prefix = "GRCh37"; promUP = 2500; promDown = 0
	#bed = trLocs;      tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0; org.keytype = "ENTREZID"
	#bed = rnar[95254]; tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0

	#Columns to add to bed Granges
	annotCols = paste(
					prefix,
					c(
						"gene",
						"promoter",
						"exon",
						"intron",
						"5utr",
						"3utr",
						"ts",
						"tsKg",
						"tsDist",
						"tsStart",
						"tsEnd",
						"tsStrand",
						"tsRelPos",
						"tss",
						"tssKg",
						"tssDist",
						"tssRelPos",
						"tssStart",
						"tssEnd",
						"tssStrand"
						),
					sep = "."
					)
	bed@elementMetadata[annotCols] = rep(NA, dim(bed@elementMetadata)[1])

	#get Ts database
	ts = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")) #, ...) # get transcripts

	#Determine the the closest transcript
	closeTs = distanceToNearest(bed, ts, ignore.strand = TRUE)
	
	#Assign ts name to annotation column (AnnotCol) "tx"
	bed@elementMetadata[[paste0(prefix, ".ts")]][queryHits(closeTs)]  = ts@elementMetadata$TXNAME[subjectHits(closeTs)]
 	bed@elementMetadata[paste0(prefix, ".tsKg")]                      = select(
 																				tx,
 																				keytype = "TXNAME",
 																				keys    = as.character(bed@elementMetadata[[paste0(prefix, ".ts")]]),
 																				columns = "GENEID"
 																				)$GENEID

	bed@elementMetadata[[paste0(prefix, ".tsStart")]][ queryHits(closeTs)] = start(ts)[subjectHits(closeTs)]
	bed@elementMetadata[[paste0(prefix, ".tsEnd")]][   queryHits(closeTs)] = end(ts)[subjectHits(closeTs)]
	bed@elementMetadata[[paste0(prefix, ".tsStrand")]][queryHits(closeTs)] = as.character(strand(ts)[subjectHits(closeTs)])

	#Determine the orientation
	sign                                                        =  sign(start(bed[queryHits(closeTs)]) - start(ts[subjectHits(closeTs)]))
	sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"] = -sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"]            #adjust for strand
	
	#Assign ts dist	
	bed@elementMetadata[[paste(prefix, "tsDist", sep = ".")]][queryHits(closeTs)] = closeTs@elementMetadata$distance * sign

	#Determine ts rel pos
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)]        = closeTs@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx = bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)] == 0
	
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx] = 
			(start(bed[queryHits(closeTs)][inTx]) - start(ts[subjectHits(closeTs)])[inTx]) / (end(ts[subjectHits(closeTs)])[inTx] - start(ts[subjectHits(closeTs)])[inTx])
	
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx & as.character(strand(ts[subjectHits(closeTs)])) == "-"] =
			(1 - bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx & as.character(strand(ts[subjectHits(closeTs)])) == "-"])

	tss = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")) #, ...)
	#make strand-specific ts range object
	end(tss)[as.character(strand(tss))   == "+"] = start(tss)[as.character(strand(tss)) == "+"]
	start(tss)[as.character(strand(tss)) == "-"] = end(tss)[as.character(strand(tss))   == "-"]

	#Determine the relative location
	closeTss = distanceToNearest(bed, tss, ignore.strand = TRUE)
	#closeTss = closeTss[!is.na(subjectHits(closeTss)), ]
	
	#Assign ts name to annotation column (AnnotCol) "ts"
	bed@elementMetadata[[paste0(prefix, ".tss")]][queryHits(closeTs)] = tss@elementMetadata$TXNAME[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssKg")]]                   = select(tx, keytype = "TXNAME", keys = bed@elementMetadata[[paste0(prefix, ".tss")]], columns = "GENEID")$GENEID

	bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)]  = start(ts)[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)]    = end(ts)[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssStrand")]][queryHits(closeTss)] = as.character(strand(ts)[subjectHits(closeTss)])

	#Determine the orientation
	sign                                                          =  sign(start(bed[queryHits(closeTss)]) - start(tss[subjectHits(closeTss)]))
	sign[as.character(strand(tss[subjectHits(closeTss)])) == "-"] = -sign[as.character(strand(tss[subjectHits(closeTss)])) == "-"]              #adjust for strand

	#Assign tss dist
	bed@elementMetadata[[paste0(prefix, ".tssDist")]][queryHits(closeTss)]   = closeTss@elementMetadata$distance * sign
	#Determine the position relative to the closest tss
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] = closeTss@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx    = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] >= 0 & bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] < (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)])
	afterTx = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] >= (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)])
	
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][inTx]    = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][inTx]    / (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)][inTx]    - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)][inTx])
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][afterTx] = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][afterTx] -  bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)][afterTx] + bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)][afterTx]

	#Assign gene overlap
	genes                                                                 = suppressWarnings(genes(tx, columns = c("GENEID")))
	overlaps                                                              = findOverlaps(bed, genes)
	overlaps                                                              = cbind(as.data.frame(overlaps), id = as.character(genes@elementMetadata$GENEID[subjectHits(overlaps)]))
	overlaps                                                              = data.table(unique(overlaps[, c(1, 3)]))
	overlaps                                                              = overlaps[!is.na(overlaps$id), ]
	ag                                                                    = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "gene", sep = ".")]][ag$queryHits] = ag$id

	#Assign promoter overlap
	proms                                                                     = suppressWarnings(promoters(tx, upstream = promUp, downstream = promDown, columns = c("GENEID")))
	overlaps                                                                  = findOverlaps(bed, proms)
	overlaps                                                                  = cbind(as.data.frame(overlaps), id = as.character(proms@elementMetadata$GENEID[subjectHits(overlaps)]))
	overlaps                                                                  = data.table(unique(overlaps[, c(1, 3)]))
	overlaps                                                                  = overlaps[!is.na(overlaps$id), ]
	ag                                                                        = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "promoter", sep = ".")]][ag$queryHits] = ag$id

	#Assign exon overlap
	exons                                                                 = exonsBy(tx, by = "tx")
	if (exists("vals") && all(!is.na(vals$tx_name))) {
		exons = exons[names(exons) %in% vals$tx_name]
	} else if (exists("vals")) {
		stop("Need vals$tx_name to subset Exon")
	} else {
	}
	overlaps                                                              = findOverlaps(bed, exons)
	overlaps                                                              = cbind(as.data.frame(overlaps), id = names(exons)[subjectHits(overlaps)])
	overlaps$id                                                           = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXID)]) #switch transcript id with gene id
	overlaps                                                              = data.table(unique(overlaps[, c(1, 3)]))                                               #Remove duplicates
	ag                                                                    = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "exon", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign intron overlap
	introns                                                                 = intronsByTranscript(tx, use.name = TRUE)
	if (exists("vals") && all(!is.na(vals$tx_name))) {
		introns = introns[names(introns) %in% vals$tx_name]
	} else if {(exists("vals")) {
		stop("Need vals$tx_name to subset Intron")
	} else {
	}
	overlaps                                                                = findOverlaps(bed, introns)
	overlaps                                                                = cbind(as.data.frame(overlaps), id = names(introns)[subjectHits(overlaps)])
	overlaps$id                                                             = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]) #switch transcript id with gene id
	overlaps                                                                = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)]))                             #Remove duplicates
	ag                                                                      = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "intron", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 5'utr classification
	fiveUtr                                                               = fiveUTRsByTranscript(tx, use.name = TRUE)
	if (exists("vals") && all(!is.na(vals$tx_name))) {
		fiveUtr = fiveUtr[names(fiveUtr) %in% vals$tx_name]
	} else if (exists("vals")) {
		stop("Need vals$tx_name to subset 5' UTR")
	} else {
	}
	overlaps                                                              = findOverlaps(bed, fiveUtr)
	overlaps                                                              = cbind(as.data.frame(overlaps), id = names(fiveUtr)[subjectHits(overlaps)])
	overlaps$id                                                           = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]) #switch transcript id with gene id
	overlaps                                                              = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)]))                             #Remove duplicates
	ag                                                                    = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "5utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 3'utr classification
	threeUtr                                                              = threeUTRsByTranscript(tx, use.name = TRUE)
	if (exists("vals") && all(!is.na(vals$tx_name))) {
		threeUtr = threeUtr[names(threeUtr) %in% vals$tx_name]
	} else if (exists("vals")) {
		stop("Need vals$tx_name to subset 3' UTR")
	} else {
	}
	overlaps                                                              = findOverlaps(bed, threeUtr)
	overlaps                                                              = cbind(as.data.frame(overlaps), id = names(threeUtr)[subjectHits(overlaps)])
	overlaps$id                                                           = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]) #switch transcript id with gene id
	overlaps                                                              = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)]))                             #Remove duplicates
	ag                                                                    = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "3utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	if (!is.null(org)) {
		bed@elementMetadata[paste0(prefix, ".tsSymbol")]    = rep(NA, length(bed))
		sym                                                 = unique(select(org, keytype = org.keytype, keys = as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]]), columns = c(org.keytype, "SYMBOL")))
		sym                                                 = aggregate(sym$SYMBOL, by = list(sym[[org.keytype]]), FUN = paste, collapse = "|")
		bed@elementMetadata[[paste0(prefix, ".tsSymbol")]]  = sym$x[match(bed@elementMetadata[[paste0(prefix, ".tsKg")]], sym$Group.1)]

		bed@elementMetadata[paste0(prefix, ".tssSymbol")]   = rep(NA, length(bed))
		sym                                                 = unique(select(org, keytype = org.keytype, keys = as.character(bed@elementMetadata[[paste0(prefix, ".tssKg")]]), columns = c(org.keytype, "SYMBOL")))
		sym                                                 = aggregate(sym$SYMBOL, by = list(sym[[org.keytype]]), FUN = paste, collapse = "|")
		bed@elementMetadata[[paste0(prefix, ".tssSymbol")]] = sym$x[match(bed@elementMetadata[[paste0(prefix, ".tssKg")]], sym$Group.1)]
	}		
	return(bed)
}

annotateGtf = function (gtfFile, outFile = c(NA)) {

	#Debug
	#gtfFile = "/media/bbarwick/RAID/tools/Genome/GRCh37/Homo_sapiens.GRCh37.74.gtf"

	gtf = rtracklayer::import(gtfFile)
	gtf = as.data.table(as.data.frame(gtf))

	#summarize by id
	genes = gtf[
		,
		.(
			chr          = unique(seqnames),
			start        = min(start),
			end          = max(end),
			strand       = unique(strand),
			symbol       = unique(gene_name),
			gene_biotype = unique(gene_biotype),
			source       = paste0(unique(source), collapse = "|")
			),
		by = gene_id
		]

	#summarize protein ids that are not NA
	pids = gtf[
		!is.na(protein_id),
		.(protein_id = paste0(unique(protein_id), collapse = "|")),
		by           = gene_id
		]

	#summarize coding start and stop
	cds = gtf[
		type == "CDS",
		.(
			cds.start = min(start),
			cds.end   = max(end)
			),
		by = gene_id
		]

	#merge
	genes = merge(genes, cds,  by = "gene_id", all.x = TRUE)
	genes = merge(genes, pids, by = "gene_id", all.x = TRUE)

	if (!is.na(outFile)) fwrite(genes, outFile, sep = "\t", row.names = FALSE)

	invisible(genes)
}
