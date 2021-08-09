#!/bin/bash
#SBATCH -p standard -o CnT.mapping.bowtie2.log -t 16:00:00
#SBATCH -c 12 --mem=144G

#################################################################################################################
# input files:
#                 1. fastq/fastq.gz
#                 2. the index reference genomes, including hg38, E. coli, and phiX
#                 3. samples.txt: contains samples' names, one per each line
#
# required programs:
#                 1. fastqc
#                 2. fastp
#                 3. bowtie2
#                 4. samtools
#                 5. picard
#                 6. deeptools
#
# note:
#                 1. fastp and bamqc are needed to be installed.
#                 2. check the names of fastq/fastq.gz files; the file names in this scripts is the output 
#                    from URMC genome sequencing core.
#                 3. deeptools commands in this pipeline are based on v2.5.3. Normalized bigWig commands need
#                    to be modified if you use higher version.
#################################################################################################################

#########################################
## loading modules: required for slurm ##
#########################################

# mapping and bam file processing
module load bowtie2
module load samtools
module load picard

# QC and visualization
module load fastp
module load fastqc
module load deeptools

# misc
module load java
module load perl

#######################
## parameter setting ##
#######################
# library type
libraryType="PE"
# SE: single-end libraries, and PE: paired-end libraries; cannot be used for mixed types of libraries

# environment parameters
numberOfProcessors=12

# program parameters
bwBinSize=10
cutWindowSize=4
lengthRequired=35
mappingQuality=10
maxFragmentLength=800
maxin=700

# references
reference="/scratch/cshih8/references/hg38/bowtie2_index/default/hg38"
reference_ecoli="/home/cshih8/references/ecoli_bowtie2/ecoli"
reference_phix="/home/cshih8/references/phix_bowtie2/phix"

# import the sample list
awk '{print $0=$0".bam"}'         < samples.txt > bamfiles.txt
awk '{print $0=$0".dupMark.bam"}' < samples.txt > dupMarkBamfiles.txt

sampleFiles=$(<samples.txt)
dupMarkBamfiles=$(<dupMarkBamfiles.txt)
dedupBamfiles=$(<dedupBamfiles.txt)

echo "Samples are:"
echo "${sampleFiles}"
echo "libraryType: ${libraryType}"
echo "bwBinSize = ${bwBinSize}"
echo "maximun fragment length: ${maxin}"
echo "mapping quality cut-off: ${mappingQuality}"
echo ""
echo ""

# mkdir
mkdir -p QC/bamqc
mkdir -p QC/bowtie2_summary
mkdir -p QC/fastqc
mkdir -p QC/fastp
mkdir -p QC/fragmentSize
mkdir -p QC/misc
mkdir -p QC/multiBamSummary

mkdir -p bamFiles
mkdir -p bamFiles.ecoli
mkdir -p bamFiles.phix

mkdir -p bigWig/dupMark/rpkm
mkdir -p bigWig/dupMark/none
mkdir -p bigWig/dedup/rpkm
mkdir -p bigWig/dedup/none

########################
## QC - raw sequences ##
########################

for FILE in *.gz
do
	echo "fastqc: ${FILE}"
	fastqc "${FILE}" \
		--threads ${numberOfProcessors} \
		--outdir QC/fastqc \
		--quiet
	echo ""
	#--threads $(echo ${numberOfProcessors}/2 | bc) \
done

for sample in ${sampleFiles[*]}
do
	echo "QC and trimming ${sample} by fastp"
	if [ "${libraryType}" == "SE" ]; then
		fastp \
			-i ${sample}.R1.fastq.gz \
			-o ${sample}.R1_trim.fastq.gz \
			--trim_poly_g \
			-x \
			--cut_window_size ${cutWindowSize} \
			--cut_tail \
			--length_required ${lengthRequired}

	elif [ "$libraryType" == "PE" ]; then
		~/tools/fastp/fastp \
			--in1 ${sample}.R1.fastq.gz \
			--in2 ${sample}.R2.fastq.gz \
			--out1 ${sample}.R1_trim.fastq.gz \
			--out2 ${sample}.R2_trim.fastq.gz \
			--trim_poly_g \
			-x \
			--cut_window_size ${cutWindowSize} \
			--cut_tail \
			--length_required ${lengthRequired}
	else
		echo "wrong library type"
	fi

	mv fastp.html QC/fastp/${sample}.fastp.html
	mv fastp.json QC/fastp/${sample}.fastp.json
	echo ""
done

#######################
## mapping - bowtie2 ##
#######################

for sample in ${sampleFiles[*]}
do
	echo "mapping ${sample} with bowtie2"

	if [ "$libraryType" == "SE" ]; then
		echo "mapping single-end library of ${sample} to human genome"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference} \
			--no-unal \
			--local \
			--very-sensitive \
			-U ${sample}.R1_trim.fastq.gz \
			-S ${sample}.fastp_trim.sam &> QC/bowtie2_summary/${sample}.bowtie2.txt
		echo ""

		echo "mapping single-end library of ${sample} to E. coli"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_ecoli} \
			--no-unal \
			--local \
			--very-sensitive-local \
			-U ${sample}.trim.fastq.gz \
			-S ${sample}.ecoli.sam &> QC/bowtie2_summary/${sample}_ecoli.bowtie2.txt
		echo ""

		echo "mapping single-end library of ${sample} to phiX"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_phix} \
			--no-unal \
			--local \
			--very-sensitive-local \
			-U ${sample}.trim.fastq.gz \
			-S ${sample}.phix.sam  &> QC/bowtie2_summary/${sample}_phix.bowtie2.txt
		echo ""

	elif [ "${libraryType}" == "PE" ]; then
		echo "mapping paired-end library of ${sample} to human genome"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference} \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin} \
			-1 ${sample}.R1_trim.fastq.gz \
			-2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.fastp_trim.sam &> QC/bowtie2_summary/${sample}.bowtie2.txt
		echo ""

		echo "mapping paired-end library of ${sample} to E. coli"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_ecoli} \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin} \
			-1 ${sample}.R1_trim.fastq.gz \
			-2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.ecoli.sam  &> QC/bowtie2_summary/${sample}_ecoli.bowtie2.txt
		echo ""

		echo "mapping paired-end library of ${sample} to phiX"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_phix} \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin} \
			-1 ${sample}.R1_trim.fastq.gz \
			-2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.phix.sam &> QC/bowtie2_summary/${sample}_phix.bowtie2.txt
		echo ""
	else
		echo "wrong library type"
	fi

	echo "filtering ${sample} with the mapping quality < ${mappingQuality}, remove unassembled contigs, mitochrondria, and chrY"
  
	sed '/random/d;/chrUn/d;/chrEBV/d;/chrM/d;/chrY/d' < ${sample}.fastp_trim.sam > ${sample}.filtered.sam
#	samtools view -bSq ${mappingQuality} ${sample}.fastp_trim.sam chrM > ${sample}.chrM.bam
  	samtools view -bSq ${mappingQuality} ${sample}.filtered.sam        > ${sample}.filtered.bam
	samtools view -bSq ${mappingQuality} ${sample}.ecoli.sam           > ${sample}.ecoli10.bam
	samtools view -bSq ${mappingQuality} ${sample}.phix.sam            > ${sample}.phix10.bam
	echo ""

	echo "samtools flagstat ${sample}.filtered.bam"
	samtools flagstat ${sample}.filtered.bam
	echo ""

	rm -f ${sample}.*_trim.fastq.gz
	rm -f ${sample}.*.sam
	mv ${sample}.ecoli10.bam bamFiles.ecoli
	mv ${sample}.phix10.bam  bamFiles.phix
#	mv ${sample}.chrM.bam    bamFiles.chrM
done

###################################################################
## bam files processing: sorting and mark duplicates with picard ##
###################################################################

for sample in ${sampleFiles[*]}
do
	echo "sorting ${sample} by picard"
	java -jar ${PICARD} SortSam \
		I=${sample}.filtered.bam \
		O=${sample}.sorted.bam \
		TMP_DIR=tmp \
		SORT_ORDER=coordinate
	rm -f ${sample}.filtered.bam
	echo ""

	echo "marking duplicates in ${sample} and indexing"
	java -jar ${PICARD} MarkDuplicates \
		I=${sample}.sorted.bam \
		O=${sample}.dupMark.bam \
		TMP_DIR=tmp \
		METRICS_FILE=${sample}.dupMark.metrics
	samtools index ${sample}.dupMark.bam
	rm -f ${sample}.sorted.bam
	echo ""

	echo "check bam files of ${sample} by bamqc"
	~/tools/BamQC/bin/bamqc ${sample}.dupMark.bam
	echo ""

	mv *.metrics QC/misc/
	mv *_bamqc   QC/bamqc/
	echo ""
done

######################################################################################
## Quality controls of aligned reads (BAM files)                                    ##
## multiBamSummary: similarity of read distribution in different sequencing samples ##
######################################################################################

echo "multiBamSummary"
multiBamSummary bins \
	--numberOfProcessors ${numberOfProcessors} \
	--bamfiles ${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--minMappingQuality ${mappingQuality} \
	-o multiBamSummary.results.npz
echo ""

# it might be better to plot PCA in R 
# --transpose can be used after version 3
#echo "plot PCA"
#plotPCA \
#	--corData multiBamSummary.results.npz \
#       --labels ${sampleFiles} \
#	--transpose \
#	-o PCA_readCounts.png
#echo ""

echo "plot heatmap for pairwise Spearman correlation coefficiencies"
plotCorrelation \
	-in multiBamSummary.results.npz \
	--corMethod spearman \
	--skipZeros \
	--plotTitle "Spearman Correlation of Read Counts" \
	--whatToPlot heatmap \
	--colorMap RdYlBu \
	--plotNumbers \
	-o heatmap_SpearmanCorr_readCounts.png \
	--outFileCorMatrix SpearmanCorr_readCounts.tab
echo ""

mv multiBamSummary.results.npz         QC/multiBamSummary/
mv SpearmanCorr_readCounts.tab         QC/multiBamSummary/
mv heatmap_SpearmanCorr_readCounts.png QC/multiBamSummary/

##############################################
## plotCoverage: coverage of the sequencing ##
##############################################
echo "plotCoverage, w/ duplicates"
plotCoverage \
	--bamfiles ${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--skipZeros \
	--numberOfSamples 1000000 \
	--numberOfProcessors ${numberOfProcessors} \
	--plotFile coverage.dupMark.png \
	--outRawCounts rawCounts.coverage.dupMark.txt
echo ""

echo "plotCoverage, w/o duplicates"
plotCoverage \
	--bamfiles ${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--skipZeros \
	--numberOfSamples 1000000 \
	--numberOfProcessors ${numberOfProcessors} \
	--plotFile coverage.dedup.png \
	--outRawCounts rawCounts.coverage.dedup.txt \
	--ignoreDuplicates
echo ""

mv coverage.*.png  QC/misc/
mv rawCounts.*.txt QC/misc/

# How well did your ChIP experiment work?
#	- deepTools: plotFingerprint
#	Processing/normalizing of alignment files
#	to remove a GC-bias from aligned reads (BAM file)
#	- deepTools: correctGCbias

##########################################################################
## FragmentSize :the fragment size distribution of paired-end data      ##
##########################################################################

if [ "${libraryType}" == "PE" ]; then
	for sample in ${sampleFiles[*]}
	do
		bamPEFragmentSize \
			-hist fragmentSize.${sample}.png \
			--plotTitle "Fragment Size" \
			--numberOfProcessors ${numberOfProcessors} \
			--maxFragmentLength ${maxFragmentLength} \
			--bamfiles ${sample}.dupMark.bam \
			--samplesLabel ${sample}
		echo ""
	done
	mv fragmentSize.*.png QC/fragmentSize
fi

# plotFingerprint
plotFingerprint \
	--bamfiles ${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--minMappingQuality ${mappingQuality} \
	--skipZeros \
	--numberOfSamples 50000 \
	--numberOfProcessors ${numberOfProcessors} \
	--plotTitle "Fingerprints of samples" \
	--plotFile fingerprints.png \
	--outRawCounts fingerprints.tab \
	--outQualityMetrics qualityMetrics.tab

mv fingerprints.*     QC/misc/
mv qualityMetrics.tab QC/misc/

# Should you be worried about a GC bias of your sample?
# - deepTools: computeGCbias, correctGCbias
# - if the data is normalized for GC bias, --ignoreDuplicates should absolutely not be used.

###############################################################################################################
## coverage: to generate a sequencing-depth-normalized continuous profile of read coverages (BAM --> bigWig) ##
###############################################################################################################

echo "making bigWig files"
echo "binSize: ${bwBinSize} and ignore chrX, chrY, and chrM for normalization"
echo ""

for sample in ${sampleFiles[*]}
do
	echo "generating normalized (RPKM) bigwig files of ${sample} w/ duplicates"
	bamCoverage \
		--bam ${sample}.dupMark.bam \
		--outFileName ${sample}.dupMark.rpkm.bw \
		--binSize ${bwBinSize} \
		--extendReads \
		--numberOfProcessors ${numberOfProcessors} \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig files of ${sample} w/ duplicates"
	bamCoverage \
		--bam ${sample}.dupMark.bam \
		--outFileName ${sample}.dupMark.bw \
		--binSize ${bwBinSize} \
		--extendReads \
		--numberOfProcessors ${numberOfProcessors}
	echo ""

	echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
	bamCoverage \
		--bam ${sample}.dupMark.bam \
		--outFileName ${sample}.dedup.rpkm.bw \
		--binSize ${bwBinSize} \
		--ignoreDuplicates \
		--extendReads \
		--numberOfProcessors ${numberOfProcessors} \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
	bamCoverage \
		--bam ${sample}.dupMark.bam \
		--outFileName ${sample}.dedup.bw \
		--binSize ${bwBinSize} \
		--ignoreDuplicates \
		--extendReads \
		--numberOfProcessors ${numberOfProcessors}
	echo ""

	mv ${sample}.dupMark.rpkm.bw bigWig/dupMark/rpkm
        mv ${sample}.dupMark.bw      bigWig/dupMark/none
        mv ${sample}.dedup.rpkm.bw   bigWig/dedup/rpkm
        mv ${sample}.dedup.bw        bigWig/dedup/none
done
mkdir -p raw
mv *gz        raw
mv *files.txt QC/misc/
mv samples.tx QC/misc/
mv *.bam*     bamFiles/

rm -rf tmp

echo "end of mapping"
