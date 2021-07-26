#!/bin/bash
#SBATCH -p standard -o chipseq.mapping.bowtie2.log -t 16:00:00
#SBATCH -c 12 --mem=144G

#################################################################################################################
# required files:
#                 1. fastq/fastq.gz
#                 2. the reference genome
#                 3. annotation files: gtf and bed formats
#                 4. samples.txt: contains samples' names, one per each line
# Notice:
#                 1. fastp is needed to be installed.
#                 2. check the names of fastq/fastq.gz files; the file names in this scripts is the output 
#                    from URMC genome sequencing core.
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
# library
libraryType="PE"
# SE: single-end libraries, and PE: paired-end libraries; cannot be used for mixed types of libraries

# references
reference="/scratch/cshih8/references/hg38/bowtie2_index/default/hg38"
reference_ecoli="/home/cshih8/references/ecoli_bowtie2/ecoli"
reference_phix="/home/cshih8/references/phix_bowtie2/phix"

# programs
bwBinSize=10
mappingQuality=10
maxin=700
numberOfProcessors=12

# import the sample list
awk '{print $0=$0".bam"}'         < samples.txt > bamfiles.txt
awk '{print $0=$0".dupMark.bam"}' < samples.txt > dupMarkBamfiles.txt
awk '{print $0=$0".dedup.bam"}'   < samples.txt > dedupBamfiles.txt

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

#########################
## QC1 - raw sequences ##
#########################

mkdir -p QC/fastqc
mkdir -p QC/fastp
mkdir -p QC/bamqc
mkdir -p QC/misc

for sample in ${sampleFiles[*]}
do
	echo "trimming $sample by fastp"
	if [ "${libraryType}" == "SE" ]; then
		fastp \
			-i ${sample}.R1.fastq.gz \
			-o ${sample}.R1_trim.fastq.gz \
			--trim_poly_g \
			-x \
			--cut_window_size 4 \
			--cut_tail \
			--length_required 35

	elif [ "$libraryType" == "PE" ]; then
		~/tools/fastp/fastp \
			--in1 ${sample}.R1.fastq.gz \
			--in2 ${sample}.R2.fastq.gz \
			--out1 ${sample}.R1_trim.fastq.gz \
			--out2 ${sample}.R2_trim.fastq.gz \
			--trim_poly_g \
			-x \
			--cut_window_size 4 \
			--cut_tail \
			--length_required 35
	else
		echo "wrong library type"
	fi

	mv fastp.html ./QC/fastp/${sample}.fastp.html
	mv fastp.json ./QC/fastp/${sample}.fastp.json
	echo ""
done

#######################
## mapping - bowtie2 ##
#######################

# to do: option: merge replicates: cat

#mkdir -p alignment/bed
#mkdir -p alignment/bedgraph

mkdir -p bamFiles
mkdir -p bamFiles.ecoli
mkdir -p bamFiles.phix
mkdir -p bamFiles.chrM
mkdir -p QC/bowtie2_summary

for sample in ${sampleFiles[*]}
do
	echo "mapping ${sample}"

	if [ "$libraryType" == "SE" ]; then
		echo "mapping ${sample} to human genome"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference} \
			--no-unal \
			--local \
			--very-sensitive \
			-U ${sample}.R1_trim.fastq.gz \
			-S ${sample}.fastp_trim.sam &> ./QC/bowtie2_summary/${sample}_bowtie2.txt
			
		echo ""

		echo "mapping ${sample} to E. coli"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_ecoli} \
			--no-unal \
			--local \
			--very-sensitive-local \
			-U ${sample}.trim.fastq.gz \
			-S ${sample}.ecoli.sam &> ./QC/bowtie2_summary/${sample}_bowtie2.txt
		echo ""

		echo "mapping ${sample} to phix"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_phix} \
			-U ${sample}.trim.fastq.gz \
			--no-unal \
			--local \
			--very-sensitive-local \
			-S ${sample}.phix.sam  &> ./QC/bowtie2_summary/${sample}_bowtie2.txt
		echo ""

	elif [ "${libraryType}" == "PE" ]; then
		echo "mapping ${sample} to human genome"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference} \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin} \			
			-1 ${sample}.R1_trim.fastq.gz -2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.fastp_trim.sam &> ./QC/bowtie2_summary/${sample}_bowtie2.txt
		echo ""



=check
		echo "mapping $sample to E. coli"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_ecoli} \
			-1 ${sample}.R1_trim.fastq.gz -2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.ecoli.sam \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin}
		echo ""

		echo "mapping $sample to phiX"
		bowtie2 -t \
			-p ${numberOfProcessors} \
			-x ${reference_phix} \
			-1 ${sample}.R1_trim.fastq.gz -2 ${sample}.R2_trim.fastq.gz \
			-S ${sample}.phix.sam \
			--no-unal \
			--local \
			--very-sensitive-local \
			--no-mixed \
			--no-discordant \
			--maxins ${maxin}
		echo ""
	else
		echo "wrong library type"
	fi

	echo "filtering ${sample} with mapping quality < ${mappingQuality}, remove mitochrondria and unassembled contigs"
  
  	# Need to add an option for filtering: unconventional chromosomes, chrY, and chrM
	sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' < $sample.fastp_trim.sam > $sample.filtered.sam
	# sed '/chrM/d;/random/d;/chrUn/d' < $sample.fastp_trim.sam > $sample.filtered.sam
#######  
  	samtools view -bSq $mappingQuality $sample.filtered.sam > $sample.filtered.bam
	samtools view -bSq $mappingQuality $sample.ecoli.sam > $sample.ecoli10.bam
	samtools view -bSq $mappingQuality $sample.phix.sam > $sample.phix10.bam
	echo ""

	echo "samtools flagstat $sample.filtered.bam"
	samtools flagstat $sample.filtered.bam
	echo ""

	rm -f $sample.*.sam
done

##########################
## bam files processing ##
##########################

for sample in ${sampleFiles[*]}
do
	echo "sorting $sample by picard"
	java -Xmx10g -Xmx16G -jar /software/picard/2.12.0/picard.jar SortSam \
		I=$sample.filtered.bam \
		O=$sample.sorted.bam \
		TMP_DIR=tmp \
		SORT_ORDER=coordinate
	rm -f $sample.filtered.bam
	echo ""

	echo "marking duplicates in $sample and indexing"
	java -Xmx10g -Xmx16G -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
		I=$sample.sorted.bam \
		O=$sample.dupMark.bam \
		TMP_DIR=tmp \
		METRICS_FILE=$sample.dupMark.metrics
	samtools index $sample.dupMark.bam
	rm -f $sample.sorted.bam
	echo ""

	echo "check bam files of $sample by bamqc"
	~/tools/BamQC/bin/bamqc $sample.fastp_trim.bam
	~/tools/BamQC/bin/bamqc $sample.dupMark.bam
	echo ""
done

###################################################
## Quality controls of aligned reads (BAM files) ##
###################################################

######################################################################################
## multiBamSummary: similarity of read distribution in different sequencing samples ##
######################################################################################
echo "multiBamSummary"
multiBamSummary bins \
	--numberOfProcessors $numberOfProcessors \
	--bamfiles $dupMarkBamfiles \
	--labels $sampleFiles \
	--minMappingQuality $mappingQuality \
	-o multiBamSummary.results.npz
echo ""

# it might be better to plot PCA in R 
# --transpose can be used after version 3

#echo "plot PCA"
#plotPCA \
#	--corData multiBamSummary.results.npz \
#	--transpose \
#	-o PCA_readCounts.png
#echo ""

echo "plot heatmap"
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

##############################################
## plotCoverage: coverage of the sequencing ##
##############################################
echo "plotCoverage, w/ duplicates"
plotCoverage \
	--bamfiles $dupMarkBamfiles \
	--labels $sampleFiles \
	--skipZeros \
	--numberOfSamples 1000000 \
	--numberOfProcessors $numberOfProcessors \
	--plotFile coverage.dupMark.png \
	--outRawCounts rawCounts.coverage.dupMark.txt
echo ""

echo "plotCoverage, w/o duplicates"
plotCoverage \
	--bamfiles $dupMarkBamfiles \
	--labels $sampleFiles \
	--skipZeros \
	--numberOfSamples 1000000 \
	--numberOfProcessors $numberOfProcessors \
	--plotFile coverage.dedup.png \
	--outRawCounts rawCounts.coverage.dedup.txt \
	--ignoreDuplicates
echo ""

##########################################################################
## FragmentSize :the fragment size distribution of your paired-end data ##
##########################################################################

if [ "$libraryType" == "PE" ]; then
	bamPEFragmentSize \
		-hist fragmentSize.png \
		--plotTitle "Fragment Size" \
		--numberOfProcessors $numberOfProcessors \
		--maxFragmentLength 800 \
		--bamfiles $dupMarkBamfiles \
		--samplesLabel $sampleFiles

	for sample in ${sampleFiles[*]}
	do
		bamPEFragmentSize \
			-hist fragmentSize.$sample.png \
			--plotTitle "Fragment Size" \
			--numberOfProcessors $numberOfProcessors \
			--maxFragmentLength 800 \
			--bamfiles $sample.dupMark.bam \
			--samplesLabel $sample
	echo ""
	done
fi


# How well did your ChIP experiment work?
#	- deepTools: plotFingerprint
#	Processing/normalizing of alignment files
#	to remove a GC-bias from aligned reads (BAM file)
#	- deepTools: correctGCbias

# plotFingerprint
plotFingerprint \
	--bamfiles $dupMarkBamfiles \
	--labels $sampleFiles \
	--minMappingQuality $mappingQuality \
	--skipZeros \
	--numberOfSamples 50000 \
	--numberOfProcessors $numberOfProcessors \
	--plotTitle "Fingerprints of samples" \
	--plotFile fingerprints.png \
	--outRawCounts fingerprints.tab \
	--outQualityMetrics qualityMetrics.tab

# Should you be worried about a GC bias of your sample?
# - deepTools: computeGCbias, correctGCbias
# - if the data is normalized for GC bias, --ignoreDuplicates should absolutely not be used.

###############################################################################################################
## coverage: to generate a sequencing-depth-normalized continuous profile of read coverages (BAM --> bigWig) ##
###############################################################################################################
echo "binSize: $bwBinSize and ignore chrX, chrY, and chrM for normalization"
echo ""

if [ "$libraryType" == "SE" ]; then
	for sample in ${sampleFiles[*]}
	do
		echo "split $sample dupMark.bam file (SE) according to the strand"
		# forward strand
		bamCoverage -b $sample.dupMark.bam -o $sample.dupMark.fwd.bw --samFlagExclude 16
		# reverse strand
		bamCoverage -b $sample.dupMark.bam -o $sample.dupMark.rev.bw --samFlagInclude 16

		# index the filtered BAM file
		samtools index $sample.dupMark.fwd.bam
		samtools index $sample.dupMark.rev.bam
		
		echo ""
	done
elif [ "$libraryType" == "PE" ]; then
	for sample in ${sampleFiles[*]}
	do
		echo "split $sample dupMark.bam file (PE) according to the strand"
		# forward strand
		# include reads that are 2nd in a pair (128);
		# exclude reads that are mapped to the reverse strand (16)
		samtools view -b -f 128 -F 16 $sample.dupMark.bam > a.fwd1.bam

		# exclude reads that are mapped to the reverse strand (16) and
		# first in a pair (64): 64 + 16 = 80
		samtools view -b -f 80 $sample.dupMark.bam > a.fwd2.bam

		# reverse strand
		# include reads that map to the reverse strand (128)
		# and are second in a pair (16): 128 + 16 = 144
		samtools view -b -f 144 $sample.dupMark.bam > a.rev1.bam

		# include reads that are first in a pair (64), but
		# exclude those ones that map to the reverse strand (16)
		samtools view -b -f 64 -F 16 $sample.dupMark.bam > a.rev2.bam

		# merge the temporary files
		samtools merge -f $sample.dupMark.fwd.bam a.fwd1.bam a.fwd2.bam
		samtools merge -f $sample.dupMark.rev.bam a.rev1.bam a.rev2.bam

		# index the filtered BAM file
		samtools index $sample.dupMark.fwd.bam
		samtools index $sample.dupMark.rev.bam

		# remove temporary files
		rm -f a.fwd*.bam
		rm -f a.rev*.bam
		
		echo ""
	done
else
	echo "wrong library type"
fi


for sample in ${sampleFiles[*]}
do
	echo "generating normalized (RPKM) bigwig file of $sample w/ duplicate"
	bamCoverage \
		--bam $sample.dupMark.bam \
		--outFileName $sample.dupMark.rpkm.bw \
		--binSize $bwBinSize \
		--extendReads \
		--numberOfProcessors $numberOfProcessors \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig file of $sample w/ duplicate"
	bamCoverage \
		--bam $sample.dupMark.bam \
		--outFileName $sample.dupMark.bw \
		--binSize $bwBinSize \
		--extendReads \
		--numberOfProcessors $numberOfProcessors \
    --ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating normalized (RPKM) bigwig file of $sample w/o duplicate"
	bamCoverage \
		--bam $sample.dupMark.bam \
		--outFileName $sample.dedupMark.rpkm.bw \
		--binSize $bwBinSize \
		--ignoreDuplicates \
		--extendReads \
		--numberOfProcessors $numberOfProcessors \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig file of $sample w/o duplicate"
	bamCoverage \
		--bam $sample.dupMark.bam \
		--outFileName $sample.dedupMark.rpkm.bw \
		--binSize $bwBinSize \
		--ignoreDuplicates \
		--extendReads \
		--numberOfProcessors $numberOfProcessors \
    --ignoreForNormalization chrM chrX chrY
	echo ""

done

# other commands need to be added.
