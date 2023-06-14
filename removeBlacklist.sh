!/bin/bash
#SBATCH -p standard -o CnT.mapping.bowtie2.log -t 2:00:00
#SBATCH -c 16 --mem=128G

module load bedtools
module load samtools

blacklist="/scratch/cshih8/references/hg38-blacklist.v2.bed"
sampleFiles=$(<samples.txt)

# library type
libraryType="PE"
# SE: single-end libraries, and PE: paired-end libraries; cannot be used for mixed types of libraries

extendReads=200
# for SE libraries to make bigWig files. No more than 4X read length.

# program parameters
bwBinSize=10
numberOfProcessors=8

dupType="dupMark"

mkdir -p dedup.rmBlacklist
mkdir -p bigWig/${dupType}.rmBlacklist/rpkm
mkdir -p bigWig/${dupType}.rmBlacklist/none

for sample in ${sampleFiles[*]}
do
	bedtools intersect \
	-abam bamFiles/${dupType}/${sample}.${dupType}.bam \
	-b ${blacklist} \
	-v \
	> bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam

	samtools index bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam
done

for sample in ${sampleFiles[*]}
do
	if [ "${libraryType}" == "SE" ]; then

		echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam \
			--outFileName ${sample}.${dupType}.rpkm.bw \
			--binSize ${bwBinSize} \
			--extendReads ${extendReads} \
			--numberOfProcessors ${numberOfProcessors} \
			--normalizeUsing RPKM \
			--ignoreForNormalization chrX
		echo ""

		echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam \
			--outFileName ${sample}.${dupType}.bw \
			--binSize ${bwBinSize} \
			--extendReads ${extendReads} \
			--numberOfProcessors ${numberOfProcessors}
   		echo ""

		mv ${sample}.${dupType}.rpkm.bw bigWig/${dupType}.rmBlacklist/rpkm
  		mv ${sample}.${dupType}.bw      bigWig/${dupType}.rmBlacklist/none

	elif [ "${libraryType}" == "PE" ]; then
		echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam \
			--outFileName ${sample}.${dupType}.rpkm.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors} \
			--normalizeUsing RPKM \
			--ignoreForNormalization chrX
		echo ""

		echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/${dupType}.rmBlacklist/${sample}.${dupType}.bam \
			--outFileName ${sample}.${dupType}.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors}
		echo ""

		mv ${sample}.${dupType}.rpkm.bw bigWig/${dupType}.rmBlacklist/rpkm
		mv ${sample}.${dupType}.bw      bigWig/${dupType}.rmBlacklist/none
	else
		echo "wrong library type"
  	fi
done
