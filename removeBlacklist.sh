!/bin/bash
#SBATCH -p standard -o CnT.mapping.bowtie2.log -t 2:00:00
#SBATCH -c 16 --mem=128G

module load bedtools
module load samtools

blacklist="/scratch/cshih8/references/hg38-blacklist.v2.bed"
sampleFiles=$(<samples.txt)

mkdir -p dedup.rmBlacklist

for sample in ${sampleFiles[*]}
do
	bedtools intersect \
	-abam dedup/${sample}.dedup.bam \
	-b ${blacklist} \
	-v \
	> dedup.rmBlacklist/${sample}.dedup.bam

	samtools index dedup.rmBlacklist/${sample}.dedup.bam
done
