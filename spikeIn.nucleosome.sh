#!/bin/bash
#SBATCH -p standard -o CnT.nucleosomeSpikeIn.counts.log -t 4:00:00
#SBATCH -c 2 --mem=32G

# Notes:
# EpiCypher considers an antibody with <20% binding to all off-target PTMs specific and suitable for downstream data analysis. 
# For IgG, data is normalized to the sum of total barcode reads.

module load R/4.1.1

sampleFiles=$(<samples.txt)

mkdir -p spikeIn

echo -e \
	sample"\t"\
	raw_counts"\t"\
	barcode"\t"\
	R1.count"\t"\
	R2.count > spikeIn/spikeIn.counts.txt

for sample in ${sampleFiles[*]}
do
	counts=$(zcat ${sample}.R1.fastq.gz | wc -l)
	counts=$(echo ${counts}"/4" | bc -l)

	for barcode in \
		TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG \
		CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG \
		CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA \
		TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG \
		TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA \
		CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA \
		ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG \
		CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG \
		CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA \
		ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA \
		CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC \
		GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG \
		CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG \
		GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA \
		CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG \
		CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
	do

		echo "${sample}	${barcode}"

		R1count=$(zcat ${sample}.R1.fastq.gz | grep -c ${barcode})
		R2count=$(zcat ${sample}.R2.fastq.gz | grep -c ${barcode})

		echo -e \
			${sample}"\t"\
			${counts}"\t"\
			${barcode}"\t"\
			${R1count}"\t"\
			${R2count} >> spikeIn/spikeIn.counts.txt
	done
done

Rscript --vanilla spikeIn.nucleosome.R

echo "end of spikeIn counting"
