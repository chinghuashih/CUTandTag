#!/bin/bash
#SBATCH -p standard  -o peaksCalling.log -t 08:00:00
#SBATCH -c 8 --mem=128G

module load bedtools
module load R

# peak calling methods: SEACR or MACS2
peakCalling="SEACR"

dupType="dedup"
# dedup or dupMark

# SEACR parameters
AUC=0.01

# macs2 parameters
inputFormat="auto"
# bam for BAMPE; bed for BEDPE; auto if not sure
genome="hs"
qvalue=0.05

echo "peak calling using ${peakCalling}"
echo ""

if [ "${peakCalling}" == "SEACR" ]; then

	# SEACR requires paired-end sequencing as inputs, both for experiments and the control(IgG) if having IgG as the control.

	mkdir -p peaks_SEACR/bedgraph
	mkdir -p peaks_SEACR/bed

	while read exp
	do
		echo ${exp}
		IFS=' ' read -r -a array <<< "${exp}"

		sample=${array[0]}
		control=${array[1]}

		echo "treatment: ${sample}"
		echo "control: ${control}"

		echo "Preparing bedgraph files: ${sample} and ${control}"
		# experiments
		bedtools genomecov -ibam ${sample}.${dupType}.bam -bg > ${sample}.bedgraph

		if [ "${control}" != "NA" ]; then
			# control
			bedtools genomecov -ibam ${control}.${dupType}.bam -bg > ${control}.bedgraph
		fi
        
        if [ "${control}" != "NA" ]; then
            ##############################################
            ## peakcalling with the stringent threshold ##
            ##############################################
            echo "peakCalls for ${sample} using normalized ${control} as the control track with stringent threshold"
            bash SEACR_1.3.sh ${sample}.bedgraph ${control}.bedgraph norm stringent ${sample}.stringent_norm
 
            echo "peakCalls for ${sample} using unnormalized ${control} as the control track with stringent threshold"
            bash SEACR_1.3.sh ${sample}.bedgraph ${control}.bedgraph non stringent ${sample}.stringent_non
            
            echo ""

            ############################################
            ## peakcalling with the relaxed threshold ##
            ############################################
            echo "peakCalls for ${sample} using normalized ${control} as the control track with relaxed threshold"
            bash SEACR_1.3.sh ${sample}.bedgraph ${control}.bedgraph norm relaxed ${sample}.relaxed_norm

            echo "peakCalls for ${sample} using unnormalized ${control} as the control track with relaxed threshold"
            bash SEACR_1.3.sh ${sample}.bedgraph ${control}.bedgraph non relaxed ${sample}.relaxed_non
                
            echo ""

        elif [ "${control}" == "NA" ]; then
        
            echo "peakCalls for ${sample} by selecting the top ${AUC} of regions by AUC"
            bash SEACR_1.3.sh ${sample}.bedgraph ${AUC} norm stringent ${sample}.AUC
            
            echo ""
            
        else
            echo "wrong library type"
        fi

    done <samplesPeaks.txt

    gzip *bedgraph
    mv *bedgraph.gz peaks_SEACR/bedgraph
    mv *bed         peaks_SEACR/bed

elif [ "${peakCalling}" == "MACS2" ]; then

    module load macs2
    
    outdir="peaks_macs2"
    
    mkdir -p ${outdir}/bedgraph
    mkdir -p ${outdir}/broadPeak
    mkdir -p ${outdir}/gappedPeak
    mkdir -p ${outdir}/narrowPeak
    mkdir -p ${outdir}/xls
    mkdir -p ${outdir}/summits
    
    echo "Samples for peak-calling are:"
    echo ${sampleFiles}
    echo "libraryType: ${inputFormat}"
    echo "output files are saved into ${outdir}"

    if [ "${inputFormat}" == "bam" ]; then
        format="BAMPE"
    elif [ "${inputFormat}" == "bed" ]; then
        format="BEDPE"
    elif [ "${inputFormat}" == "auto" ]; then
        format="AUTO"
    else
        echo "wrong input format"
    fi

    while read exp
    do
        echo ${exp}
        IFS=' ' read -r -a array <<< "${exp}"

        sample=${array[0]}
        control=${array[1]}

        echo "treatment: ${sample}"
        echo "control: ${control}"

        if [ "${control}" != "NA" ]; then

                # calling peaks
                echo "calling broad peaks: ${sample} with ${control}"
                macs2 callpeak \
                        --format ${format} \
                        --treatment ${sample}.${dupType}.bam \
                        --control ${control}.${dupType}.bam \
                        --name ${sample} \
                        --outdir ${outdir} \
                        --broad \
                        -g ${g} \
                        -q ${qvalue} \
                        -B \
                        --SPMR \
                        --nomodel
                echo ""               
                
                echo "calling narrow peaks: ${sample} with ${control}"
                macs2 callpeak \
                        --format ${format} \
                        --treatment ${sample}.${dupType}.bam \
                        --control ${control}.${dupType}.bam \
                        --name ${sample} \
                        --outdir ${outdir} \
                        -g ${g} \
                        -q ${qvalue} \
                        -B \
                        --SPMR \
                        --nomodel
                echo ""
                
        elif [ "${control}" == "NA" ]; then

            # calling peaks without control
            echo "calling broad peaks: ${sample} without control"
            macs2 callpeak \
                --format ${format} \
                --treatment ${sample}.${dupType}.bam \
                --name ${sample}.noctrl \
                --outdir ${outdir} \
                --broad \
                -g ${g} \
                -q ${qvalue} \
                -B \
                --SPMR \
                --nomodel
            echo ""

            echo "calling narrow peaks: ${sample} without control"
            macs2 callpeak \
                --format ${format} \
                --treatment ${sample}.${dupType}.bam \
                --name ${sample}.noctrl \
                --outdir ${outdir} \
                -g ${g} \
                -q ${qvalue} \
                -B \
                --SPMR \
                --nomodel
            echo ""
        fi
    done <samplesPeaks.txt

    cd ${dir}
    gzip *bdg

    mv *.bdg.gz     bedgraph
    mv *.broadPeak  broadPeak
    mv *.gappedPeak gappedPeak
    mv *.narrowPeak narrowPeak
    mv *.xls        xls
    mv *.summits*   summits

    cd ..

else
    echo "wrong peak calling algorithm"
fi

echo "end of peakCalling"
echo ""
