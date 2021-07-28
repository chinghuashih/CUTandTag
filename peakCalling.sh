#!/bin/bash
#SBATCH -p standard  -o peaksCalling.log -t 06:00:00
#SBATCH -c 8 --mem=96G

# SEACR requires paired-end sequencing as inputs, both for experiments and the control(IgG) if having IgG as the control.

# format:
# the name of experimental sample, and the name of control sample; separated by a space
#
# If there is no control sample, leave it with NA

module load bedtools
module load R

# peak calling methods: SEACR or MACS2
peakCalling="SEACR"


echo "peak calling using ${peakCalling}"
echo ""
echo ""
echo ""

if [ "${peakCalling}" == "SEACR" ]; then

    mkdir peaks_SEACR
    while read exp
    do
        echo $exp
        IFS=' ' read -r -a array <<< "$exp"

        sample=${array[0]}
        control=${array[1]}

        echo "treatment: $sample"
        echo "control: $control"

        echo "Preparing bedgraph files: $sample and $control"
        # control
        bedtools genomecov -ibam $control.dupMark.bam -bg > $control.bedgraph
        # experiments
        bedtools genomecov -ibam $sample.dupMark.bam -bg > $sample.bedgraph
        # IgG

#       bedtools bamtobed -bedpe -i $sample.bam > $sample.bed
#       awk '$1==$4 && $6-$2 < 1000 {print $0}' $sample.bed > $sample.clean.bed
#       cut -f 1,2,6 $sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > $sample.fragments.bed
#       bedtools genomecov -bg -i $sample.fragments.bed -g my.genome > $sample.fragments.bedgraph


        if [ "$control" != "NA" ]; then
                # calling peaks with the stringent threshold
                echo "Calls enriched regions in $sample using normalized $control as the control track with stringent threshold"
                bash SEACR_1.2.sh $sample.bedgraph $control.bedgraph norm stringent $sample.stringent_norm
                echo "Calls enriched regions in $sample using unnormalized No_Antibody as the control track with stringent threshold"
                bash SEACR_1.2.sh $sample.bedgraph $control.bedgraph non stringent $sample.stringent_non
                echo ""

                # calling peaks with the relaxed threshold
                echo "Calls enriched regions in $sample using normalized $control as the control track with relaxed threshold"
                bash SEACR_1.2.sh $sample.bedgraph $control.bedgraph norm relaxed $sample.relaxed_norm
                echo "Calls enriched regions in $sample using unnormalized $control as the control track with relaxed threshold"
                bash SEACR_1.2.sh $sample.bedgraph $control.bedgraph non relaxed $sample.relaxed_non
                echo ""
        fi

        echo "Calls enriched regions in $sample by selecting the top 1% of regions by AUC"
        bash SEACR_1.2.sh $sample.bedgraph 0.01 norm stringent $sample.AUC_norm01
        echo ""
        echo ""
    done <samplesPeaks.txt

    mv *bedgraph peaks_SEACR
    mv *bed peaks_SEACR
    #mv *stringent.bed peaks_SEACR
    #mv *relaxed.bed peaks_SEACR

    ## description
    #  Field 1: Target data bedgraph file in UCSC bedgraph format
    #  Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling.
    #           Alternatively, a numeric threshold *n* between 0 and 1 returns the top *n* fraction of peaks
    #           based on total signal within peaks.
    #  Field 3: "norm" denotes normalization of control to target data,
    #           "non" skips this behavior. 
    #           "norm" is recommended unless experimental and control data are already rigorously normalized
    #           to each other (e.g. via spike-in).
    #  Field 4: "relaxed" uses a total signal threshold between the knee and peak of the total signal curve, 
    #           and corresponds to the “relaxed” mode described in the text
    #           "stringent" uses the peak of the curve, and corresponds to “stringent” mode.
    #  Field 5: Output prefix

elif [ "${peakCalling}" == "MACS2" ]; then
    # format:
    # the name of experimental sample, and the name of control sample; separated by a space
    #
    # If there is no control sample, leave it with NA

    inputFormat="bam"
    # bam for BAMPE
    # bed for BEDPE
    # auto if not sure
    qvalue=0.05
    outdir="peaks_macs2"

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

    mkdir ${outdir}

    while read exp
    do
        echo ${exp}
        IFS=' ' read -r -a array <<< "${exp}"

        sample=${array[0]}
        control=${array[1]}

        echo "treatment: ${sample}"
        echo "control: ${control}"

        if [ "$control" != "NA" ]; then

                # calling peaks
                echo "calling broad peaks: ${sample} with ${control}"
                macs2 callpeak \
                        --format ${format} \
                        --treatment ${sample}.dupMark.bam \
                        --control ${control}.dupMark.bam \
                        --name ${sample}.broad \
                        --outdir ${outdir} \
                        --broad \
                        -g hs \
                        -q ${qvalue} \
                        -B \
                        --SPMR \
                        --nomodel
                echo ""               
                
                echo "calling narrow peaks: $sample with $control"
                macs2 callpeak \
                        --format $format \
                        --treatment $sample.dupMark.bam \
                        --control $control.dupMark.bam \
                        --name $sample.narrow \
                        --outdir $outdir \
                        -g hs \
                        -q $qvalue \
                        -B \
                        --SPMR \
                        --nomodel
                echo ""
        fi
        # calling peaks without control
        echo "calling broad peaks: $sample without control"
        macs2 callpeak \
                --format $format \
                --treatment $sample.dupMark.bam \
                --name $sample.broad_noctrl \
                --outdir $outdir \
                --broad \
                -g hs \
                -q $qvalue \
                -B \
                --SPMR \
                --nomodel
        echo ""

        echo "calling narrow peaks: $sample without control"
        macs2 callpeak \
                --format $format \
                --treatment $sample.dupMark.bam \
                --name $sample.narrow_noctrl \
                --outdir $outdir \
                -g hs \
                -q $qvalue \
                -B \
                --SPMR \
                --nomodel
        echo ""
    done <samplesPeaks.txt

    cd $dir
    gzip *bdg

else
    echo "wrong peak calling algorithm"
    
fi
