#!/bin/bash

################################################################################
# Pre requisities
################################################################################
# 
# If installing GATK, this must be downloaded manually from
# https://www.broadinstitute.org/gatk/download/auth?package=GATK and put in home
# directory



################################################################################
# Set up envirnoment variables
################################################################################


# catch errors
set -e
set -o pipefail

# directories
export ORIGINAL_DIR=`pwd`
export PROCESSED_DATA_DIR="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/assembled_samples"
export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt"
export CROSSES_DIR='/nfs/team112_internal/oxford_mirror/data/plasmodium/pfalciparum/pf-crosses/data/public/1.0'

# parameters
export MAX_ALTERNATE_ALLELES=2
export GVCF_BATCH_SIZE=3
export VQSR_ANNOTATIONS_SNP="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
export VQSR_ANNOTATIONS_INDEL="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
export MAX_GAUSSIANS_SNP=8
export MAX_GAUSSIANS_INDEL=4

# manifests
export SAMPLE_MANIFEST="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/meta/validation_sample_bams.txt"
export SHUFFLED_SAMPLE_MANIFEST="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/meta/validation_sample_bams_shuffled.txt"

# data files
export REF_GENOME='/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/roamato/Pf3D7_v3/3D7_sorted.fa'
export REF_GENOME_INDEX='/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/roamato/Pf3D7_v3/3D7_sorted.fa.fai'

# software versions
export PICARD_VERSION="1.136"
export BWA_VERSION="0.7.12"
export SAMTOOLS_VERSION="1.2"
export GATK_VERSION="3.4-46"

# executables
export JAVA7_EXE="/software/jre1.7.0_25/bin/java"
export SAMTOOLS_EXE="${OPT_DIR}/samtools/samtools-${SAMTOOLS_VERSION}/samtools"
export PICARD_EXE="${JAVA7_EXE} -jar ${OPT_DIR}/picard/picard-tools-${PICARD_VERSION}/picard.jar"
export SNPEFF_DIR="${OPT_DIR}/snpeff/snpEff"
export SNPEFF_EXE="${JAVA7_EXE} -jar ${SNPEFF_DIR}/snpEff.jar"
export BWA_EXE="${OPT_DIR}/bwa/bwa-${BWA_VERSION}/bwa"
export GATK_EXE="${JAVA7_EXE} -jar ${OPT_DIR}/gatk/GenomeAnalysisTK.jar"
export FIRST_LAST_100BP_EXE="python $HOME/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/first_last_100bp.py"



################################################################################
# Functions
################################################################################

get_RG () {
    $SAMTOOLS_EXE view -H "$1" | grep '@RG'
}

# get_RG /nfs/team112_internal/production_files/Pf/4_0/PFprog1/PG0008_CW/PG0008_CW.bam



################################################################################
# Install software
################################################################################

# Install picard
if [ ! -s ${OPT_DIR}/picard/picard-tools-${PICARD_VERSION}/picard.jar ]; then
    mkdir -p ${OPT_DIR}/picard
    cd ${OPT_DIR}/picard
    wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip
    unzip picard-tools-${PICARD_VERSION}.zip
    cd ${ORIGINAL_DIR}
fi

# Install SnpEff
if [ ! -s ${SNPEFF_DIR}/snpEff.jar ]; then
    mkdir -p ${OPT_DIR}/snpeff
    cd ${OPT_DIR}/snpeff
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd ${ORIGINAL_DIR}
fi

# Install bwa
if [ ! -s ${BWA_EXE} ]; then
    mkdir -p ${OPT_DIR}/bwa
    cd ${OPT_DIR}/bwa
    wget http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2
    tar -xf bwa-${BWA_VERSION}.tar.bz2
    cd bwa-${BWA_VERSION}
    make 2> /dev/null
    cd ${ORIGINAL_DIR}
fi

# Install samtools
if [ ! -s ${SAMTOOLS_EXE} ]; then
    mkdir -p ${OPT_DIR}/samtools
    cd ${OPT_DIR}/samtools
    wget http://downloads.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    cd samtools-${SAMTOOLS_VERSION}
    make 2> /dev/null
    cd ${ORIGINAL_DIR}
fi

# Install GATK
if [ ! -s ${OPT_DIR}/gatk/GenomeAnalysisTK.jar ]; then
    mkdir -p ${OPT_DIR}/gatk
    cp ~/GenomeAnalysisTK-${GATK_VERSION}.tar.bz2 ${OPT_DIR}/gatk
    cd ${OPT_DIR}/gatk
    tar xjf GenomeAnalysisTK-${GATK_VERSION}.tar.bz2
    cd ${ORIGINAL_DIR}
fi



################################################################################
# Create data directories
################################################################################

if [ ! -d "${PROCESSED_DATA_DIR}/bams/bwa_mem" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/bams/bwa_mem"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/vcfs/gvcf/samples" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/vcfs/gvcf/samples"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/vcfs/gvcf/lists" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/vcfs/gvcf/lists"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/vcfs/gvcf/combined_gvcf" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/vcfs/gvcf/combined_gvcf"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/vcfs/vcf/recal" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/vcfs/vcf/recal"
fi



################################################################################
# Miscellaneous setup
################################################################################

# Shuffle sample manifest
if [ ! -s ${SHUFFLED_SAMPLE_MANIFEST}]; then
    shuf ${SAMPLE_MANIFEST} > ${SHUFFLED_SAMPLE_MANIFEST}
fi



################################################################################
# Run pipeline
################################################################################

number_of_samples=$(wc -l < ${SAMPLE_MANIFEST})
number_of_chromosomes=$(wc -l < ${REF_GENOME_INDEX})


# Remap with bwa mem, including subsetting 250bp reads to 100bp and converting cram
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    original_bam_fn=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f2`
    bwa_mem_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.sam"
    read_group_info=`get_RG "${original_bam_fn}" | sed 's/\t/\\\\t/g'`
    if [ ! -s ${bwa_mem_fn} ]; then
        if [[ "${original_bam_fn}" == *.bam ]]; then
            ${SAMTOOLS_EXE} bamshuf -uOn 128 "${original_bam_fn}" tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${BWA_EXE} mem -M -R ${read_group_info} -p ${REF_GENOME} - > ${bwa_mem_fn} 2> /dev/null
        elif [[ "${original_bam_fn}" == *.cram ]]; then
            ${SAMTOOLS_EXE} view -b "${original_bam_fn}" | \
            ${SAMTOOLS_EXE} bamshuf -uOn 128 - tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${FIRST_LAST_100BP_EXE} - | \
            ${BWA_EXE} mem -M -R ${read_group_info} -p ${REF_GENOME} - > ${bwa_mem_fn}
            # 2> /dev/null
        fi
    fi
done


# Sort and mark duplicates
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    bwa_mem_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.sam"
    sorted_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.sorted.bam"
    dedup_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.dedup.bam"
    dedup_index_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.dedup.bai"
    dedup_metrics_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.dedup.metrics"
    if [ ! -s ${sorted_fn} ]; then
        ${PICARD_EXE} SortSam \
            INPUT=${bwa_mem_fn} \
            OUTPUT=${sorted_fn} \
            SORT_ORDER=coordinate 2> /dev/null
    fi
    if [ ! -s ${dedup_fn} ]; then
        ${PICARD_EXE} MarkDuplicates \
            INPUT=${sorted_fn} \
            OUTPUT=${dedup_fn} \
            METRICS_FILE=${dedup_metrics_fn} 2> /dev/null
    fi
    if [ ! -s ${dedup_index_fn} ]; then
        ${PICARD_EXE} BuildBamIndex \
            INPUT=${dedup_fn} 2> /dev/null
    fi
done


# Indel realignment
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    dedup_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.dedup.bam"
    targets_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.dedup.targets.list"
    realigned_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.realigned.bam"
    if [ ! -s ${targets_fn} ]; then
        ${GATK_EXE} \
            -T RealignerTargetCreator \
            -R ${REF_GENOME} \
            -I ${dedup_fn} \
            -o ${targets_fn} 2> /dev/null
    fi
    if [ ! -s ${targets_fn} ]; then
        ${GATK_EXE} \
            -T IndelRealigner \
            -R ${REF_GENOME} \
            -I ${dedup_fn} \
            -targetIntervals ${targets_fn} \
            -o ${realigned_fn} 2> /dev/null
    fi
done


# Base quality score recalibration (BQSR)
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    realigned_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.realigned.bam"
    recal_table_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.recal.table"
    post_recal_table_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.post_recal.table"
    recal_plots_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.recal.pdf"
    recal_bam_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.recal.bam"
    if [ ! -s ${recal_table_fn} ]; then
        ${GATK_EXE} \
            -T BaseRecalibrator \
            -R ${REF_GENOME} \
            -I ${realigned_fn} \
            -knownSites ${CROSSES_DIR}/7g8_gb4.combined.final.vcf.gz \
            -knownSites ${CROSSES_DIR}/hb3_dd2.combined.final.vcf.gz \
            -knownSites ${CROSSES_DIR}/3d7_hb3.combined.final.vcf.gz \
            -o ${recal_table_fn} 2> /dev/null

    fi
    if [ ! -s ${post_recal_table_fn} ]; then
        ${GATK_EXE} \
            -T BaseRecalibrator \
            -R ${REF_GENOME} \
            -I ${realigned_fn} \
            -knownSites ${CROSSES_DIR}/7g8_gb4.combined.final.vcf.gz \
            -knownSites ${CROSSES_DIR}/hb3_dd2.combined.final.vcf.gz \
            -knownSites ${CROSSES_DIR}/3d7_hb3.combined.final.vcf.gz \
            -BQSR ${recal_table_fn} \
            -o ${post_recal_table_fn} 2> /dev/null
    fi
    if [ ! -s ${recal_plots_fn} ]; then
        ${GATK_EXE} \
            -T AnalyzeCovariates \
            -R ${REF_GENOME} \
            -before ${recal_table_fn} \
            -after ${post_recal_table_fn} \
            -plots ${recal_plots_fn} 2> /dev/null
    fi
    if [ ! -s ${recal_bam_fn} ]; then
        ${GATK_EXE} \
            -T PrintReads \
            -R ${REF_GENOME} \
            -I ${realigned_fn} \
            -BQSR ${recal_table_fn} \
            -o ${recal_bam_fn} 2> /dev/null
    fi
done


# Haplotype Caller
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    recal_bam_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.recal.bam"
    gvcf_filestem="${PROCESSED_DATA_DIR}/vcfs/gvcf/samples/${sample_name}.raw.snps.indels"
    for (( j=1; j<=${number_of_chromosomes}; j++ ));
    do
        chromosome=`awk "NR==$j" ${REF_GENOME_INDEX} | cut -f1`
        gvcf_fn=${gvcf_filestem}.${chromosome}.g.vcf
        if [ ! -s ${gvcf_fn} ]; then
            ${GATK_EXE} \
                -T HaplotypeCaller \
                -R ${REF_GENOME} \
                -I ${recal_bam_fn} \
                -L ${chromosome} \
                --emitRefConfidence GVCF \
                --variant_index_type LINEAR \
                --variant_index_parameter 128000 \
                --max_alternate_alleles ${MAX_ALTERNATE_ALLELES} \
                --annotation HomopolymerRun \
                --annotation VariantType \
                -o ${gvcf_fn} 2> /dev/null
        fi
done


# Combine and genotype GVCFs
for (( chromosome_index=1; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    combined_gvcf_list_filename="${PROCESSED_DATA_DIR}/vcfs/gvcf/lists/gvcfs/${chromosome}.all.list"
    if [ -f ${combined_gvcf_list_filename} ]; then
        rm ${combined_gvcf_list_filename}
    fi
    genotyped_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/raw.snps.indels.${chromosome}.vcf"
    batch_index=1
    position_in_batch=0
    for (( sample_index=1; sample_index<=${number_of_samples}; sample_index++ ));
    do
        position_in_batch=position_in_batch+1
        if [ ${position_in_batch} -gt ${GVCF_BATCH_SIZE}]; then
            position_in_batch=1
            batch_index=batch_index+1
        fi
        if [ ${position_in_batch} == 1 ]; then
            gvcf_list_filename = "${PROCESSED_DATA_DIR}/vcfs/gvcf/lists/gvcfs/${chromosome}.${batch_index}.list"
            if [ -f ${gvcf_list_filename} ]; then
                rm ${gvcf_list_filename}
            fi
            combined_gvcf_fn = "${PROCESSED_DATA_DIR}/gvcf/combined_gvcf/combined.${chromosome}.${batch_index}.g.vcf"
            echo ${combined_gvcf_fn} >> ${combined_gvcf_list_filename}
        fi
        sample_name=`awk "NR==$sample_index" ${SHUFFLED_SAMPLE_MANIFEST} | cut -f1`
        gvcf_fn="${PROCESSED_DATA_DIR}/vcfs/gvcf/samples/${sample_name}.raw.snps.indels.${chromosome}.g.vcf"
        echo ${gvcf_fn} >> gvcf_list_filename
        if [ ${position_in_batch} == ${GVCF_BATCH_SIZE} ] || [ ${sample_index} == ${number_of_samples} ]; then
            if [ ! -s ${combined_gvcf_fn}]; then
                ${GATK_EXE} \
                    -T CombineGVCFs \
                    -R ${REF_GENOME} \
                    --variant ${gvcf_list_filename} \
                    -o ${combined_gvcf_fn} 2> /dev/null
            fi
        fi
    done
    if [ ! -s ${genotyped_vcf_fn} ]; then
        ${GATK_EXE} \
            -T GenotypeGVCFs \
            -R ${REF_GENOME} \
            --variant ${combined_gvcf_list_filename} \
            -o ${genotyped_vcf_fn} 2> /dev/null
    fi
done


# Train VQSR
chromosome=`awk "NR==1" ${REF_GENOME_INDEX} | cut -f1`
genotyped_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/raw.snps.indels.${chromosome}.vcf"
input_line="-input ${genotyped_vcf_fn}"
for (( chromosome_index=2; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    genotyped_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/raw.snps.indels.${chromosome}.vcf"
    input_line="${input_line} -input ${genotyped_vcf_fn}"
done

recal_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_SNP.recal"
tranches_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_SNP.tranches"
if [ ! -s ${recal_fn} ]; then
    ${GATK_EXE} \
        -T VariantRecalibrator \
        -R ${REF_GENOME} \
        ${input_line} \
        -resource\:7g8_gb4,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/7g8_gb4.combined.final.vcf.gz \
        -resource\:hb3_dd2,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/hb3_dd2.combined.final.vcf.gz \
        -resource\:3d7_hb3,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/3d7_hb3.combined.final.vcf.gz \
        ${VQSR_ANNOTATIONS_SNP} \
        --maxGaussians ${MAX_GAUSSIANS_SNP} \
        -mode SNP \
        -recalFile ${recal_fn} \
        -tranchesFile ${tranches_fn} 2> /dev/null
fi

recal_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_INDEL.recal"
tranches_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_INDEL.tranches"
if [ ! -s ${recal_fn} ]; then
    ${GATK_EXE} \
        -T VariantRecalibrator \
        -R ${REF_GENOME} \
        ${input_line} \
        -resource\:7g8_gb4,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/7g8_gb4.combined.final.vcf.gz \
        -resource\:hb3_dd2,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/hb3_dd2.combined.final.vcf.gz \
        -resource\:3d7_hb3,known=false,training=true,truth=true,prior=15.0 ${CROSSES_DIR}/3d7_hb3.combined.final.vcf.gz \
        ${VQSR_ANNOTATIONS_INDEL} \
        --maxGaussians ${MAX_GAUSSIANS_INDEL} \
        -mode INDEL \
        -recalFile ${recal_fn} \
        -tranchesFile ${tranches_fn} 2> /dev/null
fi



