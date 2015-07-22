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

# directories
export ORIGINAL_DIR=`pwd`
export PROCESSED_DATA_DIR="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/assembled_samples"
export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt"

# manifests
export SAMPLE_MANIFEST="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/meta/validation_sample_bams.txt"

# data files
export REF_GENOME='/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/roamato/Pf3D7_v3/3D7_sorted.fa'

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

# parameters


################################################################################
# Functions
################################################################################

get_RG () {
    $SAMTOOLS_EXE view -H "$1" | grep '@RG'
}

get_RG /nfs/team112_internal/production_files/Pf/4_0/PFprog1/PG0008_CW/PG0008_CW.bam


################################################################################
# Install software
################################################################################

# Install picard
if [ ! -f ${OPT_DIR}/picard/picard-tools-${PICARD_VERSION}/picard.jar ]; then
    mkdir -p ${OPT_DIR}/picard
    cd ${OPT_DIR}/picard
    wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip
    unzip picard-tools-${PICARD_VERSION}.zip
    cd ${ORIGINAL_DIR}
fi

# Install SnpEff
if [ ! -f ${SNPEFF_DIR}/snpEff.jar ]; then
    mkdir -p ${OPT_DIR}/snpeff
    cd ${OPT_DIR}/snpeff
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd ${ORIGINAL_DIR}
fi

# Install bwa
if [ ! -f ${BWA_EXE} ]; then
    mkdir -p ${OPT_DIR}/bwa
    cd ${OPT_DIR}/bwa
    wget http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2
    tar -xf bwa-${BWA_VERSION}.tar.bz2
    cd bwa-${BWA_VERSION}
    make 2> /dev/null
    cd ${ORIGINAL_DIR}
fi

# Install samtools
if [ ! -f ${SAMTOOLS_EXE} ]; then
    mkdir -p ${OPT_DIR}/samtools
    cd ${OPT_DIR}/samtools
    wget http://downloads.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    cd samtools-${SAMTOOLS_VERSION}
    make 2> /dev/null
    cd ${ORIGINAL_DIR}
fi

# Install GATK
if [ ! -f ${OPT_DIR}/gatk/GenomeAnalysisTK.jar ]; then
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
# Run pipeline
################################################################################

number_of_samples=$(wc -l < ${SAMPLE_MANIFEST})

# Remap with bwa mem, including subsetting 250bp reads to 100bp and converting cram
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    original_bam_fn=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f2`
    bwa_mem_fn="${PROCESSED_DATA_DIR}/bams/bwa_mem/${sample_name}.bwa_mem.sam"
    read_group_info=`get_RG ${original_bam_fn}`
    echo sample_name ${read_group_info}
    if [ ! -f ${bwa_mem_fn} ]; then
        if [[ "${original_bam_fn}" == *.bam ]]; then
            ${SAMTOOLS_EXE} bamshuf -uOn 128 ${original_bam_fn} tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${BWA_EXE} mem -M -R '${read_group_info}' -p ${REF_GENOME} - > ${bwa_mem_fn} 2> /dev/null
        elif [[ "${original_bam_fn}" == *.cram ]]; then
            ${SAMTOOLS_EXE} view -b ${original_bam_fn} | \
            ${SAMTOOLS_EXE} bamshuf -uOn 128 - tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${FIRST_LAST_100BP_EXE} - | \
            ${BWA_EXE} mem -M -R '${read_group_info}' -p ${REF_GENOME} - > ${bwa_mem_fn} 2> /dev/null
        fi
    fi
    
done
