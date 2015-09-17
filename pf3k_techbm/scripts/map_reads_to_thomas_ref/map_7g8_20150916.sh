#!/bin/bash

################################################################################
# Introduction
################################################################################
# This should do a "manual" run (as opposed to a vr-pipe run) of the Pf3k
# pipeline on the 5 validation samples

# To run (ideally using screen):
# cd $HOME/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline
# $HOME/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/gatk_pipeline_20150911.sh >> $HOME/src/github/malariagen/methods-dev/pf3k_techbm/log/gatk_pipeline_20150911.log


################################################################################
# Pre requisities
################################################################################
# 
# If installing GATK, this must be downloaded manually from
# https://www.broadinstitute.org/gatk/download/auth?package=GATK and put in home
# directory. This can't be downloaded automatically using wget.



################################################################################
# Set up envirnoment variables
################################################################################

# catch errors
# set -e
# set -o pipefail

# directories
export ORIGINAL_DIR=`pwd`
export PROCESSED_DATA_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/map_to_thomas_ref"
export MAPPING_DIR="bwa_mem"
export BAMS_DIR="sort_markdup"
export UNFILTERED_DIR="UnifiedGenotyper"
export TEMP_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/tmp"
export REF_GENOME_DIR="${PROCESSED_DATA_DIR}/reference_genomes"
export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4"
# export PROCESSED_DATA_DIR="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/assembled_samples"
# export REF_GENOME_DIR="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/Pf3D7_GeneDB"
# export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt"
export SNPEFF_DIR="${OPT_DIR}/snpeff/snpEff"
export CROSSES_DIR="/nfs/team112_internal/oxford_mirror/data/plasmodium/pfalciparum/pf-crosses/data/public/1.0"

mkdir -p ${REF_GENOME_DIR}
cp /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/Pf7G8.Jun2015.fasta ${REF_GENOME_DIR}


# parameters
export MAX_ALTERNATE_ALLELES=6
export GVCF_BATCH_SIZE=3
# export VQSR_ANNOTATIONS_SNP="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
# export VQSR_ANNOTATIONS_INDEL="-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
# export MAX_GAUSSIANS_SNP=8
# export MAX_GAUSSIANS_INDEL=4

# manifests
export SAMPLE_MANIFEST="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/meta/validation_sample_bams.txt"
# export SHUFFLED_SAMPLE_MANIFEST="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/meta/validation_sample_bams_shuffled.txt"

# reference genome
export REF_GENOME_DATE="2015-08"
# export REF_GENOME_FASTA_URL="http://www.plasmodb.org/common/downloads/release-25/Pfalciparum3D7/fasta/data/PlasmoDB-25_Pfalciparum3D7_Genome.fasta"
export REF_GENOME_FASTA_URL="ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/${REF_GENOME_DATE}/Pfalciparum.genome.fasta.gz"
export REF_GENOME_GFF_URL="ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/${REF_GENOME_DATE}/Pfalciparum.gff3.gz"
# export APICOPLAST_REF_GENOME_FASTA="/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011/Pf3D7_v3.fasta"
# export APICOPLAST_SEQUENCE_NAME="PF_apicoplast_genome_1"
export SNPEFF_DB="Pfalciparum_GeneDB_Aug2015"
export SNPEFF_CONFIG_FN="${SNPEFF_DIR}/snpEff.config"

# data files
export REF_GENOME="${REF_GENOME_DIR}/Pf7G8.Jun2015.fasta"
export REF_GENOME_INDEX="${REF_GENOME_DIR}/Pf7G8.Jun2015.fasta.fai"
export REF_GENOME_DICTIONARY="${REF_GENOME_DIR}/Pf7G8.Jun2015.dict"
# export REF_GENOME="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/roamato/Pf3D7_v3/3D7_sorted.fa"
# export REF_GENOME_INDEX="/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/roamato/Pf3D7_v3/3D7_sorted.fa.fai"
export REGIONS_FN="$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz"

# software versions
export PICARD_VERSION="1.137"
export BWA_VERSION="0.7.12"
export SAMTOOLS_VERSION="1.2"
export GATK_VERSION="3.4-46"
export VCFTOOLS_VERSION="0.1.10"
export SNPEFF_VERSION="4_1i"
# export VCFTOOLS_VERSION="0.1.12b"

# executables
export JAVA7_EXE="/software/jre1.7.0_25/bin/java"
export SAMTOOLS_EXE="${OPT_DIR}/samtools/samtools-${SAMTOOLS_VERSION}/samtools"
export PICARD_EXE="${JAVA7_EXE} -jar ${OPT_DIR}/picard/picard-tools-${PICARD_VERSION}/picard.jar"
export SNPEFF_EXE="${JAVA7_EXE} -jar ${SNPEFF_DIR}/snpEff.jar"
# export SNPEFF_EXE="${JAVA7_EXE} -jar /nfs/team112_internal/production/tools/bin/snpEff_4_1/snpEff.jar"
export BWA_EXE="${OPT_DIR}/bwa/bwa-${BWA_VERSION}/bwa"
export GATK_EXE="${JAVA7_EXE} -jar ${OPT_DIR}/gatk/GenomeAnalysisTK.jar"
export FIRST_LAST_100BP_EXE="python $HOME/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/first_last_100bp.py"
export VCF_ANNOTATE_EXE="${OPT_DIR}/vcftools/vcftools_${VCFTOOLS_VERSION}/perl/vcf-annotate"


################################################################################
# Functions
################################################################################

get_RG () {
    $SAMTOOLS_EXE view -H "$1" | grep '^@RG'
}



# Generate BWA index
if [ ! -s ${REF_GENOME}.bwt ] || [ ! -s ${REF_GENOME}.pac ] || [ ! -s ${REF_GENOME}.ann ] || [ ! -s ${REF_GENOME}.amb ] || [ ! -s ${REF_GENOME}.sa ]; then
    $BWA_EXE index ${REF_GENOME}
fi

# Generate samtools index
if [ ! -s ${REF_GENOME}.fai ]; then
    $SAMTOOLS_EXE faidx ${REF_GENOME}
fi

# Generate sequence dictionary
if [ ! -s ${REF_GENOME_DICTIONARY} ]; then
    $PICARD_EXE CreateSequenceDictionary \
    REFERENCE=${REF_GENOME} \
    OUTPUT=${REF_GENOME_DICTIONARY}
fi

number_of_samples=$(wc -l < ${SAMPLE_MANIFEST})
number_of_chromosomes=$(wc -l < ${REF_GENOME_INDEX})

################################################################################
# Create data directories
################################################################################

if [ ! -d "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/_sams" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/_sams"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files"
fi
if [ ! -d "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/${UNFILTERED_DIR}/_vcfs" ]; then
    mkdir -p "${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/${UNFILTERED_DIR}/_vcfs"
fi
if [ ! -d "${TEMP_DIR}" ]; then
    mkdir -p "${TEMP_DIR}"
fi



################################################################################
# Miscellaneous setup
################################################################################

# # Shuffle sample manifest
# if [ ! -s ${SHUFFLED_SAMPLE_MANIFEST} ]; then
#     shuf ${SAMPLE_MANIFEST} > ${SHUFFLED_SAMPLE_MANIFEST}
# fi



################################################################################
# Run pipeline
################################################################################

# Remap with bwa mem, including subsetting 250bp reads to 100bp and converting cram
# for (( i=1; i<=${number_of_samples}; i++ ));
for (( i=1; i<=1; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    original_bam_fn=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f2`
    bwa_mem_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/_sams/${sample_name}.bwa_mem.sam"
    read_group_info=`get_RG "${original_bam_fn}" | sed 's/\t/\\\\t/g'`
    if [ ! -s ${bwa_mem_fn} ]; then
        if [[ "${original_bam_fn}" == *.bam ]]; then
            ${SAMTOOLS_EXE} bamshuf -uOn 128 "${original_bam_fn}" tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${BWA_EXE} mem -M -R "${read_group_info}" -p ${REF_GENOME} - > ${bwa_mem_fn} 2> /dev/null
        elif [[ "${original_bam_fn}" == *.cram ]]; then
            ${SAMTOOLS_EXE} view -b "${original_bam_fn}" | \
            ${SAMTOOLS_EXE} bamshuf -uOn 128 - tmp | \
            ${SAMTOOLS_EXE} bam2fq - | \
            ${FIRST_LAST_100BP_EXE} - | \
            ${BWA_EXE} mem -M -R "${read_group_info}" -p ${REF_GENOME} - > ${bwa_mem_fn} 2> /dev/null
        fi
    fi
done


# Sort and mark duplicates
# for (( i=1; i<=${number_of_samples}; i++ ));
for (( i=1; i<=1; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    bwa_mem_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/_sams/${sample_name}.bwa_mem.sam"
    sorted_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.sorted.bam"
    dedup_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.bam"
    dedup_index_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.bai"
    dedup_metrics_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.metrics"
    if [ ! -s ${sorted_fn} ]; then
        ${PICARD_EXE} SortSam \
            INPUT=${bwa_mem_fn} \
            OUTPUT=${sorted_fn} \
            SORT_ORDER=coordinate \
            TMP_DIR=${TEMP_DIR} #2> /dev/null
    fi
    if [ ! -s ${dedup_fn} ]; then
        ${PICARD_EXE} MarkDuplicates \
            INPUT=${sorted_fn} \
            OUTPUT=${dedup_fn} \
            METRICS_FILE=${dedup_metrics_fn} \
            TMP_DIR=${TEMP_DIR} #2> /dev/null#2> /dev/null
    fi
    if [ ! -s ${dedup_index_fn} ]; then
        ${PICARD_EXE} BuildBamIndex \
            INPUT=${dedup_fn} \
            TMP_DIR=${TEMP_DIR} #2> /dev/null#2> /dev/null
    fi
done


# Indel realignment
for (( i=1; i<=${number_of_samples}; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    dedup_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.bam"
    targets_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.targets.list"
    realigned_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.realigned.bam"
    if [ ! -s ${targets_fn} ]; then
        ${GATK_EXE} \
            -T RealignerTargetCreator \
            -R ${REF_GENOME} \
            -I ${dedup_fn} \
            -o ${targets_fn} 2> /dev/null
    fi
    if [ ! -s ${realigned_fn} ]; then
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
    realigned_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.realigned.bam"
    recal_table_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_final_bams/${sample_name}.bwa_mem.recal.table"
    post_recal_table_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_final_bams/${sample_name}.bwa_mem.post_recal.table"
    recal_plots_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_final_bams/${sample_name}.bwa_mem.recal.pdf"
    recal_bam_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_final_bams/${sample_name}.bwa_mem.recal.bam"
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
    # Note the following is in recommended best practices, but makes no difference to 
    # if [ ! -s ${post_recal_table_fn} ]; then
    #     ${GATK_EXE} \
    #         -T BaseRecalibrator \
    #         -R ${REF_GENOME} \
    #         -I ${realigned_fn} \
    #         -knownSites ${CROSSES_DIR}/7g8_gb4.combined.final.vcf.gz \
    #         -knownSites ${CROSSES_DIR}/hb3_dd2.combined.final.vcf.gz \
    #         -knownSites ${CROSSES_DIR}/3d7_hb3.combined.final.vcf.gz \
    #         -BQSR ${recal_table_fn} \
    #         -o ${post_recal_table_fn} 2> /dev/null
    # fi
    # if [ ! -s ${recal_plots_fn} ]; then
    #     ${GATK_EXE} \
    #         -T AnalyzeCovariates \
    #         -R ${REF_GENOME} \
    #         -before ${recal_table_fn} \
    #         -after ${post_recal_table_fn} \
    #         -plots ${recal_plots_fn} 2> /dev/null
    # fi
    if [ ! -s ${recal_bam_fn} ]; then
        ${GATK_EXE} \
            -T PrintReads \
            -R ${REF_GENOME} \
            -I ${realigned_fn} \
            -BQSR ${recal_table_fn} \
            -o ${recal_bam_fn} 2> /dev/null
    fi
done


# Unified genotyper
# for (( i=1; i<=${number_of_samples}; i++ ));
for (( i=1; i<=1; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    bam_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.bam"
    vcf_filestem="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/${UNFILTERED_DIR}/_vcfs/${sample_name}.raw.snps.indels"
    for (( j=1; j<=${number_of_chromosomes}; j++ ));
    do
        chromosome=`awk "NR==$j" ${REF_GENOME_INDEX} | cut -f1`
        vcf_fn=${vcf_filestem}.${chromosome}.vcf
        if [ ! -s ${vcf_fn} ]; then
            ${GATK_EXE} \
                -T UnifiedGenotyper \
                -R ${REF_GENOME} \
                -I ${bam_fn} \
                -L ${chromosome} \
                -o ${vcf_fn}
             # --annotation HomopolymerRun \
             # --annotation VariantType \
        fi
    done
done

# Unified genotyper
# for (( i=1; i<=${number_of_samples}; i++ ));
for (( i=1; i<=1; i++ ));
do
    sample_name=`awk "NR==$i" ${SAMPLE_MANIFEST} | cut -f1`
    bam_fn="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/_intermediate_files/${sample_name}.bwa_mem.dedup.bam"
    vcf_filestem="${PROCESSED_DATA_DIR}/${MAPPING_DIR}/${BAMS_DIR}/${UNFILTERED_DIR}/_vcfs/${sample_name}.raw.haploid.hrun.snps.indels"
    for (( j=1; j<=${number_of_chromosomes}; j++ ));
    do
        chromosome=`awk "NR==$j" ${REF_GENOME_INDEX} | cut -f1`
        vcf_fn=${vcf_filestem}.${chromosome}.vcf
        if [ ! -s ${vcf_fn} ]; then
            ${GATK_EXE} \
                -T UnifiedGenotyper \
                -R ${REF_GENOME} \
                -I ${bam_fn} \
                -L ${chromosome} \
                --sample_ploidy 1 \
                --genotype_likelihoods_model BOTH \
                -o ${vcf_fn} \
                --annotation HomopolymerRun
             # --annotation VariantType \
        fi
    done
done


# Combine and genotype GVCFs
for (( chromosome_index=1; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    combined_gvcf_list_filename="${PROCESSED_DATA_DIR}/vcfs/gvcf/lists/gvcfs.${chromosome}.all.list"
    if [ -f ${combined_gvcf_list_filename} ]; then
        rm ${combined_gvcf_list_filename}
    fi
    genotyped_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/raw.snps.indels.${chromosome}.vcf"
    batch_index=1
    position_in_batch=0
    for (( sample_index=1; sample_index<=${number_of_samples}; sample_index++ ));
    do
        position_in_batch=$((${position_in_batch}+1))
        if [ ${position_in_batch} -gt ${GVCF_BATCH_SIZE} ]; then
            position_in_batch=1
            batch_index=$((${batch_index}+1))
        fi
        if [ ${position_in_batch} == 1 ]; then
            gvcf_list_filename="${PROCESSED_DATA_DIR}/vcfs/gvcf/lists/gvcfs.${chromosome}.${batch_index}.list"
            if [ -f ${gvcf_list_filename} ]; then
                rm ${gvcf_list_filename}
            fi
            combined_gvcf_fn="${PROCESSED_DATA_DIR}/vcfs/gvcf/combined_gvcf/combined.${chromosome}.${batch_index}.g.vcf"
            touch ${combined_gvcf_list_filename}
            echo ${combined_gvcf_fn} >> ${combined_gvcf_list_filename}
        fi
        sample_name=`awk "NR==$sample_index" ${SAMPLE_MANIFEST} | cut -f1`
        gvcf_fn="${PROCESSED_DATA_DIR}/vcfs/gvcf/samples/${sample_name}.raw.snps.indels.${chromosome}.g.vcf"
        touch ${gvcf_list_filename}
        echo ${gvcf_fn} >> ${gvcf_list_filename}
        # echo ${position_in_batch} ${GVCF_BATCH_SIZE} ${sample_index} ${number_of_samples}
        if [ ${position_in_batch} == ${GVCF_BATCH_SIZE} ] || [ ${sample_index} == ${number_of_samples} ]; then
            if [ ! -s ${combined_gvcf_fn} ]; then
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
            --max_alternate_alleles ${MAX_ALTERNATE_ALLELES} \
            --annotation QualByDepth \
            --annotation FisherStrand \
            --annotation StrandOddsRatio \
            --annotation VariantType \
            -o ${genotyped_vcf_fn}
         # 2> /dev/null
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


# Apply VQSR
for (( chromosome_index=1; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    genotyped_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/raw.snps.indels.${chromosome}.vcf"
    unfiltered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/unfiltered.SNP.${chromosome}.vcf"
    filtered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.SNP.${chromosome}.vcf"
    recal_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_SNP.recal"
    tranches_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_SNP.tranches"
    if [ ! -s ${unfiltered_vcf_fn} ]; then
        ${GATK_EXE} \
            -T SelectVariants \
            -R ${REF_GENOME} \
            -V ${genotyped_vcf_fn} \
            -selectType SNP \
            -o ${unfiltered_vcf_fn} 2> /dev/null
    fi
    if [ ! -s ${filtered_vcf_fn}.gz ]; then
        ${GATK_EXE} \
            -T ApplyRecalibration \
            -R ${REF_GENOME} \
            -input ${unfiltered_vcf_fn} \
            -tranchesFile ${tranches_fn} \
            -recalFile ${recal_fn} \
            --ts_filter_level 99.5 \
            -mode SNP \
            -o ${filtered_vcf_fn} 2> /dev/null
        bgzip -f ${filtered_vcf_fn}
        tabix -p vcf -f ${filtered_vcf_fn}.gz
    fi
    unfiltered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/unfiltered.INDEL.${chromosome}.vcf"
    filtered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.INDEL.${chromosome}.vcf"
    recal_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_INDEL.recal"
    tranches_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/recal/recalibrate_INDEL.tranches"
    if [ ! -s ${unfiltered_vcf_fn} ]; then
        ${GATK_EXE} \
            -T SelectVariants \
            -R ${REF_GENOME} \
            -V ${genotyped_vcf_fn} \
            -xlSelectType SNP  \
            -o ${unfiltered_vcf_fn} 2> /dev/null
    fi
    if [ ! -s ${filtered_vcf_fn}.gz ]; then
        ${GATK_EXE} \
            -T ApplyRecalibration \
            -R ${REF_GENOME} \
            -input ${unfiltered_vcf_fn} \
            -tranchesFile ${tranches_fn} \
            -recalFile ${recal_fn} \
            --ts_filter_level 99.0 \
            -mode INDEL \
            -o ${filtered_vcf_fn} 2> /dev/null
        bgzip -f ${filtered_vcf_fn}
        tabix -p vcf -f ${filtered_vcf_fn}.gz
    fi
done


# Apply SnpEff
for (( chromosome_index=1; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    filtered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.SNP.${chromosome}.vcf.gz"
    snpeff_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.snpeff.SNP.${chromosome}.vcf"
    snpeff_annotated_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.snpeff_annotated.SNP.${chromosome}.vcf"
    annotated_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.annotated.SNP.${chromosome}.vcf"    
    if [ ! -s ${snpeff_vcf_fn} ]; then
        ${SNPEFF_EXE} \
            -v -o gatk ${SNPEFF_DB} \
            ${filtered_vcf_fn} \
            -no-downstream \
            -no-upstream \
            -onlyProtein \
            -c /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_3/snpeff/snpEff/snpEff.config \
            -dataDir /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_3/snpeff/snpEff/data/ \
            > ${snpeff_vcf_fn} \
            # 2> /dev/null
    fi
    if [ ! -s ${snpeff_annotated_vcf_fn} ]; then
        ${GATK_EXE} \
            -T VariantAnnotator \
            -R ${REF_GENOME} \
            -A SnpEff \
            --variant ${filtered_vcf_fn} \
            --snpEffFile ${snpeff_vcf_fn} \
            -o ${snpeff_annotated_vcf_fn} \
            2> /dev/null
    fi
    if [ ! -s ${annotated_vcf_fn} ]; then
        cat ${snpeff_annotated_vcf_fn} \
        | ${VCF_ANNOTATE_EXE} -a ${REGIONS_FN} \
           -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
           -c CHROM,FROM,TO,INFO/RegionType \
        > ${annotated_vcf_fn}
    fi
    filtered_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.INDEL.${chromosome}.vcf.gz"
    snpeff_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.snpeff.INDEL.${chromosome}.vcf"
    snpeff_annotated_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.snpeff_annotated.INDEL.${chromosome}.vcf"
    annotated_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.annotated.INDEL.${chromosome}.vcf"
    if [ ! -s ${snpeff_vcf_fn} ]; then
        ${SNPEFF_EXE} \
            -v -o gatk ${SNPEFF_DB} \
            ${filtered_vcf_fn} \
            -no-downstream \
            -no-upstream \
            -onlyProtein \
            -c /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_3/snpeff/snpEff/snpEff.config \
            -dataDir /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_3/snpeff/snpEff/data/ \
            > ${snpeff_vcf_fn} \
            # 2> /dev/null
    fi
    if [ ! -s ${snpeff_annotated_vcf_fn} ]; then
        ${GATK_EXE} \
            -T VariantAnnotator \
            -R ${REF_GENOME} \
            -A SnpEff \
            --variant ${filtered_vcf_fn} \
            --snpEffFile ${snpeff_vcf_fn} \
            -o ${snpeff_annotated_vcf_fn} \
            2> /dev/null
    fi
    if [ ! -s ${annotated_vcf_fn} ]; then
        cat ${snpeff_annotated_vcf_fn} \
        | ${VCF_ANNOTATE_EXE} -a ${REGIONS_FN} \
           -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
           -c CHROM,FROM,TO,INFO/RegionType \
        > ${annotated_vcf_fn}
    fi
done


# Combine SNP and INDEL
for (( chromosome_index=1; chromosome_index<=${number_of_chromosomes}; chromosome_index++ ));
do
    chromosome=`awk "NR==$chromosome_index" ${REF_GENOME_INDEX} | cut -f1`
    annotated_snp_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.annotated.SNP.${chromosome}.vcf"
    annotated_indel_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.annotated.INDEL.${chromosome}.vcf"
    annotated_combined_vcf_fn="${PROCESSED_DATA_DIR}/vcfs/vcf/filtered.annotated.SNP_INDEL.${chromosome}.vcf"
    if [ ! -s ${annotated_combined_vcf_fn}.gz ]; then
        ${GATK_EXE} \
            -T CombineVariants \
            -R ${REF_GENOME} \
            --variant:snp ${annotated_snp_vcf_fn} \
            --variant:indel ${annotated_indel_vcf_fn} \
            -o ${annotated_combined_vcf_fn} \
            -genotypeMergeOptions PRIORITIZE \
            -priority snp,indel \
            2> /dev/null
        bgzip -f ${annotated_combined_vcf_fn}
        tabix -p vcf -f ${annotated_combined_vcf_fn}.gz
    fi
done

