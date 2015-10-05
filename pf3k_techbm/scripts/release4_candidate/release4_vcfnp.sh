export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4"
export JAVA7_EXE="/software/jre1.7.0_25/bin/java"
export GATK_EXE="${JAVA7_EXE} -jar ${OPT_DIR}/gatk/GenomeAnalysisTK.jar"

export REF_GENOME_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015"
export REF_GENOME="${REF_GENOME_DIR}/Pfalciparum.genome.fasta"

export VCF_FN="/lustre/scratch109/malaria/pf3k_methods/output/f/d/2/8/336049/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz"
export SAMPLE_FILE="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_samples.txt"
export OUTPUT_FN="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/Pf3D7_01_v3.5validation.vcf.gz"

${GATK_EXE} \
    -T SelectVariants \
    -R ${REF_GENOME} \
    --variant ${VCF_FN} \
    --out /dev/stdout 
    --sample_file ${SAMPLE_FILE} \
    --excludeNonVariants \
    | bgzip -c > ${OUTPUT_FN}

tabix -p vcf ${OUTPUT_FN}
