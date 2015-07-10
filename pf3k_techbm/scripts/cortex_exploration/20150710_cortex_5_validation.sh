bash

# Downlaod Cortex from github
cd /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex
git clone https://github.com/iqbal-lab/cortex.git
cd cortex
bash install.sh 2> /dev/null
cd

export DATA_DIR=/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/cortex_exploration/data
# export CORTEX_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21
export CORTEX_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/cortex
export CORTEX_EXE=$CORTEX_DIR/bin/cortex_var_31_c1
export REF_GENOME=$DATA_DIR/roamato/Pf3D7_v3/3D7_sorted.fa
export STAMPY_EXE="/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/stampy/stampy-1.0.27/stampy.py"
export SAMTOOLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/samtools/samtools-1.2/samtools
# export RUN_CALLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21/scripts/calling/run_calls.pl
export RUN_CALLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/cortex/scripts/calling/run_calls.pl
# export VCFTOOLS_DIR=/nfs/team112_internal/production/tools/bin/
export VCFTOOLS_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/vcftools/vcftools_0.1.12b
export SAMPLE_NAMES=(7G8 GB4 KE01 KH02 GN01)
export BAM_CRAM_FNS=(\
/nfs/team112_internal/production_files/Pf/4_0/PFprog3/PG0083_C/PG0083_C.bam \
/nfs/team112_internal/production_files/Pf/4_0/PFprog3/PG0084_C/PG0084_C.bam \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#3.cram \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#2.cram \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#6.cram\
)
export FASTQ_DIR=$DATA_DIR/fastq

cd $CORTEX_DIR
make cortex_var
make NUM_COLS=2 cortex_var
make NUM_COLS=3 cortex_var
cd -

mkdir -p $DATA_DIR/binaries
mkdir -p $DATA_DIR/fofns

for sample_name in $SAMPLE_NAMES
do
    mkdir -p $DATA_DIR/cortex_output/$sample_name/log
done

echo $REF_GENOME > $DATA_DIR/fofns/ref.fofn

$CORTEX_EXE --kmer_size 31 \
--mem_height 20 --mem_width 100 \
--se_list $DATA_DIR/fofns/ref.fofn \
--max_read_len 10000 \
--dump_binary $DATA_DIR/binaries/ref.k31.ctx --sample_id REF

cd $DATA_DIR/roamato/Pf3D7_v3
$STAMPY_EXE -G 3D7_sorted $REF_GENOME
$STAMPY_EXE -g 3D7_sorted -H 3D7_sorted
cd -

# Convert bams/crams to fastq
number_of_samples=${#SAMPLE_NAMES[@]}
for (( i=0; i<${number_of_samples}; i++ ));
do
    echo ${BAM_CRAM_FNS[$i]}
    if [[ "${BAM_CRAM_FNS[$i]}" == *.bam ]]
        then
            $SAMTOOLS_EXE bamshuf -uOn 128 ${BAM_CRAM_FNS[$i]} tmp | $SAMTOOLS_EXE bam2fq - > ${FASTQ_DIR}/${SAMPLE_NAMES[$i]}.fastq
            # echo bam ${BAM_CRAM_FNS[$i]}
        elif [[ "${BAM_CRAM_FNS[$i]}" == *.cram ]]
        then
            $SAMTOOLS_EXE view -b ${BAM_CRAM_FNS[$i]} | $SAMTOOLS_EXE bamshuf -uOn 128 - tmp | $SAMTOOLS_EXE bam2fq - > ${FASTQ_DIR}/${SAMPLE_NAMES[$i]}.fastq
            # echo cram ${BAM_CRAM_FNS[$i]}
        fi
    echo
done

# Create index files
number_of_samples=${#SAMPLE_NAMES[@]}
for (( i=0; i<${number_of_samples}; i++ ));
do
    echo ${FASTQ_DIR}/${SAMPLE_NAMES[$i]}.fastq > $DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.se.fofn
    printf "%s\t%s\t.\t.\n" "${SAMPLE_NAMES[$i]}" "$DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.se.fofn" > $DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.fastaq.index
    # echo "${SAMPLE_NAMES[$i]} \t $DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.se.fofn\\t.\\t." > $DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.fastaq.index
done



# Run Cortex
number_of_samples=${#SAMPLE_NAMES[@]}
for (( i=0; i<${number_of_samples}; i++ ));
do
    echo ${SAMPLE_NAMES[$i]}
	mkdir -p $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/log
    perl $RUN_CALLS_EXE --first_kmer 31 \
    --fastaq_index $DATA_DIR/fofns/${SAMPLE_NAMES[$i]}.fastaq.index \
    --auto_cleaning yes \
    --bc yes \
    --pd no \
    --outdir $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]} \
    --outvcf ${SAMPLE_NAMES[$i]}.vcf \
    --ploidy 2 \
    --stampy_hash $DATA_DIR/roamato/Pf3D7_v3/3D7_sorted \
    --stampy_bin $STAMPY_EXE \
    --list_ref_fasta $DATA_DIR/fofns/ref.fofn \
    --refbindir $DATA_DIR/binaries/ \
    --genome_size 24000000 \
    --qthresh 15 \
    --mem_height 20 \
    --mem_width 150 \
    --vcftools_dir $VCFTOOLS_DIR \
    --do_union yes \
    --ref CoordinatesAndInCalling \
    --workflow independent \
    --logfile $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/log/${SAMPLE_NAMES[$i]}.log.txt
done
# --apply_pop_classifier

# Create VCFs fofn
number_of_samples=${#SAMPLE_NAMES[@]}
for (( i=0; i<${number_of_samples}; i++ ));
do
    printf "%s\n" "$DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/vcfs/${SAMPLE_NAMES[$i]}.vcf_union_BC_calls_k31.decomp.vcf" >> $DATA_DIR/fofns/all_decomp_vcfs.fofn
done

cortex_output/GN01/vcfs/GN01.vcf_union_BC_calls_k31.decomp.vcf



# Combine VCFs
mkdir -p $DATA_DIR/cortex_output/combined

perl $CORTEX_DIR/scripts/analyse_variants/combine/combine_vcfs.pl \
$DATA_DIR/fofns/all_decomp_vcfs.fofn \
$CORTEX_DIR \
$VCFTOOLS_DIR \
$CORTEX_DIR/scripts/analyse_variants/bioinf-perl/vcf_scripts \
$CORTEX_DIR/scripts/analyse_variants \
$CORTEX_DIR/scripts/analyse_variants/combine \
$DATA_DIR/cortex_output/combined \
combined \
31 \
Pf3D7_v3 \
$REF_GENOME \
$DATA_DIR/binaries/ref.k31.ctx \
20 \
150 \
$DATA_DIR/cortex_output



# genotype each samples against combined sites
mkdir -p $DATA_DIR/cortex_output/7G8/genotyped_vcf

perl $CORTEX_DIR/scripts/calling/genotype_1sample_against_sites.pl \
    --cortex_dir $CORTEX_DIR \
    --invcf $DATA_DIR/cortex_output/combined/combined.sites_vcf \
    --bubble_graph $DATA_DIR/cortex_output/combined/combined.pseudo_callfile.branches.k31.ctx \
    --bubble_callfile  $DATA_DIR/cortex_output/combined/combined.pseudo_callfile \
    --outdir $DATA_DIR/cortex_output/7G8/genotyped_vcf \
    --sample 7G8 \
    --ref_overlap_bubble_graph $DATA_DIR/cortex_output/combined/ref.k31.ctx.list \
    --ref_overlap_bubble_graph $DATA_DIR/cortex_output/combined/ref.k31.ctx.list_intersect_bubbles.ctx \
    --genome_size 24000000 \
    --sample_graph $DATA_DIR/cortex_output/7G8/binaries/cleaned/k31/7G8.kmer31.q15cleaned_7.ctx \
    --overlap_log $DATA_DIR/cortex_output/7G8/genotyped_vcf/log.txt \
    --mem_height 18 \
    --mem_width 100

