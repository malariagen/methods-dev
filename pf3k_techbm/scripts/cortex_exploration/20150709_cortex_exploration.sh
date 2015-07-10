export DATA_DIR=/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/cortex_exploration/data
export CORTEX_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21
export CORTEX_EXE=$CORTEX_DIR/bin/cortex_var_31_c1
export REF_GENOME=$DATA_DIR/roamato/Pf3D7_v3/3D7_sorted.fa
export STAMPY_EXE="/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/stampy/stampy-1.0.27/stampy.py"
export SAMTOOLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/samtools/samtools-1.2/samtools
export RUN_CALLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21/scripts/calling/run_calls.pl
export VCFTOOLS_DIR=/nfs/team112_internal/production/tools/bin/

cd $CORTEX_DIR
make cortex_var
make NUM_COLS=2 cortex_var
make NUM_COLS=3 cortex_var
cd -

mkdir -p $DATA_DIR/binaries
mkdir -p $DATA_DIR/fofns
mkdir -p $DATA_DIR/cortex_output/dd2_ga01
mkdir -p $DATA_DIR/cortex_output/dd2_ga01/log

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
$SAMTOOLS_EXE bamshuf -uOn 128 /nfs/team112_internal/production_files/Pf/4_0/PFprog1/PG0008_CW/PG0008_CW.bam tmp | \
$SAMTOOLS_EXE bam2fq - > /nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/cortex_exploration/data/fastq/PG0008_CW.fastq

$SAMTOOLS_EXE view -b /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#5.cram |
$SAMTOOLS_EXE bamshuf -uOn 128 - tmp | \
$SAMTOOLS_EXE bam2fq - > /nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/cortex_exploration/data/fastq/16503_1#5.fastq

# Run Cortex
perl $RUN_CALLS_EXE --first_kmer 31 \
--fastaq_index $DATA_DIR/fofns/dd2_ga01.fofn \
--auto_cleaning yes \
--bc yes \
--pd no \
--outdir $DATA_DIR/cortex_output/dd2_ga01 \
--outvcf dd2_ga01.vcf \
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
--logfile $DATA_DIR/cortex_output/dd2_ga01/log/dd2_ga01.log.txt
# --apply_pop_classifier

