bash

# Downlaod Cortex from github
cd /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex
git clone https://github.com/iqbal-lab/cortex.git
cd cortex
bash install.sh 2> /dev/null
cd

# Install picard
mkdir -p /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/picard
cd /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/picard
wget https://github.com/broadinstitute/picard/releases/download/1.135/picard-tools-1.135.zip
unzip picard-tools-1.135.zip
cd -

mkdir -p /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/snpeff
cd /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/snpeff
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd -

export DATA_DIR=/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/cortex_exploration/data
# export CORTEX_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21
export CORTEX_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/cortex
export CORTEX_EXE=$CORTEX_DIR/bin/cortex_var_31_c1
export REF_GENOME=$DATA_DIR/roamato/Pf3D7_v3/3D7_sorted.fa
export STAMPY_EXE="/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/stampy/stampy-1.0.27/stampy.py"
export SAMTOOLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/samtools/samtools-1.2/samtools
export GATK_EXE="/software/jre1.7.0_25/bin/java -jar /nfs/team112_internal/production/tools/bin/gatk/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -et NO_ET -K /nfs/team112_internal/production/tools/bin/gatk/dj6_sanger.ac.uk.key"
export PICARD_EXE="/software/jre1.7.0_25/bin/java -jar /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/picard/picard-tools-1.135/picard.jar"
export SNPEFF_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/snpeff/snpEff
export SNPEFF_EXE="/software/jre1.7.0_25/bin/java -jar $SNPEFF_DIR/snpEff.jar"
export VCF_ANNOTATE_EXE=/nfs/team112_internal/production/tools/bin/vcftools_0.1.10/bin/vcf-annotate
# export RUN_CALLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/CORTEX_release_v1.0.5.21/scripts/calling/run_calls.pl
export RUN_CALLS_EXE=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/cortex/cortex/scripts/calling/run_calls.pl
# export VCFTOOLS_DIR=/nfs/team112_internal/production/tools/bin/
export VCFTOOLS_DIR=/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/vcftools/vcftools_0.1.12b
export SAMPLE_NAMES=(7G8 GB4 KE01 KH02 GN01)
number_of_samples=${#SAMPLE_NAMES[@]}
export BAM_CRAM_FNS=(\
/nfs/team112_internal/production_files/Pf/4_0/PFprog3/PG0083_C/PG0083_C.bam \
/nfs/team112_internal/production_files/Pf/4_0/PFprog3/PG0084_C/PG0084_C.bam \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#3.cram \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#2.cram \
/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#6.cram\
)
export FASTQ_DIR=$DATA_DIR/fastq
export REGIONS_FN=~/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz

mkdir -p $SNPEFF_DIR/data/Pf3D7july2015
cd $SNPEFF_DIR/data/Pf3D7july2015
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz
gunzip Pfalciparum.gff3.gz
mv Pfalciparum.gff3 genes.gff
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.genome.fasta.gz
gunzip Pfalciparum.genome.fasta.gz
mv Pfalciparum.genome.fasta sequences.fa
cd
echo "
# Pf3D7 GeneDB July 2015. Added by Richard Pearson. Downloaded July 13, 2015 from ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz 
Pf3D7july2015.genome: Plasmodium_falciparum
        Pf3D7july2015.Pf_M76611.codonTable: Protozoan_Mitochondrial
        Pf3D7july2015.reference: ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz
        Pf3D7july2015.mal_mito_1.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_2.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_3.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_10\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_15\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_16\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_17\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_RNA19\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_1\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_20\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_21\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_9\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUC\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUG\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA11\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA12\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA14\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA18\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA22\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA4\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA5\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA6\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA7\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_SSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_SSUB\:rRNA.codonTable : Vertebrate_Mitochondrial
" >> $SNPEFF_DIR/snpEff.config

cd $SNPEFF_DIR
$SNPEFF_EXE build -gff3 -v Pf3D7july2015
cd

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
for (( i=0; i<${number_of_samples}; i++ ));
do
    printf "%s\n" "$DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/vcfs/${SAMPLE_NAMES[$i]}.vcf_union_BC_calls_k31.decomp.vcf" >> $DATA_DIR/fofns/all_decomp_vcfs.fofn
done

for (( i=0; i<${number_of_samples}; i++ ));
do
    printf "%s\n" "$DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/vcfs/${SAMPLE_NAMES[$i]}.vcf_union_BC_calls_k31.raw.vcf" >> $DATA_DIR/fofns/all_raw_vcfs.fofn
done



# Combine VCFs
mkdir -p $DATA_DIR/cortex_output/combined

perl $CORTEX_DIR/scripts/analyse_variants/combine/combine_vcfs.pl \
$DATA_DIR/fofns/all_raw_vcfs.fofn \
$CORTEX_DIR \
$VCFTOOLS_DIR \
$CORTEX_DIR/scripts/analyse_variants/bioinf-perl/vcf_scripts \
$CORTEX_DIR/scripts/analyse_variants \
$CORTEX_DIR/scripts/analyse_variants/combine \
$DATA_DIR/cortex_output/combined \
combined_raw \
31 \
Pf3D7_v3 \
$REF_GENOME \
$DATA_DIR/binaries/ref.k31.ctx \
20 \
150 \
$DATA_DIR/cortex_output



# genotype each samples against combined sites
for (( i=0; i<${number_of_samples}; i++ ));
do
    echo ${SAMPLE_NAMES[$i]}
    mkdir -p $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/genotyped_vcf
    sample_graph_fn=`ls $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/binaries/cleaned/k31/*.ctx`
    perl $CORTEX_DIR/scripts/calling/genotype_1sample_against_sites.pl \
        --cortex_dir $CORTEX_DIR \
        --invcf $DATA_DIR/cortex_output/combined/combined_raw.sites_vcf \
        --bubble_graph $DATA_DIR/cortex_output/combined/combined_raw.pseudo_callfile.branches.k31.ctx \
        --bubble_callfile  $DATA_DIR/cortex_output/combined/combined_raw.pseudo_callfile \
        --outdir $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/genotyped_vcf \
        --sample ${SAMPLE_NAMES[$i]} \
        --ref_overlap_bubble_graph $DATA_DIR/cortex_output/combined/ref.k31.ctx.list_intersect_bubbles.ctx \
        --genome_size 24000000 \
        --sample_graph $sample_graph_fn \
        --overlap_log $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/genotyped_vcf/log.txt \
        --mem_height 18 \
        --mem_width 100
done
# --sample_graph $DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/binaries/cleaned/k31/${SAMPLE_NAMES[$i]}.kmer31.q15cleaned_7.ctx \

# Create VCFs fofn
for (( i=0; i<${number_of_samples}; i++ ));
do
    printf "%s\n" "$DATA_DIR/cortex_output/${SAMPLE_NAMES[$i]}/genotyped_vcf/${SAMPLE_NAMES[$i]}.vcf" >> $DATA_DIR/fofns/all_genotyped_vcfs.list
done


$GATK_EXE \
    -T CombineVariants \
    -R $REF_GENOME \
    --variant $DATA_DIR/fofns/all_genotyped_vcfs.list \
    -o $DATA_DIR/cortex_output/combined/combined_genotyped.vcf \
    --unsafe

# Install picard


$PICARD_EXE SortVcf \
    INPUT=$DATA_DIR/cortex_output/combined/combined_genotyped.vcf \
    OUTPUT=$DATA_DIR/cortex_output/combined/combined_genotyped_sorted.vcf \
    SEQUENCE_DICTIONARY=$DATA_DIR/roamato/Pf3D7_v3/3D7_sorted.dict

# Annotate with SnpEff
$SNPEFF_EXE \
-v -o gatk Pf3D7july2015 \
$DATA_DIR/cortex_output/combined/combined_genotyped.vcf \
-no-downstream \
-no-upstream \
> $DATA_DIR/cortex_output/combined/combined_genotyped.snpeff.vcf \
2> /dev/null

$GATK_EXE \
-T VariantAnnotator \
-R $REF_GENOME \
-A SnpEff \
--variant $DATA_DIR/cortex_output/combined/combined_genotyped.vcf  \
--snpEffFile $DATA_DIR/cortex_output/combined/combined_genotyped.snpeff.vcf \
-o $DATA_DIR/cortex_output/combined/combined_genotyped.snpeff_annotated.vcf \
--unsafe \
2> /dev/null

cat $DATA_DIR/cortex_output/combined/combined_genotyped.snpeff_annotated.vcf \
| $VCF_ANNOTATE_EXE -a $REGIONS_FN \
   -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
   -c CHROM,FROM,TO,INFO/RegionType \
> $DATA_DIR/cortex_output/combined/combined_genotyped.annotated.vcf
