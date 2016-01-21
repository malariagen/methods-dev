#!/bin/bash


export ORIGINAL_DIR=`pwd`
export PROCESSED_DATA_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/exploration"
if [ ! -d "$PROCESSED_DATA_DIR" ]; then
    mkdir -p "$PROCESSED_DATA_DIR"
fi
export OPT_DIR="$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4"

export REF_GENOME_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015"
export REF_GENOME="${REF_GENOME_DIR}/Pfalciparum.genome.fasta"

export VCFS_DIR="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/validation_results/bwa_mem/gatk_rec/HaplotypeCaller/alistair_ann/SnpEff_region/_final_vcfs"

export MUMMER_VERSION="3.23"

export NUCMER_EXE="${OPT_DIR}/mummer/MUMmer${MUMMER_VERSION}/nucmer"
export VCF_CONSENSUS_EXE="/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt/vcftools/vcftools_0.1.12b/bin/vcf-consensus"
export NUCMER_SCRIPT="/nfs/users/nfs_t/tdo/Bin/little.donucmer.sh"



# Install mummer
if [ ! -s ${NUCMER_EXE} ]; then
    mkdir -p ${OPT_DIR}/mummer
    cd ${OPT_DIR}/mummer
    wget http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz
    tar -xvzf MUMmer${MUMMER_VERSION}.tar.gz
    cd MUMmer${MUMMER_VERSION}
    make check
    make install
    cd ${ORIGINAL_DIR}
fi

# Copy Thomas's script
mkdir -p ${OPT_DIR}/tdo
cp ${NUCMER_SCRIPT} ${OPT_DIR}/tdo
export NUCMER_SCRIPT="${OPT_DIR}/tdo/little.donucmer.sh"


samtools faidx $REF_GENOME Pf3D7_01_v3 > $PROCESSED_DATA_DIR/3D7_Pf3D7_01_v3.fa
samtools faidx $REF_GENOME Pf3D7_01_v3 | $VCF_CONSENSUS_EXE --sample PG0083-C $VCFS_DIR/final.SNP_INDEL.Pf3D7_01_v3.vcf.gz > $PROCESSED_DATA_DIR/GATK_7G8_Pf3D7_01_v3.fa
samtools faidx $REF_GENOME Pf3D7_01_v3 | $VCF_CONSENSUS_EXE --sample PG0084-C $VCFS_DIR/final.SNP_INDEL.Pf3D7_01_v3.vcf.gz > $PROCESSED_DATA_DIR/GATK_GB4_Pf3D7_01_v3.fa
cd $PROCESSED_DATA_DIR
export ref3D7=$REF_GENOME
$NUCMER_SCRIPT 3D7 3D7_Pf3D7_01_v3.fa /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta
$NUCMER_SCRIPT GATK_7G8 GATK_7G8_Pf3D7_01_v3.fa /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta
$NUCMER_SCRIPT GATK_GB4 GATK_GB4_Pf3D7_01_v3.fa /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta

cd ..
mkdir nucmer
cd nucmer
$NUCMER_SCRIPT 7G8 /lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/fasta/7G8.fasta /nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta


