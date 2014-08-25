#!/bin/bash

# log all commands to stderr
set -x

# catch errors
set -e
set -o pipefail

# check all variables are set
set -u

# handle command line args
while getopts "c:j:C:v:a:z:b:y:1:d:" opt; do
	case $opt in
	c)
		CHUNKS_MANIFEST=$OPTARG ;;
	j)
		CHUNKED_GENOTYPES_DIR=$OPTARG ;;
	C)
	    CFSCRIPT=$OPTARG ;;
	v)
		ORIGINAL_VCF=$OPTARG ;;
	a)
	    MIN_CODING_COVERAGE=$OPTARG ;;
	z)
	    MAX_CODING_COVERAGE=$OPTARG ;;
	1)
	    REGIONS=$OPTARG ;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1 ;;
	:)
		echo "Option -$OPTARG requires an argument" >&2
		exit 1 ;;
	esac
done


# dereference task ID to chunk strings
CHUNK_SAMTOOLS=`awk "NR==${SGE_TASK_ID}" $CHUNKS_MANIFEST | cut -f1`
CHROMOSOME=`awk "NR==${SGE_TASK_ID}" $CHUNKS_MANIFEST | cut -f3`
START=`awk "NR==${SGE_TASK_ID}" $CHUNKS_MANIFEST | cut -f4`
END=`awk "NR==${SGE_TASK_ID}" $CHUNKS_MANIFEST | cut -f5`

# check output directories exist
mkdir -p $CHUNKED_GENOTYPES_DIR


# determine output file paths
CF_VCF=${CHUNKED_GENOTYPES_DIR}/${CHUNK_SAMTOOLS}.coverageFiltered.vcf
ANNOTATED_VCF=${CHUNKED_GENOTYPES_DIR}/${CHUNK_SAMTOOLS}.annotated.vcf

MD5=$ANNOTATED_VCF.gz.md5

# here we use an md5 checksum file to indicate that the task was
# run successfully to completion - if the md5 file is present and
# the checksum verifies, we can assume the task was previously run
# and can be skipped this time

if ( md5sum --check $MD5 )
then

	echo "skipping ${CHUNK_SAMTOOLS}"

else

    echo "input verified"


    ###########################################################################
    # 1) change coverage filters
    ###########################################################################

    $CFSCRIPT --input_vcf_fn $ORIGINAL_VCF \
        --output_vcf_fn $CF_VCF \
        --min_coding_coverage $MIN_CODING_COVERAGE \
        --max_coding_coverage $MAX_CODING_COVERAGE \
        --chromosome $CHROMOSOME \
        --start $START \
        --end $END

    bgzip -f $CF_VCF

    tabix -f -p vcf $CF_VCF.gz


    ###########################################################################
    # 2) merge annotations to create annotated vcfs
    ###########################################################################

    zcat $CF_VCF.gz \
    | vcf-annotate -a $REGIONS \
                 -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomereHypervariable: region of poor conservation between the 3D7 reference genome and other samples, where the region is adjacent to a telomere. InternalHypervariable: region of poor conservation between the 3D7 reference genome and other samples, where the region is not adjacent to a telomere. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
                 -c CHROM,FROM,TO,INFO/RegionType \
    > $ANNOTATED_VCF
    # | vcf-annotate -a $PRESOM_DEPTH \
    #              -d key=INFO,ID=BWA_DP,Number=1,Type=Integer,Description='Total pre-SNP-o-matic depth' \
    #              -c CHROM,FROM,INFO/BWA_DP \
    
	bgzip -f $ANNOTATED_VCF
	tabix -f -p vcf $ANNOTATED_VCF.gz
	

	# compute MD5
	md5sum $ANNOTATED_VCF.gz > $MD5    

fi

