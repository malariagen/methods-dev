#!/bin/bash

# log all commands to stderr
set -x

# catch errors
set -e
set -o pipefail

# check all variables are set
set -u

# handle command line args
while getopts "n:7:V:H:C:" opt; do
    case $opt in
    n)
        SETUPNAME=$OPTARG ;;
    7)
        FINAL_GENOTYPES_DIR=$OPTARG ;;
    V)
        FILTERED_VCF_FOFN=$OPTARG ;;
    H)
        FINAL_VCF=$OPTARG ;;
    C)
        VCFCONCAT=$OPTARG ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1 ;;
    :)
        echo "Option -$OPTARG requires an argument" >&2
        exit 1 ;;
    esac
done


# determine file paths
# FINAL_VCF=$FINAL_GENOTYPES_DIR/${SETUPNAME}_final.vcf
MD5=${FINAL_VCF}.gz.md5

# here we use an md5 checksum file to indicate that the task was
# run successfully to completion - if the md5 file is present and
# the checksum verifies, we can assume the task was previously run
# and can be skipped this time

if ( md5sum --check $MD5 )
then

    echo "skipping"

else
    
    echo "input verified"

    $VCFCONCAT -f $FILTERED_VCF_FOFN > $FINAL_VCF
    
    bgzip -f $FINAL_VCF
    
    tabix -f -p vcf $FINAL_VCF.gz

    # compute MD5
    md5sum $FINAL_VCF.gz > $MD5    

fi
