#!/bin/bash

# 
# One step pipeline to call genotypes.
# 

# bail out on first error
set -e

# make sure variables are set
set -u

# debug
set -x

# external variables expected (source these from a setup)
echo "setup name: $SETUPNAME"
echo "input vcf filename: $VCF_FN"
echo "input samples filename: $SAMPLES_FN"
echo "calldata file suffix: $CALLDATA_FILE_SUFFIX"

##########
# step 1 #
##########

# job name
STEP1_JOBNAME="pf_build_numpy_arrays_$SETUPNAME"

# one task per chromosome
STEP1_NUMTASKS=14

# submit job
qsub -S /usr/bin/python \
-N $STEP1_JOBNAME \
-cwd \
-t 1-$STEP1_NUMTASKS \
-l h_vmem=4000m \
-o log \
-j y \
./scripts/cluster/steps/build_numpy_arrays_4_0.py -v $VCF_FN -s $SAMPLES_FN -c $CALLDATA_FILE_SUFFIX

# end of pipeline
