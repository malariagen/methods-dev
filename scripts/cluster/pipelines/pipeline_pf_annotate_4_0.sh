#!/bin/bash

# 
# pipeline to align fastq with bwa
# 

# bail out on first error
set -e

# make sure variables are set
set -u

# debug
set -x

# external variables expected (source these from a setup)
echo "setup name: $SETUPNAME"

##################################
# step 0 - create manifest files #
##################################

# job name
STEP0_JOBNAME="pf_manifests_$SETUPNAME"

# submit job
qsub -S /bin/bash \
-N $STEP0_JOBNAME \
-cwd \
-l h_vmem=1000M \
-o log \
-j y \
./scripts/cluster/steps/manifests.sh \
    -g $GENOME \
    -h $CHUNK_SIZE \
    -c $CHUNKS_MANIFEST \
    -j $CHUNKED_GENOTYPES_DIR \
    -4 $ANNOTATED_VCF_SUFFIX \
    -5 $ANNOTATED_VCF_FOFN \
    -R $REGIONCHUNKS


##############################################
# step 1 - set coverage filters and annotate #
##############################################

# job name
STEP1_JOBNAME="pf_annotate_and_set_coverage_filters_$SETUPNAME"

# how many tasks?
STEP1_NUMTASKS=`cat $CHUNKS_MANIFEST | wc -l`

# submit job
qsub -S /bin/bash \
-N $STEP1_JOBNAME \
-hold_jid $STEP0_JOBNAME \
-cwd \
-t 1-$STEP1_NUMTASKS \
-l h_vmem=1000M \
-o log \
-j y \
./scripts/cluster/steps/step_annotate_and_set_coverage_filters_4_0.sh \
    -c $CHUNKS_MANIFEST \
    -j $CHUNKED_GENOTYPES_DIR \
    -C $CFSCRIPT \
    -v $ORIGINAL_VCF \
    -a $MIN_CODING_COVERAGE \
    -z $MAX_CODING_COVERAGE \
    -1 $REGIONS
    # -b $MIN_NONCODING_COVERAGE \
    # -y $MAX_NONCODING_COVERAGE \
    # -d $PRESOM_DEPTH

#######################################
# step 2 - merge all genotyped chunks #
#######################################

# job name
STEP2_JOBNAME="pf_merge_final_genotyped_$SETUPNAME"

# submit job
qsub -S /bin/bash \
-N $STEP2_JOBNAME \
-hold_jid $STEP1_JOBNAME \
-cwd \
-l h_vmem=4000m \
-o log \
-j y \
./scripts/cluster/steps/merge_vcf_fofn.sh \
    -n $SETUPNAME \
    -7 $FINAL_GENOTYPES_DIR \
    -V $ANNOTATED_VCF_FOFN \
    -H $FINAL_GENOTYPED_VCF_FN \
    -C $VCFCONCAT
