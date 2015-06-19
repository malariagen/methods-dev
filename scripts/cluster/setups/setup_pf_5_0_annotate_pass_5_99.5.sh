# To kick off the full pipeline, cd to base directory (e.g. ~/src/github/malariagen/methods-dev), then:
# ./scripts/cluster/setups/setup_pf_5_0_annotate_pass_5_99.5.sh

export SETUPNAME=pf_5_0_pass_5_99.5 #n

export ORIGINAL_VCF=/data/mirror/nfs/team112_internal/production_files/Pv/5_0/pf_50_vfp1.vcf.gz
# export ORIGINAL_VCF=/data/plasmodium/pfalciparum/pf_community/data/4_0/vcf/pf_4_0_20140712_vfp1.vcf.gz #v
export REGIONS=/data/plasmodium/pfalciparum/pf_community/data/3_0/pre_release/misc/regions-20130225.onebased.txt.gz #1
# export PRESOM_DEPTH=/data/plasmodium/pfalciparum/methods_dev/Pf/misc/joined.depth.gz #d

export CHUNK_SIZE=10000 #h
export CHUNKS_MANIFEST=./meta/${SETUPNAME}_regions_chunks_${CHUNK_SIZE}bp.tab #c - note this will get created as part of pipeline

# export MIN_CODING_COVERAGE=209186 #a
# export MAX_CODING_COVERAGE=322716 #z
# export MIN_CODING_COVERAGE=158141 #a
# export MAX_CODING_COVERAGE=358381 #z
export MIN_CODING_COVERAGE=190053 #a
export MAX_CODING_COVERAGE=711257 #z
# export MIN_NONCODING_COVERAGE=62150 #b
# export MAX_NONCODING_COVERAGE=233117 #y

export CHUNKED_GENOTYPES_DIR=data/${SETUPNAME}/chunked_genotypes/snps #j
mkdir -p ${CHUNKED_GENOTYPES_DIR}
export MERGED_GENOTYPES_DIR=data/${SETUPNAME}/merged_genotypes/snps #6
mkdir -p ${MERGED_GENOTYPES_DIR}
export ANNOTATED_VCF_SUFFIX=.annotated.vcf.gz #4
export ANNOTATED_VCF_FOFN=${MERGED_GENOTYPES_DIR}/${SETUPNAME}_annotated_vcfs.list #5 - note this will get created as part of pipeline
export FINAL_GENOTYPES_DIR=data/${SETUPNAME}/final_genotypes/snps #7
mkdir -p ${FINAL_GENOTYPES_DIR}
export FINAL_GENOTYPED_VCF_FN=${FINAL_GENOTYPES_DIR}/pf_50_vfp1.newCoverageFilters_pass_5_99.5.vcf #L

export CFSCRIPT=./scripts/python/coverageFilters_4_0.py #C
export REGIONCHUNKS=./scripts/python/regionChunks.py #R
export VCFCONCAT=vcf-concat #C

export GENOME=/data/mirror/nfs/team112_internal/production_files/Pf/3_0/3D7_V3.fasta #g

scripts/cluster/pipelines/pipeline_pf_annotate_4_0.sh

