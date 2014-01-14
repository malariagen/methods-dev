# To kick off the full pipeline, cd to base directory (e.g. /data/plasmodium/pfalciparum/3_0/methods-dev), then:
# ./scripts/cluster/setups/pf_3_0_annotate_25_99.sh

export SETUPNAME=pf_3_0_25_99 #n

export ORIGINAL_VCF=/data/plasmodium/pfalciparum/methods_dev/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.gz #v - note that this is not the same as the version at the Sanger, which was zipped with gzip and not bgzip
export REGIONS=/data/plasmodium/pfalciparum/methods_dev/Pf/misc/regions-20130225.onebased.txt.gz #1
export PRESOM_DEPTH=/data/plasmodium/pfalciparum/methods_dev/Pf/misc/joined.depth.gz #d

export CHUNK_SIZE=10000 #h
export CHUNKS_MANIFEST=./meta/${SETUPNAME}_regions_chunks_${CHUNK_SIZE}bp.tab #c - note this will get created as part of pipeline

# export MIN_CODING_COVERAGE=209186 #a
# export MAX_CODING_COVERAGE=322716 #z
export MIN_CODING_COVERAGE=158141 #a
export MAX_CODING_COVERAGE=358381 #z
export MIN_NONCODING_COVERAGE=62150 #b
export MAX_NONCODING_COVERAGE=233117 #y

export CHUNKED_GENOTYPES_DIR=data/${SETUPNAME}/chunked_genotypes/snps #j
mkdir -p ${CHUNKED_GENOTYPES_DIR}
export MERGED_GENOTYPES_DIR=data/${SETUPNAME}/merged_genotypes/snps #6
mkdir -p ${MERGED_GENOTYPES_DIR}
export ANNOTATED_VCF_SUFFIX=.annotated.vcf.gz #4
export ANNOTATED_VCF_FOFN=${MERGED_GENOTYPES_DIR}/${SETUPNAME}_annotated_vcfs.list #5 - note this will get created as part of pipeline
export FINAL_GENOTYPES_DIR=data/${SETUPNAME}/final_genotypes/snps #7
mkdir -p ${FINAL_GENOTYPES_DIR}
export FINAL_GENOTYPED_VCF_FN=${FINAL_GENOTYPES_DIR}/${SETUPNAME}_merged_hetuniq_newbiallelic_20130809.newCoverageFilters_25_99.vcf #L

export CFSCRIPT=./scripts/python/coverageFilters.py #C
export REGIONCHUNKS=./scripts/python/regionChunks.py #R
export VCFCONCAT=vcf-concat #C

export GENOME=/data/mirror/nfs/team112_internal/production_files/Pf/3_0/3D7_V3.fasta #g

scripts/cluster/pipelines/pf_annotate.sh

