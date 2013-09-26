# To kick off, cd to base directory (e.g. /data/plasmodium/pfalciparum/3_0/methods-dev/), then:
# . ./scripts/cluster/setups/pf_3_0_build_numpy_arrays.sh

export SETUPNAME=3_0
export VCF_FN=/data/plasmodium/pfalciparum/methods_dev/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.bgz
export SAMPLES_FN=./meta/evaluation_samples.txt

. ./scripts/cluster/pipelines/pf_build_numpy_arrays.sh
