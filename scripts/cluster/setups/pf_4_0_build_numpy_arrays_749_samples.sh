# To kick off, cd to base directory (e.g. ~/src/github/malariagen/methods-dev/), then:
# ./scripts/cluster/setups/pf_4_0_build_numpy_arrays_749_samples.sh

export SETUPNAME=4_0
export VCF_FN=/data/plasmodium/pfalciparum/pf_community/data/4_0/vcf/pf_4_0_20140712_vfp1.vcf.gz
export SAMPLES_FN=/data/plasmodium/pfalciparum/pf_community/data/4_0/meta/evaluation_samples_749.txt
export CALLDATA_FILE_SUFFIX=.calldata_749samples.npy

./scripts/cluster/pipelines/pf_build_numpy_arrays_4_0.sh
