# To kick off, cd to base directory (e.g. ~/src/github/malariagen/methods-dev/), then:
# ./scripts/cluster/setups/pf_ppq_build_numpy_arrays_25_samples.sh

export SETUPNAME=ppq
export VCF_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/ppq.40.vfp1.vcf.gz
export SAMPLES_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/meta/evaluation_samples_25.txt
export CALLDATA_FILE_SUFFIX=.calldata_25samples.npy

./scripts/cluster/pipelines/pf_build_numpy_arrays_4_0.sh
