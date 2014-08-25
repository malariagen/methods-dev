# To kick off, cd to base directory (e.g. ~/src/github/malariagen/methods-dev/), then:
# ./scripts/cluster/setups/setup_build_numpy_arrays_ppq_pf_4_0_common_samples.sh

export SETUPNAME=build_numpy_arrays_ppq_pf_4_0_common_samples
export VCF_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/ppq.40.vfp1.newCoverageFilters_pass_25_99.5.HyperHet.vcf.gz
export SAMPLES_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/meta/common_samples_ppq_4_0.txt
export CALLDATA_FILE_SUFFIX=.calldata_common_samples.npy

./scripts/cluster/pipelines/pf_build_numpy_arrays_4_0.sh
