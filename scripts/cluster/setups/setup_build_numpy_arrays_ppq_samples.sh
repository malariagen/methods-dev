# To kick off, cd to base directory (e.g. ~/src/github/malariagen/methods-dev/), then:
# ./scripts/cluster/setups/setup_build_numpy_arrays_ppq_samples.sh

export SETUPNAME=build_numpy_arrays_ppq_samples
export VCF_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/ppq.40.vfp1.newCoverageFilters_pass_25_99.5.HyperHet.vcf.gz
export SAMPLES_FN=/data/plasmodium/pfalciparum/pf_community/data/ppq/meta/ppq_samples.txt
export CALLDATA_FILE_SUFFIX=.calldata_ppq_samples.npy

./scripts/cluster/pipelines/pf_build_numpy_arrays_4_0.sh
