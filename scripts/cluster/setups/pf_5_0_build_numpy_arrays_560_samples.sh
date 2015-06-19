# To kick off, cd to base directory (e.g. ~/src/github/malariagen/methods-dev/), then:
# ./scripts/cluster/setups/pf_5_0_build_numpy_arrays_560_samples.sh

export SETUPNAME=5_0
export VCF_FN=/data/mirror/nfs/team112_internal/production_files/Pv/5_0/pf_50_vfp1.vcf.gz
export SAMPLES_FN=/data/plasmodium/pfalciparum/pf_community/data/5_0/meta/evaluation_samples_560.txt
export CALLDATA_FILE_SUFFIX=.calldata_560samples.npy

./scripts/cluster/pipelines/pf_build_numpy_arrays_4_0.sh
