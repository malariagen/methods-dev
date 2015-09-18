rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/variant_annotation_using_snpeff_gatk_vcf_annotate_snps \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150904_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150907_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants2 \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  >> ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150907_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants2 \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150908_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/assembled_samples_2/vcfs/vcf \
  /lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/assembled_samples_2/vcfs \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150909_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch109/malaria/pfalciparum/resources \
  /lustre/scratch109/malaria/pfalciparum \
  >> ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150909_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/resources \
  /lustre/scratch110/malaria/dj6/pf3k \
  >> ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150909_rsync_2.log

rsync -avL /nfs/team112_internal/production_files/Pf3k/methods/assembled_samples/truth_vcfs \
  rp7@malsrv2:/nfs/team112_internal/rp7/Pf3k/methods/assembled_samples \
  >> ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150910_rsync.log


rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/assembled_samples_2/vcfs/vcf \
  /lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/assembled_samples_2/vcfs \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150910_rsync.log
