rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/variant_annotation_using_snpeff_gatk_vcf_annotate_snps \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150904_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  > ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150907_rsync.log

rsync -avL rp7@malsrv2:/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants2 \
  /lustre/scratch110/malaria/dj6/pf3k/output \
  >> ~/src/github/malariagen/methods-dev/pf3k_techbm/log/20150907_rsync.log

