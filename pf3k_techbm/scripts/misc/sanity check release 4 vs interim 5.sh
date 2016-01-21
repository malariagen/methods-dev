cat /lustre/scratch109/malaria/pf3k_methods/input/output_fofn/pf3kgatk_variant_filtration_ps583for_2640samples.output

VCF_ORIGINAL=/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz
VCF_NEW=/lustre/scratch109/malaria/pf3k_methods/output/1/f/5/c/394734/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz

zgrep PASS $VCF_ORIGINAL | less -S
zgrep PASS $VCF_NEW | less -S

tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep PASS | less -S
tabix $VCF_NEW Pf3D7_01_v3:92914-93374 | grep PASS | less -S
tabix $VCF_NEW Pf3D7_01_v3:92914-93374 | less -S

/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS" && VQSLOD>6.0' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_ORIGINAL | wc -l

/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS" && VQSLOD>6.0' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_NEW | wc -l


/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS"' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_ORIGINAL | wc -l

/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS"' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_NEW | wc -l

VCF_API_OLD=/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/SNP_INDEL_Pf3D7_API_v3.combined.filtered.vcf.gz
VCF_API_NEW=/lustre/scratch109/malaria/pf3k_methods/output/3/b/1/d/394690/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_API_v3.combined.filtered.vcf.gz

/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS"' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_API_OLD | wc -l

/nfs/team112_internal/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/bcftools/bcftools view \
--include 'FILTER="PASS"' \
--min-alleles 2 \
--max-alleles 2 \
--types snps $VCF_API_NEW | wc -l

l -L $VCF_ORIGINAL
l -L $VCF_NEW
l -L $VCF_API_OLD
l -L $VCF_API_NEW


# Check metrics for two variants, 1 in training set, and one not
zgrep PASS /nfs/team112_internal/oxford_mirror/data/plasmodium/pfalciparum/pf-crosses/data/public/1.0/7g8_gb4.combined.final.vcf.gz | head -n 1

tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep 93157 | head -n 1 | cut -f 1-10
tabix $VCF_NEW Pf3D7_01_v3:92914-93374 | grep 93157 | head -n 1 | cut -f 1-10

tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep 92914 | head -n 1 | cut -f 1-10
tabix $VCF_NEW Pf3D7_01_v3:92914-93374 | grep 92914 | head -n 1 | cut -f 1-10

# Why do some variants not have ReadPosRankSum in original VCF? Because there are no het calls!
tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep ReadPosRankSum | less -S
tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep ReadPosRankSum | grep '0/1' | less -S
tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep ReadPosRankSum | grep '0/2' | less -S
tabix $VCF_ORIGINAL Pf3D7_01_v3:92914-93374 | grep ReadPosRankSum | grep -v '0/1' | less -S

