import vcfnp
import collections

final_vcfs = collections.OrderedDict()
final_vcfs['Pf3D7_01_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/f/d/2/8/336049/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_02_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/b/b/c/e/336050/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_02_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_03_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/3/f/0/6/336051/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_03_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_04_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/2/9/4/0/336052/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_05_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/2/9/e/5/336053/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_05_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_06_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/3/7/b/6/336054/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_06_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_07_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/5/c/c/8/336055/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_07_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_08_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/6/9/1/9/336056/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_08_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_09_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/f/f/a/3/336057/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_09_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_10_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/d/b/5/0/336058/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_10_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_11_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/e/4/a/7/336059/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_11_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_12_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/f/a/0/3/336060/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_12_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_13_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/7/5/6/d/336061/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_13_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_14_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/b/c/4/5/336062/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_14_v3.combined.filtered.vcf.gz'
final_vcfs['Pf3D7_API_v3'] = '/lustre/scratch109/malaria/pf3k_methods/output/9/7/b/4/336063/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_API_v3.combined.filtered.vcf.gz'
final_vcfs['Pf_M76611'] = '/lustre/scratch109/malaria/pf3k_methods/output/9/4/3/9/336064/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf_M76611.combined.filtered.vcf.gz'

output_dir = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate'

def run_vcfnp(chromosome='Pf3D7_01_v3'):
    vcf_fn = final_vcfs[chromosome]
    v = vcfnp.variants(
        vcf_fn=vcf_fn,
        count=100,
        progress=10,
        cache=True,
        cachedir=output_dir+'/vcfnp'
    )
    return(v)

run_vcfnp()

