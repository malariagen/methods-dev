
# coding: utf-8

# # Setup

# In[2]:

get_ipython().magic('run _shared_setup.ipynb')


# In[3]:

GATK_EXE="/software/jre1.7.0_25/bin/java -jar /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar"
REF_GENOME="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta"
FINAL_VCF_FORMAT="/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/SNP_INDEL_%s.combined.filtered.vcf.gz"
OUTGROUP_VCF_FORMAT="/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/with_outgroup_alleles/SNP_INDEL_%s.combined.filtered.vcf.gz"
SAMPLE_FILE="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_samples.txt"
VALIDATION_SAMPLES_VCF_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/%s.5validation.vcf.gz"
INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/%s.%s.vcf.gz"
FULL_NPY_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/with_outgroup_%s"


# In[4]:

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))
ref_dict.keys()


# In[5]:

validation_vcfs = collections.OrderedDict()
validation_vcfs['Pf3D7_01_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/e/b/f/a/15179/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_02_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/4/8/3/9/15180/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_02_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_03_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/4/2/4/0/15181/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_03_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_04_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/5/c/c/9/15182/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_05_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/6/1/c/9/15183/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_05_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_06_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/f/4/2/9/15184/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_06_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_07_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/7/c/e/d/15185/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_07_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_08_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/0/5/1/0/15186/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_08_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_09_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/c/4/c/f/15187/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_09_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_10_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/a/2/4/2/15188/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_10_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_11_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/e/1/d/15189/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_11_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_12_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/6/6/c/8/15190/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_12_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_13_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/b/6/7/4/15191/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_13_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_14_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/e/2/d/15192/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_14_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf3D7_API_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/0/3/b/e/15193/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_API_v3.combined.filtered.vcf.gz'
validation_vcfs['Pf_M76611'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/7/6/3/15194/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf_M76611.combined.filtered.vcf.gz'


# # Process data

# In[ ]:

for chrom in ref_dict.keys():
    print(chrom)
    get_ipython().system('{GATK_EXE}         -T SelectVariants         -R {REF_GENOME}         --variant {FINAL_VCF_FORMAT % chrom}         --out /dev/stdout         --sample_file {SAMPLE_FILE}         --excludeNonVariants         | bgzip -c > {VALIDATION_SAMPLES_VCF_FORMAT % chrom}')
        
    get_ipython().system('tabix -p vcf {VALIDATION_SAMPLES_VCF_FORMAT % chrom}')

    


# In[36]:

# !{GATK_EXE} \
#     org.broadinstitute.gatk.tools.CatVariants \

# !/software/jre1.7.0_25/bin/java -cp \
#     /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar \
#     org.broadinstitute.gatk.tools.CatVariants \

get_ipython().system("{GATK_EXE}     -T CombineVariants     -R {REF_GENOME}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_01_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_02_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_03_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_04_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_05_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_06_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_07_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_08_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_09_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_10_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_11_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_12_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_13_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_14_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf3D7_API_v3'}     -V {VALIDATION_SAMPLES_VCF_FORMAT % 'Pf_M76611'}     -genotypeMergeOptions UNSORTED     --out /dev/stdout     | bgzip -c > {VALIDATION_SAMPLES_VCF_FORMAT % 'WG'}")

get_ipython().system("tabix -p vcf {VALIDATION_SAMPLES_VCF_FORMAT % 'WG'}")

#     --excludeNonVariants \
#     --removeUnusedAlternates \


# In[5]:

get_ipython().system("{GATK_EXE}     -T CombineVariants     -R {REF_GENOME}     -V {validation_vcfs['Pf3D7_01_v3']}     -V {validation_vcfs['Pf3D7_02_v3']}     -V {validation_vcfs['Pf3D7_03_v3']}     -V {validation_vcfs['Pf3D7_04_v3']}     -V {validation_vcfs['Pf3D7_05_v3']}     -V {validation_vcfs['Pf3D7_06_v3']}     -V {validation_vcfs['Pf3D7_07_v3']}     -V {validation_vcfs['Pf3D7_08_v3']}     -V {validation_vcfs['Pf3D7_09_v3']}     -V {validation_vcfs['Pf3D7_10_v3']}     -V {validation_vcfs['Pf3D7_11_v3']}     -V {validation_vcfs['Pf3D7_12_v3']}     -V {validation_vcfs['Pf3D7_13_v3']}     -V {validation_vcfs['Pf3D7_14_v3']}     -V {validation_vcfs['Pf3D7_API_v3']}     -V {validation_vcfs['Pf_M76611']}     -genotypeMergeOptions UNSORTED     --out /dev/stdout     | bgzip -c > {VALIDATION_SAMPLES_VCF_FORMAT % 'validation_WG'}")

get_ipython().system("tabix -p vcf {VALIDATION_SAMPLES_VCF_FORMAT % 'validation_WG'}")


# In[32]:

get_ipython().system("{GATK_EXE}     -T SelectVariants     -R ${REF_GENOME}     --variant {VALIDATION_SAMPLES_VCF_FORMAT % 'WG'}     --sample_name {sample}     --excludeNonVariants     --removeUnusedAlternates     --out /dev/stdout     | bgzip -c > {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % (chrom, sample)}")

    get_ipython().system('tabix -p vcf {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % (chrom, sample)}')



# In[6]:

for sample in ['7G8', 'ERS740936', 'ERS740937', 'ERS740940', 'GB4']:
    print(sample)
    get_ipython().system("{GATK_EXE}         -T SelectVariants         -R {REF_GENOME}         --variant {VALIDATION_SAMPLES_VCF_FORMAT % 'WG'}         --sample_name {sample}         --excludeNonVariants         --removeUnusedAlternates         --out /dev/stdout         | bgzip -c > {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % ('WG', sample)}")

    get_ipython().system("tabix -p vcf {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % ('WG', sample)}")

    


# In[6]:

for sample in ['7G8', 'GB4', 'GN01', 'KE01', 'KH02']:
    print(sample)
    get_ipython().system("{GATK_EXE}         -T SelectVariants         -R {REF_GENOME}         --variant {VALIDATION_SAMPLES_VCF_FORMAT % 'validation_WG'}         --sample_name {sample}         --excludeNonVariants         --removeUnusedAlternates         --out /dev/stdout         | bgzip -c > {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % ('validation_WG', sample)}")

    get_ipython().system("tabix -p vcf {INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % ('validation_WG', sample)}")

    


# In[8]:

import vcfnp

def run_vcfnp(sample='7G8', file_prefix='WG'):
    vcf_fn = INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT % (file_prefix, sample)
    v = vcfnp.variants(
        vcf_fn=vcf_fn,
        progress=10000,
        arities={
            'ALT': 2,
            'AF': 2,
            'AC': 2,
            'MLEAF': 2,
            'MLEAC': 2,
            'RPA': 3
        },
        dtypes={
            'REF': 'a400', 
            'ALT': 'a400',
            'RegionType': 'a25', 
            'VariantType': 'a40',
            'RU': 'a40',
            'set': 'a40',
            'SNPEFF_AMINO_ACID_CHANGE':'a20',
            'SNPEFF_CODON_CHANGE':'a20',
            'SNPEFF_EFFECT':'a33',
            'SNPEFF_EXON_ID':'a2',
            'SNPEFF_FUNCTIONAL_CLASS':'a8',
            'SNPEFF_GENE_BIOTYPE':'a14',
            'SNPEFF_GENE_NAME':'a20',
            'SNPEFF_IMPACT':'a8',
            'SNPEFF_TRANSCRIPT_ID':'a20',
            'VariantType':'a60',
            'culprit':'a14',
        },
        cache=True
    )
    c = vcfnp.calldata_2d(
        vcf_fn=vcf_fn,
        progress=10000,
        fields=['GT', 'AD'],
        arities={'AD': 3},
        cache=True,
    )
    print(sample, max(v['num_alleles']), max([len(x) for x in v['REF']]))
    return(v, c)


# In[9]:

for sample in ['7G8', 'ERS740936', 'ERS740937', 'ERS740940', 'GB4']:
    print(sample)
    _ = run_vcfnp(sample)


# In[8]:

for sample in ['7G8', 'GB4', 'GN01', 'KE01', 'KH02']:
    print(sample)
    _ = run_vcfnp(sample, 'validation_WG')


# In[11]:

def run_vcfnp_all(chrom='Pf3D7_01_v3'):
    vcf_fn = OUTGROUP_VCF_FORMAT % chrom
    v = vcfnp.variants(
        vcf_fn=vcf_fn,
        progress=10000,
        arities={
            'ALT': 6,
            'AF': 6,
            'AC': 6,
            'MLEAF': 6,
            'MLEAC': 6,
            'RPA': 7
        },
        dtypes={
#             'REF': 'a100', 
#             'ALT': 'a100',
            'RegionType': 'a25', 
            'VariantType': 'a40',
            'RU': 'a40',
            'set': 'a40',
            'SNPEFF_AMINO_ACID_CHANGE':'a20',
            'SNPEFF_CODON_CHANGE':'a20',
            'SNPEFF_EFFECT':'a33',
            'SNPEFF_EXON_ID':'a2',
            'SNPEFF_FUNCTIONAL_CLASS':'a8',
            'SNPEFF_GENE_BIOTYPE':'a14',
            'SNPEFF_GENE_NAME':'a20',
            'SNPEFF_IMPACT':'a8',
            'SNPEFF_TRANSCRIPT_ID':'a20',
            'VariantType':'a60',
            'culprit':'a14',
        },
        cache=True,
        cachedir=FULL_NPY_FORMAT % chrom
    )
    return(0)


# In[12]:

for chrom in ref_dict.keys():
    print(chrom)
    run_vcfnp_all(chrom)


# In[10]:

def run_vcfnp_calls_all(chrom='Pf3D7_01_v3'):
    vcf_fn = OUTGROUP_VCF_FORMAT % chrom
    c = vcfnp.calldata_2d(
        vcf_fn=vcf_fn,
        progress=10000,
#         fields=['GT', 'AD'],
        fields=['AD'],
        arities={'AD': 7},
        cache=True,
        cachedir=FULL_NPY_FORMAT % chrom
    )
    return(0)


# In[ ]:

for chrom in ref_dict.keys():
    print(chrom)
    run_vcfnp_calls_all(chrom)


# In[ ]:



