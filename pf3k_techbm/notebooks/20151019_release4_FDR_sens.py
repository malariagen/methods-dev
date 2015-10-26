
# coding: utf-8

# In this workbook, three different methods are used for determining FDR and sensitivity:
# 1. "Align ref": Align 3D7 and PacBio assembly with nucmer, and use differences as truth set
# 1. "Align consensus": Apply all variants to 3D7, align to PacBio assembly with nucmer, and use no diff as TP, diffs as TN
# 1. "Map ref, GATK": Map PacBio assembly to 3D7 and call variants with GATK
# 
# Possible future methods include:
# 1. "Map ref, samtools": Map PacBio assembly to 3D7 and call variants with samtools
# 1. "Synthetic reads around variants": Create synthetic reads around variants and see if they map perfectly to PacBio assembly
# 1. "Synthetic reads across genome": Create synthetic reads for whole genome

# # Plan
# Inputs:
# - Original GATK vcf file format (e.g. full build or 5 samples)
# - Names of samples in these vcfs
# - Ref genome used by GATK
# - Chromosomes to analyse (default is all)
# - Output directory
# - Fasta file format for assembled samples
# - VCF file format for assembled samples
# - Threshold variable (VQSLOD)
# - Threshold values (dict for each method)
# - Distance to consider clustered
# 
# Outputs:
# - FDR and sensitivity at each threshold value, broken down by coding/non-coding, SNP/non-STR/STR, clustered/non-clustered, Core/non-core
# - Plot of above
# - Maybe some plots of FDR/sensitivity/heterozygosity/TiTv/inframe_indels broken down by variant-level metrics

# # Setup

# In[1]:

get_ipython().magic('run _shared_setup.ipynb')


# In[2]:

FINAL_VCF_FORMAT="/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/SNP_INDEL_%s.combined.filtered.vcf.gz"


# In[29]:

vcfs_to_evaluate = collections.OrderedDict()
vcfs_to_evaluate['validation5'] = collections.OrderedDict()
vcfs_to_evaluate['release4'] = collections.OrderedDict()
vcfs_to_evaluate['test'] = collections.OrderedDict()

vcfs_to_evaluate['validation5']['Pf3D7_01_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/e/b/f/a/15179/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_02_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/4/8/3/9/15180/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_02_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_03_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/4/2/4/0/15181/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_03_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_04_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/5/c/c/9/15182/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_05_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/6/1/c/9/15183/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_05_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_06_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/f/4/2/9/15184/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_06_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_07_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/7/c/e/d/15185/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_07_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_08_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/0/5/1/0/15186/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_08_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_09_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/c/4/c/f/15187/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_09_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_10_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/a/2/4/2/15188/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_10_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_11_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/e/1/d/15189/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_11_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_12_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/6/6/c/8/15190/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_12_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_13_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/b/6/7/4/15191/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_13_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_14_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/e/2/d/15192/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_14_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf3D7_API_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/0/3/b/e/15193/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_API_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['validation5']['Pf_M76611'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/3/7/6/3/15194/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf_M76611.combined.filtered.vcf.gz'

vcfs_to_evaluate['test']['Pf3D7_04_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/5/c/c/9/15182/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz'
vcfs_to_evaluate['test']['Pf3D7_05_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/6/1/c/9/15183/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_05_v3.combined.filtered.vcf.gz'

for chrom in vcfs_to_evaluate['validation5'].keys():
    vcfs_to_evaluate['release4'][chrom] = FINAL_VCF_FORMAT % chrom


# In[27]:

sample_ids = collections.OrderedDict()
sample_ids['validation5'] = collections.OrderedDict()
sample_ids['release4'] = collections.OrderedDict()
sample_ids['test'] = collections.OrderedDict()

sample_ids['validation5']['7G8'] = '7G8'
sample_ids['validation5']['GB4'] = 'GB4'
sample_ids['validation5']['KH02'] = 'KH02'
sample_ids['validation5']['KE01'] = 'KE01'
sample_ids['validation5']['GN01'] = 'GN01'

sample_ids['test']['7G8'] = '7G8'
sample_ids['test']['GB4'] = 'GB4'
sample_ids['test']['KH02'] = 'KH02'
sample_ids['test']['KE01'] = 'KE01'
sample_ids['test']['GN01'] = 'GN01'

sample_ids['release4']['7G8'] = '7G8'
sample_ids['release4']['GB4'] = 'GB4'
sample_ids['release4']['KH02'] = 'ERS740936'
sample_ids['release4']['KE01'] = 'ERS740937'
sample_ids['release4']['GN01'] = 'ERS740940'


# In[2]:

release4_vcfnp_dir = collections.OrderedDict()
release4_vcfnp_dir['7G8'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['GB4'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.GB4.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['KH02'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740936.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['KE01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740937.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['GN01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740940.vcf.gz.vcfnp_cache'

validation5_vcfnp_dir = collections.OrderedDict()
validation5_vcfnp_dir['7G8'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.7G8.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['GB4'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.GB4.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['KH02'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.KH02.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['KE01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.KE01.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['GN01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.GN01.vcf.gz.vcfnp_cache'


# In[ ]:

GATK_EXE="/software/jre1.7.0_25/bin/java -jar /nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar"
REF_GENOME="/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta"
FINAL_VCF_FORMAT="/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/SNP_INDEL_%s.combined.filtered.vcf.gz"
OUTGROUP_VCF_FORMAT="/nfs/team112_internal/production/release_build/Pf3K/pilot_4_0/with_outgroup_alleles/SNP_INDEL_%s.combined.filtered.vcf.gz"
SAMPLE_FILE="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_samples.txt"
VALIDATION_SAMPLES_VCF_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/%s.5validation.vcf.gz"
INDIVIDAL_VALIDATION_SAMPLES_VCF_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/%s.%s.vcf.gz"
FULL_NPY_FORMAT="/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/with_outgroup_%s"


# In[5]:

non_coding_effects = np.array([b'INTERGENIC', b'INTRON'])
snp_types = np.array([b'SNP', b'MULTIALLELIC_SNP'])


# # Functions

# In[11]:

def subset_vcf_to_samples(
    vcf_to_evaluate='/lustre/scratch110/malaria/dj6/pf3k/output/variant_filtration/5/c/c/9/15182/1_gatk_variant_filter_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz',
    samples_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/vcfnp/validation_samples.txt',
    sample_subset_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/SNP_INDEL_Pf3D7_04_v3.combined.filtered.sample_subset.vcf.gz',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
):
    get_ipython().system('{gatk_exe}         -T SelectVariants         -R {ref_genome_fn}         --variant {vcf_to_evaluate}         --out /dev/stdout         --sample_file {samples_fn}         --excludeNonVariants         | bgzip -c > {sample_subset_vcf_fn}')

    get_ipython().system('tabix -p vcf {sample_subset_vcf_fn}')


# In[37]:

def subset_vcfs_to_samples(
    vcfs_to_evaluate_dict=vcfs_to_evaluate['validation5'],
    sample_ids_dict=sample_ids['validation5'],
    sample_subset_vcf_dir='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    rewrite=False
):
    samples_fn = "%s/validation_samples.txt" % sample_subset_vcf_dir
    if not os.path.exists(sample_subset_vcf_dir):
        os.makedirs(sample_subset_vcf_dir)
    if not os.path.exists(samples_fn) or rewrite:
        samples_fo = open(samples_fn, 'w')
        for sample_id in sample_ids_dict.values():
            print(sample_id, file=samples_fo)
        samples_fo.close()

    for vcf_to_evaluate in vcfs_to_evaluate_dict.values():
        sample_subset_vcf_fn = "%s/%s" % (sample_subset_vcf_dir, os.path.basename(vcf_to_evaluate).replace(".vcf", ".sample_subset.vcf"))
        if not os.path.exists(sample_subset_vcf_fn) or rewrite:
            subset_vcf_to_samples(
                vcf_to_evaluate=vcf_to_evaluate,
                samples_fn=samples_fn,
                sample_subset_vcf_fn=sample_subset_vcf_fn,
                ref_genome_fn=ref_genome_fn,
                gatk_exe=gatk_exe,
            )                                  


# In[34]:

def combine_chromosome_vcfs(
    chromosomes=vcfs_to_evaluate['validation5'].keys(),
    sample_subset_vcf_format='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/SNP_INDEL_%s.combined.filtered.sample_subset.vcf.gz',
    combined_sample_subset_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/SNP_INDEL.all.combined.filtered.sample_subset.vcf.gz',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    rewrite=False
):
    if not os.path.exists(combined_sample_subset_vcf_fn) or rewrite:
        variants_string = ''
        for chromosome in chromosomes:
            variants_string=variants_string+"-V %s " % (sample_subset_vcf_format % chromosome)
        get_ipython().system('{gatk_exe}             -T CombineVariants             -R {ref_genome_fn}             {variants_string}             -genotypeMergeOptions UNSORTED             --out /dev/stdout             | bgzip -c > {combined_sample_subset_vcf_fn}')

        get_ipython().system('tabix -p vcf {combined_sample_subset_vcf_fn}')
    


# In[39]:

def subset_vcf_to_sample(
    combined_sample_subset_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/SNP_INDEL.all.combined.filtered.sample_subset.vcf.gz',
    sample='7G8',
    sample_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/7G8.vcf.gz',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    rewrite=False
):
    if not os.path.exists(sample_vcf_fn) or rewrite:
        get_ipython().system('{gatk_exe}             -T SelectVariants             -R {ref_genome_fn}             --variant {combined_sample_subset_vcf_fn}             --out /dev/stdout             --sample_name {sample}             --excludeNonVariants             --removeUnusedAlternates             | bgzip -c > {sample_vcf_fn}')

        get_ipython().system('tabix -p vcf {sample_vcf_fn}')


# In[48]:

import vcfnp

def run_vcfnp(
    sample_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/7G8.vcf.gz',
    sample_vcfnp_dir='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/sample_subset_vcfs/vcfnp/7G8'
):
    v = vcfnp.variants(
        vcf_fn=sample_vcf_fn,
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
        cache=True,
        cachedir=sample_vcfnp_dir
    )
    c = vcfnp.calldata_2d(
        vcf_fn=sample_vcf_fn,
        progress=10000,
        fields=['GT', 'AD'],
        arities={'AD': 3},
        cache=True,
        cachedir=sample_vcfnp_dir
    )
    return(v, c)


# In[49]:

def full_pipeline(
    callset='validation5',
    threshold_variable='VQSLOD',
    threshold_logic='>',
    align_consensus_threshold_values=[-2, 0, 2, 4, 6],
    delta_filter_is=[97],
    output_dir='/lustre/scratch110/malaria/rp7/Pf3k/fdr',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23',
    regions_fn = '$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz',
    vcfstreamsort_exe='/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    snpeff_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/snpEff.jar',
    snpeff_db='Pfalciparum_GeneDB_Aug2015',
    vcfannotate_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcftools/vcftools_0.1.10/bin/vcf-annotate',
    rewrite=False
):
    vcfs_to_evaluate_dict=vcfs_to_evaluate[callset]
    sample_ids_dict=sample_ids[callset]
    sample_subset_vcf_dir="%s/%s/sample_subset_vcfs" % (output_dir, callset)
    sample_subset_vcf_format="%s/SNP_INDEL_%%s.combined.filtered.sample_subset.vcf.gz" % sample_subset_vcf_dir
    combined_sample_subset_vcf_fn=sample_subset_vcf_format % 'all'
    sample_vcf_format='%s/%%s.vcf.gz' % sample_subset_vcf_dir
    sample_vcfnp_dir="%s/%s/sample_subsets_vcfnp" % (output_dir, callset)
    sample_vcfnp_format='%s/%%s' % sample_vcfnp_dir
    consensus_dir="%s/%s/consensus_alignment" % (output_dir, callset)

    print("Subsetting vcfs to samples %s" % callset)
    subset_vcfs_to_samples(vcfs_to_evaluate_dict, sample_ids_dict, sample_subset_vcf_dir, ref_genome_fn, gatk_exe, rewrite)

    print("Combining chromosomal subset vcfs %s" % callset)
    combine_chromosome_vcfs(vcfs_to_evaluate_dict.keys(), sample_subset_vcf_format, combined_sample_subset_vcf_fn, ref_genome_fn, gatk_exe, rewrite)

    align_consensus_evaluation_metrics = collections.OrderedDict()
    for sample in sample_ids_dict.keys():
        print("Creating subset vcf %s" % sample)
        sample_vcf_fn = sample_vcf_format % sample
        subset_vcf_to_sample(
            combined_sample_subset_vcf_fn,
            sample_ids_dict[sample],
            sample_vcf_fn,
            ref_genome_fn,
            gatk_exe,
            rewrite
        )
        print("Creating vcfnp %s" % sample)
        v, c = run_vcfnp(sample_vcf_fn, sample_vcfnp_format % sample)
        
        align_consensus_evaluation_metrics[sample] = collections.OrderedDict()
        for align_consensus_threshold_value in align_consensus_threshold_values:
            align_consensus_evaluation_metrics[sample][str(align_consensus_threshold_value)] = collections.OrderedDict()
            for delta_filter_i in delta_filter_is:
                print("Running align consensus pipeline %s %s %s" % sample, align_consensus_threshold_value, delta_filter_i)
                align_consensus_evaluation_metrics[sample][str(align_consensus_threshold_value)][str(delta_filter_i)],
                _, _ = align_consensus_pipeline(
                    v,
                    c,
                    sample,
                    threshold_variable=threshold_variable,
                    threshold_logic=threshold_logic,
                    align_consensus_threshold_value=align_consensus_threshold_value,
                    delta_filter_i=delta_filter_i,
                    ref_genome_fn=ref_genome_fn,
                    consensus_dir=consensus_dir,
                    nucmer_dir=nucmer_dir,
                    regions_fn=regions_fn,
                    vcfstreamsort_exe=vcfstreamsort_exe,
                    gatk_exe=gatk_exe,
                    snpeff_exe=snpeff_exe,
                    snpeff_db=snpeff_db,
                    vcfannotate_exe=vcfannotate_exe,
                    rewrite=False
                )
    
    return(align_consensus_evaluation_metrics)
    


# In[274]:

def align_ref_pipeline(
    v,
    c,
    sample='7G8',
    threshold_variable='VQSLOD',
    threshold_logic='>',
    align_consensus_threshold_value=0,
#     vqslod_threshold=0.0,
    delta_filter_i=97,
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    align_ref_dir='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/align_ref',
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23',
    regions_fn = '$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz',
    vcfstreamsort_exe='/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    snpeff_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/snpEff.jar',
    snpeff_db='Pfalciparum_GeneDB_Aug2015',
    vcfannotate_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcftools/vcftools_0.1.10/bin/vcf-annotate',
    rewrite=False
):

    this_align_ref_dir = "%s/%s/%s/%s" % (align_ref_dir, threshold_variable, align_consensus_threshold_value, delta_filter_i)
    pos_fn = '%s/pos/%s.pos' % (this_align_ref_dir, sample)
    consensus_fasta_fn = '%s/fasta/%s.fasta' % (this_align_ref_dir, sample)
    assembly_fasta_fn = determine_assembly_fasta(sample)
    nucmer_filestem = '%s/nucmer/%s' % (this_align_ref_dir, sample)
    
    if not os.path.exists(os.path.dirname(nucmer_filestem)):
        os.makedirs(os.path.dirname(nucmer_filestem))

    print("Creating array %s" % sample)
    variants_array = create_variants_array(v, c)
    if threshold_logic == '>':
        flt_pass = variants_array[threshold_variable] > align_consensus_threshold_value
    elif threshold_logic == '<':
        flt_pass = variants_array[threshold_variable] < align_consensus_threshold_value
    else:
        stop("Error: unknown threshold_logic %s" % threshold_logic)
        
    flt_core = variants_array['RegionType'] == b'Core'
    flt_core_pass = (flt_pass & flt_core)

    print("Running nucmer %s" % sample)
    run_nucmer(ref_genome_fn, assembly_fasta_fn, nucmer_filestem, delta_filter_i, nucmer_dir, rewrite)

    print("Coverting nucmer output to vcf %s" % sample)
    convert_Csnp(sample, consensus_fasta_fn, nucmer_filestem, rewrite)

    print("Annotate nucmer vcf %s" % sample)
    annotate_vcf(unannotated_vcf_fn, ref_genome_fn, regions_fn, vcfstreamsort_exe, gatk_exe, snpeff_exe, snpeff_db, vcfannotate_exe, rewrite)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_array = run_vcfnp(annotated_vcf_fn)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_closest_array = distance_to_closest(nucmer_array)
    
    print("Determining closest nucmer calls %s" % sample)
    variants_nearest_array = distance_to_closest(distance_to_nearest(variants_array[flt_core_pass], nucmer_array))

    print("Determing accessibility %s" % sample)
    variants_accessible_array = within_accessible(variants_nearest_array)
    
    print("Characterising variants %s" % sample)
    plot_histograms(variants_accessible_array)
    
    print("Plotting histogram of distances %s" % sample)
    align_consensus_evaluation_metrics = characterise_nucmer(nucmer_closest_array, variants_accessible_array)
    
    return(align_consensus_evaluation_metrics, variants_accessible_array, nucmer_closest_array)
        


# In[274]:

def align_consensus_pipeline(
    v,
    c,
    sample='7G8',
    threshold_variable='VQSLOD',
    threshold_logic='>',
    align_consensus_threshold_value=0,
#     vqslod_threshold=0.0,
    delta_filter_i=97,
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    consensus_dir='/lustre/scratch110/malaria/rp7/Pf3k/fdr/validation5/consensus_alignment',
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23',
    regions_fn = '$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz',
    vcfstreamsort_exe='/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    snpeff_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/snpEff.jar',
    snpeff_db='Pfalciparum_GeneDB_Aug2015',
    vcfannotate_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcftools/vcftools_0.1.10/bin/vcf-annotate',
    rewrite=False
):

    this_consensus_dir = "%s/%s/%s/%s" % (consensus_dir, threshold_variable, align_consensus_threshold_value, delta_filter_i)
    pos_fn = '%s/pos/%s.pos' % (this_consensus_dir, sample)
    consensus_fasta_fn = '%s/fasta/%s.fasta' % (this_consensus_dir, sample)
    assembly_fasta_fn = determine_assembly_fasta(sample)
    nucmer_filestem = '%s/nucmer/%s' % (this_consensus_dir, sample)
    unannotated_vcf_fn = "%s.3d7coordinates.vcf" % nucmer_filestem
    annotated_vcf_fn = "%s.3d7coordinates.annotated.vcf.gz" % nucmer_filestem
    
    if not os.path.exists(os.path.dirname(nucmer_filestem)):
        os.makedirs(os.path.dirname(nucmer_filestem))

    print("Creating array %s" % sample)
    variants_array = create_variants_array(v, c)
    if threshold_logic == '>':
        flt_pass = variants_array[threshold_variable] > align_consensus_threshold_value
    elif threshold_logic == '<':
        flt_pass = variants_array[threshold_variable] < align_consensus_threshold_value
    else:
        stop("Error: unknown threshold_logic %s" % threshold_logic)
        
    flt_core = variants_array['RegionType'] == b'Core'
    flt_core_pass = (flt_pass & flt_core)

    print("Creating consensus %s" % sample)
    create_consensus(variants_array[['CHROM', 'POS', 'REF', 'alt']][flt_core_pass], ref_genome_fn, pos_fn, consensus_fasta_fn, rewrite)
    
    print("Running nucmer %s" % sample)
    run_nucmer(consensus_fasta_fn, assembly_fasta_fn, nucmer_filestem, delta_filter_i, nucmer_dir, rewrite)

    print("Coverting nucmer output to vcf %s" % sample)
    convert_Csnp(sample, consensus_fasta_fn, nucmer_filestem, rewrite)

    print("Converting nucmer vcf coordinates %s" % sample, rewrite)
    _ = convert_coords(pos_fn, nucmer_filestem, rewrite)

    print("Annotate nucmer vcf %s" % sample)
    annotate_vcf(unannotated_vcf_fn, ref_genome_fn, regions_fn, vcfstreamsort_exe, gatk_exe, snpeff_exe, snpeff_db, vcfannotate_exe, rewrite)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_array = run_vcfnp(annotated_vcf_fn)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_closest_array = distance_to_closest(nucmer_array)
    
    print("Determining closest nucmer calls %s" % sample)
    variants_nearest_array = distance_to_closest(distance_to_nearest(variants_array[flt_core_pass], nucmer_array))

    print("Determing accessibility %s" % sample)
    variants_accessible_array = within_accessible(variants_nearest_array)
    
    print("Characterising variants %s" % sample)
    plot_histograms(variants_accessible_array)
    
    print("Plotting histogram of distances %s" % sample)
    align_consensus_evaluation_metrics = characterise_nucmer(nucmer_closest_array, variants_accessible_array)
    
    return(align_consensus_evaluation_metrics, variants_accessible_array, nucmer_closest_array)
        


# In[50]:

full_pipeline('test')


# In[160]:

def create_variants_array(
    variants_array,
    calls_array
):
    alt_allele_num = np.array([int(x.astype(str)[0]) for x in calls_array['GT'][:,0]])
    is_het = np.array([x[0] != x[1] for x in calls_array['GT'][:,0]])
    highest_cov_allele = np.argmax(calls_array['AD'][:,0,:], axis=1)
    alt_allele_num[is_het] = highest_cov_allele[is_het]
    flt_minor_het = (alt_allele_num > 0)
    alt_allele_num = alt_allele_num-1
    alt_allele_num[alt_allele_num == -1] = 0
    alts = variants_array['ALT'][np.arange(len(alt_allele_num)), alt_allele_num]
    flt_spanning_del = (alts != b'*')
    variants_array = np.lib.recfunctions.append_fields(
        variants_array,
        'alt',
        alts
    ).data
    flt_all = (
        flt_spanning_del &
        flt_minor_het
    )
    return(variants_array[flt_all])


# In[160]:

def create_variants_array_from_files(
    variants_npy_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache/variants.npy',
    calls_npy_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache/calldata_2d.npy'
):
    variants_array = np.load(variants_npy_fn)
    calls_array = np.load(calls_npy_fn)
#     flt_spanning_del = np.array([x[0] != b'*' for x in variants_array['ALT']])
#     print(np.unique(flt_spanning_del, return_counts=True))
    alt_allele_num = np.array([int(x.astype(str)[0]) for x in calls_array['GT'][:,0]])
    is_het = np.array([x[0] != x[1] for x in calls_array['GT'][:,0]])
    highest_cov_allele = np.argmax(calls_array['AD'][:,0,:], axis=1)
    alt_allele_num[is_het] = highest_cov_allele[is_het]
    flt_minor_het = (alt_allele_num > 0)
#     print(np.unique(flt_minor_het, return_counts=True))
    alt_allele_num = alt_allele_num-1
    alt_allele_num[alt_allele_num == -1] = 0
    alts = variants_array['ALT'][np.arange(len(alt_allele_num)), alt_allele_num]
    flt_spanning_del = (alts != b'*')
    variants_array = np.lib.recfunctions.append_fields(
        variants_array,
        'alt',
        alts
    ).data
    flt_all = (
        flt_spanning_del &
        flt_minor_het
    )
    return(variants_array[flt_all])


# In[43]:

import copy
def create_consensus(
    variants_array,
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    pos_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/pos/7G8.pos',
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
    rewrite=False,
):
    if not os.path.exists(pos_fn) or not os.path.exists(consensus_fasta_fn) or rewrite:
        ref_dict=SeqIO.to_dict(SeqIO.parse(open(ref_genome_fn), "fasta"))
        pos_dir = os.path.dirname(pos_fn)
        if not os.path.exists(pos_dir):
            os.makedirs(pos_dir)
        fasta_dir = os.path.dirname(consensus_fasta_fn)
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)
        pos_fo = open(pos_fn, 'w')
        fasta_fo = open(consensus_fasta_fn, 'w')
        mutable_ref_dict = copy.deepcopy(ref_dict)
        for chrom in mutable_ref_dict:
            mutable_ref_dict[chrom].seq = mutable_ref_dict[chrom].seq.tomutable()
        previous_new_pos = 1
        current_offset = 0
        current_chrom = 'Pf3D7_01_v3'
        for i, rec in enumerate(etl.fromarray(variants_array).data()):
            if i > 0 and current_chrom != rec[0].decode("ascii"):
                print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  
                current_chrom = rec[0].decode("ascii")
                current_offset = 0
                previous_new_pos = 1
            reflen = len(rec[2])
            altlen = len(rec[3])
            pos = int(rec[1])
            new_pos = pos + current_offset
            print("%s\t%d\t%d\t%d\t%d" % (rec[0].decode("ascii"), pos, new_pos, previous_new_pos, current_offset), file=pos_fo)  
            previous_new_pos = new_pos
            startpos = pos + current_offset - 1
            endpos = pos + current_offset + reflen - 1   
            mutable_ref_dict[rec[0].decode("ascii")].seq[startpos:endpos] = rec[3].decode("ascii")
            current_offset = current_offset + altlen - reflen
    #         if i%1000 == 0:
    #             print(current_chrom, current_offset)

        print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  

        for chrom in mutable_ref_dict:
            SeqIO.write(mutable_ref_dict[chrom], fasta_fo, "fasta")
    


# In[146]:

def run_nucmer(
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
    assembly_fasta_fn='/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
    delta_filter_i=98,
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23',
    rewrite=False,
):
    
    if not os.path.exists(nucmer_filestem + '.delta') or rewrite:
        get_ipython().system('{nucmer_dir}/nucmer -p {nucmer_filestem} {consensus_fasta_fn} {assembly_fasta_fn}')
    if not os.path.exists(nucmer_filestem + '.filter.delta') or rewrite:
        get_ipython().system('{nucmer_dir}/delta-filter -m -i {delta_filter_i} -l 1000 {nucmer_filestem}.delta > {nucmer_filestem}.filter.delta')
    if not os.path.exists(nucmer_filestem + '.filter.coords') or rewrite:
        get_ipython().system('{nucmer_dir}/show-coords -THqcl -o {nucmer_filestem}.filter.delta > {nucmer_filestem}.filter.coords')
    if not os.path.exists(nucmer_filestem + '.Csnp') or rewrite:
        get_ipython().system('{nucmer_dir}/show-snps -CHTr {nucmer_filestem}.filter.delta > {nucmer_filestem}.Csnp')


# In[39]:

def convert_Csnp(
    sample='7G8',
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
    rewrite=False,
):

    Csnp_fn=nucmer_filestem + ".Csnp"
    nucmer_vcf_fn=nucmer_filestem + ".vcf"
    if not os.path.exists(nucmer_vcf_fn) or rewrite:
        ref_dict=SeqIO.to_dict(SeqIO.parse(open(consensus_fasta_fn), "fasta"))

    #     Write VCF header
        fo = open(nucmer_vcf_fn, 'w')
        fo.write("##fileformat=VCFv4.1\n")
        fo.write("##description=This file created with convert_Csnp function in 20150929_consensus_alignment.ipynb\n")
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    #     Write data
        ref_chroms = ['Pf_M76611', 'Pf3D7_API_v3'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
    #     ref_chroms = ['M76611', 'PFC10_API_IRAB'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
        query_chroms = ["Pf%s_MT" % sample, "Pf%s_API" % sample] + ["Pf%s_%02d" % (sample, i) for i in range(1, 15)]
        chrom_dict = dict(zip(ref_chroms, query_chroms))

        tbl_csnp = (etl.fromtsv(Csnp_fn).select(lambda rec: rec[9] == chrom_dict[rec[8]]))

        current_state = "SNP"
        ins_sequence = ''
        previous_chrom = ''
        for rec in tbl_csnp:
            (pos, ref, alt, pos2, buff, dist, frm, frm2, chrom, query_chrom) = rec
            pos = int(pos)
            pos2 = int(pos2)
            if (
                frm != '1' or
                not(ref in ['A', 'C', 'T', 'G', '.']) or
                not(alt in ['A', 'C', 'T', 'G', '.'])
            ):
                return("error")
            if alt == '.':
                variant_type = 'Del'
            elif ref == '.':
                variant_type = 'Ins'
            else:
                variant_type = 'SNP'
            if variant_type == 'SNP':
                if current_state == 'Del':
                    del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                    del_alt = ref_dict[previous_chrom][variant_start-1]
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
                if current_state == 'Ins':
                    ins_ref = ref_dict[previous_chrom][variant_start-1]
                    ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                    ins_sequence = ''
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, pos, ref, alt))
                current_state = "SNP"
            elif variant_type == 'Del':
                if current_state == 'Del':
                    if pos > (variant_end+1) or chrom != previous_chrom: # i.e. we have moved into a different del so need to print out previous
                        if variant_start < 1 or variant_start > len(ref_dict[chrom]):
                            print(chrom, variant_start, ref_dict[chrom], len(ref_dict[chrom]))
                        del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                        del_alt = ref_dict[previous_chrom][variant_start-1]
                        fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
                        variant_start = pos-1
                if current_state == 'Ins':
                    ins_ref = ref_dict[previous_chrom][variant_start-1]
                    ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                    variant_start = pos-1
                if current_state == 'SNP':
                    variant_start = pos-1
                variant_end = pos
                current_state = "Del"
            elif variant_type == 'Ins':
                if current_state == 'Del':
                    del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                    del_alt = ref_dict[previous_chrom][variant_start-1]
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
                    variant_start = pos
                if current_state == 'Ins':
                    if pos > (variant_end+1) or chrom != previous_chrom: # i.e. we have moved into a different del so need to print out previous
                        ins_ref = ref_dict[previous_chrom][variant_start-1]
                        ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                        fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                        variant_start = pos
                        ins_sequence = ''
                if current_state == 'SNP':
                    variant_start = pos
                ins_sequence = ins_sequence + alt
                variant_end = pos
                current_state = "Ins"
            else:
                return("error")
            previous_chrom = chrom
        fo.close()
    
    return(nucmer_vcf_fn)


# In[40]:

def convert_coords(
    pos_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/pos/7G8.pos',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
    rewrite=False,
):
    
    nucmer_vcf_fn = nucmer_filestem + ".vcf"
    converted_vcf_fn = nucmer_filestem + ".3d7coordinates.vcf"
    
    tbl_vcf = etl.fromvcf(nucmer_vcf_fn).convert('POS', int).addfield('POS1', lambda rec: rec['POS']+1)
    tbl_pos = (etl.fromtsv(pos_fn)
        .setheader(('CHROM', 'var_pos', 'end', 'start', 'diff'))
        .convert('var_pos', int)
        .convert('start', int)
        .convert('end', int)
        .convert('diff', int)
        .select(lambda rec: rec['end'] > rec['start'])
    )
    
    tbl_vcf_converted = (tbl_vcf
        .intervalleftjoin(tbl_pos, lstart='POS', lstop='POS1', lkey='CHROM',
                          rstart='start', rstop='end', rkey='CHROM')
        .addfield('NEW_POS', lambda rec: rec['POS'] if rec['diff'] is None else rec['POS'] - rec['diff'])
    )

    if not os.path.exists(converted_vcf_fn) or rewrite:
    #     Write VCF header
        fo = open(converted_vcf_fn, 'w')
        fo.write("##fileformat=VCFv4.1\n")
        fo.write("##description=This file created with tbl_converted_coords function in 20150929_consensus_alignment.ipynb\n")
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    #     Write VCF rows
        for rec in tbl_vcf_converted.data():
            fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (rec[0], rec[13], rec[3], rec[4][0]))
        fo.close()
        
    return tbl_vcf_converted
    
    


# In[101]:

def annotate_vcf(
    unannotated_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.3d7coordinates.vcf',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    regions_fn = '$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz',
#     vcfstreamsort_exe='/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    vcfstreamsort_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    snpeff_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/snpEff.jar',
    snpeff_db='Pfalciparum_GeneDB_Aug2015',
    vcfannotate_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcftools/vcftools_0.1.10/bin/vcf-annotate',
    rewrite=False
):
    
    sorted_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".sorted.vcf")
    leftaligned_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".leftaligned.vcf")
    snpeff_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".snpeff.vcf")
    snpeff_annotated_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".snpeff_annotated.vcf")
    annotated_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".annotated.vcf")
    
    if not os.path.isfile(sorted_vcf_fn) or rewrite:
        get_ipython().system('{vcfstreamsort_exe} < {unannotated_vcf_fn} > {sorted_vcf_fn}')

    if (not os.path.isfile(leftaligned_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{gatk_exe} -T LeftAlignAndTrimVariants         -R {ref_genome_fn}         -V {sorted_vcf_fn}         -o {leftaligned_vcf_fn} #         2> /dev/null')

    if (not os.path.isfile(snpeff_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{snpeff_exe}         -v -o gatk {snpeff_db}         {leftaligned_vcf_fn}         -no-downstream         -no-upstream         -onlyProtein         > {snpeff_vcf_fn} #         2> /dev/null')

    if (not os.path.isfile(snpeff_annotated_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{gatk_exe}         -T VariantAnnotator         -R {ref_genome_fn}         -A TandemRepeatAnnotator         -A SnpEff         --variant {leftaligned_vcf_fn}         --snpEffFile {snpeff_vcf_fn}         -o {snpeff_annotated_vcf_fn} #         2> /dev/null')

    if not os.path.isfile(annotated_vcf_fn+'.gz') or rewrite:
        get_ipython().system("cat {snpeff_annotated_vcf_fn}         | {vcfannotate_exe} -a {regions_fn}            -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.'            -c CHROM,FROM,TO,INFO/RegionType         > {annotated_vcf_fn}")

        get_ipython().system('bgzip -f {annotated_vcf_fn}')
        get_ipython().system('tabix -p vcf -f {annotated_vcf_fn}.gz')

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.SNP.vcf')+'.gz') or rewrite:
        get_ipython().system("{gatk_exe}         -T SelectVariants         -R {ref_genome_fn}         -V {annotated_vcf_fn}.gz         -selectType SNP         -o {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')} 2> /dev/null")

        get_ipython().system("bgzip -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}")
        get_ipython().system("tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}.gz")

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')+'.gz') or rewrite:
        get_ipython().system("{gatk_exe}         -T SelectVariants         -R {ref_genome_fn}         -V {annotated_vcf_fn}.gz         -xlSelectType SNP         -o {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')} 2> /dev/null")

        get_ipython().system("bgzip -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}")
        get_ipython().system("tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}.gz")


# In[96]:

annotate_vcf()


# In[190]:

import vcfnp

def run_vcfnp(
    annotated_vcf_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.3d7coordinates.annotated.vcf.gz',
    rewrite=False
):
    output_fn = annotated_vcf_fn + ".vcfnp_cache/variants.npy"
    if rewrite and os.path.exists(output_fn):
        os.remove(output_fn)
    v = vcfnp.variants(
        vcf_fn=annotated_vcf_fn,
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
    return(v)
#             'VariantType': 'a40',


# In[19]:

def find_nearest(a, b):
    nearest_before = abs(a - b[np.searchsorted(b, a)-1])
    nearest_after = abs(a - b[np.searchsorted(b, a) % len(b)])
    return(np.minimum(nearest_before, nearest_after))
find_nearest(
    np.array([-1, 1, 2, 5, 9, 12, 22]),
    np.array([0, 1, 4, 6, 12])
)


# In[251]:

def find_closest(a):
    nearest_before = abs(a - np.roll(a, 1))
    nearest_after = abs(np.roll(a, -1) - a)
    return(np.minimum(nearest_before, nearest_after))

find_closest(
    np.array([0, 1, 4, 6, 12])
)


# In[195]:

def distance_to_nearest(
    variants_array,
    nucmer_array
):
    distance_array = np.zeros(len(variants_array), dtype='u4')
    for chrom in np.unique(variants_array['CHROM']):
#         print(chrom)
        variants_pos_this_chrom = variants_array['POS'][variants_array['CHROM'] == chrom]
        nucmer_pos_this_chrom = nucmer_array['POS'][nucmer_array['CHROM'] == chrom]
        nearest_this_chrom = find_nearest(variants_pos_this_chrom, nucmer_pos_this_chrom)
        distance_array[variants_array['CHROM'] == chrom] = nearest_this_chrom
    new_array = np.lib.recfunctions.append_fields(
        variants_array,
        'nearest_nucmer',
        distance_array
    ).data
#     new_array = variants_array.addfield('nearest_nucmer', distance_array)
    return(new_array)
        


# In[250]:

def distance_to_closest(
    variants_array
):
    closest_distance = np.zeros(len(variants_array), dtype='u4')
    for chrom in np.unique(variants_array['CHROM']):
        variants_pos_this_chrom = variants_array['POS'][variants_array['CHROM'] == chrom]
        closest_this_chrom = find_closest(variants_pos_this_chrom)
        closest_distance[variants_array['CHROM'] == chrom] = closest_this_chrom
    new_array = np.lib.recfunctions.append_fields(
        variants_array,
        'distance_to_closest',
        closest_distance
    ).data
    return(new_array)
        


# In[163]:

def within_accessible(
    variants_array,
    assembly_sample='7G8',
    filter_coords_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.filter.coords'
):
#     Figure out how to do this in numpy. Numpy intervals? Or just use in1d?
    is_accessible = np.zeros(len(variants_array), dtype='bool')
    array_coords = (
        etl.fromtsv(filter_coords_fn)
        .pushheader(['ref_start', 'ref_end', 'assembly_start', 'assembly_end', 'ref_length', 'assembly_length', 'identity',
             'ref_chrom_length', 'assembly_chrom_length', 'ref_pct', 'assembly_pct', 'ref_chrom', 'assembly_chrom',
             'ignore'])
        .cut(['ref_chrom', 'assembly_chrom', 'ref_start', 'ref_end'])
        .convertnumbers()
        .toarray()
    )
    for chrom in np.unique(variants_array['CHROM']):
        assembly_chrom = "Pf%s_%s" % (assembly_sample, chrom.decode("ascii")[6:8])
        starts = array_coords['ref_start'][
            (array_coords['ref_chrom'] == chrom.decode("utf-8")) &
            (array_coords['assembly_chrom'] == assembly_chrom)
        ]
        ends = array_coords['ref_end'][
            (array_coords['ref_chrom'] == chrom.decode("utf-8")) &
            (array_coords['assembly_chrom'] == assembly_chrom)
        ]
        for i, start in enumerate(starts):
            accessible_pos = np.arange(start, ends[i]+1)
            flt_accessible = ((variants_array['CHROM'] == chrom) & (np.in1d(variants_array['POS'], accessible_pos)))
            is_accessible[flt_accessible] = True
    new_array = np.lib.recfunctions.append_fields(
        variants_array,
        'is_accessible',
        is_accessible
    ).data
#     new_array = variants_array.addfield('is_accessible', is_accessible)
    return(new_array)
        


# In[244]:

def plot_histograms(
    variants_accessible_array,
):
    flt_types = collections.OrderedDict()
    flt_types['Coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['Non-coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['Coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['Non-coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    fig = plt.figure(figsize=(8,6))
    for i, flt_type in enumerate(flt_types):
        num_variants = np.count_nonzero(flt_types[flt_type])
        nearest_nucmers = variants_accessible_array['nearest_nucmer'][
            (flt_types[flt_type]) &
            (variants_accessible_array['is_accessible'])
        ]
        num_accessible = len(nearest_nucmers)        
        num_0bp = np.count_nonzero(nearest_nucmers == 0)
        num_10bp = np.count_nonzero(nearest_nucmers <= 10)
#         print("%s: #variants=%d, #accessible=%d(%4.2f%%), #0bp=%d(%4.2f%%), #10bp=%d(%4.2f%%)" % (
#                 flt_type,
#                 num_variants,
#                 num_accessible,
#                 (num_accessible/num_variants) * 100,
#                 num_0bp,
#                 (num_0bp/num_accessible) * 100,
#                 num_10bp,
#                 (num_10bp/num_accessible) * 100
#             )
#         )
        ax = fig.add_subplot(4, 1, i+1)
        ax.hist(
            variants_accessible_array['nearest_nucmer'][
                (flt_types[flt_type]) &
                (variants_accessible_array['is_accessible'])
            ],
            bins=np.arange(0, 100)
        )
        ax.set_title(flt_type)
    


# In[260]:

def characterise_nucmer(
    nucmer_array,
    variants_accessible_array,
    clustered_distance_threshold=20
):
    nucmer_is_snp = np.array([(len(x['REF'])-len(x['ALT'][0])) == 0 for x in nucmer_array])
    nucmer_is_clustered = (nucmer_array['distance_to_closest'] <= clustered_distance_threshold)
    gatk_is_clustered = (variants_accessible_array['distance_to_closest'] <= clustered_distance_threshold)
    flt_types = collections.OrderedDict()
    flt_types['nucmer'] = collections.OrderedDict()
    flt_types['gatk'] = collections.OrderedDict()
    flt_types['nucmer']['Unclustered Coding SNP'] = (
        (nucmer_is_snp) &
        np.logical_not(nucmer_is_clustered) &
        np.logical_not(np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Unclustered Non-coding SNP'] = (
        (nucmer_is_snp) &
        np.logical_not(nucmer_is_clustered) &
        (np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Unclustered Coding indel'] = (
        np.logical_not(nucmer_is_snp) &
        np.logical_not(nucmer_is_clustered) &
        np.logical_not(np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Unclustered Non-coding indel'] = (
        np.logical_not(nucmer_is_snp) &
        np.logical_not(nucmer_is_clustered) &
        (np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Clustered Coding SNP'] = (
        (nucmer_is_snp) &
        (nucmer_is_clustered) &
        np.logical_not(np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Clustered Non-coding SNP'] = (
        (nucmer_is_snp) &
        (nucmer_is_clustered) &
        (np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Clustered Coding indel'] = (
        np.logical_not(nucmer_is_snp) &
        (nucmer_is_clustered) &
        np.logical_not(np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['nucmer']['Clustered Non-coding indel'] = (
        np.logical_not(nucmer_is_snp) &
        (nucmer_is_clustered) &
        (np.in1d(nucmer_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Unclustered Coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(gatk_is_clustered) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Unclustered Non-coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(gatk_is_clustered) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Unclustered Coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(gatk_is_clustered) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Unclustered Non-coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        np.logical_not(gatk_is_clustered) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Clustered Coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (gatk_is_clustered) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Clustered Non-coding SNP'] = (
        (np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (gatk_is_clustered) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Clustered Coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (gatk_is_clustered) &
        np.logical_not(np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    flt_types['gatk']['Clustered Non-coding indel'] = (
        np.logical_not(np.in1d(variants_accessible_array['VariantType'], snp_types)) &
        (gatk_is_clustered) &
        (np.in1d(variants_accessible_array['SNPEFF_EFFECT'], non_coding_effects))
    )
    fig = plt.figure(figsize=(8,6))
    evaluation_metrics = collections.OrderedDict()
    for i, flt_type in enumerate(flt_types['nucmer']):
        evaluation_metrics[flt_type] = collections.OrderedDict()
        evaluation_metrics[flt_type]['num_variants'] = np.count_nonzero(flt_types['gatk'][flt_type])
        nearest_nucmers = variants_accessible_array['nearest_nucmer'][
            (flt_types['gatk'][flt_type]) &
            (variants_accessible_array['is_accessible'])
        ]
        evaluation_metrics[flt_type]['num_accessible'] = len(nearest_nucmers)
        evaluation_metrics[flt_type]['num_FP'] = np.count_nonzero(nearest_nucmers == 0)
        evaluation_metrics[flt_type]['num_TP'] = np.count_nonzero(nearest_nucmers > 0)
        evaluation_metrics[flt_type]['num_FN'] = np.count_nonzero(flt_types['nucmer'][flt_type])
        print("%s: #variants=%d, #accessible=%d(%4.2f%%), #TP=%d(%4.2f%%), #FP=%d(FDR=%4.2f%%), #FN=%d(sensitivity=%4.2f%%)" % (
                flt_type,
                evaluation_metrics[flt_type]['num_variants'],
                evaluation_metrics[flt_type]['num_accessible'],
                (evaluation_metrics[flt_type]['num_accessible']/evaluation_metrics[flt_type]['num_variants']) * 100,
                evaluation_metrics[flt_type]['num_TP'],
                (evaluation_metrics[flt_type]['num_TP']/evaluation_metrics[flt_type]['num_accessible']) * 100,
                evaluation_metrics[flt_type]['num_FP'],
                (evaluation_metrics[flt_type]['num_FP']/evaluation_metrics[flt_type]['num_accessible']) * 100,
                evaluation_metrics[flt_type]['num_FN'],
                (evaluation_metrics[flt_type]['num_TP']/(evaluation_metrics[flt_type]['num_TP']+evaluation_metrics[flt_type]['num_FN'])) * 100
            )
        )
    
        ax = fig.add_subplot(8, 2, (i*2)+1)
        ax.hist(
            variants_accessible_array['distance_to_closest'][
                (flt_types['gatk'][flt_type])
            ],
            bins=np.arange(0, 100)
        )
        ax.set_title('GATK ' + flt_type)
    
        ax = fig.add_subplot(8, 2, (i*2)+2)
        ax.hist(
            nucmer_array['distance_to_closest'][
                (flt_types['nucmer'][flt_type])
            ],
            bins=np.arange(0, 100)
        )
        ax.set_title('nucmer ' + flt_type)
    
    return(evaluation_metrics)


# In[22]:

def determine_assembly_fasta(sample='7G8'):
    if sample in ('7G8', 'GB4'):
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf%s.Jul2015.fasta" % sample)
    else:
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/Pf%s.Jun2015.fasta" % sample)


# In[274]:

def run_fdr_pipeline(
    sample='7G8',
    vqslod_threshold=0.0,
    delta_filter_i=97,
    vcfnp_dir_dict=release4_vcfnp_dir,
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    consensus_dir='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment',
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23',
    regions_fn = '$HOME/src/github/malariagen/pf-crosses/meta/regions-20130225.bed.gz',
    vcfstreamsort_exe='/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcflib/vcflib/bin/vcfstreamsort',
    gatk_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/gatk/GenomeAnalysisTK.jar',
    snpeff_exe='/software/jre1.7.0_25/bin/java -jar $HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/snpEff.jar',
    snpeff_db='Pfalciparum_GeneDB_Aug2015',
    vcfannotate_exe='$HOME/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/vcftools/vcftools_0.1.10/bin/vcf-annotate',
    rewrite=False
):

    variants_npy_fn = vcfnp_dir_dict[sample] + '/variants.npy'
    calls_npy_fn = vcfnp_dir_dict[sample] + '/calldata_2d.npy'
    pos_fn = '%s/pos/%s.pos' % (consensus_dir, sample)
    consensus_fasta_fn = '%s/fasta/%s.fasta' % (consensus_dir, sample)
    assembly_fasta_fn = determine_assembly_fasta(sample)
    nucmer_filestem = '%s/nucmer/%s' % (consensus_dir, sample)
    unannotated_vcf_fn = "%s.3d7coordinates.vcf" % nucmer_filestem
    annotated_vcf_fn = "%s.3d7coordinates.annotated.vcf.gz" % nucmer_filestem
    
    if not os.path.exists(os.path.dirname(nucmer_filestem)):
        os.makedirs(os.path.dirname(nucmer_filestem))

    print("Creating array %s" % sample)
    variants_array = create_variants_array_from_files(variants_npy_fn, calls_npy_fn)
    flt_pass = variants_array['VQSLOD'] > vqslod_threshold
    flt_core = variants_array['RegionType'] == b'Core'
    flt_core_pass = (flt_pass & flt_core)

    print("Creating consensus %s" % sample)
    create_consensus(variants_array[['CHROM', 'POS', 'REF', 'alt']][flt_core_pass], ref_genome_fn, pos_fn, consensus_fasta_fn, rewrite)
    
    print("Running nucmer %s" % sample)
    run_nucmer(consensus_fasta_fn, assembly_fasta_fn, nucmer_filestem, delta_filter_i, nucmer_dir, rewrite)

    print("Coverting nucmer output to vcf %s" % sample)
    convert_Csnp(sample, consensus_fasta_fn, nucmer_filestem, rewrite)

    print("Converting nucmer vcf coordinates %s" % sample, rewrite)
    _ = convert_coords(pos_fn, nucmer_filestem, rewrite)

    print("Annotate nucmer vcf %s" % sample)
    annotate_vcf(unannotated_vcf_fn, ref_genome_fn, regions_fn, vcfstreamsort_exe, gatk_exe, snpeff_exe, snpeff_db, vcfannotate_exe, rewrite)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_array = run_vcfnp(annotated_vcf_fn)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_closest_array = distance_to_closest(nucmer_array)
    
    print("Determining closest nucmer calls %s" % sample)
    variants_nearest_array = distance_to_closest(distance_to_nearest(variants_array[flt_core_pass], nucmer_array))

    print("Determing accessibility %s" % sample)
    variants_accessible_array = within_accessible(variants_nearest_array)
    
    print("Characterising variants %s" % sample)
    plot_histograms(variants_accessible_array)
    
    print("Plotting histogram of distances %s" % sample)
    evaluation_metrics = characterise_nucmer(nucmer_closest_array, variants_accessible_array)
    
    return(evaluation_metrics, variants_accessible_array, nucmer_closest_array)
        


# In[ ]:




# # Analysis

# In[270]:

eval_metrics = collections.OrderedDict()
for sample in release4_vcfnp_dir:
    print(sample)
    eval_metrics[sample], _, _ = run_fdr_pipeline(sample, rewrite=True)


# In[283]:

for flt_type in eval_metrics['7G8']:
    FDR = np.zeros(len(eval_metrics))
    sensitivity = np.zeros(len(eval_metrics))
    for i, sample in enumerate(eval_metrics):
        FDR[i] = (eval_metrics[sample][flt_type]['num_FP']/eval_metrics[sample][flt_type]['num_accessible'])
        sensitivity[i] = (eval_metrics[sample][flt_type]['num_TP']/(eval_metrics[sample][flt_type]['num_TP']+eval_metrics[sample][flt_type]['num_FN']))
    print("%s: mean FDR = %6.4f, mean sensitivity = %6.4f" % (flt_type, np.mean(FDR), np.mean(sensitivity)))


# In[272]:

for sample in eval_metrics:
    flt_type = 'Unclustered Coding SNP'
    print("%s %s: #variants=%d, #accessible=%d(%4.2f%%), #TP=%d(%4.2f%%), #FP=%d(FDR=%4.2f%%), #FN=%d(sensitivity=%4.2f%%)" % (
            sample,
            flt_type,
            eval_metrics[sample][flt_type]['num_variants'],
            eval_metrics[sample][flt_type]['num_accessible'],
            (eval_metrics[sample][flt_type]['num_accessible']/eval_metrics[sample][flt_type]['num_variants']) * 100,
            eval_metrics[sample][flt_type]['num_TP'],
            (eval_metrics[sample][flt_type]['num_TP']/eval_metrics[sample][flt_type]['num_accessible']) * 100,
            eval_metrics[sample][flt_type]['num_FP'],
            (eval_metrics[sample][flt_type]['num_FP']/eval_metrics[sample][flt_type]['num_accessible']) * 100,
            eval_metrics[sample][flt_type]['num_FN'],
            (eval_metrics[sample][flt_type]['num_TP']/(eval_metrics[sample][flt_type]['num_TP']+eval_metrics[sample][flt_type]['num_FN'])) * 100
        )
    )


# In[273]:

for sample in eval_metrics:
    flt_type = 'Unclustered Non-coding indel'
    print("%s %s: #variants=%d, #accessible=%d(%4.2f%%), #TP=%d(%4.2f%%), #FP=%d(FDR=%4.2f%%), #FN=%d(sensitivity=%4.2f%%)" % (
            sample,
            flt_type,
            eval_metrics[sample][flt_type]['num_variants'],
            eval_metrics[sample][flt_type]['num_accessible'],
            (eval_metrics[sample][flt_type]['num_accessible']/eval_metrics[sample][flt_type]['num_variants']) * 100,
            eval_metrics[sample][flt_type]['num_TP'],
            (eval_metrics[sample][flt_type]['num_TP']/eval_metrics[sample][flt_type]['num_accessible']) * 100,
            eval_metrics[sample][flt_type]['num_FP'],
            (eval_metrics[sample][flt_type]['num_FP']/eval_metrics[sample][flt_type]['num_accessible']) * 100,
            eval_metrics[sample][flt_type]['num_FN'],
            (eval_metrics[sample][flt_type]['num_TP']/(eval_metrics[sample][flt_type]['num_TP']+eval_metrics[sample][flt_type]['num_FN'])) * 100
        )
    )


# In[275]:

eval_metrics_validation5 = collections.OrderedDict()
for sample in validation5_vcfnp_dir:
    print(sample)
    eval_metrics_validation5[sample], _, _ = run_fdr_pipeline(
        sample,
        vcfnp_dir_dict=validation5_vcfnp_dir,
        consensus_dir='/lustre/scratch110/malaria/rp7/Pf3k/validation5/consensus_alignment',
        rewrite=False
    )
    


# In[276]:

for sample in eval_metrics_validation5:
    flt_type = 'Unclustered Coding SNP'
    print("%s %s: #variants=%d, #accessible=%d(%4.2f%%), #TP=%d(%4.2f%%), #FP=%d(FDR=%4.2f%%), #FN=%d(sensitivity=%4.2f%%)" % (
            sample,
            flt_type,
            eval_metrics_validation5[sample][flt_type]['num_variants'],
            eval_metrics_validation5[sample][flt_type]['num_accessible'],
            (eval_metrics_validation5[sample][flt_type]['num_accessible']/eval_metrics_validation5[sample][flt_type]['num_variants']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_TP'],
            (eval_metrics_validation5[sample][flt_type]['num_TP']/eval_metrics_validation5[sample][flt_type]['num_accessible']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_FP'],
            (eval_metrics_validation5[sample][flt_type]['num_FP']/eval_metrics_validation5[sample][flt_type]['num_accessible']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_FN'],
            (eval_metrics_validation5[sample][flt_type]['num_TP']/(eval_metrics_validation5[sample][flt_type]['num_TP']+eval_metrics_validation5[sample][flt_type]['num_FN'])) * 100
        )
    )


# In[277]:

for sample in eval_metrics_validation5:
    flt_type = 'Unclustered Non-coding indel'
    print("%s %s: #variants=%d, #accessible=%d(%4.2f%%), #TP=%d(%4.2f%%), #FP=%d(FDR=%4.2f%%), #FN=%d(sensitivity=%4.2f%%)" % (
            sample,
            flt_type,
            eval_metrics_validation5[sample][flt_type]['num_variants'],
            eval_metrics_validation5[sample][flt_type]['num_accessible'],
            (eval_metrics_validation5[sample][flt_type]['num_accessible']/eval_metrics_validation5[sample][flt_type]['num_variants']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_TP'],
            (eval_metrics_validation5[sample][flt_type]['num_TP']/eval_metrics_validation5[sample][flt_type]['num_accessible']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_FP'],
            (eval_metrics_validation5[sample][flt_type]['num_FP']/eval_metrics_validation5[sample][flt_type]['num_accessible']) * 100,
            eval_metrics_validation5[sample][flt_type]['num_FN'],
            (eval_metrics_validation5[sample][flt_type]['num_TP']/(eval_metrics_validation5[sample][flt_type]['num_TP']+eval_metrics_validation5[sample][flt_type]['num_FN'])) * 100
        )
    )


# In[284]:

for flt_type in eval_metrics_validation5['7G8']:
    FDR = np.zeros(len(eval_metrics_validation5))
    sensitivity = np.zeros(len(eval_metrics_validation5))
    for i, sample in enumerate(eval_metrics_validation5):
        FDR[i] = (eval_metrics_validation5[sample][flt_type]['num_FP']/eval_metrics_validation5[sample][flt_type]['num_accessible'])
        sensitivity[i] = (eval_metrics_validation5[sample][flt_type]['num_TP']/(eval_metrics_validation5[sample][flt_type]['num_TP']+eval_metrics_validation5[sample][flt_type]['num_FN']))
    print("%s: mean FDR = %6.4f, mean sensitivity = %6.4f" % (flt_type, np.mean(FDR), np.mean(sensitivity)))


# In[262]:

eval_metrics, temp, temp2 = run_fdr_pipeline(rewrite=False)


# In[265]:

eval_metrics['Unclustered Coding SNP']


# In[259]:

temp, temp2 = run_fdr_pipeline(rewrite=True)


# In[230]:

temp_alt = np.array([len(x['REF'])-len(x['ALT'][0]) for x in temp2])
temp_is_snp = np.array([(len(x['REF'])-len(x['ALT'][0])) == 0 for x in temp2])
np.unique(temp_alt==0, return_counts=True)
np.unique(temp_is_snp, return_counts=True)


# In[205]:

temp = run_fdr_pipeline(vqslod_threshold=-10.0, rewrite=True)


# In[206]:

temp = run_fdr_pipeline(vqslod_threshold=6.0, rewrite=True)


# In[238]:

temp, temp2 = run_fdr_pipeline(delta_filter_i=98, rewrite=True)


# In[239]:

temp, temp2 = run_fdr_pipeline(delta_filter_i=97, rewrite=True)


# In[175]:

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.hist(temp['nearest_nucmer'], bins=np.arange(0, 100))


# In[179]:

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.hist(temp['nearest_nucmer'][temp['is_accessible']], bins=np.arange(0, 100))


# In[176]:

len(temp)


# In[181]:

print(np.unique(temp['nearest_nucmer'] == 0, return_counts=True))
print(np.unique(temp['nearest_nucmer'] < 0, return_counts=True))
print(np.unique(temp['nearest_nucmer'] > 0, return_counts=True))
print(np.max(temp['nearest_nucmer']))
print(np.max(temp['nearest_nucmer'][temp['is_accessible']]))
print(np.median(temp['nearest_nucmer']))
print(np.median(temp['nearest_nucmer'][temp['is_accessible']]))
print(np.where(temp['nearest_nucmer'] > 50000))


# In[183]:

print(np.unique(temp['is_accessible'], return_counts=True))


# In[120]:

len(temp)


# In[121]:

831/37416


# In[124]:

36585 + 831 


# In[184]:

temp


# # Sandbox

# In[30]:

sample = '7G8'
ref_chroms = ['M76611', 'PFC10_API_IRAB'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
query_chroms = ["Pf%s_MT" % sample, "Pf%s_API" % sample] + ["Pf%s_%02d" % (sample, i) for i in range(1, 15)]
chrom_dict = dict(zip(ref_chroms, query_chroms))


# In[31]:

chrom_dict.keys()


# In[83]:

temp=SeqIO.to_dict(SeqIO.parse(open('/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta'), "fasta"))


# In[86]:

len(temp['Pf3D7_05_v3'])


# In[91]:

temp=SeqIO.to_dict(SeqIO.parse(open('/nfs/users/nfs_r/rp7/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/snpeff/snpEff/data/Pfalciparum_GeneDB_Aug2015/sequences.fa'), "fasta"))
len(temp['Pf3D7_05_v3'])


# In[ ]:




# In[149]:

variants_array = np.load('/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache/variants.npy')


# In[150]:

variants_array[0]


# In[ ]:



