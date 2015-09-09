# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

!mkdir {os.path.join(CACHE_DIR, 'tbl_gatk2')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_vrpipe')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_vrpipe_comparison')}

# <codecell>

vrpipe_vcfs = collections.OrderedDict()
vrpipe_vcfs['PF_apicoplast_genome_1'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/6/d/c/0/14285/1_gatk_combine_variants_gatk3/SNP_INDEL_PF_apicoplast_genome_1.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_01_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/b/1/6/9/14286/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_01_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_02_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/b/b/a/1/14287/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_02_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_03_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/7/4/f/f/14288/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_03_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_04_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/b/f/0/3/14289/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_04_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_05_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/5/4/c/3/14290/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_05_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_06_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/c/8/0/a/14291/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_06_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_07_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/0/4/6/0/14292/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_07_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_08_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/5/d/a/2/14293/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_08_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_09_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/6/4/0/9/14294/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_09_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_10_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/3/5/8/b/14295/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_10_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_11_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/f/c/a/1/14296/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_11_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_12_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/8/6/6/b/14297/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_12_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_13_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/9/6/7/b/14298/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_13_v3.combined.vcf.gz'
vrpipe_vcfs['Pf3D7_14_v3'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/9/3/d/7/14299/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf3D7_14_v3.combined.vcf.gz'
vrpipe_vcfs['Pf_M76611'] = '/lustre/scratch110/malaria/dj6/pf3k/output/combine_variants/9/d/8/9/14300/1_gatk_combine_variants_gatk3/SNP_INDEL_Pf_M76611.combined.vcf.gz'

# <codecell>

PLOTS_DIR = '/Users/rpearson/Documents/projects/Pf3k_techbm/slides/20150904_vrpipe_assessment/plots'
!mkdir -p {PLOTS_DIR}

# <codecell>

# CODING_EFFECTS = ['NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED']
NONCODING_EFFECTS = ['INTERGENIC', 'INTRON', 'TRANSCRIPT']
 

# <codecell>

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))

# <headingcell level=1>

# Determine samples

# <codecell>

tbl_samples_to_process = (tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
    .selecteq('To be used for', 'Validation')
    .addfield('truth_vcf_filestem', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'truth_vcfs',
        "truth_%s" % rec['Isolate code']
    ))
    .cutout('To be used for')
    .convert('ox_code', lambda v: 'ERS740937', where=lambda rec: rec['Isolate code'] == 'KE01')
    .convert('ox_code', lambda v: 'ERS740936', where=lambda rec: rec['Isolate code'] == 'KH02')
    .convert('ox_code', lambda v: 'ERS740940', where=lambda rec: rec['Isolate code'] == 'GN01')
)
tbl_samples_to_process.displayall(index_header=True)

# <headingcell level=1>

# Functions

# <codecell>

# def determine_consensus_region(rec, column1=4, column2=11):
#     if rec[column1] is None:
#         return rec[column2]
#     elif rec[column2] is None:
#         return rec[column1]
#     elif rec[column1] == rec[column2]:
#         return rec[column1]
#     else:
#         return "ERROR__%s__%s" % (rec[column1], rec[column2])
    

# <codecell>

# def upstream(prv, cur, nxt):
#     if prv is None:
#         return None
#     else:
#         return cur.POS - prv.POS

# def downstream(prv, cur, nxt):
#     if nxt is None:
#         return None
#     else:
#         return nxt.POS - cur.POS

# def nearest(distance_previous, distance_next):
#     if distance_previous is None:
#         return(distance_next)
#     elif distance_next is None:
#         return(distance_previous)
#     else:
#         return(min(distance_previous, distance_next))

# <codecell>

def determine_mode(rec):
    if len(rec['REF']) == len(str(rec['ALT'])):
        return('SNP')
    else:
        return('INDEL')

# <codecell>

# def trim_variants(ref, alt):
# #     print(ref, alt)
#     while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
#         ref = ref[0:-1]
#         alt = alt[0:-1]
#     return((ref, alt))
 
# trim_variants('AATC', 'TC')[1]

# <codecell>

def determine_filter(rec):
    if len(rec['all_FILTERS']) == 0:
        return('')
    else:
        return(str(rec['all_FILTERS'][0]))

# <codecell>

fields_to_compare = ['QUAL', 'FILTER', 'AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ',
                     'MQRankSum', 'QD', 'ReadPosRankSum', 'RegionType', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS',
                     'SNPEFF_IMPACT', 'SOR', 'VQSLOD', 'culprit', 'set', 'AD', 'GT']
vrpipe_renames = dict(zip(fields_to_compare, ['vrpipe_%s' % field_name for field_name in fields_to_compare]))
gatk_renames = dict(zip(fields_to_compare, ['gatk_%s' % field_name for field_name in fields_to_compare]))

# <codecell>

def assess_vrpipe(isolate_code='7G8', chromosome='Pf3D7_01_v3', rewrite=False):
    ox_code = tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('ox_code')[0]
    print(isolate_code, ox_code)
    
    tbl_vrpipe_cache_fn = os.path.join(CACHE_DIR, 'tbl_vrpipe', "tbl_vrpipe_%s.%s" % (isolate_code, chromosome))
    tbl_gatk2_cache_fn = os.path.join(CACHE_DIR, 'tbl_gatk2', "tbl_gatk2_%s.%s" % (isolate_code, chromosome))
    tbl_vrpipe_comparison_cache_fn = os.path.join(CACHE_DIR, 'tbl_vrpipe_comparison', "tbl_vrpipe_comparison_%s.%s" % (isolate_code, chromosome))
    
    if not os.path.exists(tbl_vrpipe_cache_fn) or rewrite:
        vrpipe_vcf_fn = vrpipe_vcfs[chromosome]
        
        tbl_vrpipe = (etl.fromvcf(vrpipe_vcf_fn)
            .sort(['CHROM', 'POS'])
            .vcfmeltsamples()
            .vcfunpackcall()
            .unpackdict('INFO', samplesize=1000000)
            .selecteq('SAMPLE', isolate_code)
            .selectnotnone('GT')
            .selectne('GT', '0/0')
            .rename('FILTER', 'all_FILTERS')
            .addfield('FILTER', lambda rec: determine_filter(rec))
#             .rename('ALT', 'all_ALTs')
#             .addfield('ALT', lambda rec: str(rec[4][0]))
            .addfield('Truth_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
#             .cut(['CHROM', 'POS', 'REF', 'ALT', 'RegionType', 'Truth_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
#                   'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID'])
            .rename(vrpipe_renames)
#             .addfieldusingcontext('Truth_distance_previous', upstream)
#             .addfieldusingcontext('Truth_distance_next', downstream)
#             .addfield('Truth_distance_nearest', lambda rec: nearest(rec['Truth_distance_previous'], rec['Truth_distance_next']))
            .addfield('Truth_mode', determine_mode)
        )
        etl.topickle(tbl_vrpipe, tbl_vrpipe_cache_fn)
    else:
        tbl_vrpipe = etl.frompickle(tbl_vrpipe_cache_fn)
        
    tbl_vrpipe.display()
        
    if not os.path.exists(tbl_gatk2_cache_fn) or rewrite:
        gatk_vcf_fn = "%s.%s.SNP_INDEL.%s.vcf.gz" % (
            os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
            'Pf3D7july2015',
            chromosome
        )
        
        tbl_gatk2 = (etl.fromvcf(gatk_vcf_fn)
            .sort(['CHROM', 'POS'])
            .vcfmeltsamples()
            .vcfunpackcall()
            .unpackdict('INFO', samplesize=1000000)
            .selecteq('SAMPLE', ox_code)
            .selectnotnone('GT')
            .selectne('GT', '0/0')
            .rename('FILTER', 'all_FILTERS')
            .addfield('FILTER', lambda rec: determine_filter(rec))
#             .rename('ALT', 'all_ALTs')
#             .addfield('ALT', lambda rec: str(rec[4][int(rec[11][2])-1]))
            .addfield('gatk_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
#             .cut(['CHROM', 'POS', 'REF', 'ALT', 'RegionType', 'GATK_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
#                   'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID',
#                   'SAMPLE', 'AD', 'DP', 'GQ', 'GT', 'VQSLOD', 'culprit'])
            .rename(gatk_renames)
#             .select('ALT', lambda rec: str(rec) != '<*:DEL>')
#             .addfield('temp_REF', lambda rec: rec['REF'])
#             .addfield('temp_ALT', lambda rec: rec['ALT'])
#             .convert('REF', lambda ref, rec: trim_variants(rec['temp_REF'], rec['temp_ALT'])[0], pass_row=True)
#             .convert('ALT', lambda ref, rec: trim_variants(rec['temp_REF'], rec['temp_ALT'])[1], pass_row=True)
#             .addfieldusingcontext('GATK_distance_previous', upstream)
#             .addfieldusingcontext('GATK_distance_next', downstream)
#             .addfield('GATK_distance_nearest', lambda rec: nearest(rec['GATK_distance_previous'], rec['GATK_distance_next']))
            .addfield('GATK_mode', determine_mode)
        )
        etl.topickle(tbl_gatk2, tbl_gatk2_cache_fn)
    else:
        tbl_gatk2 = etl.frompickle(tbl_gatk2_cache_fn)

    if not os.path.exists(tbl_vrpipe_comparison_cache_fn) or rewrite:
#         truth_coords_fn = "%s/Pf%s.filter.coords" % (
#             '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/SNP/',
#             isolate_code
#         )
#         assembly_chrom = "Pf%s_%s" % (isolate_code, chromosome[6:8])
#         tbl_coords = (etl.fromtsv(truth_coords_fn)
#             .pushheader(['ref_start', 'ref_end', 'assembly_start', 'assembly_end', 'ref_length', 'assembly_length', 'identity',
#                  'ref_chrom_length', 'assembly_chrom_length', 'ref_pct', 'assembly_pct', 'ref_chrom', 'assembly_chrom',
#                  'ignore'])
#             .selecteq('ref_chrom', chromosome)
#             .selecteq('assembly_chrom', assembly_chrom)
#             .cut(['ref_start', 'ref_end', 'assembly_chrom'])
#             .convertnumbers()
#         )
#         lkp = etl.intervallookup(tbl_coords, 'ref_start', 'ref_end', include_stop=True, value='assembly_chrom')
        
        tbl_vrpipe_comparison = (tbl_vrpipe
#             .outerjoin(tbl_gatk2, key=['CHROM', 'POS', 'REF', 'ALT'])
            .outerjoin(tbl_gatk2, key=['CHROM', 'POS'])
            .addfield('IsInVRPIPE', lambda rec: rec['vrpipe_RegionType'] is not None)
            .addfield('IsInGATK', lambda rec: rec['gatk_RegionType'] is not None)
#             .addfield('RegionType', lambda rec: determine_consensus_region(rec, 4, 15))
#             .addfield('IsCoding', lambda rec: determine_consensus_region(rec, 5, 16))
#             .addfield('Effect', lambda rec: determine_consensus_region(rec, 7, 18))
#             .addfield('AAChange', lambda rec: determine_consensus_region(rec, 6, 17))
#             .addfield('Gene', lambda rec: determine_consensus_region(rec, 8, 19))
#             .addfield('Impact', lambda rec: determine_consensus_region(rec, 9, 20))
#             .addfield('Transcript', lambda rec: determine_consensus_region(rec, 10, 21))
#             .addfield('mode', lambda rec: determine_consensus_region(rec, 14, 34))
#             .addfield('IsAccessible', lambda rec: len(lkp.search(rec['POS'])) > 0)
#             .cut(['CHROM', 'POS', 'REF', 'ALT', 'IsInTruth', 'IsInGATK', 'VQSLOD', 'RegionType', 'IsCoding', 'Effect', 'AAChange', 'Gene', 'Impact', 'Transcript',
#                   'SAMPLE', 'AD', 'DP', 'GQ', 'GT', 'culprit',
#                   'Truth_distance_nearest', 'GATK_distance_nearest', 'IsAccessible', 'mode'])
#             .sort('VQSLOD', reverse=True)
        )
        etl.topickle(tbl_vrpipe_comparison, tbl_vrpipe_comparison_cache_fn)
    else:
        tbl_vrpipe_comparison = etl.frompickle(tbl_vrpipe_comparison_cache_fn)
        
    return(tbl_vrpipe_comparison)

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3 = assess_vrpipe(rewrite=True)

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3.valuecounts('IsInVRPIPE', 'IsInGATK')

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3.valuecounts('IsInVRPIPE', 'IsInGATK')

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3.valuecounts('IsInVRPIPE', 'IsInGATK')

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3.select(lambda rec: rec['IsInVRPIPE']==True and rec['IsInGATK']==False)

# <codecell>

tbl_vrpipe_comparison_7g8_Pf3D7_01_v3.select(lambda rec: rec['IsInVRPIPE']==False and rec['IsInGATK']==True)

# <codecell>

(tbl_vrpipe_comparison_7g8_Pf3D7_01_v3
    .select(lambda rec: rec['IsInVRPIPE']==False and rec['IsInGATK']==True)
    .valuecounts('gatk_FILTER')
)

# <codecell>

(tbl_vrpipe_comparison_7g8_Pf3D7_01_v3
    .select(lambda rec: rec['IsInVRPIPE']==False and rec['IsInGATK']==True)
    .valuecounts('gatk_GT')
    .displayall()
)

# <codecell>

(tbl_vrpipe_comparison_7g8_Pf3D7_01_v3
    .select(lambda rec: rec['IsInVRPIPE']==True and rec['IsInGATK']==False)
    .valuecounts('vrpipe_FILTER')
)

# <codecell>

def plotROC(isolate_code='7G8', RegionType='Core', IsCoding=True, mode='SNP',
            min_Truth_distance_nearest=0, min_GATK_distance_nearest=0, sort_by=None, sort_reverse=True,
            annotate_values=[0, 2, 4, 6, 8], show_annotations=False, ymax=1.0, xmax=None, ax=None):
    
    # set up axes
    if ax is None:
        fig, ax = plt.subplots()

#     # set up figure
#     fig = figure(figsize=(8, 8))

#     # set up roc plot
#     ax = fig.add_subplot(1, 1, 1)
    
    if(sort_by is None):
        array_to_plot = (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('IsInGATK', True)
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
            .select(lambda rec: rec['GATK_distance_nearest'] is None or rec['GATK_distance_nearest'] >= min_GATK_distance_nearest)
            .cut(['IsInGATK', 'IsInTruth', 'VQSLOD'])
        ).toarray()
        annotate_indices = np.array([np.sum(array_to_plot['VQSLOD'] >= value) for value in annotate_values])
        
    else:
        array_to_plot = (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('IsInGATK', True)
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
            .select(lambda rec: rec['GATK_distance_nearest'] is None or rec['GATK_distance_nearest'] >= min_GATK_distance_nearest)
            .sort(sort_by, reverse=sort_reverse)
            .cut(['IsInGATK', 'IsInTruth', sort_by])
        ).toarray()
        annotate_indices = np.array([np.sum(array_to_plot[sort_by] >= value) for value in annotate_values])

    TPs = np.cumsum(array_to_plot['IsInTruth'])
    FPs = np.cumsum(np.logical_not(array_to_plot['IsInTruth']))
    
    number_of_positives = len(array_to_plot)
    number_of_true_variants = len(
        (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('IsInTruth', True)
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
            .select(lambda rec: rec['GATK_distance_nearest'] is None or rec['GATK_distance_nearest'] >= min_GATK_distance_nearest)
#             .selectge('GATK_distance_nearest', min_GATK_distance_nearest)
            .cut(['IsInGATK', 'IsInTruth', 'VQSLOD'])
            .data()
        )
    )
    
    print(isolate_code, RegionType, IsCoding, mode, number_of_positives, number_of_true_variants)
    
    sensitivity = TPs / number_of_true_variants
    FDR = FPs / number_of_positives
    
    
    ax.plot(FDR, sensitivity, label=isolate_code)
    if show_annotations:
        for i, annotate_index in enumerate(annotate_indices):
            ax.annotate('%s (%4.1f%% sens, %4.1f%% FDR)' % (annotate_values[i], sensitivity[annotate_index]*100, FDR[annotate_index]*100),
                        (FDR[annotate_index], sensitivity[annotate_index]), xycoords='data', xytext=(10, -10),
                        textcoords='offset points')
    if ymax is not None:
        ax.set_ylim([0, ymax])
    if xmax is not None:
        ax.set_xlim([0, xmax])
            
    return(ax)


# <headingcell level=1>

# Process all data

# <codecell>

rewrite=False

isolate_codes = list(tbl_samples_to_process.values('Isolate code'))
# isolate_codes = ['7G8']
chromosomes = ['Pf3D7_%02d_v3' % chrom for chrom in range(1,15)]

tbl_comparisons = collections.OrderedDict()
tbl_whole_genome_comparisons = collections.OrderedDict()

tbl_comparison_all_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_comparison_all")

# if not os.path.exists(tbl_comparison_all_cache_fn) or rewrite:
for isolate_code in isolate_codes:
    tbl_comparisons[isolate_code] = collections.OrderedDict()
    print(isolate_code)
    for chromosome in chromosomes:
        tbl_comparisons[isolate_code][chromosome] = assess_validation(isolate_code, chromosome, rewrite=rewrite)

    tbl_comparison_sample_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_comparison_%s" % (isolate_code))
    if not os.path.exists(tbl_comparison_sample_cache_fn) or rewrite:
        tbl_whole_genome_comparisons[isolate_code] = (
            tbl_comparisons[isolate_code]['Pf3D7_01_v3']
            .cat(tbl_comparisons[isolate_code]['Pf3D7_02_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_03_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_04_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_05_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_06_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_07_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_08_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_09_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_10_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_11_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_12_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_13_v3'])
            .cat(tbl_comparisons[isolate_code]['Pf3D7_14_v3'])
            .sort('VQSLOD', reverse=True)
        )
        etl.topickle(tbl_whole_genome_comparisons[isolate_code], tbl_comparison_sample_cache_fn)
    else:
        tbl_whole_genome_comparisons[isolate_code] = etl.frompickle(tbl_comparison_sample_cache_fn)

#     etl.topickle(tbl_whole_genome_comparisons, tbl_comparison_all_cache_fn)
# else:
#     tbl_whole_genome_comparisons = etl.frompickle(tbl_comparison_all_cache_fn)


# <headingcell level=1>

# Analysis - cross sample check

# <codecell>

for isolate_code in isolate_codes:
#     if isolate_code != 'GN01':
    print(isolate_code)
    (tbl_whole_genome_comparisons[isolate_code]
        .selecteq('IsInGATK', True)
        .selecteq('mode', 'SNP')
        .selecteq('RegionType', 'Core')
        .selecteq('IsCoding', True)
        .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
        .displayall()
    )

# <headingcell level=1>

# Analysis - variant types check

# <codecell>

for RegionType in ['Core', 'SubtelomericHypervariable', 'SubtelomericRepeat', 'Centromere']:
    for mode in ['SNP', 'INDEL']:
        for IsCoding in [True, False]:
            print(RegionType, mode, IsCoding)
            (tbl_whole_genome_comparisons['7G8']
#                 .selecteq('IsInGATK', True)
                .selecteq('mode', mode)
                .selecteq('RegionType', RegionType)
                .selecteq('IsCoding', IsCoding)
                .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
                .displayall()
            )

# <codecell>

aggregation = collections.OrderedDict()
aggregation['#Truth'] = 'IsInTruth', sum
aggregation['#GATK'] = 'IsInGATK', sum
aggregation['#GATK accessible'] = ('IsInGATK', 'IsAccessible'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
aggregation['#TP'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
aggregation['#FP'] = ('IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([x[0] and not x[1] and x[2] for x in list(rec)])
aggregation['#FN'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([not x[0] and x[1] for x in list(rec)])
aggregation['#TPmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and x[3] for x in list(rec)])
aggregation['#FPmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and not x[3] and x[4] for x in list(rec)])
aggregation['#FNmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and not x[2] and x[3] for x in list(rec)])

(tbl_whole_genome_comparisons['7G8']
    .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation)
    .convertnumbers()
    .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
    .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
    .addfield('%TPmod3', lambda rec: 0.0 if rec['#TP'] == 0 else round(rec['#TPmod3'] / rec['#TP'], 3))
    .addfield('%FPmod3', lambda rec: 0.0 if rec['#FP'] == 0 else round(rec['#FPmod3'] / rec['#FP'], 3))
    .addfield('%FNmod3', lambda rec: 0.0 if rec['#FN'] == 0 else round(rec['#FNmod3'] / rec['#FN'], 3))
).displayall()

# <codecell>

tbl_whole_genome_comparisons['7G8'].selecteq('mode', 'INDEL').selecteq('IsCoding', True)

# <codecell>

def is_transition(rec):
    return (((rec[0] == 'A') and (rec[1] == 'G')) 
             or ((rec[0] == 'G') and (rec[1] == 'A')) 
             or ((rec[0] == 'C') and (rec[1] == 'T')) 
             or ((rec[0] == 'T') and (rec[1] == 'C')))

def is_transversion(rec):
    return (((rec[0] == 'A') and ((rec[1] == 'T') or (rec[1] == 'C')))
             or ((rec[0] == 'G') and ((rec[1] == 'T') or (rec[1] == 'C')))
             or ((rec[0] == 'C') and ((rec[1] == 'A') or (rec[1] == 'G')))
             or ((rec[0] == 'T') and ((rec[1] == 'A') or (rec[1] == 'G'))))

def calc_titv(rec_list, variant_type='TP'):
    if variant_type=='TP':
        new_list = [x for x in rec_list if x[2] and x[3]]
    elif variant_type=='FP':
        new_list = [x for x in rec_list if x[2] and not x[3] and x[4]]
    if variant_type=='FN':
        new_list = [x for x in rec_list if not x[2] and x[3]]
    ti = np.count_nonzero([is_transition(rec) for rec in new_list])
    tv = np.count_nonzero([is_transversion(rec) for rec in new_list])
    if tv == 0:
        titv = 0.0
    else:
        titv = 1.0*ti/tv
#     print(ti, tv, titv)
    return(titv)

# <codecell>

def calc_variant_summary_table(isolate_code='7G8', min_Truth_distance_nearest=0, min_VQSLOD=-99999):
    aggregation = collections.OrderedDict()
    aggregation['#Truth'] = 'IsInTruth', sum
    aggregation['#GATK'] = 'IsInGATK', sum
    aggregation['#GATK accessible'] = ('IsInGATK', 'IsAccessible'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#TP'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#FP'] = ('IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([x[0] and not x[1] and x[2] for x in list(rec)])
    aggregation['#FN'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([not x[0] and x[1] for x in list(rec)])
    aggregation['TiTv TP'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'TP')
    aggregation['TiTv FP'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FP')
    aggregation['TiTv FN'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FN')
    aggregation['#TPmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and x[3] for x in list(rec)])
    aggregation['#FPmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and not x[3] and x[4] for x in list(rec)])
    aggregation['#FNmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and not x[2] and x[3] for x in list(rec)])

    tbl_variant_summary = (tbl_whole_genome_comparisons[isolate_code]
        .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
        .select(lambda rec: rec['VQSLOD'] is None or rec['VQSLOD'] >= min_VQSLOD)
        .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation)
#         .convertnumbers()
        .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
        .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
        .addfield('%TPmod3', lambda rec: 0.0 if rec['#TP'] == 0 else round(rec['#TPmod3'] / rec['#TP'], 3))
        .addfield('%FPmod3', lambda rec: 0.0 if rec['#FP'] == 0 else round(rec['#FPmod3'] / rec['#FP'], 3))
        .addfield('%FNmod3', lambda rec: 0.0 if rec['#FN'] == 0 else round(rec['#FNmod3'] / rec['#FN'], 3))
    )
    
    return(tbl_variant_summary)

# <codecell>

calc_variant_summary_table().displayall()

# <codecell>

rewrite=False

tbl_variant_summary = collections.OrderedDict()
for isolate_code in isolate_codes:
    tbl_variant_summary[isolate_code] = collections.OrderedDict()
    for min_Truth_distance_nearest in [0, 10, 100]:
#     for min_Truth_distance_nearest in [0]:
        tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)] = collections.OrderedDict()
        for min_VQSLOD in [-99999, 0, 2, 6]:
            print(isolate_code, min_Truth_distance_nearest, min_VQSLOD)
            tbl_variant_summary_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_variant_summary_%s_%d_%d" % (isolate_code, min_Truth_distance_nearest, min_VQSLOD))
            if not os.path.exists(tbl_variant_summary_cache_fn) or rewrite:
                tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)][str(min_VQSLOD)] = calc_variant_summary_table(isolate_code, min_Truth_distance_nearest, min_VQSLOD)
                etl.topickle(tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)][str(min_VQSLOD)], tbl_variant_summary_cache_fn)
            else:
                tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)][str(min_VQSLOD)] = etl.frompickle(tbl_variant_summary_cache_fn)

# <codecell>

rewrite = False

aggregation_allsamples = collections.OrderedDict()
aggregation_allsamples['#Truth'] = '#Truth', sum
aggregation_allsamples['#GATK'] = '#GATK', sum
aggregation_allsamples['#GATK accessible'] = '#GATK accessible', sum
aggregation_allsamples['#TP'] = '#TP', sum
aggregation_allsamples['#FP'] = '#FP', sum
aggregation_allsamples['#FN'] = '#FN', sum
aggregation_allsamples['#TPmod3'] = '#TPmod3', sum
aggregation_allsamples['#FPmod3'] = '#FPmod3', sum
aggregation_allsamples['#FNmod3'] = '#FNmod3', sum

tbl_allsamples_summary = collections.OrderedDict()
for min_Truth_distance_nearest in [0, 10, 100]:
# for min_Truth_distance_nearest in [0]:
    print(min_Truth_distance_nearest)
    tbl_allsamples_summary_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_allsamples_summary_%d" % (min_Truth_distance_nearest))
    if not os.path.exists(tbl_variant_summary_cache_fn) or rewrite:
        tbl_allsamples_summary[str(min_Truth_distance_nearest)] = (
            tbl_variant_summary['7G8'][str(min_Truth_distance_nearest)]
            .cat(tbl_variant_summary['GB4'][str(min_Truth_distance_nearest)])
            .cat(tbl_variant_summary['KE01'][str(min_Truth_distance_nearest)])
            .cat(tbl_variant_summary['KH02'][str(min_Truth_distance_nearest)])
            .cat(tbl_variant_summary['GN01'][str(min_Truth_distance_nearest)])
            .update('IsCoding', 'Coding', where=lambda rec: rec['IsCoding'] == 1)
            .update('IsCoding', 'Non-coding', where=lambda rec: rec['IsCoding'] != 'Coding')
            .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation_allsamples)
            .convertnumbers()
            .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
            .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
            .addfield('%TPmod3', lambda rec: 0.0 if rec['#TP'] == 0 else round(rec['#TPmod3'] / rec['#TP'], 3))
            .addfield('%FPmod3', lambda rec: 0.0 if rec['#FP'] == 0 else round(rec['#FPmod3'] / rec['#FP'], 3))
            .addfield('%FNmod3', lambda rec: 0.0 if rec['#FN'] == 0 else round(rec['#FNmod3'] / rec['#FN'], 3))
            .sort('Unfilt sensitivity', reverse=True)
        )
        etl.topickle(tbl_allsamples_summary[str(min_Truth_distance_nearest)], tbl_allsamples_summary_cache_fn)
    else:
        tbl_allsamples_summary[str(min_Truth_distance_nearest)] = etl.frompickle(tbl_allsamples_summary_cache_fn)

# <codecell>

rewrite = True

aggregation_allsamples = collections.OrderedDict()
aggregation_allsamples['#Truth'] = '#Truth', sum
aggregation_allsamples['#GATK'] = '#GATK', sum
aggregation_allsamples['#GATK accessible'] = '#GATK accessible', sum
aggregation_allsamples['#TP'] = '#TP', sum
aggregation_allsamples['#FP'] = '#FP', sum
aggregation_allsamples['#FN'] = '#FN', sum
# aggregation_allsamples['TiTv TP'] = 'TiTv TP', np.nanmean
# aggregation_allsamples['TiTv FP'] = 'TiTv FP', np.nanmean
# aggregation_allsamples['TiTv FN'] = 'TiTv FN', np.nanmean
aggregation_allsamples['#TPmod3'] = '#TPmod3', sum
aggregation_allsamples['#FPmod3'] = '#FPmod3', sum
aggregation_allsamples['#FNmod3'] = '#FNmod3', sum

tbl_allsamples_filtered_summary = collections.OrderedDict()
# for min_Truth_distance_nearest in [0, 10, 100]:
for min_Truth_distance_nearest in [0]:
    tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)] = collections.OrderedDict()
    for min_VQSLOD in [-99999, 0, 2, 6]:
        print(min_Truth_distance_nearest, min_VQSLOD)
#         tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)] = collections.OrderedDict()
        tbl_allsamples_filtered_summary_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_allsamples_filtered_summary_%d_%d" % (min_Truth_distance_nearest, min_VQSLOD))
        if not os.path.exists(tbl_allsamples_filtered_summary_cache_fn) or rewrite:
            tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)][str(min_VQSLOD)] = (
                tbl_variant_summary['7G8'][str(min_Truth_distance_nearest)][str(min_VQSLOD)]
                .cat(tbl_variant_summary['GB4'][str(min_Truth_distance_nearest)][str(min_VQSLOD)])
                .cat(tbl_variant_summary['KE01'][str(min_Truth_distance_nearest)][str(min_VQSLOD)])
                .cat(tbl_variant_summary['KH02'][str(min_Truth_distance_nearest)][str(min_VQSLOD)])
                .cat(tbl_variant_summary['GN01'][str(min_Truth_distance_nearest)][str(min_VQSLOD)])
#                 .select(lambda rec: rec['VQSLOD'] is None or rec['VQSLOD'] >= VQSLOD_threshold)
                .update('IsCoding', 'Coding', where=lambda rec: rec['IsCoding'] == 1)
                .update('IsCoding', 'Non-coding', where=lambda rec: rec['IsCoding'] != 'Coding')
                .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation_allsamples)
                .convertnumbers()
                .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
                .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
                .addfield('%TPmod3', lambda rec: 0.0 if rec['#TP'] == 0 else round(rec['#TPmod3'] / rec['#TP'], 3))
                .addfield('%FPmod3', lambda rec: 0.0 if rec['#FP'] == 0 else round(rec['#FPmod3'] / rec['#FP'], 3))
                .addfield('%FNmod3', lambda rec: 0.0 if rec['#FN'] == 0 else round(rec['#FNmod3'] / rec['#FN'], 3))
                .sort('Unfilt sensitivity', reverse=True)
            )
            etl.topickle(tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)][str(min_VQSLOD)], tbl_allsamples_filtered_summary_cache_fn)
        else:
            tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)][str(min_VQSLOD)] = etl.frompickle(tbl_allsamples_filtered_summary_cache_fn)

# <codecell>

tbl_allsamples_filtered_summary['0']['-99999'].selecteq('RegionType', 'Core')

# <codecell>

tbl_allsamples_filtered_summary['0']['2'].selecteq('RegionType', 'Core')

# <codecell>

list(tbl_allsamples_filtered_summary['0'].keys())

# <codecell>

# tbl_variant_summary['7G8']['0']

# <codecell>

tbl_variant_summary['7G8']['0']['-99999'].selecteq('RegionType', 'Core').displayall()

# <codecell>

tbl_variant_summary['7G8']['0']['0'].selecteq('RegionType', 'Core').displayall()

# <codecell>

tbl_variant_summary['7G8']['0']['2'].selecteq('RegionType', 'Core').displayall()

# <codecell>

tbl_variant_summary['7G8']['0']['6'].selecteq('RegionType', 'Core').displayall()

# <codecell>

tbl_allsamples_summary['0'].selecteq('RegionType', 'Core').selecteq('mode', 'INDEL')

# <codecell>

for min_Truth_distance_nearest in [0, 10, 100]:
    (tbl_allsamples_summary[str(min_Truth_distance_nearest)]
        .selecteq('RegionType', 'Core')
        .selecteq('mode', 'INDEL')
        .cut(['IsCoding', '#TP', '#FP', '#FN'])
        .melt('IsCoding')
        .rename('variable', 'Variant set')
        .convert('Variant set', lambda val: val[1:])
        .rename('value', '# variants')
    ).annex(tbl_allsamples_summary[str(min_Truth_distance_nearest)]
        .selecteq('RegionType', 'Core')
        .selecteq('mode', 'INDEL')
        .cut(['IsCoding', '#TPmod3', '#FPmod3', '#FNmod3'])
        .melt('IsCoding')
        .cut('value')
        .rename('value', '# mod 3')
    ).annex(tbl_allsamples_summary[str(min_Truth_distance_nearest)]
        .selecteq('RegionType', 'Core')
        .selecteq('mode', 'INDEL')
        .cut(['IsCoding', '%TPmod3', '%FPmod3', '%FNmod3'])
        .melt('IsCoding')
        .cut('value')
        .rename('value', '% mod 3')
    ).sort('IsCoding').toxlsx("%s/mod3_indels_Core_%d.xlsx" % (PLOTS_DIR, min_Truth_distance_nearest))


# <codecell>

df_variant_summary = collections.OrderedDict()
for min_Truth_distance_nearest in [0, 10, 100]:
    df_variant_summary[str(min_Truth_distance_nearest)] = (
        tbl_variant_summary['7G8'][str(min_Truth_distance_nearest)].addfield('Sample', '7G8')
        .cat(tbl_variant_summary['GB4'][str(min_Truth_distance_nearest)].addfield('Sample', 'GB4'))
        .cat(tbl_variant_summary['KE01'][str(min_Truth_distance_nearest)].addfield('Sample', 'KE01'))
        .cat(tbl_variant_summary['KH02'][str(min_Truth_distance_nearest)].addfield('Sample', 'KH02'))
        .cat(tbl_variant_summary['GN01'][str(min_Truth_distance_nearest)].addfield('Sample', 'GN01'))
        .selecteq('RegionType', 'Core')
        .selectin('IsCoding', [0, 1])
        .update('IsCoding', 'Coding', where=lambda rec: rec['IsCoding'] == 1)
        .update('IsCoding', 'Non-coding', where=lambda rec: rec['IsCoding'] != 'Coding')
        .addfield('Variant type', lambda rec: "%s %s" % (rec['IsCoding'], rec['mode']))
        .sort('IsCoding')
        .sort('mode', reverse=True)
        .todataframe()
    )

# <codecell>

for min_Truth_distance_nearest in [0, 10, 100]:
    fig = figure(figsize=(8, 6))
    ax = fig.add_subplot(2, 1, 1)
    g = sns.barplot(x="Variant type", y="Unfilt sensitivity", hue="Sample", ci=None,
                    data=df_variant_summary[str(min_Truth_distance_nearest)])
    g.legend(loc="lower center")
    g.set_ylabel("Unfiltered sensitivity")
    g.set_xlabel("")
    ax = fig.add_subplot(2, 1, 2)
    g = sns.barplot(x="Variant type", y="Unfilt FDR", hue="Sample", ci=None,
                    data=df_variant_summary[str(min_Truth_distance_nearest)])
    g.legend(loc="lower center")
    g.set_ylabel("Unfiltered FDR")
    g.set_xlabel("")
    fig.tight_layout()

# <codecell>

tbl_variant_summary['7G8']['0'].selecteq('RegionType', 'Core')

# <codecell>

tbl_variant_summary['7G8']['0']

# <codecell>

for min_Truth_distance_nearest in [0, 10, 100]:
    (tbl_allsamples_summary[str(min_Truth_distance_nearest)]
        .selecteq('RegionType', 'Core')
        .cutout('RegionType')
        .sort('IsCoding')
        .sort('mode', reverse=True)
        .cut([0, 1, 2, 3, 4, 5, 6, 7, 11, 12])
    #     .displayall()
        .toxlsx("%s/all_samples_summary_Core_%d.xlsx" % (PLOTS_DIR, min_Truth_distance_nearest))
    )

# <codecell>

# for min_Truth_distance_nearest in [0, 10, 100]:
for min_Truth_distance_nearest in [0]:
    for min_VQSLOD in [-99999, 0, 2, 6]:
        (tbl_allsamples_filtered_summary[str(min_Truth_distance_nearest)][str(min_VQSLOD)]
            .selecteq('RegionType', 'Core')
            .cutout('RegionType')
            .sort('IsCoding')
            .sort('mode', reverse=True)
            .cut([0, 1, 2, 3, 4, 5, 6, 7, 11, 12])
        #     .displayall()
            .toxlsx("%s/tbl_allsamples_filtered_summary_Core_%d_%d.xlsx" % (PLOTS_DIR, min_Truth_distance_nearest, min_VQSLOD))
        )

# <codecell>

(tbl_allsamples_summary['10']
    .selecteq('RegionType', 'Core')
    .cutout('RegionType')
    .sort('IsCoding')
    .sort('mode', reverse=True)
    .displayall()
)

# <codecell>

(tbl_allsamples_summary['100']
    .selecteq('RegionType', 'Core')
    .cutout('RegionType')
    .sort('IsCoding')
    .sort('mode', reverse=True)
    .displayall()
)

# <codecell>

(tbl_variant_summary['7G8']['0']
    .selecteq('RegionType', 'Core')
    .addfield('Coding', lambda rec: 'Coding' if rec['IsCoding'] == 1 else 'Non-coding')
    .cutout('RegionType')
    .cutout('IsCoding')
    .sort('Coding')
    .sort('mode', reverse = True)
    .displayall()
)

(tbl_variant_summary['7G8']['10']
    .selecteq('RegionType', 'Core')
    .addfield('Coding', lambda rec: 'Coding' if rec['IsCoding'] == 1 else 'Non-coding')
    .cutout('RegionType')
    .cutout('IsCoding')
    .sort('Coding')
    .sort('mode', reverse = True)
    .displayall()
)

(tbl_variant_summary['7G8']['100']
    .selecteq('RegionType', 'Core')
    .addfield('Coding', lambda rec: 'Coding' if rec['IsCoding'] == 1 else 'Non-coding')
    .cutout('RegionType')
    .cutout('IsCoding')
    .sort('Coding')
    .sort('mode', reverse = True)
    .displayall()
)

# <codecell>


# <codecell>

isolate_code = '7G8'
mode = 'SNP'
RegionType = 'Core'
IsCoding = True

df_vqslod_gt = collections.OrderedDict()
for isolate_code in isolate_codes:
    df_vqslod_gt[isolate_code] = collections.OrderedDict()
    for mode in ['SNP', 'INDEL']:
        df_vqslod_gt[isolate_code][mode] = collections.OrderedDict()
        for IsCoding in [True, False]:
            df_vqslod_gt[isolate_code][mode][IsCoding] = (tbl_whole_genome_comparisons[isolate_code]
                .selecteq('IsInGATK', True)
                .selecteq('mode', mode)
                .selecteq('RegionType', RegionType)
                .selecteq('IsCoding', IsCoding)
                .selecteq('IsAccessible', True)
                .addfield('call_type', lambda rec: 'Hom' if rec['GT'][0] == rec['GT'][2] else 'Het')
                .selectge('VQSLOD', -20)
                .selectle('VQSLOD', 20)
                .cut(['VQSLOD', 'call_type', 'GT'])
            ).todataframe()

# <codecell>

(tbl_whole_genome_comparisons[isolate_code]
    .selecteq('IsInGATK', True)
    .selecteq('mode', mode)
    .selecteq('RegionType', RegionType)
    .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .addfield('call_type', lambda rec: 'Hom' if rec['GT'][0] == rec['GT'][2] else 'Het')
    .valuecounts('GT')
    .displayall()
)
(tbl_whole_genome_comparisons[isolate_code]
    .selecteq('IsInGATK', True)
    .selecteq('mode', mode)
    .selecteq('RegionType', RegionType)
    .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .addfield('call_type', lambda rec: 'Hom' if rec['GT'][0] == rec['GT'][2] else 'Het')
    .valuecounts('call_type')
    .displayall()
)
(tbl_whole_genome_comparisons[isolate_code]
    .selecteq('IsInGATK', True)
    .selecteq('mode', mode)
    .selecteq('RegionType', RegionType)
    .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .addfield('VQSLODgt8', lambda rec: rec['VQSLOD'] >= 8.0)
    .valuecounts('VQSLODgt8')
    .displayall()
)

# <codecell>

upper_limit={'SNP':8, 'INDEL':4}
lower_limit={'SNP':2, 'INDEL':2}

for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate([True, False]):
        fig = figure(figsize=(8, 4))
        for k, isolate_code in enumerate(isolate_codes):
            ax = fig.add_subplot(1, 5, k+1)
            g = sns.violinplot(x="call_type", y="VQSLOD", data=df_vqslod_gt[isolate_code][mode][IsCoding], ax=ax)
            ax.set_ylim([-15, 15])
            ax.set_title(isolate_code)
            ax.axhline(upper_limit[mode], color='green')
            ax.axhline(lower_limit[mode], color='red')
        fig.tight_layout()


# <headingcell level=1>

# Plot ROC curves

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, xmax=0.15, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, min_Truth_distance_nearest=10, min_GATK_distance_nearest=10,
                 xmax=0.15, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, min_Truth_distance_nearest=100, min_GATK_distance_nearest=100,
                 xmax=0.15, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, IsCoding=False, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations,
                 min_Truth_distance_nearest=100, min_GATK_distance_nearest=100,
                 xmax=0.45, IsCoding=False, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, mode='INDEL', IsCoding=True, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations,
                 min_Truth_distance_nearest=100, min_GATK_distance_nearest=100,
                 xmax=0.40, mode='INDEL', IsCoding=True, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations, mode='INDEL', IsCoding=False, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotROC(isolate_code, show_annotations=show_annotations,
                 min_Truth_distance_nearest=100, min_GATK_distance_nearest=100,
                 xmax=0.35, mode='INDEL', IsCoding=False, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

for isolate_code in isolate_codes:
    plotROC(isolate_code, sort_by='GQ')

# <codecell>

for isolate_code in isolate_codes:
    plotROC(isolate_code, sort_by='GATK_distance_nearest')

# <codecell>

for RegionType in ['Core', 'SubtelomericHypervariable']:
    for mode in ['SNP', 'INDEL']:
        for IsCoding in [True, False]:
            print(RegionType, mode, IsCoding)
            plotROC('7G8', RegionType, IsCoding, mode)

# <codecell>

for RegionType in ['Core', 'SubtelomericHypervariable']:
    for mode in ['SNP']:
        for IsCoding in [True, False]:
#             print(RegionType, mode, IsCoding)
            plotROC('7G8', RegionType, IsCoding, mode)

# <codecell>

for RegionType in ['Core', 'SubtelomericHypervariable']:
    for mode in ['SNP']:
        for IsCoding in [True, False]:
            plotROC('7G8', RegionType, IsCoding, mode, min_Truth_distance_nearest=100)

# <codecell>

for RegionType in ['Core']:
    for mode in ['INDEL']:
        for IsCoding in [True, False]:
            print(RegionType, mode, IsCoding)
            plotROC('7G8', RegionType, IsCoding, mode)

# <codecell>

for RegionType in ['SubtelomericHypervariable']:
    for mode in ['INDEL']:
        for IsCoding in [True, False]:
            print(RegionType, mode, IsCoding)
            plotROC('7G8', RegionType, IsCoding, mode)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .selecteq('IsAccessible', True)
    .selecteq('IsInTruth', False)    
).display(20)

# <codecell>

tbl_comparisons['7G8']['Pf3D7_03_v3'].sort('POS').selectgt('POS', 540000).display(4)

# <codecell>

ref_dict['Pf3D7_03_v3'].seq[540020:540032]

# <codecell>

tbl_comparisons['7G8']['Pf3D7_13_v3'].sort('POS').selectgt('POS', 987000)

# <codecell>

tbl_comparisons['7G8']['Pf3D7_08_v3'].sort('POS').selectgt('POS', 242200).display(3)

# <codecell>

tbl_comparisons['7G8']['Pf3D7_14_v3'].sort('POS').selectgt('POS', 2296600).display(2)

# <codecell>

ref_dict['Pf3D7_14_v3'].seq[2296690:2296700]

# <codecell>

tbl_comparisons['7G8']['Pf3D7_14_v3'].sort('POS').selectgt('POS', 1031200).display(20)

# <codecell>

tbl_comparisons['7G8']['Pf3D7_03_v3'].sort('POS').selectgt('POS', 540000)

# <codecell>


