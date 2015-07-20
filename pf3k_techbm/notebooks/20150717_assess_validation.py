# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

!mkdir {os.path.join(CACHE_DIR, 'tbl_truth')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_gatk')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_comparison')}

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

def determine_consensus_region(rec, column1=4, column2=11):
    if rec[column1] is None:
        return rec[column2]
    elif rec[column2] is None:
        return rec[column1]
    elif rec[column1] == rec[column2]:
        return rec[column1]
    else:
        return "ERROR__%s__%s" % (rec[column1], rec[column2])
    

# <codecell>

def upstream(prv, cur, nxt):
    if prv is None:
        return None
    else:
        return cur.POS - prv.POS

def downstream(prv, cur, nxt):
    if nxt is None:
        return None
    else:
        return nxt.POS - cur.POS

def nearest(distance_previous, distance_next):
    if distance_previous is None:
        return(distance_next)
    elif distance_next is None:
        return(distance_previous)
    else:
        return(min(distance_previous, distance_next))

# <codecell>

def determine_mode(rec):
    if len(rec['REF']) == len(str(rec['ALT'])):
        return('SNP')
    else:
        return('INDEL')

# <codecell>

def trim_variants(ref, alt):
#     print(ref, alt)
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[0:-1]
        alt = alt[0:-1]
    return((ref, alt))
 
trim_variants('AATC', 'TC')[1]

# <codecell>

def assess_validation(isolate_code='7G8', chromosome='Pf3D7_01_v3', rewrite=False):
    ox_code = tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('ox_code')[0]
    
    tbl_truth_cache_fn = os.path.join(CACHE_DIR, 'tbl_truth', "tbl_truth_%s.%s" % (isolate_code, chromosome))
    tbl_gatk_cache_fn = os.path.join(CACHE_DIR, 'tbl_gatk', "tbl_gatk_%s.%s" % (isolate_code, chromosome))
    tbl_comparison_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_comparison_%s.%s" % (isolate_code, chromosome))
    
    if not os.path.exists(tbl_truth_cache_fn) or rewrite:
        truth_vcf_fn = "%s.%s.annotated.vcf.gz" % (
            tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('truth_vcf_filestem')[0],
            chromosome
        )
        
        tbl_truth = (etl.fromvcf(truth_vcf_fn)
            .unpackdict('INFO', samplesize=1000000)
            .rename('ALT', 'all_ALTs')
            .addfield('ALT', lambda rec: str(rec[4][0]))
            .addfield('Truth_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
            .cut(['CHROM', 'POS', 'REF', 'ALT', 'RegionType', 'Truth_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
                  'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID'])
            .rename('RegionType', 'Truth_RegionType')
            .rename('SNPEFF_AMINO_ACID_CHANGE', 'Truth_AAChange')
            .rename('SNPEFF_EFFECT', 'Truth_Effect')
            .rename('SNPEFF_GENE_NAME', 'Truth_Gene')
            .rename('SNPEFF_IMPACT', 'Truth_Impact')
            .rename('SNPEFF_TRANSCRIPT_ID', 'Truth_Transcript')
            .addfieldusingcontext('Truth_distance_previous', upstream)
            .addfieldusingcontext('Truth_distance_next', downstream)
            .addfield('Truth_distance_nearest', lambda rec: nearest(rec['Truth_distance_previous'], rec['Truth_distance_next']))
            .addfield('Truth_mode', determine_mode)
        )
        etl.topickle(tbl_truth, tbl_truth_cache_fn)
    else:
        tbl_truth = etl.frompickle(tbl_truth_cache_fn)
        
    if not os.path.exists(tbl_gatk_cache_fn) or rewrite:
        gatk_vcf_fn = "%s.%s.SNP_INDEL.%s.vcf.gz" % (
            os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
            'Pf3D7july2015',
            chromosome
        )
        
        tbl_gatk = (etl.fromvcf(gatk_vcf_fn)
            .sort(['CHROM', 'POS'])
            .vcfmeltsamples()
            .vcfunpackcall()
            .unpackdict('INFO', samplesize=1000000)
            .selecteq('SAMPLE', ox_code)
            .selectnotnone('GT')
            .selectne('GT', '0/0')
            .rename('ALT', 'all_ALTs')
            .addfield('ALT', lambda rec: str(rec[4][int(rec[11][2])-1]))
            .addfield('GATK_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
            .cut(['CHROM', 'POS', 'REF', 'ALT', 'RegionType', 'GATK_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
                  'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID',
                  'SAMPLE', 'AD', 'DP', 'GQ', 'GT', 'VQSLOD', 'culprit'])
            .rename('RegionType', 'GATK_RegionType')
            .rename('SNPEFF_AMINO_ACID_CHANGE', 'GATK_AAChange')
            .rename('SNPEFF_EFFECT', 'GATK_Effect')
            .rename('SNPEFF_GENE_NAME', 'GATK_Gene')
            .rename('SNPEFF_IMPACT', 'GATK_Impact')
            .rename('SNPEFF_TRANSCRIPT_ID', 'GATK_Transcript')
            .select('ALT', lambda rec: str(rec) != '<*:DEL>')
            .addfield('temp_REF', lambda rec: rec['REF'])
            .addfield('temp_ALT', lambda rec: rec['ALT'])
            .convert('REF', lambda ref, rec: trim_variants(rec['temp_REF'], rec['temp_ALT'])[0], pass_row=True)
            .convert('ALT', lambda ref, rec: trim_variants(rec['temp_REF'], rec['temp_ALT'])[1], pass_row=True)
            .addfieldusingcontext('GATK_distance_previous', upstream)
            .addfieldusingcontext('GATK_distance_next', downstream)
            .addfield('GATK_distance_nearest', lambda rec: nearest(rec['GATK_distance_previous'], rec['GATK_distance_next']))
            .addfield('GATK_mode', determine_mode)
        )
        etl.topickle(tbl_gatk, tbl_gatk_cache_fn)
    else:
        tbl_gatk = etl.frompickle(tbl_gatk_cache_fn)

    if not os.path.exists(tbl_comparison_cache_fn) or rewrite:
        truth_coords_fn = "%s/Pf%s.filter.coords" % (
            '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/SNP/',
            isolate_code
        )
        assembly_chrom = "Pf%s_%s" % (isolate_code, chromosome[6:8])
        tbl_coords = (etl.fromtsv(truth_coords_fn)
            .pushheader(['ref_start', 'ref_end', 'assembly_start', 'assembly_end', 'ref_length', 'assembly_length', 'identity',
                 'ref_chrom_length', 'assembly_chrom_length', 'ref_pct', 'assembly_pct', 'ref_chrom', 'assembly_chrom',
                 'ignore'])
            .selecteq('ref_chrom', chromosome)
            .selecteq('assembly_chrom', assembly_chrom)
            .cut(['ref_start', 'ref_end', 'assembly_chrom'])
            .convertnumbers()
        )
        lkp = etl.intervallookup(tbl_coords, 'ref_start', 'ref_end', include_stop=True, value='assembly_chrom')
        
        tbl_comparison = (tbl_truth
            .outerjoin(tbl_gatk, key=['CHROM', 'POS', 'REF', 'ALT'])
            .addfield('IsInTruth', lambda rec: rec['Truth_RegionType'] is not None)
            .addfield('IsInGATK', lambda rec: rec['GATK_RegionType'] is not None)
            .addfield('RegionType', lambda rec: determine_consensus_region(rec, 4, 15))
            .addfield('IsCoding', lambda rec: determine_consensus_region(rec, 5, 16))
            .addfield('Effect', lambda rec: determine_consensus_region(rec, 7, 18))
            .addfield('AAChange', lambda rec: determine_consensus_region(rec, 6, 17))
            .addfield('Gene', lambda rec: determine_consensus_region(rec, 8, 19))
            .addfield('Impact', lambda rec: determine_consensus_region(rec, 9, 20))
            .addfield('Transcript', lambda rec: determine_consensus_region(rec, 10, 21))
            .addfield('mode', lambda rec: determine_consensus_region(rec, 14, 34))
            .addfield('IsAccessible', lambda rec: len(lkp.search(rec['POS'])) > 0)
            .cut(['CHROM', 'POS', 'REF', 'ALT', 'IsInTruth', 'IsInGATK', 'VQSLOD', 'RegionType', 'IsCoding', 'Effect', 'AAChange', 'Gene', 'Impact', 'Transcript',
                  'SAMPLE', 'AD', 'DP', 'GQ', 'GT', 'culprit',
                  'Truth_distance_nearest', 'GATK_distance_nearest', 'IsAccessible', 'mode'])
            .sort('VQSLOD', reverse=True)
        )
        etl.topickle(tbl_comparison, tbl_comparison_cache_fn)
    else:
        tbl_comparison = etl.frompickle(tbl_comparison_cache_fn)
        
    return(tbl_comparison)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3 = assess_validation(rewrite=True)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3

# <codecell>

def plotROC(isolate_code='7G8', RegionType='Core', IsCoding=True, mode='SNP',
            min_Truth_distance_nearest=0, min_GATK_distance_nearest=0, sort_by=None, sort_reverse=True):
    
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
    
    plot(FDR, sensitivity)


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

(tbl_whole_genome_comparisons['7G8']
    .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation)
    .convertnumbers()
    .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
    .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
).displayall()

# <codecell>

def calc_variant_summary_table(isolate_code='7G8', min_Truth_distance_nearest=0):
    aggregation = collections.OrderedDict()
    aggregation['#Truth'] = 'IsInTruth', sum
    aggregation['#GATK'] = 'IsInGATK', sum
    aggregation['#GATK accessible'] = ('IsInGATK', 'IsAccessible'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#TP'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#FP'] = ('IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([x[0] and not x[1] and x[2] for x in list(rec)])
    aggregation['#FN'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([not x[0] and x[1] for x in list(rec)])

    tbl_variant_summary = (tbl_whole_genome_comparisons[isolate_code]
        .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
        .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation)
        .convertnumbers()
        .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
        .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
    )
    
    return(tbl_variant_summary)

# <codecell>

rewrite=True

tbl_variant_summary = collections.OrderedDict()
for isolate_code in isolate_codes:
    tbl_variant_summary[isolate_code] = collections.OrderedDict()
    for min_Truth_distance_nearest in [0, 10, 100]:
        print(isolate_code, min_Truth_distance_nearest)
        tbl_variant_summary_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_variant_summary_%s_%d" % (isolate_code, min_Truth_distance_nearest))
        if not os.path.exists(tbl_variant_summary_cache_fn) or rewrite:
            tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)] = calc_variant_summary_table(isolate_code, min_Truth_distance_nearest)
            etl.topickle(tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)], tbl_variant_summary_cache_fn)
        else:
            tbl_variant_summary[isolate_code][str(min_Truth_distance_nearest)] = etl.frompickle(tbl_variant_summary_cache_fn)

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

# <headingcell level=1>

# Plot ROC curves

# <codecell>

for isolate_code in isolate_codes:
    plotROC(isolate_code)

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


