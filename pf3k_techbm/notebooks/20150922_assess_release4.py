# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb
CACHE_DIR = os.path.join(CACHE_DIR, '20150918_assess_indel_filtering')
PROCESSED_ASSEMBLED_SAMPLES_DIR = '/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/validation_results'

# <codecell>

PIPELINES_FN='%s/src/github/malariagen/methods-dev/pf3k_techbm/scripts/gatk_pipeline/pipeline_parameters_20150918.xlsx' % os.environ['HOME']
PLOTS_DIR = '/Users/rpearson/Documents/projects/Pf3k_techbm/slides/20150915_assess_release4/plots'
!mkdir -p {PLOTS_DIR}

# <codecell>

# CODING_EFFECTS = ['NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED']
NONCODING_EFFECTS = ['INTERGENIC', 'INTRON', 'TRANSCRIPT']
 

# <codecell>

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))

# <headingcell level=1>

# Copy data from malsrv2

# <codecell>

tbl_pipelines = etl.fromxlsx(PIPELINES_FN)
tbl_pipelines.display(index_header=True)

# <codecell>


# <headingcell level=1>

# Determine samples

# <codecell>

tbl_samples_to_process = (tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
    .selecteq('To be used for', 'Validation')
    .addfield('truth_vcf_filestem', lambda rec: os.path.join(
        '/nfs/team112_internal/production_files/Pf3k/methods/assembled_samples',
        'truth_vcfs_2',
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

def determine_consensus_region2(rec, column1=4, column2=41, column3=79):
    vals = []
    if(rec[column1] is not None):
        vals.append(str(rec[column1]))
    if(rec[column2] is not None):
        vals.append(str(rec[column2]))
    if(rec[column3] is not None):
        vals.append(str(rec[column3]))
    unique_vals = np.unique(vals)
    return('__'.join(unique_vals))
    

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

fields_to_compare = ['QUAL', 'FILTER', 'AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ',
                     'MQRankSum', 'QD', 'ReadPosRankSum', 'RegionType', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS',
                     'SNPEFF_IMPACT', 'SOR', 'VQSLOD', 'culprit', 'set', 'AD', 'GT']
gatk_renames = dict(zip(fields_to_compare, ['gatk_%s' % field_name for field_name in fields_to_compare]))

# <codecell>

def assess_validation(isolate_code='7G8', chromosome='Pf3D7_01_v3', MAPPING_DIR='bwa_mem', BAMS_DIR='gatk_rec',
                      UNFILTERED_DIR='HaplotypeCaller', FILTERED_DIR='alistair_ann', ANNOTATED_DIR='SnpEff_region',
                      rewrite=False):
    ox_code = tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('ox_code')[0]
    
    tbl_truth_cache_dir = os.path.join(CACHE_DIR, 'tbl_truth', MAPPING_DIR, BAMS_DIR, UNFILTERED_DIR, FILTERED_DIR, ANNOTATED_DIR)
    tbl_gatk_cache_dir = os.path.join(CACHE_DIR, 'tbl_gatk', MAPPING_DIR, BAMS_DIR, UNFILTERED_DIR, FILTERED_DIR, ANNOTATED_DIR)
    tbl_comparison_cache_dir = os.path.join(CACHE_DIR, 'tbl_comparison', MAPPING_DIR, BAMS_DIR, UNFILTERED_DIR, FILTERED_DIR, ANNOTATED_DIR)
    if not os.path.exists(tbl_truth_cache_dir):
        os.makedirs(tbl_truth_cache_dir)
    if not os.path.exists(tbl_gatk_cache_dir):
        os.makedirs(tbl_gatk_cache_dir)
    if not os.path.exists(tbl_comparison_cache_dir):
        os.makedirs(tbl_comparison_cache_dir)
    tbl_truth_cache_fn = os.path.join(tbl_truth_cache_dir, "tbl_truth_%s.%s" % (isolate_code, chromosome))
    tbl_gatk_cache_fn = os.path.join(tbl_gatk_cache_dir, "tbl_gatk_%s.%s" % (isolate_code, chromosome))
    tbl_comparison_cache_fn = os.path.join(tbl_comparison_cache_dir, "tbl_comparison_%s.%s" % (isolate_code, chromosome))
    
#     print(tbl_truth_cache_fn)
#     print(tbl_gatk_cache_fn)
#     print(tbl_comparison_cache_fn)
    
    if not os.path.exists(tbl_truth_cache_fn) or rewrite:
        truth_vcf_fn = "%s.%s.annotated.vcf.gz" % (
            tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('truth_vcf_filestem')[0],
            chromosome
        )
        
        tbl_truth = (etl.fromvcf(truth_vcf_fn)
            .unpackdict('INFO', samplesize=1000000)
            .rename('ALT', 'truth_ALTs')
            .addfield('ALT', lambda rec: str(rec[4][0]))
            .addfield('Truth_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
            .cut(['CHROM', 'POS', 'REF', 'ALT', 'truth_ALTs', 'RegionType', 'Truth_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
                  'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID', 'RPA', 'RU', 'STR'])
            .addfield('Truth_RPA_REF', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][0]))
            .addfield('Truth_RPA_ALT', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][1]))
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
        gatk_vcf_fn = "%s.SNP_INDEL.%s.vcf.gz" % (
            os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, MAPPING_DIR, BAMS_DIR, UNFILTERED_DIR, FILTERED_DIR, ANNOTATED_DIR,
                         '_annotated_vcfs', 'filtered.annotated'),
            chromosome
        )
        print(gatk_vcf_fn)
        
        tbl_gatk = (etl.fromvcf(gatk_vcf_fn)
            .sort(['CHROM', 'POS'])
            .vcfmeltsamples()
            .vcfunpackcall()
            .unpackdict('INFO', samplesize=1000000)
            .selecteq('SAMPLE', ox_code)
            .selectnotnone('GT')
            .selectne('GT', '0/0')
            .rename('ALT', 'gatk_ALTs')
            .addfield('ALT', lambda rec: str(rec['gatk_ALTs'][int(rec['GT'][2])-1]))
            .addfield('gatk_RPA_REF', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][0]))
            .addfield('gatk_RPA_ALT', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][int(rec['GT'][2])]))
            .rename('RPA', 'gatk_RPA')
            .rename('STR', 'gatk_STR')
            .rename('RU', 'gatk_RU')
            .addfield('GATK_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
            .addfield('is_het', lambda rec: None if rec['GT'] is None else rec['GT'][0] != rec['GT'][2])
            .rename(gatk_renames)
#             .select('ALT', lambda rec: str(rec) != '<*:DEL>')
            .select('ALT', lambda rec: str(rec) != '*')
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
            .addfield('IsInGATK', lambda rec: rec['gatk_RegionType'] is not None)
            .addfield('RegionType', lambda rec: determine_consensus_region(rec, 'Truth_RegionType', 'gatk_RegionType'))
            .addfield('IsCoding', lambda rec: determine_consensus_region(rec, 'Truth_is_coding', 'GATK_is_coding'))
            .addfield('Effect', lambda rec: determine_consensus_region(rec, 'Truth_Effect', 'gatk_SNPEFF_EFFECT'))
            .addfield('mode', lambda rec: determine_consensus_region(rec, 'Truth_mode', 'GATK_mode'))
            .addfield('IsAccessible', lambda rec: len(lkp.search(rec['POS'])) > 0)
            .sort('gatk_VQSLOD', reverse=True)
        )
        etl.topickle(tbl_comparison, tbl_comparison_cache_fn)
    else:
        tbl_comparison = etl.frompickle(tbl_comparison_cache_fn)
        
    return(tbl_comparison)

# <codecell>

pipelines_complete = ['default', 'alistair_ann']

# <codecell>

rewrite=False

for rec in tbl_pipelines.data():
#     if True:
    if rec[0] in pipelines_complete:
        annotated_vcfs_dir = "%s/%s/%s/%s/%s/%s/" % (PROCESSED_ASSEMBLED_SAMPLES_DIR, rec[1], rec[2], rec[3], rec[4], rec[5])

        isolate_codes = list(tbl_samples_to_process.values('Isolate code'))
        # isolate_codes = ['7G8']
        chromosomes = ['Pf3D7_%02d_v3' % chrom for chrom in range(1,15)]

        tbl_comparisons = collections.OrderedDict()
        tbl_whole_genome_comparisons = collections.OrderedDict()

    #     tbl_comparison_all_cache_dir = os.path.join(CACHE_DIR, 'tbl_comparison', rec[1], rec[2], rec[3], rec[4], rec[5])
    #     tbl_comparison_all_cache_fn = os.path.join(tbl_comparison_all_cache_dir, "tbl_comparison_all")

        # if not os.path.exists(tbl_comparison_all_cache_fn) or rewrite:
        for isolate_code in isolate_codes:
#             print(rec[0], isolate_code)
            tbl_comparisons[isolate_code] = collections.OrderedDict()
            for chromosome in chromosomes:
                print(rec[0], isolate_code, chromosome)
                tbl_comparisons[isolate_code][chromosome] = assess_validation(isolate_code, chromosome, MAPPING_DIR=rec[1],
                                                                              BAMS_DIR=rec[2], UNFILTERED_DIR=rec[3],
                                                                              FILTERED_DIR=rec[4], ANNOTATED_DIR=rec[5],
                                                                              rewrite=rewrite)

            tbl_comparison_cache_dir = os.path.join(CACHE_DIR, 'tbl_comparison', rec[1], rec[2], rec[3], rec[4], rec[5])
            tbl_comparison_sample_cache_fn = os.path.join(tbl_comparison_cache_dir, "tbl_comparison_%s" % (isolate_code))
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
                    .sort('gatk_VQSLOD', reverse=True)
                )
                etl.topickle(tbl_whole_genome_comparisons[isolate_code], tbl_comparison_sample_cache_fn)
            else:
                tbl_whole_genome_comparisons[isolate_code] = etl.frompickle(tbl_comparison_sample_cache_fn)
                
        tbl_comparison_sample_cache_fn = os.path.join(tbl_comparison_cache_dir, "tbl_comparison_All")
        if not os.path.exists(tbl_comparison_sample_cache_fn) or rewrite:
            tbl_whole_genome_comparisons['All'] = (
                tbl_whole_genome_comparisons['7G8'].addfield('Sample', '7G8')
                .cat(tbl_whole_genome_comparisons['GB4'].addfield('Sample', 'GB4'))
                .cat(tbl_whole_genome_comparisons['KE01'].addfield('Sample', 'KE01'))
                .cat(tbl_whole_genome_comparisons['KH02'].addfield('Sample', 'KH02'))
                .cat(tbl_whole_genome_comparisons['GN01'].addfield('Sample', 'GN01'))
                .sort('gatk_VQSLOD', reverse=True)
            )
            etl.topickle(tbl_whole_genome_comparisons['All'], tbl_comparison_sample_cache_fn)
        else:
            tbl_whole_genome_comparisons['All'] = etl.frompickle(tbl_comparison_sample_cache_fn)


        #     etl.topickle(tbl_whole_genome_comparisons, tbl_comparison_all_cache_fn)
        # else:
        #     tbl_whole_genome_comparisons = etl.frompickle(tbl_comparison_all_cache_fn)


# <codecell>

rewrite=False

tbl_whole_genome_comparisons = collections.OrderedDict()
comparison_arrays = collections.OrderedDict()

for rec in tbl_pipelines.sort('PIPELINE_NAME').data():
#     print(rec[0])
#     if not rec[0] in ['hrun']:
    if rec[0] in pipelines_complete:
#     if True:
        print(rec[0])
        tbl_comparison_cache_dir = os.path.join(CACHE_DIR, 'tbl_comparison', rec[1], rec[2], rec[3], rec[4], rec[5])
        tbl_comparison_sample_cache_fn = os.path.join(tbl_comparison_cache_dir, "tbl_comparison_All")
        tbl_whole_genome_comparisons[rec[0]] = etl.frompickle(tbl_comparison_sample_cache_fn)
        comparison_array_cache_fn = os.path.join(tbl_comparison_cache_dir, "comparison_array.npy")
        if not os.path.exists(comparison_array_cache_fn) or rewrite:
            comparison_array = etl.toarray(
                tbl_whole_genome_comparisons[rec[0]]
                .addfield('num_alts', lambda rec: 0 if rec['gatk_ALTs'] is None else len(rec['gatk_ALTs']))
                .addfield('num_filters', lambda rec: 0 if rec['gatk_FILTER'] is None else len(rec['gatk_FILTER']))
                .replace('Truth_distance_nearest', None, 0)
                .replace('STR', None, False)
                .replace('Effect', None, '')
                .replace('gatk_VQSLOD', None, np.nan)
#                 .replace('gatk_ALTs', None, np.nan)
#                 .replace('gatk_FILTER', None, '')
                .replace('gatk_QUAL', None, np.nan)
#                 .replace('gatk_AC', None, np.nan)
#                 .replace('gatk_AF', None, np.nan)
#                 .replace('gatk_AN', None, np.nan)
                .replace('gatk_BaseQRankSum', None, np.nan)
                .replace('gatk_ClippingRankSum', None, np.nan)
                .replace('GC', None, np.nan)
                .replace('gatk_DP', None, np.nan)
                .replace('gatk_FS', None, np.nan)
                .replace('gatk_MQ', None, np.nan)
                .replace('gatk_MQRankSum', None, np.nan)
                .replace('gatk_QD', None, np.nan)
                .replace('gatk_ReadPosRankSum', None, np.nan)
                .replace('gatk_SOR', None, np.nan)
                .replace('GQ', None, 0)
#                 .replace('is_het', None, False)
#                 .replace('NEGATIVE_TRAIN_SITE', None, False)
#                 .replace('POSITIVE_TRAIN_SITE', None, False)
#                 .replace('gatk_RPA_REF', None, )
                .convert('gatk_RPA_REF', lambda rec: np.nan if rec is None else int(rec))
                .convert('gatk_RPA_ALT', lambda rec: np.nan if rec is None else int(rec))
#                 .replace('gatk_RPA_ALT', None, '')
                .replace('gatk_RU', None, '')
                .replace('gatk_SNPEFF_IMPACT', None, '')
                .replace('gatk_STR', None, False)
                .replace('gatk_culprit', None, '')
                .replace('gatk_set', None, '')
                .cut(['CHROM', 'POS', 'IsInTruth', 'IsInGATK', 'RegionType', 'IsCoding', 'Effect', 'mode', 'IsAccessible',
                    'Truth_distance_nearest', 'STR',
                    'gatk_VQSLOD',
                    'num_alts',
                    'num_filters',
                    'gatk_QUAL',
# #                     #'gatk_AC', 'gatk_AF', 'gatk_AN',
                    'gatk_BaseQRankSum', 'gatk_ClippingRankSum',
                    'GQ', 'GC', 'gatk_MQ', 'gatk_MQRankSum', 'gatk_QD',
                    'gatk_ReadPosRankSum', 'gatk_SOR',
# 'GQ',
#                     'NEGATIVE_TRAIN_SITE', 'POSITIVE_TRAIN_SITE',
                    'gatk_RPA_REF',
                    'gatk_RPA_ALT', 'gatk_RU',
                    'gatk_SNPEFF_IMPACT', 'gatk_STR', 'gatk_culprit', 'gatk_set'
                ])
            )
            np.save(comparison_array_cache_fn, comparison_array)
        else:
            comparison_array = np.load(comparison_array_cache_fn)
        comparison_arrays[rec[0]] = comparison_array

# <codecell>

print(len(tbl_whole_genome_comparisons['default'].data()))
print(len(tbl_whole_genome_comparisons['alistair_ann'].data()))

# <codecell>

print(shape(comparison_arrays['default']))
print(shape(comparison_arrays['alistair_ann']))

# <codecell>

test_filters = collections.OrderedDict()
test_variants = collections.OrderedDict()
# for pipeline in pipelines_complete:
for pipeline in pipelines_complete:
    print(pipeline)
    test_filters[pipeline] = (
        (comparison_arrays[pipeline]['IsInGATK']) &
        (comparison_arrays[pipeline]['RegionType'] == 'Core') &
        (comparison_arrays[pipeline]['IsCoding'] == 'True') &
        (comparison_arrays[pipeline]['IsAccessible'] == True) &
        (np.in1d(comparison_arrays[pipeline]['mode'], ['IND', 'False__INDEL', 'False__True__IND', 'INDEL__True']))
    #     (comparison_arrays['alistair_ann']['STR'] == True)
    )
    test_variants[pipeline] = comparison_arrays[pipeline][test_filters[pipeline]]

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>

def calc_indel_vqslod_summary_table(isolate_code='All', min_Truth_distance_nearest=0, min_VQSLOD=-99999):
    aggregation = collections.OrderedDict()
    aggregation['# Variants'] = 'IsInGATK', sum
    aggregation['# TP'] = 'IsInTruth', sum
    aggregation['# Coding'] = ('IsCoding'), lambda rec: np.sum([x=='True' for x in list(rec)])
    aggregation['# Coding indels'] = ('REF', 'ALT', 'IsCoding'), lambda rec: np.sum([x[2]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) for x in list(rec)])
#     aggregation['Ti'] = ('REF', 'ALT'), lambda rec: calc_ti(list(rec))
#     aggregation['Tv'] = ('REF', 'ALT'), lambda rec: calc_tv(list(rec))
#     aggregation['TiTv'] = ('REF', 'ALT'), lambda rec: calc_titv(list(rec))
    aggregation['# Coding indels mod3'] = ('REF', 'ALT', 'IsCoding'), lambda rec: np.sum([x[2]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) and (len(x[0]) - len(x[1])) % 3 == 0 for x in list(rec)])
    aggregation['NUM_MULTIALLELIC_COMPLEX'] = 'vt1', lambda rec: np.sum([x == 'MULTIALLELIC_COMPLEX' for x in list(rec)])
    aggregation['NUM_MULTIALLELIC_MIXED'] = 'vt1', lambda rec: np.sum([x == 'MULTIALLELIC_MIXED' for x in list(rec)])
    aggregation['NUM_INSERTION'] = 'vt1', lambda rec: np.sum([x == 'INSERTION' for x in list(rec)])
    aggregation['NUM_DELETION'] = 'vt1', lambda rec: np.sum([x == 'DELETION' for x in list(rec)])
    aggregation['NUM_1rep'] = 'event_len', lambda rec: np.sum([x == 1 for x in list(rec)])
    aggregation['NUM_1rep_DEL'] = ('event_len', 'vt1'), lambda rec: np.sum([x[0] == 1 and x[1] == 'DELETION' for x in list(rec)])
    aggregation['# het'] = 'is_het', sum
    aggregation['# <10bp nearest'] = 'GATK_distance_nearest', lambda rec: np.sum([x < 10 for x in list(rec)])
    aggregation['# >=100bp nearest'] = 'GATK_distance_nearest', lambda rec: np.sum([x >= 100 for x in list(rec)])
    aggregation['median nearest'] = 'GATK_distance_nearest', lambda rec: np.nanmedian([x for x in list(rec)])
    aggregation['median QUAL'] = 'gatk_QUAL', lambda rec: np.nanmedian([x for x in list(rec)])
    aggregation['median DP'] = 'gatk_DP', lambda rec: np.nanmedian([x for x in list(rec)])
    aggregation['median GQ'] = 'GQ', lambda rec: np.nanmedian([x for x in list(rec)])
    aggregation['median gatk_BaseQRankSum'] = 'gatk_BaseQRankSum', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['# gatk_BaseQRankSum'] = 'gatk_BaseQRankSum', lambda rec: np.sum([x is not None for x in list(rec) if x is not None])
    aggregation['# gatk_FS'] = 'gatk_FS', lambda rec: np.sum([x is not None for x in list(rec) if x is not None])
    aggregation['median gatk_ClippingRankSum'] = 'gatk_ClippingRankSum', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['median gatk_FS'] = 'gatk_FS', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['median gatk_MQ'] = 'gatk_MQ', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['median gatk_MQRankSum'] = 'gatk_MQRankSum', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['median gatk_QD'] = 'gatk_QD', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['median gatk_SOR'] = 'gatk_SOR', lambda rec: np.nanmedian([x for x in list(rec) if x is not None])
    aggregation['mean #ALTs'] = 'gatk_ALTs', lambda rec: np.mean([len(x) for x in list(rec) if x is not None])
    aggregation['# multiallelic'] = 'gatk_ALTs', lambda rec: np.sum([len(x) > 1 for x in list(rec) if x is not None])
      
    
    tbl_variant_summary = (tbl_whole_genome_comparisons[isolate_code]
        .selecteq('mode', 'INDEL')
        .selecteq('RegionType', 'Core')
        .selecteq('IsCoding', 'True')
        .selecteq('IsAccessible', True)
        .selectnotnone('gatk_VQSLOD')
        .addfield('int_VQSLOD', lambda rec: None if rec['gatk_VQSLOD'] is None else 0 if rec['gatk_VQSLOD'] <= 0 else 6 if rec['gatk_VQSLOD'] >= 6 else int(rec['gatk_VQSLOD']))
        .addfield('vt1', lambda rec: None if rec['VariantType'] is None else rec['VariantType'].split('.')[0])
        .addfield('num_reps', lambda rec: None if re.match(r'.*NumRepetitions_([0-9]+)\..*', rec['VariantType']) is None else int(re.sub(r'.*NumRepetitions_([0-9]+)\..*', r'\1', rec['VariantType'])))
        .addfield('event_len', lambda rec: None if re.match(r'.*EventLength_([0-9]+).*', rec['VariantType']) is None else int(re.sub(r'.*EventLength_([0-9]+).*', r'\1', rec['VariantType'])))
        .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#         .selectge('Truth_distance_nearest', 100)
#         .addfield('Truth_distance_log10', lambda rec: 2 if rec['Truth_distance_nearest'] > 100 else 0 if rec['Truth_distance_nearest'] < 1 else int(log10(rec['Truth_distance_nearest'])))
#         .aggregate(('CHROM'), aggregation)
        .aggregate(('int_VQSLOD'), aggregation) 
        .addfield('% Frameshift indels', lambda rec: 0.0 if rec['# Coding indels'] == 0 else round(1.0 - (rec['# Coding indels mod3'] / rec['# Coding indels']), 3))
        .addfield('% MULTIALLELIC_COMPLEX', lambda rec: round(rec['NUM_MULTIALLELIC_COMPLEX'] / rec['# Variants'], 3))
        .addfield('% MULTIALLELIC_MIXED', lambda rec: round(rec['NUM_MULTIALLELIC_MIXED'] / rec['# Variants'], 3))
        .addfield('% INSERTION', lambda rec: round(rec['NUM_INSERTION'] / rec['# Variants'], 3))
        .addfield('% DELETION', lambda rec: round(rec['NUM_DELETION'] / rec['# Variants'], 3))
        .addfield('% NUM_1rep', lambda rec: round(rec['NUM_1rep'] / rec['# Variants'], 3))
        .addfield('% NUM_1rep_DEL', lambda rec: round(rec['NUM_1rep_DEL'] / rec['# Variants'], 3))
        .addfield('% het', lambda rec: round(rec['# het'] / rec['# Variants'], 3))
        .addfield('% <10bp nearest', lambda rec: round(rec['# <10bp nearest'] / rec['# Variants'], 3))
        .addfield('% >=100bp nearest', lambda rec: round(rec['# >=100bp nearest'] / rec['# Variants'], 3))
        .addfield('% gatk_BaseQRankSum', lambda rec: round(rec['# gatk_BaseQRankSum'] / rec['# Variants'], 3))
        .addfield('% multiallelic', lambda rec: round(rec['# multiallelic'] / rec['# Variants'], 3))
    )
    
    return(tbl_variant_summary)

# <codecell>

def plotFDRsens(pipeline='default', RegionType='Core', IsCoding=True, mode='SNP', x_var='FDR',
            min_Truth_distance_nearest=100, min_GATK_distance_nearest=0, sort_by=None, sort_reverse=True,
            annotate_values=None, selected_threshold=0, show_annotations=False, ymax=1.0, xmax=1.0, ax=None):
    
    if annotate_values is None:
        if mode=='SNP':
            annotate_values = [-999, 0, 4, 6, 8]
        else:
            annotate_values = [-999, 0, 2, 4, 6]
    # set up axes
    if ax is None:
        fig, ax = plt.subplots()

    if(sort_by is None):
        array_to_plot = (tbl_whole_genome_comparisons[pipeline]
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .addfield('is_frameshift', lambda rec: (len(rec['REF']) - len(rec['ALT'])) % 3 != 0)
            .addfield('TruthUnclustered', lambda rec: False if rec['Truth_distance_nearest'] is None else rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
            .cut(['IsInGATK', 'IsInTruth', 'TruthUnclustered', 'gatk_VQSLOD'])
        ).toarray()
        annotate_indices = np.array([np.sum(array_to_plot['gatk_VQSLOD'] >= value) for value in annotate_values])
        
    else:
        array_to_plot = (tbl_whole_genome_comparisons[pipeline]
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .addfield('TruthUnclustered', lambda rec: False if rec['Truth_distance_nearest'] is None else rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
            .sort(sort_by, reverse=sort_reverse)
            .cut(['IsInGATK', 'IsInTruth', 'TruthUnclustered', 'gatk_VQSLOD', sort_by])
        ).toarray()
        annotate_indices = np.array([np.sum(array_to_plot[sort_by] >= value) for value in annotate_values])

    TPs_unclustered = np.cumsum(array_to_plot['IsInGATK'] & array_to_plot['TruthUnclustered'])
    Truth_unclustered = np.cumsum(array_to_plot['TruthUnclustered'])
    TPs_all = np.cumsum(array_to_plot['IsInGATK'] & array_to_plot['IsInTruth'])
    FPs = np.cumsum(array_to_plot['IsInGATK'] & np.logical_not(array_to_plot['IsInTruth']))
    
#     print(TPs_unclustered, len(TPs_unclustered))
#     print(Truth_unclustered, len(Truth_unclustered))
#     print(TPs_all, len(TPs_all))
#     print(FPs, len(FPs))
#     print(FPs[1000:1010])
    
    number_of_positives = len(array_to_plot)
    number_of_true_variants = len(
        (tbl_whole_genome_comparisons[pipeline]
            .selecteq('IsInTruth', True)
            .selecteq('mode', mode)
            .selecteq('RegionType', RegionType)
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
#             .select(lambda rec: rec['GATK_distance_nearest'] is None or rec['GATK_distance_nearest'] >= min_GATK_distance_nearest)
#             .selectge('GATK_distance_nearest', min_GATK_distance_nearest)
            .cut(['IsInGATK', 'IsInTruth', 'gatk_VQSLOD'])
            .data()
        )
    )
    
    print(pipeline, RegionType, IsCoding, mode, number_of_positives, number_of_true_variants)
    
    sensitivity = TPs_unclustered / number_of_true_variants
    if x_var == 'FDR':
        FDR = FPs / (TPs_all + FPs)
    else:
        FDR = FPs / number_of_positives
        
#     print(FDR)
#     print(sensitivity)
    
    
    ax.plot(FDR, sensitivity, label=pipeline)
    if show_annotations:
        for i, annotate_index in enumerate(annotate_indices):
            point_type = 'bo'
            if annotate_values[i] == selected_threshold:
                point_type = 'ro'
            ax.plot(FDR[annotate_index], sensitivity[annotate_index], point_type)
            annotate_text = '%s (%4.1f%% sens, %4.1f%% FDR)' % (annotate_values[i], sensitivity[annotate_index]*100, FDR[annotate_index]*100)
#             print(annotate_text)
            ax.annotate(annotate_text,
                        (FDR[annotate_index], sensitivity[annotate_index]), xycoords='data', xytext=(10, -10),
                        textcoords='offset points', rotation=-45)
    if ymax is not None:
        ax.set_ylim([0, ymax])
    if xmax is not None:
        ax.set_xlim([0, xmax])
            
    return(ax)


# <codecell>

plotFDRsens(show_annotations=True, min_Truth_distance_nearest=0)

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate([True, False]):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        plotFDRsens(pipeline='default', show_annotations=True, mode=mode, IsCoding=IsCoding, min_Truth_distance_nearest=100, ax=ax)
        ax.set_title("%s %s" % (mode, 'Coding' if IsCoding==True else 'Non-coding'))
# fig.tight_layout()
fig.savefig('%s/default_unclustered_sensitivity_vs_FDR.pdf' % PLOTS_DIR, format='pdf')

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate([True, False]):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        plotFDRsens(pipeline='alistair_ann', show_annotations=True, mode=mode, IsCoding=IsCoding, min_Truth_distance_nearest=100, ax=ax)
        ax.set_title("%s %s" % (mode, 'Coding' if IsCoding==True else 'Non-coding'))
# fig.tight_layout()
fig.savefig('%s/improved_unclustered_sensitivity_vs_FDR.pdf' % PLOTS_DIR, format='pdf')

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate([True, False]):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        plotFDRsens(pipeline='default', show_annotations=True, mode=mode, IsCoding=IsCoding, min_Truth_distance_nearest=0, ax=ax)
        ax.set_title("%s %s" % (mode, 'Coding' if IsCoding==True else 'Non-coding'))
# fig.tight_layout()
fig.savefig('%s/default_full_sensitivity_vs_FDR.pdf' % PLOTS_DIR, format='pdf')

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate([True, False]):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        plotFDRsens(pipeline='alistair_ann', show_annotations=True, mode=mode, IsCoding=IsCoding, min_Truth_distance_nearest=0, ax=ax)
        ax.set_title("%s %s" % (mode, 'Coding' if IsCoding==True else 'Non-coding'))
# fig.tight_layout()
fig.savefig('%s/improved_full_sensitivity_vs_FDR.pdf' % PLOTS_DIR, format='pdf')

# <codecell>


# <headingcell level=1>

# Oddness of VQSLOD plots

# <codecell>

variants_classes = collections.OrderedDict()
variants_classes['Coding indels'] = ('True', 'INDEL')
variants_classes['Non-coding indels'] = ('False', 'INDEL')
variants_classes['Coding SNPs'] = ('True', 'SNP')
variants_classes['Non-coding SNPs'] = ('False', 'SNP')
test_filters = collections.OrderedDict()
test_variants = collections.OrderedDict()
for pipeline in ['default', 'alistair_ann']:
    print(pipeline)
    test_filters[pipeline] = collections.OrderedDict()
    test_variants[pipeline] = collections.OrderedDict()
    for variants_class in variants_classes:
        test_filters[pipeline][variants_class] = (
            (comparison_arrays[pipeline]['IsInGATK']) &
            (comparison_arrays[pipeline]['RegionType'] == 'Core') &
            (comparison_arrays[pipeline]['IsCoding'] == variants_classes[variants_class][0]) &
            (comparison_arrays[pipeline]['IsAccessible'] == True) &
            (comparison_arrays[pipeline]['mode'] == variants_classes[variants_class][1])
        )
        test_variants[pipeline][variants_class] = comparison_arrays[pipeline][test_filters[pipeline][variants_class]]

# <codecell>

TP_states = (True, False)
TP_colors = ('green', 'red')
TP_labels = ('TP', 'FP')

def plot_annotation_diagnostics(variants, annotation, bins):
    
    TPs = variants['IsInTruth']
    ann = variants[annotation]
    n_variants, bins = np.histogram(ann, bins=bins)

    fig = plt.figure(figsize=(9, 3))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    pc_bottom = 0
    n_bottom = 0
    counts = dict()
    for st, lbl, clr in zip(TP_states, TP_labels, TP_colors):
        values = (TPs == st)
        n, _, _ = scipy.stats.binned_statistic(ann, values, statistic='sum', bins=bins)
        counts[st] = n
        left = bins[:-1]
        width = np.diff(bins)
        # plot absolute numbers
        ax1.bar(left, n, width=width, bottom=n_bottom, color=clr, label=lbl, lw=1)
        n_bottom += n
        # plot percentages
        pc = n * 100 / n_variants
        ax2.bar(left, pc, width=width, bottom=pc_bottom, color=clr, label=lbl, lw=0)
        pc_bottom += pc
        
    for s in 'top', 'left', 'right':
        ax1.spines[s].set_visible(False)
    ax1.xaxis.tick_bottom()
    ax1.set_yticks([])
    ax1.set_xlim(bins[0], bins[-1])
    ax1.set_ylabel('#calls')
    ax2.set_ylim(0, 100)
    ax2.set_xlim(bins[0], bins[-1])
    ax2.set_ylabel('%calls')
    fig.suptitle(annotation, fontsize=14, fontweight='bold')
    fig.tight_layout()

# <codecell>

for pipeline in ['default', 'alistair_ann']:
    print(pipeline)
    plot_annotation_diagnostics(test_variants[pipeline], 'gatk_VQSLOD', np.linspace(-2, 10, 100))
    plt.show()

# <codecell>


# <headingcell level=1>

# Variant summary tables

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

def calc_titv(rec_list, variant_type='TP', VQSLOD_threshold=None):
    if variant_type=='TP':
        new_list = [x for x in rec_list if x[2] and x[3] and (VQSLOD_threshold is None or x[5] > VQSLOD_threshold)]
    elif variant_type=='FP':
        new_list = [x for x in rec_list if x[2] and not x[3] and x[4] and (VQSLOD_threshold is None or x[5] > VQSLOD_threshold)]
    if variant_type=='FN':
        new_list = [x for x in rec_list if (not x[2] and x[3]) or (VQSLOD_threshold is None or (x[5] is not None and x[5] <= VQSLOD_threshold))]
    ti = np.count_nonzero([is_transition(rec) for rec in new_list])
    tv = np.count_nonzero([is_transversion(rec) for rec in new_list])
    if tv == 0:
        titv = 0.0
    else:
        titv = 1.0*ti/tv
#     print(ti, tv, titv)
    return(titv)

# def calc_titv(rec_list):
#     ti = np.count_nonzero([is_transition(rec) for rec in rec_list])
#     tv = np.count_nonzero([is_transversion(rec) for rec in rec_list])
#     if tv == 0:
#         titv = 0.0
#     else:
#         titv = 1.0*ti/tv
# #     print(ti, tv, titv)
#     return(titv)

# def calc_ti(rec_list):
#     return np.count_nonzero([is_transition(rec) for rec in rec_list])
# def calc_tv(rec_list):
#     return np.count_nonzero([is_transversion(rec) for rec in rec_list])

# <codecell>

def calc_variant_summary_table(pipeline='default', min_Truth_distance_nearest=100, min_VQSLOD=-99999, VQSLOD_threshold=0.0):
    aggregation = collections.OrderedDict()
    aggregation['#Variant sites discovered'] = ('IsInGATK', 'CHROM', 'POS'), lambda tbl: len(etl.wrap(tbl).selecteq(0, True).distinct([1, 2]).data())
    aggregation['#Variant sites after filtering'] = ('IsInGATK', 'CHROM', 'POS', 'gatk_VQSLOD'), lambda tbl: len(etl.wrap(tbl).selecteq(0, True).selectgt(3, VQSLOD_threshold).distinct([1, 2]).data())
    aggregation['Non-ref genotypes'] = 'IsInGATK', sum
    aggregation['#Non-ref genotypes accessible in PacBio'] = ('IsInGATK', 'IsAccessible'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#Non-ref genotypes after filtering'] = ('IsInGATK', 'gatk_VQSLOD'), lambda rec: np.sum([x[0] and x[1] > VQSLOD_threshold for x in list(rec)])
    aggregation['#Non-ref genotypes accessible after filtering'] = ('IsInGATK', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: np.sum([x[0] and x[1] and x[2] > VQSLOD_threshold for x in list(rec)])
    aggregation['#PacBio all variants sites'] = ('IsInTruth', 'CHROM', 'POS'), lambda tbl: len(etl.wrap(tbl).selecteq(0, True).distinct([1, 2]).data())
    aggregation['#PacBio all genotypes'] = 'IsInTruth', sum
    aggregation['#PacBio unclustered variants sites'] = ('IsInTruth', 'CHROM', 'POS', 'Truth_distance_nearest'), lambda tbl: len(etl.wrap(tbl).selecteq(0, True).selectge(3, min_Truth_distance_nearest).distinct([1, 2]).data())
    aggregation['#PacBio unclustered genotypes'] = ('IsInTruth', 'Truth_distance_nearest'), lambda rec: np.sum([x[0] and (x[1] is not None and x[1] >= min_Truth_distance_nearest) for x in list(rec)])
    aggregation['#TP unfiltered'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
    aggregation['#TP unclustered unfiltered'] = ('IsInGATK', 'IsInTruth', 'Truth_distance_nearest'), lambda rec: np.sum([x[0] and x[1] and (x[2] is not None and x[2] >= min_Truth_distance_nearest) for x in list(rec)])
    aggregation['#FP unfiltered'] = ('IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([x[0] and not x[1] and x[2] for x in list(rec)])
    aggregation['#FN all unfiltered'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([not x[0] and x[1] for x in list(rec)])
    aggregation['#FN unclustered unfiltered'] = ('IsInGATK', 'IsInTruth', 'Truth_distance_nearest'), lambda rec: np.sum([not x[0] and x[1] and (x[2] is not None and x[2] >= min_Truth_distance_nearest) for x in list(rec)])
    aggregation['#TP filtered'] = ('IsInGATK', 'IsInTruth', 'gatk_VQSLOD'), lambda rec: np.sum([x[0] and x[1] and x[2] > VQSLOD_threshold for x in list(rec)])
    aggregation['#TP unclustered filtered'] = ('IsInGATK', 'IsInTruth', 'gatk_VQSLOD', 'Truth_distance_nearest'), lambda rec: np.sum([x[0] and x[1] and x[2] > VQSLOD_threshold and (x[3] is not None and x[3] >= min_Truth_distance_nearest) for x in list(rec)])
    aggregation['#FP filtered'] = ('IsInGATK', 'IsInTruth', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: np.sum([x[0] and not x[1] and x[2] and x[3] > VQSLOD_threshold for x in list(rec)])
    aggregation['#FN all filtered'] = ('IsInGATK', 'IsInTruth', 'gatk_VQSLOD'), lambda rec: np.sum([(not x[0] and x[1]) or (x[2] is not None and x[2] <= VQSLOD_threshold) for x in list(rec)])
    aggregation['#FN unclustered filtered'] = ('IsInGATK', 'IsInTruth', 'gatk_VQSLOD', 'Truth_distance_nearest'), lambda rec: np.sum([((not x[0] and x[1]) or (x[2] is not None and x[2] <= VQSLOD_threshold)) and (x[3] is not None and x[3] >= min_Truth_distance_nearest) for x in list(rec)])
    aggregation['TiTv TP unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'TP')
    aggregation['TiTv FP unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FP')
    aggregation['TiTv FN unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FN')
#     aggregation['TiTv FN unclustered unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible', 'Truth_distance_nearest'), lambda rec: calc_titv(list(rec), 'FN')
    aggregation['TiTv TP filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: calc_titv(list(rec), 'TP', VQSLOD_threshold)
    aggregation['TiTv FP filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: calc_titv(list(rec), 'FP', VQSLOD_threshold)
    aggregation['TiTv FN filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: calc_titv(list(rec), 'FN', VQSLOD_threshold)
    aggregation['#TPmod3 unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and x[3] for x in list(rec)])
    aggregation['#FPmod3 unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and not x[3] and x[4] for x in list(rec)])
    aggregation['#FNmod3 unfiltered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and not x[2] and x[3] for x in list(rec)])
    aggregation['#TPmod3 filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'gatk_VQSLOD'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and x[3] and x[4] > VQSLOD_threshold for x in list(rec)])
    aggregation['#FPmod3 filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible', 'gatk_VQSLOD'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and not x[3] and x[4] and x[5] > VQSLOD_threshold for x in list(rec)])
    aggregation['#FNmod3 filtered'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'gatk_VQSLOD'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and (not x[2] and x[3]) or (x[4] is not None and x[4] <= VQSLOD_threshold) for x in list(rec)])

    tbl_variant_summary = (tbl_whole_genome_comparisons[pipeline]
        .selecteq('RegionType', 'Core')
        .select(lambda rec: not str(rec['IsCoding']).startswith('ERROR'))
#         .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
#         .select(lambda rec: rec['gatk_VQSLOD'] is None or rec['gatk_VQSLOD'] >= min_VQSLOD)
        .aggregate(('IsCoding', 'mode'), aggregation)
        .sort(('IsCoding', 'mode'), reverse=True)
        .addfield('Unfiltered all sensitivity', lambda rec: 0.0 if rec['#PacBio all genotypes'] == 0 else round(rec['#TP unfiltered'] / rec['#PacBio all genotypes'], 3))
        .addfield('Unfiltered unclustered sensitivity', lambda rec: 0.0 if rec['#PacBio unclustered genotypes'] == 0 else round(rec['#TP unclustered unfiltered'] / rec['#PacBio unclustered genotypes'], 3))
        .addfield('Unfiltered FDR', lambda rec: 0.0 if rec['#Non-ref genotypes accessible in PacBio'] == 0 else round(rec['#FP unfiltered'] / rec['#Non-ref genotypes accessible in PacBio'], 3))
        .addfield('Filtered all sensitivity', lambda rec: 0.0 if rec['#PacBio all genotypes'] == 0 else round(rec['#TP filtered'] / rec['#PacBio all genotypes'], 3))
        .addfield('Filtered unclustered sensitivity', lambda rec: 0.0 if rec['#PacBio unclustered genotypes'] == 0 else round(rec['#TP unclustered filtered'] / rec['#PacBio unclustered genotypes'], 3))
        .addfield('Filtered FDR', lambda rec: 0.0 if rec['#Non-ref genotypes accessible in PacBio'] == 0 else round(rec['#FP filtered'] / rec['#Non-ref genotypes accessible in PacBio'], 3))
        .addfield('Unfiltered %TPmod3', lambda rec: 0.0 if rec['#TP unfiltered'] == 0 else round(rec['#TPmod3 unfiltered'] / rec['#TP unfiltered'], 3))
        .addfield('Unfiltered %FPmod3', lambda rec: 0.0 if rec['#FP unfiltered'] == 0 else round(rec['#FPmod3 unfiltered'] / rec['#FP unfiltered'], 3))
        .addfield('Unfiltered %FNmod3', lambda rec: 0.0 if rec['#FN all unfiltered'] == 0 else round(rec['#FNmod3 unfiltered'] / rec['#FN all unfiltered'], 3))
        .addfield('Filtered %TPmod3', lambda rec: 0.0 if rec['#TP filtered'] == 0 else round(rec['#TPmod3 filtered'] / rec['#TP filtered'], 3))
        .addfield('Filtered %FPmod3', lambda rec: 0.0 if rec['#FP filtered'] == 0 else round(rec['#FPmod3 filtered'] / rec['#FP filtered'], 3))
        .addfield('Filtered %FNmod3', lambda rec: 0.0 if rec['#FN all filtered'] == 0 else round(rec['#FNmod3 filtered'] / rec['#FN all filtered'], 3))
        .addfield('Filtered %TP frameshift', lambda rec: 1.0 - rec['Filtered %TPmod3'])
        .addfield('Filtered %FP frameshift', lambda rec: 1.0 - rec['Filtered %FPmod3'])
        .addfield('Filtered %FN frameshift', lambda rec: 1.0 - rec['Filtered %FNmod3'])
    )
    
    return(tbl_variant_summary)

# <codecell>

calc_variant_summary_table(pipeline='alistair_ann')

# <codecell>

tbl_whole_genome_comparisons['default'].valuecounts('IsCoding', 'mode').displayall()

# <codecell>

tbl_whole_genome_comparisons['alistair_ann'].valuecounts('IsCoding', 'mode').displayall()

# <codecell>

tbl_whole_genome_comparisons['default'].valuecounts('IsInGATK').displayall()

# <codecell>

tbl_whole_genome_comparisons['alistair_ann'].valuecounts('IsInGATK').displayall()

# <codecell>

tbl_whole_genome_comparisons['default'].valuecounts('RegionType').displayall()

# <codecell>

tbl_whole_genome_comparisons['alistair_ann'].valuecounts('RegionType').displayall()

# <codecell>

rewrite=False

tbl_variant_summaries = collections.OrderedDict()

for rec in tbl_pipelines.sort('PIPELINE_NAME').data():
    if rec[0] in pipelines_complete:
#     if rec[0] == 'alistair_ann':
#     if rec[0] == 'default':
        print(rec[0])
        tbl_comparison_cache_dir = os.path.join(CACHE_DIR, 'tbl_comparison', rec[1], rec[2], rec[3], rec[4], rec[5])
        tbl_variant_summary_cache_fn = os.path.join(tbl_comparison_cache_dir, "tbl_variant_summary")
        if not os.path.exists(tbl_variant_summary_cache_fn) or rewrite:
            tbl_variant_summaries[rec[0]] = calc_variant_summary_table(rec[0])
            etl.topickle(tbl_variant_summaries[rec[0]], tbl_variant_summary_cache_fn)
        else:
            tbl_variant_summaries[rec[0]] = etl.frompickle(tbl_variant_summary_cache_fn)

# <codecell>

# tbl_variant_summary = calc_variant_summary_table()

# <codecell>

tbl_variant_summaries['default'].transpose().displayall()

# <codecell>

tbl_variant_summaries['alistair_ann'].transpose().displayall()

# <codecell>

tbl_variant_summaries['default'].totsv("%s/tbl_variant_summary_default.txt" % (PLOTS_DIR))

# <codecell>

tbl_variant_summaries['default'].toxlsx("%s/tbl_variant_summary_default.xlsx" % (PLOTS_DIR))

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>

min(tbl_whole_genome_comparisons['All'].selectnotnone('gatk_VQSLOD').values('gatk_VQSLOD'))

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate(['True', 'False']):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        plotFDRsens(show_annotations=True, mode=mode, IsCoding=IsCoding, min_Truth_distance_nearest=100, ax=ax)
        ax.set_title("%s %s" % (mode, 'Coding' if IsCoding=='True' else 'Non-coding'))
fig.tight_layout()

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate(['True', 'False']):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        array_to_plot = (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('mode', mode)
            .selecteq('RegionType', 'Core')
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .selectnotnone('gatk_VQSLOD')
            .cut(['gatk_VQSLOD'])
        ).toarray()
        sns.distplot(array_to_plot['gatk_VQSLOD'], kde=False, ax=ax, bins=np.arange(-10, 14, 0.1))
        ax.set_xlim([-10, 14])
        

# <codecell>

for gatk_var in ['gatk_QD', 'gatk_MQ', 'gatk_MQRankSum', 'gatk_ReadPosRankSum', 'gatk_FS', 'gatk_SOR', 'gatk_DP']:
    print(gatk_var)
    fig = figure(figsize=(8, 8))
    for i, mode in enumerate(['SNP', 'INDEL']):
        for j, IsCoding in enumerate(['True', 'False']):
            ax = fig.add_subplot(2, 2, (i*2)+j+1)
            array_to_plot = (tbl_whole_genome_comparisons[isolate_code]
                .selecteq('mode', mode)
                .selecteq('RegionType', 'Core')
                .selecteq('IsCoding', IsCoding)
                .selecteq('IsAccessible', True)
                .selectnotnone(gatk_var)
                .cut([gatk_var])
            ).toarray()
            sns.distplot(array_to_plot[gatk_var], kde=False, ax=ax)
            ax.set_title("%s %s %s" % (gatk_var, mode, 'Coding' if IsCoding=='True' else 'Non-coding'))
#             , bins=np.arange(-2, 14, 0.1))
#             ax.set_xlim([-2, 14])
        

# <codecell>

fig = figure(figsize=(8, 8))
for i, mode in enumerate(['SNP', 'INDEL']):
    for j, IsCoding in enumerate(['True', 'False']):
        ax = fig.add_subplot(2, 2, (i*2)+j+1)
        array_to_plot = (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('mode', mode)
            .selecteq('RegionType', 'Core')
            .selecteq('IsCoding', IsCoding)
            .selecteq('IsAccessible', True)
            .selectnotnone('gatk_QUAL')
            .cut(['gatk_QUAL'])
        ).toarray()
        sns.distplot(array_to_plot['gatk_QUAL'], kde=False, ax=ax)
#                      , bins=np.arange(-2, 14, 0.1))
#         ax.set_xlim([-2, 14])
        

# <codecell>

fig.tight_layout()

# <codecell>

isolate_code='All'
RegionType='Core'
IsCoding=True
mode='SNP'
min_Truth_distance_nearest=100
(tbl_whole_genome_comparisons[isolate_code]
            .selecteq('mode', mode)
            .valuecounts('RegionType')
#             .selecteq('RegionType', RegionType)
#             .selecteq('IsCoding', 'True')
#             .selecteq('IsAccessible', True)
#             .addfield('TruthUnclustered', lambda rec: False if rec['Truth_distance_nearest'] is None else rec['Truth_distance_nearest'] >= min_Truth_distance_nearest)
#             .cut(['IsInGATK', 'IsInTruth', 'TruthUnclustered', 'gatk_VQSLOD'])
        )

# <codecell>

tbl_whole_genome_comparisons[isolate_code].valuecounts('RegionType')

# <headingcell level=1>

# Process all data

# <codecell>

rewrite=False
tbl_comparison_sample_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_comparison_All")
if not os.path.exists(tbl_comparison_sample_cache_fn) or rewrite:
    tbl_whole_genome_comparisons['All'] = (
        tbl_whole_genome_comparisons['7G8'].addfield('Sample', '7G8')
        .cat(tbl_whole_genome_comparisons['GB4'].addfield('Sample', 'GB4'))
        .cat(tbl_whole_genome_comparisons['KE01'].addfield('Sample', 'KE01'))
        .cat(tbl_whole_genome_comparisons['KH02'].addfield('Sample', 'KH02'))
        .cat(tbl_whole_genome_comparisons['GN01'].addfield('Sample', 'GN01'))
        .sort('gatk_VQSLOD', reverse=True)
    )
    etl.topickle(tbl_whole_genome_comparisons['All'], tbl_comparison_sample_cache_fn)
else:
    tbl_whole_genome_comparisons['All'] = etl.frompickle(tbl_comparison_sample_cache_fn)


# <codecell>

tbl_whole_genome_comparisons['7G8']

# <codecell>

tbl_whole_genome_comparisons['All']

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == 2)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .addfield('indel_length', lambda rec: (len(rec['REF']) - len(rec['ALT'])))
    .valuecounts('indel_length')
    .display(20)
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .addfield('indel_length', lambda rec: (len(rec['REF']) - len(rec['ALT'])))
    .valuecounts('indel_length')
    .display(20)
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == 1)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == 1)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == -1)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == -1)
    .valuecounts('IsInGATK', 'IsInTruth')
)

# <codecell>

(tbl_whole_genome_comparisons['All']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_alt')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .selectgt('gatk_VQSLOD', 6)
    .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_alt')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .selectgt('gatk_VQSLOD', 6)
    .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#     .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_het')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .selectgt('gatk_VQSLOD', 0)
    .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#     .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_het')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .selectgt('gatk_VQSLOD', 2)
    .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#     .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_het')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .selectlt('gatk_VQSLOD', 0)
    .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#     .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_het')
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
#     .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
#     .selecteq('IsCoding', IsCoding)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
#     .selectlt('gatk_VQSLOD', 0)
    .addfield('is_het', lambda rec: None if rec['gatk_GT'] is None else rec['gatk_GT'][0] != rec['gatk_GT'][2])
#     .addfield('is_alt', lambda rec: rec['gatk_AD'][1] > rec['gatk_AD'][0])
    .valuecounts('is_het')
)

# <codecell>

tbl_whole_genome_comparisons['7G8']

# <codecell>

array_all = etl.toarray(
    tbl_whole_genome_comparisons['All']
    .replace('Truth_distance_nearest', None, 0)
    .cut(['CHROM', 'POS', 'IsInTruth', 'IsInGATK', 'RegionType', 'IsCoding', 'Effect', 'mode', 'IsAccessible', 'gatk_QUAL',
          'Truth_distance_nearest'])
)

# <codecell>

array_all['Truth_distance_nearest']

# <codecell>

1312/66111

# <codecell>

6336 / 55406

# <codecell>

for isolate_code in isolate_codes:
    print(isolate_code)
    (tbl_whole_genome_comparisons[isolate_code]
        .selecteq('mode', 'INDEL')
        .selecteq('RegionType', 'Core')
    #     .selecteq('IsCoding', IsCoding)
        .selecteq('IsAccessible', True)
        .selecteq('IsInGATK', True)
        .select(lambda rec: (len(rec['REF']) - len(rec['ALT'])) == -1)
        .valuecounts('IsInGATK', 'IsInTruth')
        .display()
    )

# <codecell>

tbl_whole_genome_comparisons['7G8'].selecteq('IsAccessible', True).valuecounts('IsInTruth')

# <codecell>

(tbl_whole_genome_comparisons['All'].valuecounts('REF').display(20)

# <headingcell level=1>

# Ti/Tv of clustered variants

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

def calc_titv(rec_list):
    ti = np.count_nonzero([is_transition(rec) for rec in rec_list])
    tv = np.count_nonzero([is_transversion(rec) for rec in rec_list])
    if tv == 0:
        titv = 0.0
    else:
        titv = 1.0*ti/tv
#     print(ti, tv, titv)
    return(titv)

def calc_ti(rec_list):
    return np.count_nonzero([is_transition(rec) for rec in rec_list])
def calc_tv(rec_list):
    return np.count_nonzero([is_transversion(rec) for rec in rec_list])

# <codecell>

tbl_whole_genome_comparisons['7G8'].valuecounts('IsCoding')

# <codecell>

def calc_variant_summary_table(isolate_code='7G8', min_Truth_distance_nearest=0, min_VQSLOD=-99999):
    aggregation = collections.OrderedDict()
    aggregation['# Variants'] = 'IsInTruth', sum
    aggregation['# Coding'] = ('IsInTruth', 'IsCoding'), lambda rec: np.sum([x[0] and x[1]=='True' for x in list(rec)])
    aggregation['# Coding indels'] = ('REF', 'ALT', 'IsInTruth', 'IsCoding'), lambda rec: np.sum([x[2] and x[3]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) for x in list(rec)])
#     aggregation['#GATK'] = 'IsInGATK', sum
#     aggregation['#GATK accessible'] = ('IsInGATK', 'IsAccessible'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
#     aggregation['#TP'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([x[0] and x[1] for x in list(rec)])
#     aggregation['#FP'] = ('IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([x[0] and not x[1] and x[2] for x in list(rec)])
#     aggregation['#FN'] = ('IsInGATK', 'IsInTruth'), lambda rec: np.sum([not x[0] and x[1] for x in list(rec)])
    aggregation['Ti'] = ('REF', 'ALT'), lambda rec: calc_ti(list(rec))
    aggregation['Tv'] = ('REF', 'ALT'), lambda rec: calc_tv(list(rec))
    aggregation['TiTv'] = ('REF', 'ALT'), lambda rec: calc_titv(list(rec))
#     aggregation['TiTv FP'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FP')
#     aggregation['TiTv FN'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: calc_titv(list(rec), 'FN')
    aggregation['# Coding indels mod3'] = ('REF', 'ALT', 'IsInTruth', 'IsCoding'), lambda rec: np.sum([x[2] and x[3]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) and (len(x[0]) - len(x[1])) % 3 == 0 for x in list(rec)])
#     aggregation['#FPmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth', 'IsAccessible'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and x[2] and not x[3] and x[4] for x in list(rec)])
#     aggregation['#FNmod3'] = ('REF', 'ALT', 'IsInGATK', 'IsInTruth'), lambda rec: np.sum([(len(x[0]) - len(x[1])) % 3 == 0 and not x[2] and x[3] for x in list(rec)])

    tbl_variant_summary = (tbl_whole_genome_comparisons[isolate_code]
        .selecteq('Truth_RegionType', 'Core')
#         .selecteq('IsInGATK', True)
#         .selecteq('Truth_mode', 'SNP')
#         .addfield('Truth_distance_100', lambda rec: 100 if rec['Truth_distance_nearest'] > 1000 else int(rec['Truth_distance_nearest']/10))
        .addfield('Truth_distance_log10', lambda rec: 2 if rec['Truth_distance_nearest'] > 100 else 0 if rec['Truth_distance_nearest'] < 1 else int(log10(rec['Truth_distance_nearest'])))
#         .addfield('Truth_distance_gt100', lambda rec: True if rec['Truth_distance_nearest'] > 100 else False)
        .aggregate(('Truth_distance_log10'), aggregation)
#         .aggregate(('VQSLOD_int'), aggregation) 
#         .aggregate(('RegionType', 'IsCoding', 'mode'), aggregation)
#         .convertnumbers()
#         .addfield('Unfilt sensitivity', lambda rec: 0.0 if rec['#Truth'] == 0 else round(rec['#TP'] / rec['#Truth'], 3))
#         .addfield('Unfilt FDR', lambda rec: 0.0 if rec['#GATK accessible'] == 0 else round(rec['#FP'] / rec['#GATK accessible'], 3))
        .addfield('% Frameshift indels', lambda rec: 0.0 if rec['# Coding indels'] == 0 else round(1.0 - (rec['# Coding indels mod3'] / rec['# Coding indels']), 3))
#         .addfield('%FPmod3', lambda rec: 0.0 if rec['#FP'] == 0 else round(rec['#FPmod3'] / rec['#FP'], 3))
#         .addfield('%FNmod3', lambda rec: 0.0 if rec['#FN'] == 0 else round(rec['#FNmod3'] / rec['#FN'], 3))
    )
    
    return(tbl_variant_summary)

# <codecell>

def calc_region_summary_table(isolate_code='All', min_Truth_distance_nearest=0, min_VQSLOD=-99999):
    aggregation = collections.OrderedDict()
    aggregation['# Variants'] = 'IsInTruth', sum
    aggregation['# Coding'] = ('IsInTruth', 'IsCoding'), lambda rec: np.sum([x[0] and x[1]=='True' for x in list(rec)])
    aggregation['# Coding indels'] = ('REF', 'ALT', 'IsInTruth', 'IsCoding'), lambda rec: np.sum([x[2] and x[3]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) for x in list(rec)])
    aggregation['Ti'] = ('REF', 'ALT'), lambda rec: calc_ti(list(rec))
    aggregation['Tv'] = ('REF', 'ALT'), lambda rec: calc_tv(list(rec))
    aggregation['TiTv'] = ('REF', 'ALT'), lambda rec: calc_titv(list(rec))
    aggregation['# Coding indels mod3'] = ('REF', 'ALT', 'IsInTruth', 'IsCoding'), lambda rec: np.sum([x[2] and x[3]=='True' and (len(x[0]) > 1 or len(x[1]) > 1 ) and (len(x[0]) - len(x[1])) % 3 == 0 for x in list(rec)])

    tbl_variant_summary = (tbl_whole_genome_comparisons[isolate_code]
        .selectnotnone('Truth_RegionType')
#         .selectge('Truth_distance_nearest', 100)
#         .addfield('Truth_distance_log10', lambda rec: 2 if rec['Truth_distance_nearest'] > 100 else 0 if rec['Truth_distance_nearest'] < 1 else int(log10(rec['Truth_distance_nearest'])))
        .aggregate(('Truth_RegionType'), aggregation)
        .addfield('% Frameshift indels', lambda rec: 0.0 if rec['# Coding indels'] == 0 else round(1.0 - (rec['# Coding indels mod3'] / rec['# Coding indels']), 3))
    )
    
    return(tbl_variant_summary)

# <codecell>

calc_region_summary_table('7G8')

# <codecell>

tbl_region_summary_all = calc_region_summary_table()

# <codecell>

tbl_region_summary_all

# <codecell>

tbl_truth_clustering_summary = collections.OrderedDict()
for isolate_code in ['All'] + isolate_codes:
    print(isolate_code)
    tbl_truth_clustering_summary[isolate_code] = calc_variant_summary_table(isolate_code)
    tbl_truth_clustering_summary[isolate_code].displayall()

# <codecell>

isolate_codes + ['All']

# <codecell>

df_truth_clustering_summary = (
    tbl_truth_clustering_summary['7G8'].addfield('Sample', '7G8')
    .cat(tbl_truth_clustering_summary['GB4'].addfield('Sample', 'GB4'))
    .cat(tbl_truth_clustering_summary['KE01'].addfield('Sample', 'KE01'))
    .cat(tbl_truth_clustering_summary['KH02'].addfield('Sample', 'KH02'))
    .cat(tbl_truth_clustering_summary['GN01'].addfield('Sample', 'GN01'))
    .cat(tbl_truth_clustering_summary['All'].addfield('Sample', 'All'))
    .addfield('Distance to nearest', lambda rec: '1-9' if rec['Truth_distance_log10'] == 0 else '10-99' if rec['Truth_distance_log10'] == 1 else '100+')
    .todataframe()
)

# <codecell>

df_truth_clustering_summary

# <codecell>

fig = figure(figsize=(8, 6))
ax = fig.add_subplot(2, 1, 1)
g = sns.barplot(x="Distance to nearest", y="TiTv", hue="Sample", ci=None,
                data=df_truth_clustering_summary)
g.legend(loc="lower center")
g.set_ylabel("Ti/Tv")
g.set_xlabel("")
ax = fig.add_subplot(2, 1, 2)
g = sns.barplot(x="Distance to nearest", y="% Frameshift indels", hue="Sample", ci=None,
                data=df_truth_clustering_summary)
g.legend(loc="lower center")
g.set_ylabel("% Frameshift indels")
g.set_xlabel("")
fig.tight_layout()

# <codecell>

fig = figure(figsize=(8, 6))
ax = fig.add_subplot(2, 1, 1)
g = sns.barplot(x="Sample", y="TiTv", hue="Distance to nearest", ci=None,
                data=df_truth_clustering_summary)
ax.set_ylim([0.7, 1.1])
g.legend(loc="center", title='Distance to nearest')
g.set_ylabel("Ti/Tv")
g.set_xlabel("")
ax = fig.add_subplot(2, 1, 2)
g = sns.barplot(x="Sample", y="% Frameshift indels", hue="Distance to nearest", ci=None,
                data=df_truth_clustering_summary)
g.legend(loc="center", title='Distance to nearest')
# g.legend(loc=None)
g.set_ylabel("% Frameshift indels")
g.set_xlabel("")
fig.tight_layout()
fig.savefig('%s/truth_clustering_summary.pdf' % PLOTS_DIR, format='pdf')

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
    ax = plotFDRsens(isolate_code, show_annotations=show_annotations, xmax=0.17, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotFDRsens(isolate_code, show_annotations=show_annotations, IsCoding='False', xmax=None, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotFDRsens(isolate_code, show_annotations=show_annotations, mode='INDEL', IsCoding='True', xmax=None, ax=ax)
    ax.legend(loc="lower right")

# <codecell>

fig = figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
for isolate_code in isolate_codes:
    if isolate_code=='KE01':
        show_annotations=True
    else:
        show_annotations=False
    ax = plotFDRsens(isolate_code, show_annotations=show_annotations, mode='INDEL', IsCoding='False', xmax=None, ax=ax)
    ax.legend(loc="lower right")

# <codecell>


# <codecell>


# <codecell>


# <codecell>


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


