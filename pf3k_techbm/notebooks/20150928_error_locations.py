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
PLOTS_DIR = '/Users/rpearson/Documents/projects/Pf3k_techbm/slides/20150928_error_locations/plots'
!mkdir -p {PLOTS_DIR}

# <codecell>

# CODING_EFFECTS = ['NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED']
NONCODING_EFFECTS = ['INTERGENIC', 'INTRON', 'TRANSCRIPT']
 

# <codecell>

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))

# <codecell>


# <codecell>


# <headingcell level=1>

# Load data

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


# <headingcell level=1>

# Plot errors along genome

# <codecell>

import allel
print(allel.__version__)

# <codecell>

def plot_variant_densities(pipeline='alistair_ann', mode=['SNP', 'IND'], IsCoding=['True', 'False'],
                           TP_lim=100, FP_lim=50, FN_lim=100):
    fig = figure(figsize=(24, 36))
    for i, chrom in enumerate(['Pf3D7_%02d_v3' % chrom for chrom in range(1,15)]):
        TPs = sort(comparison_arrays[pipeline]['POS'][
            (comparison_arrays[pipeline]['CHROM'] == chrom) &
            (comparison_arrays[pipeline]['RegionType'] == 'Core') &
            (comparison_arrays[pipeline]['IsAccessible'] == True) &
            (comparison_arrays[pipeline]['gatk_VQSLOD'] > 0.0) &
            (np.in1d(comparison_arrays[pipeline]['mode'], mode)) &
            (np.in1d(comparison_arrays[pipeline]['IsCoding'], IsCoding)) &
            (comparison_arrays[pipeline]['IsInTruth'] == True)
        ])
        var_density = allel.stats.window.windowed_count(TPs, size=1000)
        ax = fig.add_subplot(14*4, 1, (i*4)+1)
        ax.plot([(x[0]+x[1])/2 for x in var_density[1]], var_density[0], color='green')
        ax.set_xlim([0, 3.3e+6])
        ax.set_ylim([0, TP_lim])
        ax.set_title(chrom)

        FPs = sort(comparison_arrays[pipeline]['POS'][
            (comparison_arrays[pipeline]['CHROM'] == chrom) &
            (comparison_arrays[pipeline]['RegionType'] == 'Core') &
            (comparison_arrays[pipeline]['IsAccessible'] == True) &
            (comparison_arrays[pipeline]['gatk_VQSLOD'] > 0.0) &
            (np.in1d(comparison_arrays[pipeline]['mode'], mode)) &
            (np.in1d(comparison_arrays[pipeline]['IsCoding'], IsCoding)) &
            (comparison_arrays[pipeline]['IsInTruth'] == False)
        ])
        var_density = allel.stats.window.windowed_count(FPs, size=1000)
        ax = fig.add_subplot(14*4, 1, (i*4)+2)
        ax.plot([(x[0]+x[1])/2 for x in var_density[1]], var_density[0], color='red')
        ax.set_xlim([0, 3.3e+6])
        ax.set_ylim([0, FP_lim])

        FNs = sort(comparison_arrays[pipeline]['POS'][
            (comparison_arrays[pipeline]['CHROM'] == chrom) &
            (comparison_arrays[pipeline]['RegionType'] == 'Core') &
            (comparison_arrays[pipeline]['IsAccessible'] == True) &
            (comparison_arrays[pipeline]['IsInGATK'] == False) &
            (np.in1d(comparison_arrays[pipeline]['mode'], mode)) &
            (np.in1d(comparison_arrays[pipeline]['IsCoding'], IsCoding)) &
            (comparison_arrays[pipeline]['IsInTruth'] == True)
        ])
        var_density = allel.stats.window.windowed_count(FNs, size=1000)
        ax = fig.add_subplot(14*4, 1, (i*4)+3)
        ax.plot([(x[0]+x[1])/2 for x in var_density[1]], var_density[0], color='orange')
        ax.set_xlim([0, 3.3e+6])
        ax.set_ylim([0, FN_lim])



# <codecell>

plot_variant_densities()

# <codecell>

plot_variant_densities(mode=['IND'], IsCoding=['True'], TP_lim=20, FP_lim=20, FN_lim=20)

# <codecell>

plot_variant_densities(mode=['IND'], IsCoding=['False'], TP_lim=50, FP_lim=50, FN_lim=50)

# <codecell>

plot_variant_densities(mode=['SNP'], IsCoding=['True'], TP_lim=50, FP_lim=50, FN_lim=50)

# <codecell>

plot_variant_densities(mode=['SNP'], IsCoding=['False'], TP_lim=50, FP_lim=50, FN_lim=50)

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


