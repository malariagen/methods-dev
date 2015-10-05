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
OUTPUT_DIR = '/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment'
!mkdir -p {OUTPUT_DIR}
NUCMER_DIR='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23'

# <codecell>

# CODING_EFFECTS = ['NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED']
NONCODING_EFFECTS = ['INTERGENIC', 'INTRON', 'TRANSCRIPT']
 

# <codecell>

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))

# <codecell>

ref_dict['Pf3D7_10_v3']

# <codecell>

# location of CRT "CVMNK" etc haplotypes 
ref_dict['Pf3D7_07_v3'][403610:403633]

# <codecell>

NUCMER_DIR

# <headingcell level=1>

# Load data

# <codecell>

tbl_pipelines = etl.fromxlsx(PIPELINES_FN)
tbl_pipelines.display(index_header=True)

# <codecell>


# <headingcell level=1>

# Determine samples

# <codecell>

def determine_assembly_fasta(isolate_code='7G8'):
    if isolate_code in ('7G8', 'GB4'):
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf%s.Jul2015.fasta" % isolate_code)
    else:
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/Pf%s.Jun2015.fasta" % isolate_code)

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
    .addfield('assembly_fasta', lambda rec: determine_assembly_fasta(rec['Isolate code']))
)
tbl_samples_to_process.displayall(index_header=True)

# <codecell>

assembly_fastas = etl.lookupone(tbl_samples_to_process, 'Isolate code', 'assembly_fasta')
print(assembly_fastas['7G8'])
print(assembly_fastas['KH02'])

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
    if prv is None or prv.CHROM != cur.CHROM:
        return None
    else:
        return cur.POS - prv.POS

def downstream(prv, cur, nxt):
    if nxt is None or nxt.CHROM != cur.CHROM:
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

# def nearest_consensus(prv, cur, nxt):
    

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

(tbl_whole_genome_comparisons['alistair_ann']
.selecteq('Sample', '7G8')
.selecteq('IsInGATK', True)
.selecteq('CHROM', 'Pf3D7_01_v3')
.sort('POS')
.selectgt('POS', 17816)
.display(20)
)

# <codecell>

(tbl_whole_genome_comparisons['alistair_ann']
.selecteq('Sample', '7G8')
.selecteq('IsInGATK', True)
.selecteq('CHROM', 'Pf3D7_01_v3')
.sort('POS')
.selectgt('POS', 120000)
.cut(['CHROM', 'POS', 'REF', 'ALT', 'gatk_VQSLOD'])
.display(20)
)

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

# Create consensus alignments

# <codecell>

import copy
def create_consensus(sample='7G8',
                     pos_output_format = '/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/pos/%s.pos',
                     fasta_output_format = '/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/fasta/%s.fasta'):
    pos_dir = os.path.dirname(pos_output_format % sample)
    if not os.path.exists(pos_dir):
        os.makedirs(pos_dir)
    fasta_dir = os.path.dirname(fasta_output_format % sample)
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)
    pos_fo = open(pos_output_format % sample, 'w')
    fasta_fo = open(fasta_output_format % sample, 'w')
    mutable_ref_dict = copy.deepcopy(ref_dict)
    for chrom in mutable_ref_dict:
        mutable_ref_dict[chrom].seq = mutable_ref_dict[chrom].seq.tomutable()
    previous_new_pos = 1
    current_offset = 0
    current_chrom = ''
    for i, rec in enumerate(tbl_whole_genome_comparisons['alistair_ann']
        .selecteq('Sample', sample)
        .selecteq('IsInGATK', True)
        .selecteq('RegionType', 'Core')
        .selectgt('gatk_VQSLOD', 0.0)
        .sort(('CHROM', 'POS'))
        .data()
    ):
        if current_chrom != rec[0]:
            print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  
            current_chrom = rec[0]
            current_offset = 0
            previous_new_pos = 1
    #     if current_chrom == '':
    #         current_chrom = rec[0]
        reflen = len(rec[2])
        altlen = len(rec[3])
        pos = int(rec[1])
        new_pos = pos + current_offset
        print("%s\t%d\t%d\t%d\t%d" % (rec[0], pos, new_pos, previous_new_pos, current_offset), file=pos_fo)  
        previous_new_pos = new_pos
        startpos = pos + current_offset - 1
        endpos = pos + current_offset + reflen - 1   
        mutable_ref_dict[rec[0]].seq[startpos:endpos] = rec[3]
        current_offset = current_offset + altlen - reflen
#         if i%1000 == 0:
#             print(current_chrom, current_offset)

    print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  
    
    for chrom in mutable_ref_dict:
        SeqIO.write(mutable_ref_dict[chrom], fasta_fo, "fasta")
    

# <codecell>

for sample in tbl_samples_to_process['Isolate code']:
    print(sample)
    create_consensus(sample)

# <codecell>

def run_nucmer(sample='7G8',
               fasta_output_format='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/fasta/%s.fasta',
               nucmer_output_format='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/nucmer/%s',
               delta_filter_i=99):
    gatk_fasta = fasta_output_format % sample
    assembly_fasta = assembly_fastas[sample]
    output_filestem = nucmer_output_format % sample
    
    !{NUCMER_DIR}/nucmer -p {output_filestem} {gatk_fasta} {assembly_fasta}
    !{NUCMER_DIR}/delta-filter -m -i {delta_filter_i} -l 1000 {output_filestem}.delta > {output_filestem}.filter.delta
    !{NUCMER_DIR}/show-coords -THqcl -o {output_filestem}.filter.delta > {output_filestem}.filter.coords
    !{NUCMER_DIR}/show-snps -CHTr {output_filestem}.filter.delta > {output_filestem}.Csnp

# <codecell>

run_nucmer()

# <codecell>

for sample in tbl_samples_to_process['Isolate code']:
    print(sample)
    run_nucmer(sample)

# <codecell>

ref_chroms = ['M76611', 'PFC10_API_IRAB'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
query_chroms = ["Pf%s_MT" % '7G8', "Pf%s_API" % '7G8'] + ["Pf%s_%02d" % ('7G8', i) for i in range(1, 15)]
chrom_dict = dict(zip(ref_chroms, query_chroms))
print(chrom_dict['M76611'])
print(chrom_dict['Pf3D7_08_v3'])

temp = dict(zip(('M76611', 'PFC10_API_IRAB'), ('Pf7G8_MT', 'Pf7G8_API')))
temp['M76611']

# <codecell>

def convert_Csnp(sample='7G8',
                 fasta_output_format='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/fasta/%s.fasta',
                 nucmer_output_format='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/nucmer/%s'):
    gatk_fasta = fasta_output_format % sample
    assembly_fasta = assembly_fastas[sample]
    output_filestem = nucmer_output_format % sample

    Csnp_fn=output_filestem + ".Csnp"
    vcf_fn=output_filestem + ".vcf"
    ref_dict=SeqIO.to_dict(SeqIO.parse(open(gatk_fasta), "fasta"))

#     Write VCF header
    fo = open(vcf_fn, 'w')
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("##description=This file created with convert_Csnp function in 20150929_consensus_alignment.ipynb\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
#     Write data
    ref_chroms = ['M76611', 'PFC10_API_IRAB'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
    query_chroms = ["Pf%s_MT" % sample, "Pf%s_API" % sample] + ["Pf%s_%02d" % (sample, i) for i in range(1, 15)]
    chrom_dict = dict(zip(ref_chroms, query_chroms))

    tbl_csnp = (etl.fromtsv(Csnp_fn).select(lambda rec: rec[9] == chrom_dict[rec[8]]))
    
    current_state = "SNP"
    ins_sequence = ''
    previous_chrom = ''
    for rec in tbl_csnp:
#         print('.', end='')
        (pos, ref, alt, pos2, buff, dist, frm, frm2, chrom, query_chrom) = rec
#         (chrom, bba, variant_type, pos, pos2, zero, strand, dot, note) = rec
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
#             ref = ref_dict[chrom][pos-1]
#             alt = note[20]
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
#             ins_sequence = ins_sequence + note[26]
            ins_sequence = ins_sequence + alt
            variant_end = pos
            current_state = "Ins"
#         elif variant_type == 'Synteny':
#             if current_state == 'Del':
#                 ref = str(ref_dict[chrom][variant_start-1:variant_end].seq)
#                 alt = ref_dict[chrom][variant_start]
#                 fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
#             if current_state == 'Ins':
#                 ref = ref_dict[chrom][variant_start-1]
#                 alt = ref_dict[chrom][variant_start-1] + ins_sequence
#                 fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
#                 ins_sequence = ''
#             current_state = "Synteny"
        else:
            return("error")
        previous_chrom = chrom
    fo.close()
    
    return(vcf_fn)

# <codecell>

convert_Csnp()

# <codecell>

def tbl_converted_coords(sample='7G8',
        pos_output_format = '/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/pos/%s.pos',
        nucmer_output_format='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/nucmer/%s'):
    
    pos_fn = pos_output_format % sample
    output_filestem = nucmer_output_format % sample
    vcf_fn=output_filestem + ".vcf"
    output_vcf_fn=output_filestem + ".3d7coordinates.vcf"
    
    tbl_vcf = etl.fromvcf(vcf_fn).convert('POS', int).addfield('POS1', lambda rec: rec['POS']+1)
    tbl_pos = (etl.fromtsv(pos_fn)
        .setheader(('CHROM', 'var_pos', 'end', 'start', 'diff'))
        .convert('var_pos', int)
        .convert('start', int)
        .convert('end', int)
        .convert('diff', int)
        .select(lambda rec: rec['end'] > rec['start'])
    )
    
#     tbl_vcf.display()
#     tbl_pos.display()
    
    tbl_vcf_converted = (tbl_vcf
        .intervalleftjoin(tbl_pos, lstart='POS', lstop='POS1', lkey='CHROM',
                          rstart='start', rstop='end', rkey='CHROM')
        .addfield('NEW_POS', lambda rec: rec['POS'] if rec['diff'] is None else rec['POS'] - rec['diff'])
    )
#     tbl_vcf_converted.display()

#     Write VCF header
    fo = open(output_vcf_fn, 'w')
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("##description=This file created with tbl_converted_coords function in 20150929_consensus_alignment.ipynb\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for rec in tbl_vcf_converted.data():
        fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (rec[0], rec[13], rec[3], rec[4][0]))
        
    return tbl_vcf_converted
    
    

# <codecell>

tbl_vcf_converted = tbl_converted_coords()

# <codecell>

tbl_vcf_converted.display(50, index_header=True)

# <codecell>

print(ref_dict['Pf3D7_01_v3'][102910 - 1]) #A
print(ref_dict['Pf3D7_01_v3'][123620 - 1]) #C
print(ref_dict['Pf3D7_01_v3'][139447 - 1]) #C
print(ref_dict['Pf3D7_01_v3'][141192 - 1]) #C
print(ref_dict['Pf3D7_01_v3'][147835 - 1]) #T
print(ref_dict['Pf3D7_01_v3'][188251 - 1]) #T

# <headingcell level=1>

# Install vcflib

# <codecell>

# On mac, simply:
# brew install vcflib

# <codecell>


# <codecell>


# <headingcell level=1>

# Annotate consensus alignment vcf

# <codecell>

rewrite=True

for sample in tbl_samples_to_process['Isolate code']:
    print(sample)
    nucmer_output_filestem='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/nucmer/%s' % sample
    
    unannotated_vcf_fn = "%s.3d7coordinates.vcf" % (nucmer_output_filestem)
    sorted_vcf_fn = "%s.sorted.vcf" % (nucmer_output_filestem)
    left_aligned_vcf_fn = "%s.leftaligned.vcf" % (nucmer_output_filestem)
    snpeff_vcf_fn = "%s.snpeff.vcf" % (nucmer_output_filestem)
    snpeff_annotated_vcf_fn = "%s.snpeff_annotated.vcf" % (nucmer_output_filestem)
    annotated_vcf_fn = "%s.annotated.vcf" % (nucmer_output_filestem)

    convert_Csnp(sample)
    tbl_converted_coords(sample)
    
    !vcfstreamsort < {unannotated_vcf_fn} > {sorted_vcf_fn}

    if (not os.path.isfile(left_aligned_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        !{gatk_exe} -T LeftAlignAndTrimVariants \
        -R {REF_GENOME} \
        -V {sorted_vcf_fn} \
        -o {left_aligned_vcf_fn} \
#         2> /dev/null

    if (not os.path.isfile(snpeff_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        !{snpeff_exe} \
        -v -o gatk Pf3D7july2015 \
        {left_aligned_vcf_fn} \
        -no-downstream \
        -no-upstream \
        > {snpeff_vcf_fn} \
#         2> /dev/null

    if (not os.path.isfile(snpeff_annotated_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        !{gatk_exe} \
        -T VariantAnnotator \
        -R {REF_GENOME} \
        -A TandemRepeatAnnotator \
        -A SnpEff \
        --variant {left_aligned_vcf_fn} \
        --snpEffFile {snpeff_vcf_fn} \
        -o {snpeff_annotated_vcf_fn} \
#         2> /dev/null

    if not os.path.isfile(annotated_vcf_fn+'.gz') or rewrite:
        !cat {snpeff_annotated_vcf_fn} \
        | vcf-annotate -a {regions_fn} \
           -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
           -c CHROM,FROM,TO,INFO/RegionType \
        > {annotated_vcf_fn}

        !bgzip -f {annotated_vcf_fn}
        !tabix -p vcf -f {annotated_vcf_fn}.gz

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.SNP.vcf')+'.gz') or rewrite:
        !{gatk_exe} \
        -T SelectVariants \
        -R {REF_GENOME} \
        -V {annotated_vcf_fn}.gz \
        -selectType SNP \
        -o {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')} 2> /dev/null

        !bgzip -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}
        !tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}.gz

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')+'.gz') or rewrite:
        !{gatk_exe} \
        -T SelectVariants \
        -R {REF_GENOME} \
        -V {annotated_vcf_fn}.gz \
        -xlSelectType SNP \
        -o {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')} 2> /dev/null

        !bgzip -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}
        !tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}.gz

#         !rm {unannotated_vcf_fn}
#         !rm {left_aligned_vcf_fn}
#         !rm {snpeff_vcf_fn}
#         !rm {snpeff_annotated_vcf_fn}
#         !rm {unannotated_vcf_fn}.idx
#         !rm {left_aligned_vcf_fn}.idx
#         !rm {snpeff_vcf_fn}.idx
#         !rm {snpeff_annotated_vcf_fn}.idx

# <codecell>

def merge_consensus_alignment(sample='7G8',
                      rewrite=False):
    ox_code = tbl_samples_to_process.selecteq('Isolate code', sample).values('ox_code')[0]
    nucmer_output_filestem='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/consensus_alignment/nucmer/%s' % sample
    annotated_vcf_fn = "%s.annotated.vcf.gz" % (nucmer_output_filestem)

    
    tbl_consensus = (etl.fromvcf(annotated_vcf_fn)
        .unpackdict('INFO', samplesize=1000000)
        .rename('ALT', 'consensus_ALTs')
        .addfield('ALT', lambda rec: str(rec[4][0]))
        .addfield('Consensus_is_coding', lambda rec: not rec['SNPEFF_EFFECT'] is None and not rec['SNPEFF_EFFECT'] in NONCODING_EFFECTS)
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'consensus_ALTs', 'RegionType', 'Consensus_is_coding', 'SNPEFF_AMINO_ACID_CHANGE', 'SNPEFF_EFFECT',
              'SNPEFF_GENE_NAME', 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID', 'RPA', 'RU', 'STR'])
        .addfield('Consensus_RPA_REF', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][0]))
        .addfield('Consensus_RPA_ALT', lambda rec: None if rec['RPA'] is None else str(rec['RPA'][1]))
        .rename('RegionType', 'Consensus_RegionType')
        .rename('SNPEFF_AMINO_ACID_CHANGE', 'Consensus_AAChange')
        .rename('SNPEFF_EFFECT', 'Consensus_Effect')
        .rename('SNPEFF_GENE_NAME', 'Consensus_Gene')
        .rename('SNPEFF_IMPACT', 'Consensus_Impact')
        .rename('SNPEFF_TRANSCRIPT_ID', 'Consensus_Transcript')
        .addfieldusingcontext('Consensus_distance_previous', upstream)
        .addfieldusingcontext('Consensus_distance_next', downstream)
        .addfield('Consensus_distance_nearest', lambda rec: nearest(rec['Consensus_distance_previous'], rec['Consensus_distance_next']))
        .addfield('Consensus_mode', determine_mode)
    )
        
    return(tbl_consensus)

# <codecell>

tbl_consensus_7G8 = merge_consensus_alignment()

# <codecell>

tbl_consensus_7G8.display(10, index_header=True)

# <codecell>

(tbl_consensus_7G8
    .addfield('is_clustered', lambda rec: rec['Consensus_distance_nearest'] <= 10)
    .valuecounts('STR', 'is_clustered')
)

# <codecell>

(tbl_consensus_7G8
    .addfield('is_clustered', lambda rec: rec['Consensus_distance_nearest'] <= 20)
    .valuecounts('Consensus_mode', 'STR', 'is_clustered')
    .displayall()
)

# <codecell>

(tbl_whole_genome_comparisons['alistair_ann']
.selecteq('Sample', '7G8')
.selecteq('IsInGATK', False)
.addfield('is_clustered', lambda rec: rec['Truth_distance_nearest'] <= 20)
.valuecounts('Truth_mode', 'STR', 'is_clustered')
.display(20)
)

# <codecell>

tbl_whole_genome_comparisons['alistair_ann']

# <codecell>

(tbl_consensus_7G8
    .addfield('is_clustered', lambda rec: rec['Consensus_distance_nearest'] <= 100)
    .valuecounts('Consensus_mode', 'STR', 'is_clustered')
    .displayall()
)

# <codecell>

len(tbl_consensus_7G8)

# <codecell>

consensus_distances_7G8 = tbl_consensus_7G8['Consensus_distance_nearest'].array()

# <codecell>

ax = hist(consensus_distances_7G8, bins=arange(0, 100, 1))

# <codecell>

len(consensus_distances_7G8)

# <codecell>

np.sum(consensus_distances_7G8>10)

# <codecell>

arange(0, 30, 1)

# <codecell>

tbl_consensus_7G8.selecteq('Consensus_distance_nearest', 1).display(10)

# <codecell>


# <markdowncell>

#     - create new version of pos file with start, end, diff of new coordinates to be used for converting back
#     - create truth vcf from .Csnp file
#     - annotate truth vcf (convert to original coordinates then as 20150918_truth_vcfs)
# - create new table merging tbl_whole_genome_comparisons, pos file, truth vcf, accessible (overlap with filter-coords), distance of each GATK to nearest truth
# - plot histograms of distance of GATK variants to nearest truth and set sensible threshold
# - write function to do everything (including nucmer), with params for distance to assume TP, delta-filter i, 
# - rerun everything with much lower VQSLOD, to sanity check many more FPs
# - create ROC curves based on new definitions of TP and FP
# - recreate everything from scratch, doing only what is necessary to get FDR and sens for SNP/INDEL, coding/non-coding, ensuring hets appropriately dealt with, and making as fast as possible (use numpy arrays throughout?)
# - check comparable results from 5 samples vcf and full vcf (how best to read in full vcf - talk to Alistair? maybe create small vcfs just from variant positions for that sample?)
# - turn pipeline into a single script which takes VCF, sample name, assembly file, TP_distance and delta_filter_i, and returns #TP, #FP, #FN, #inaccessible broken down by SNP/INDEL, coding/non-coding. Is this really worth it?

# <codecell>


# <headingcell level=1>

# Testing stuff

# <codecell>

mutable_ref_dict = copy.deepcopy(ref_dict)
print(mutable_ref_dict['Pf3D7_01_v3'])
mutable_ref_dict['Pf3D7_01_v3'].seq.tomutable()
print(mutable_ref_dict['Pf3D7_01_v3'])
mutable_ref_dict['Pf3D7_01_v3'].seq = mutable_ref_dict['Pf3D7_01_v3'].seq.tomutable()
print(mutable_ref_dict['Pf3D7_01_v3'])

# <codecell>

print(mutable_ref_dict['Pf3D7_01_v3'].id)

# <codecell>

print(ref_dict['Pf3D7_01_v3'].seq[93150:93170])
print(mutable_ref_dict['Pf3D7_01_v3'][93150:93170])
print()
print(ref_dict['Pf3D7_01_v3'].seq[94250:94370])
print(mutable_ref_dict['Pf3D7_01_v3'][94250:94370])

# <codecell>

print(len(ref_dict['Pf3D7_01_v3']))
print(len(mutable_ref_dict['Pf3D7_01_v3']))

# <codecell>

print(ref_dict['Pf3D7_01_v3'].seq[-10:-1])

# <codecell>

print(mutable_ref_dict['Pf3D7_01_v3'][-10:-1])

# <codecell>

mutable_ref_dict['Pf3D7_01_v3']

# <codecell>

ref_dict['Pf3D7_01_v3']

# <codecell>

temp = ref_dict.copy()

# <codecell>

for chrom in temp:
    temp[chrom].seq = temp[chrom].seq.tomutable()

# <codecell>

temp['Pf3D7_01_v3']

# <codecell>


# <codecell>

def second_allele_is_major(rec):
    if (rec['gatk_AD'][int(rec['gatk_GT'][2]) - 1]) >= int(rec['gatk_AD'][int(rec['gatk_GT'][0]) -1]):
        return True
    else:
        return False

# <codecell>

(tbl_whole_genome_comparisons['alistair_ann']
.selecteq('Sample', '7G8')
.selecteq('IsInGATK', True)
.selecteq('RegionType', 'Core')
.selectgt('gatk_VQSLOD', 0.0)
.addfield('is_2nd_major', second_allele_is_major)
# .sort(('CHROM', 'POS'))
# .valuecounts('is_het', 'IsInTruth', 'is_2nd_major', 'mode')
.valuecounts('is_het', 'mode')
.displayall()
)

# <codecell>

(tbl_whole_genome_comparisons['alistair_ann']
    .selecteq('Sample', '7G8')
    .selecteq('IsInGATK', True)
    .selecteq('RegionType', 'Core')
    .selectgt('gatk_VQSLOD', 0.0)
#     .valuecounts('gatk_GT').displayall()
    .sort(('CHROM', 'POS'))
    .display(index_header=True)
)

# <codecell>

(tbl_whole_genome_comparisons['alistair_ann']
    .selecteq('Sample', '7G8')
    .selecteq('IsInGATK', True)
    .selecteq('RegionType', 'Core')
    .selectgt('gatk_VQSLOD', 0.0)
#     .valuecounts('gatk_GT').displayall()
    .sort(('CHROM', 'POS'))
    .selecteq('IsInTruth', False)
    .selecteq('mode', 'SNP')
    .selecteq('is_het', True)
    .display(50, index_header=True)
)

