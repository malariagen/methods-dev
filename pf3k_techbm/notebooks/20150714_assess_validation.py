# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

annotated_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated')
annotated_vcf_filestem

# <codecell>

!mkdir {os.path.join(CACHE_DIR, 'tbl_truth')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_gatk')}
!mkdir {os.path.join(CACHE_DIR, 'tbl_comparison')}

# <headingcell level=1>

# Analyse one sample

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
#     .cutout('bam_fn')
    .cutout('To be used for')
    .convert('ox_code', lambda v: 'ERS740937', where=lambda rec: rec['Isolate code'] == 'KE01')
    .convert('ox_code', lambda v: 'ERS740936', where=lambda rec: rec['Isolate code'] == 'KH02')
    .convert('ox_code', lambda v: 'ERS740940', where=lambda rec: rec['Isolate code'] == 'GN01')
#     .cut([0, 1, 14, 15])
#     .head(4)
)
tbl_samples_to_process.displayall(index_header=True)

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

# CODING_EFFECTS = ['NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED']
NONCODING_EFFECTS = ['INTERGENIC', 'INTRON', 'TRANSCRIPT']
 

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

def assess_validation(isolate_code='7G8', chromosome='Pf3D7_01_v3', mode='SNP', rewrite=False):
    ox_code = tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('ox_code')[0]
    
#     tbl_temp = petl.util.materialise.cache(etl.fromvcf(truth_vcf_fn))
    
    tbl_truth_cache_fn = os.path.join(CACHE_DIR, 'tbl_truth', "tbl_truth_%s.%s.%s" % (isolate_code, chromosome, mode))
    tbl_gatk_cache_fn = os.path.join(CACHE_DIR, 'tbl_gatk', "tbl_gatk_%s.%s.%s" % (isolate_code, chromosome, mode))
    tbl_comparison_cache_fn = os.path.join(CACHE_DIR, 'tbl_comparison', "tbl_comparison_%s.%s.%s" % (isolate_code, chromosome, mode))
    
    if not os.path.exists(tbl_truth_cache_fn) or rewrite:
        truth_vcf_fn = "%s.%s.annotated.%s.vcf.gz" % (
            tbl_samples_to_process.selecteq('Isolate code', isolate_code).values('truth_vcf_filestem')[0],
            chromosome,
            mode
        )
        
        tbl_truth = (etl.fromvcf(truth_vcf_fn)
            .unpackdict('INFO', samplesize=1000000)
    #         .vcfunpackinfo()
            .rename('ALT', 'all_ALTs')
            .addfield('ALT', lambda rec: rec[4][0])
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
    #         .addfield('variant_id', lambda rec: "%s__%d__%s__%s" % (rec[0], rec[1], rec[3], rec[4][0]))
        )
    
#     tbl_truth.valuecounts('Truth_Effect').displayall()
#     tbl_truth.valuecounts('is_coding').displayall()
#     tbl_truth.valuecounts('Truth_Effect', 'is_coding').displayall()

#     tbl_truth.selectgt('POS', 111500).display(index_header=True)
#     tbl_truth.selectgt('POS', 92900).display(index_header=True)
        etl.topickle(tbl_truth, tbl_truth_cache_fn)
    else:
        tbl_truth = etl.frompickle(tbl_truth_cache_fn)
        
    if not os.path.exists(tbl_gatk_cache_fn) or rewrite:
        gatk_vcf_fn = "%s.%s.%s.%s.vcf.gz" % (
            os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
            'Pf3D7july2015',
            mode,
            chromosome
        )
        
    #     print(gatk_vcf_fn)
    #     tbl_temp_gatk = petl.util.materialise.cache(etl.fromvcf(gatk_vcf_fn))
        tbl_gatk = (etl.fromvcf(gatk_vcf_fn)
            .sort(['CHROM', 'POS'])
            .vcfmeltsamples()
            .vcfunpackcall()
            .unpackdict('INFO', samplesize=1000000)
    #         .vcfunpackinfo()
            .selecteq('SAMPLE', ox_code)
            .selectnotnone('GT')
            .selectne('GT', '0/0')
            .rename('ALT', 'all_ALTs')
            .addfield('ALT', lambda rec: rec[4][int(rec[11][2])-1])
    #         .addfield('variant_id', lambda rec: "%s__%d__%s__%s" % (rec[0], rec[1], rec[3], rec[4][int(rec[11][2])-1]))
    #         .selectne('GT', '0/1')
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
            .addfieldusingcontext('GATK_distance_previous', upstream)
            .addfieldusingcontext('GATK_distance_next', downstream)
            .addfield('GATK_distance_nearest', lambda rec: nearest(rec['GATK_distance_previous'], rec['GATK_distance_next']))
    #         .selecteq('GT', '0/1')
    #         .selectgt('VQSLOD', 0)
        )
    #     tbl_gatk.selectgt('POS', 111500).display(5, index_header=True)
    #     tbl_gatk.selectgt('POS', 92900).display(5, index_header=True)
        etl.topickle(tbl_gatk, tbl_gatk_cache_fn)
    else:
        tbl_gatk = etl.frompickle(tbl_gatk_cache_fn)

    if not os.path.exists(tbl_comparison_cache_fn) or rewrite:
    #     tbl_comparison = tbl_truth.outerjoin(tbl_gatk, key='variant_id')
    
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
            .addfield('RegionType', lambda rec: determine_consensus_region(rec, 4, 14))
            .addfield('IsCoding', lambda rec: determine_consensus_region(rec, 5, 15))
            .addfield('Effect', lambda rec: determine_consensus_region(rec, 7, 17))
            .addfield('AAChange', lambda rec: determine_consensus_region(rec, 6, 16))
            .addfield('Gene', lambda rec: determine_consensus_region(rec, 8, 18))
            .addfield('Impact', lambda rec: determine_consensus_region(rec, 9, 19))
            .addfield('Transcript', lambda rec: determine_consensus_region(rec, 10, 20))
#             .intervalleftjoin(tbl_coords_subset, lstart='POS', lstop='POS', rstart='ref_start',
#                               rstop='ref_end', include_stop=True)
#             .addfield('IsAccessible', lambda rec: rec['assembly_chrom'] is not None)
            .addfield('IsAccessible', lambda rec: len(lkp.search(rec['POS'])) > 0)
            .cut(['CHROM', 'POS', 'REF', 'ALT', 'IsInTruth', 'IsInGATK', 'VQSLOD', 'RegionType', 'IsCoding', 'Effect', 'AAChange', 'Gene', 'Impact', 'Transcript',
                  'SAMPLE', 'AD', 'DP', 'GQ', 'GT', 'culprit',
                  'Truth_distance_nearest', 'GATK_distance_nearest', 'IsAccessible'])
            .sort('VQSLOD', reverse=True)
        )
        etl.topickle(tbl_comparison, tbl_comparison_cache_fn)
    else:
        tbl_comparison = etl.frompickle(tbl_comparison_cache_fn)
        
    return(tbl_comparison)
    

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
        

# <codecell>

# def trim_variants(rec):
#     ref = rec[2]
#     alt = rec[3]
#     while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
#         ref = ref[0:-1]
#         alt = alt[0:-1]
#     return((ref, alt))
        

# <codecell>

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
#             .convert(('REF', 'ALT'), lambda rec: str(rec)[0], where=lambda rec: len(rec[2]) > 1 and len(rec[2]) == len(str(rec[3])))
            .addfield('temp_REF', lambda rec: rec['REF'])
            .addfield('temp_ALT', lambda rec: rec['ALT'])
#             .convert(('REF', 'ALT'), lambda ref, rec: trim_variants(rec), pass_row=True)
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

isolate_codes = list(tbl_samples_to_process.values('Isolate code'))
# isolate_codes = ['7G8']
chromosomes = ['Pf3D7_%02d_v3' % chrom for chrom in range(1,15)]

tbl_comparisons = collections.OrderedDict()
tbl_whole_genome_comparisons = collections.OrderedDict()

for isolate_code in isolate_codes:
    tbl_comparisons[isolate_code] = collections.OrderedDict()
    for chromosome in chromosomes:
        print(isolate_code, chromosome)
        tbl_comparisons[isolate_code][chromosome] = assess_validation(isolate_code, chromosome, rewrite=False)
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

# <codecell>

for isolate_code in isolate_codes:
    if isolate_code != 'GN01':
        print(isolate_code)
        (tbl_whole_genome_comparisons[isolate_code]
            .selecteq('IsInGATK', True)
            .selecteq('mode', 'SNP')
            .selecteq('RegionType', 'Core')
            .selecteq('IsCoding', True)
            .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
            .displayall()
        )

# <codecell>

(tbl_whole_genome_comparisons['7G8']
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'SNP')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'SNP')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', False)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

(tbl_whole_genome_comparisons['7G8']
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', False)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3.sort('POS').display(index_header=True)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3.selecteq('mode', 'SNP').valuecounts('GATK_mode', 'Truth_mode', 'mode').displayall()

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3.sort('POS').select(lambda rec: len(rec[2]) > 1 and len(rec[2]) == len(str(rec[3]))).displayall()

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3.sort('POS').selectgt('POS', 100000).display(20)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3.sort('POS').selectgt('POS', 100000).display(20)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3.sort('POS')
    .select(lambda rec: len(rec[2]) > 1 and len(rec[2]) == len(str(rec[3])))
    .convert(('REF', 'ALT'), lambda rec: str(rec)[0], where=lambda rec: len(rec[2]) > 1 and len(rec[2]) == len(str(rec[3])))
).displayall()

# <codecell>

str(tbl_comparison_7g8_Pf3D7_01_v3[1][3])

# <codecell>

(etl.fromvcf("/nfs/team112_internal/production_files/Pf3k/methods/assembled_samples/vcfs/vcf/filtered.annotated.Pf3D7july2015.INDEL.Pf3D7_01_v3.vcf.gz")
    .vcfmeltsamples()
    .vcfunpackcall()
    .unpackdict('INFO', samplesize=1000000)
    .selectgt('POS', 350000)
    .selectnotnone('PGT')
    .selecteq('SAMPLE', 'PG0083-C')
).display(20)

# <codecell>

truth_coords_fn = "%s/Pf%s.filter.coords" % (
            '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/SNP/',
            isolate_code
        )

# <codecell>

tbl_coords = (etl.fromtsv(truth_coords_fn)
    .pushheader(['ref_start', 'ref_end', 'assembly_start', 'assembly_end', 'ref_length', 'assembly_length', 'identity',
                 'ref_chrom_length', 'assembly_chrom_length', 'ref_pct', 'assembly_pct', 'ref_chrom', 'assembly_chrom',
                 'ignore'])
)
tbl_coords_subset = (tbl_coords
    .selecteq('ref_chrom', 'Pf3D7_01_v3')
    .selecteq('assembly_chrom', 'Pf7G8_01')
    .cut(['ref_start', 'ref_end', 'assembly_chrom'])
    .convertnumbers()
)
(tbl_comparison_7g8_Pf3D7_01_v3_SNP
    .intervalleftjoin(tbl_coords_subset, lstart='POS', lstop='POS', rstart='ref_start', rstop='ref_end', include_stop=True)
)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_SNP.selectnone('POS').displayall()

# <codecell>

(3330 / 2038340) * 100 

# <codecell>

(3363  / 40060 ) * 100 

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_SNP = assess_validation(rewrite=True)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_SNP

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3_SNP
    .selectge('VQSLOD', 0)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .displayall()
)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3
    .selectge('VQSLOD', 0)
    .selecteq('mode', 'INDEL')
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .displayall()
)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'INDEL')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', False)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3
    .selecteq('IsInGATK', True)
    .selecteq('mode', 'SNP')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', False)
    .valuecounts('RegionType', 'IsCoding', 'IsInTruth', 'IsInGATK', 'IsAccessible')
    .displayall()
)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3_SNP
#         .selectgt('POS', 92900)
#         .select(lambda rec: rec['VQSLOD'] is None or rec['VQSLOD'] >= 4)
#         .selectnone('consensus_region')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .selecteq('IsInTruth', False)
    .display(50, index_header=True)
)

# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3_SNP
#         .selectgt('POS', 92900)
#         .select(lambda rec: rec['VQSLOD'] is None or rec['VQSLOD'] >= 4)
#         .selectnone('consensus_region')
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .selecteq('IsInGATK', False)
    .displayall(index_header=True)
)


# <codecell>

(tbl_comparison_7g8_Pf3D7_01_v3_SNP
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
    .selecteq('IsInGATK', True)
    .selecteq('IsInTruth', True)
    .displayall(index_header=True)
)


# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_indel = assess_validation(mode='INDEL', rewrite=True)

# <codecell>

tbl_gatk_snp = etl.frompickle(os.path.join(CACHE_DIR, 'tbl_gatk', "tbl_gatk_%s.%s.%s" % ('7G8', 'Pf3D7_01_v3', 'SNP')))
tbl_gatk_indel = etl.frompickle(os.path.join(CACHE_DIR, 'tbl_gatk', "tbl_gatk_%s.%s.%s" % ('7G8', 'Pf3D7_01_v3', 'INDEL')))

# <codecell>

tbl_gatk_snp.selectgt('POS', 164900)

# <codecell>

tbl_gatk_snp.selectgt('POS', 199000)

# <codecell>

tbl_gatk_indel.selectgt('POS', 197000)

# <codecell>

def number_of_TP(prv, cur, nxt):
    if prv is None:
        if cur['IsInGATK'] and cur['IsInTruth']:
            return 1
        else:
            return 0
    else:
        if cur['IsInGATK'] and cur['IsInTruth']:
            return prv['number_TP'] + 1
        else:
            return prv['number_TP']

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_indel.selecteq('IsCoding', True).sort('POS').displayall()

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_indel.sort('POS').selectgt('POS', 348000).display(25)

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_indel.sort('POS').selectgt('POS', 348000).display(25)

# <codecell>

ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))
ref_dict['Pf3D7_01_v3'].seq[349910:349945]

# <codecell>

ref_dict['Pf3D7_01_v3'].seq[350740:350749]

# <codecell>

tbl_comparison_7g8_Pf3D7_01_v3_indel.sort('POS').selectgt('POS', 613000).display(10)

# <codecell>

def plotROC(isolate_code='7G8', RegionType='Core', IsCoding=True, mode='SNP', min_GATK_distance_nearest=0):
    array_to_plot = (tbl_whole_genome_comparisons['7G8']
        .selecteq('IsInGATK', True)
        .selecteq('mode', mode)
        .selecteq('RegionType', RegionType)
        .selecteq('IsCoding', IsCoding)
        .selecteq('IsAccessible', True)
        .selectge('GATK_distance_nearest', min_GATK_distance_nearest)
        .cut(['IsInGATK', 'IsInTruth', 'VQSLOD'])
    ).toarray()
    TPs = np.cumsum(array_to_plot['IsInTruth'])
    FPs = np.cumsum(np.logical_not(array_to_plot['IsInTruth']))
    plot(FPs, TPs)


# <codecell>

tbl_whole_genome_comparisons['7G8']

# <codecell>

plotROC()

# <codecell>

plotROC(IsCoding=False)

# <codecell>

plotROC(mode='INDEL', IsCoding=True)

# <codecell>

plotROC(mode='INDEL', IsCoding=True, min_GATK_distance_nearest=10)

# <codecell>

plotROC(mode='INDEL', IsCoding=True, min_GATK_distance_nearest=100)

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

plotROC(mode='INDEL', IsCoding=False)

# <codecell>

# tbl_temp = (tbl_comparison_7g8_Pf3D7_01_v3_SNP
tbl_temp = (tbl_comparison_7g8_Pf3D7_01_v3_indel
#     .selectge('VQSLOD', 0)
    .selecteq('RegionType', 'Core')
    .selecteq('IsCoding', True)
#     .select(lambda rec: rec['Truth_distance_nearest'] is not None and rec['Truth_distance_nearest'] > 10)
    .select(lambda rec: rec['Truth_distance_nearest'] is None or rec['Truth_distance_nearest'] > 10)
    .selecteq('IsAccessible', True)
    .selecteq('IsInGATK', True)
    .cut(['IsInGATK', 'IsInTruth', 'VQSLOD'])
#     .addfieldusingcontext('number_TP', number_of_TP)
)

# <codecell>

array_temp = etl.toarray(tbl_temp)

# <codecell>

TPs = np.cumsum(array_temp['IsInTruth'])
FPs = np.cumsum(np.logical_not(array_temp['IsInTruth']))

# <codecell>

TPs

# <codecell>

FPs

# <codecell>

plot(FPs, TPs)

# <codecell>

array_temp[84]

# <codecell>

array_temp

# <codecell>

array_temp[104]

# <codecell>

array_temp[85]

# <codecell>

array_temp[138:142]

# <codecell>

tbl_temp.valuecounts('IsInTruth', 'IsInGATK')

# <codecell>

tbl_temp.valuecounts('IsInTruth', 'IsInGATK')

# <codecell>


