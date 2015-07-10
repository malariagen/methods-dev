# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <headingcell level=1>

# Check data

# <codecell>

def get_RG(bam_fn = '/nfs/team112_internal/production_files/Pf/4_0/PFprog1/PG0008_CW/PG0008_CW.bam'):
    temp = !{samtools_exe} view -H {bam_fn} | grep '@RG'
    return(temp[0])
get_RG()

# <codecell>

tbl_samples_to_process = (tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
#     .selectne('Sample name', 'Pf3D7_II')
#     .selectne('Sample name', 'ITA4?')
#     .addfield('bam_dir', lambda rec: os.path.dirname(rec['bam_fn']) if rec['bam_fn'] is not None else None)
    .addfield('remapped_sam_fn', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'bams',
        'bwa_mem',
        os.path.basename(rec['bam_fn']).replace('.bam', '.bwa_mem.sam').replace('.cram', '.bwa_mem.sam')
    ))
    .addfield('sorted_bam_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.sorted.bam'))
    .addfield('dedupped_bam_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.dedup.bam'))
    .addfield('dedup_metrics_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.dedup.metrics'))
    .addfield('targets_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.dedup.targets.list'))
    .addfield('realigned_bam_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.realigned.bam'))
    .addfield('recal_table_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.recal.table'))
    .addfield('post_recal_table_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.post_recal.table'))
    .addfield('recal_plots_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.recal.pdf'))
    .addfield('recal_bam_fn', lambda rec: rec['remapped_sam_fn'].replace('.sam', '.recal.bam'))
    .addfield('RG', lambda rec: get_RG(rec['bam_fn']).replace('\t', '\\t'))
    .addfield('gvcf_filestem', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'vcfs',
        'gvcf',
        'samples',
        os.path.basename(rec['bam_fn']).replace('.bam', '.raw.snps.indels')
    ))
#     .head(4)
)
tbl_samples_to_process.displayall(index_header=True)

# <headingcell level=1>

# Annotate variants

# <codecell>

def apply_snpeff(chromosome='Pf3D7_01_v3',
               snpeff_db = 'Pf3D7july2015',
               filtered_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered'),
               snpeff_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.snpeff'),
               snpeff_annotated_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.snpeff_annotated'),
               annotated_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
               mode='SNP',
               reference_fasta = REF_GENOME,
               rewrite=False):
    
    filtered_vcf_fn = "%s.%s.%s.vcf.gz" % (filtered_vcf_filestem, mode, chromosome)
    snpeff_vcf_fn = "%s.%s.%s.%s.vcf" % (snpeff_vcf_filestem, snpeff_db, mode, chromosome)
    snpeff_annotated_vcf_fn = "%s.%s.%s.%s.vcf" % (snpeff_annotated_vcf_filestem, snpeff_db, mode, chromosome)
    annotated_vcf_fn = "%s.%s.%s.%s.vcf" % (annotated_vcf_filestem, snpeff_db, mode, chromosome)

    if not os.path.isfile(snpeff_vcf_fn) or rewrite:
        !{snpeff_exe} \
        -v -o gatk {snpeff_db} \
        {filtered_vcf_fn} \
        -no-downstream \
        -no-upstream \
        > {snpeff_vcf_fn} \
        2> /dev/null
    else:
        print("\tSkipping %s" % snpeff_vcf_fn)
        
    if not os.path.isfile(snpeff_annotated_vcf_fn) or rewrite:
        !{gatk_exe} \
        -T VariantAnnotator \
        -R {reference_fasta} \
        -A SnpEff \
        --variant {filtered_vcf_fn} \
        --snpEffFile {snpeff_vcf_fn} \
        -o {snpeff_annotated_vcf_fn} \
        2> /dev/null
    
    if not os.path.isfile(annotated_vcf_fn+'.gz') or rewrite:
        !cat {snpeff_annotated_vcf_fn} \
        | vcf-annotate -a {regions_fn} \
           -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
           -c CHROM,FROM,TO,INFO/RegionType \
        > {annotated_vcf_fn}
        
        !bgzip -f {annotated_vcf_fn}
        !tabix -p vcf -f {annotated_vcf_fn}.gz
    else:
        print("\tSkipping %s" % annotated_vcf_fn+'.gz')

# <codecell>

for snpeff_db in ['Pf3D7v91', 'Pf3D7v24', 'Pf3D7july2015']:
    apply_snpeff(snpeff_db=snpeff_db, rewrite=True)

# <codecell>

for snpeff_db in ['Pf3D7v91', 'Pf3D7v24', 'Pf3D7july2015']:
    apply_snpeff(snpeff_db=snpeff_db, mode='INDEL', rewrite=True)

# <codecell>

(etl
    .fromvcf("%s.%s.%s.%s.vcf.gz" % (
        os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
        'Pf3D7july2015', 'INDEL', 'Pf3D7_01_v3'), samples=None)
    .vcfunpackinfo()
)

# <codecell>

tbl_snpeff = collections.OrderedDict()

for snpeff_db in ['Pf3D7v91', 'Pf3D7v24', 'Pf3D7july2015']:
    annotated_vcf_fn = "%s.%s.%s.%s.vcf.gz" % (
        os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
        snpeff_db, 'SNP', 'Pf3D7_01_v3')
    tbl_snpeff[snpeff_db] = (
        etl
        .fromvcf(annotated_vcf_fn, samples=None)
        .vcfunpackinfo()
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_IMPACT', 'RegionType', 'VQSLOD'])
        .addfield('FILTER1', lambda rec: rec[4][0] if len(rec[4]) == 1 else 'PASS')
        .cutout('FILTER')
        .rename('SNPEFF_EFFECT', 'SNPEFF_EFFECT_%s' % snpeff_db)
        .rename('SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_FUNCTIONAL_CLASS_%s' % snpeff_db)
        .rename('SNPEFF_IMPACT', 'SNPEFF_IMPACT_%s' % snpeff_db)
    )

# <codecell>

tbl_snpeff_indel = collections.OrderedDict()

for snpeff_db in ['Pf3D7v91', 'Pf3D7v24', 'Pf3D7july2015']:
    annotated_vcf_fn = "%s.%s.%s.%s.vcf.gz" % (
        os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
        snpeff_db, 'INDEL', 'Pf3D7_01_v3')
    tbl_snpeff_indel[snpeff_db] = (
        etl
        .fromvcf(annotated_vcf_fn, samples=None)
        .vcfunpackinfo()
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_IMPACT', 'RegionType', 'VQSLOD'])
        .addfield('FILTER1', lambda rec: rec[4][0] if len(rec[4]) == 1 else 'PASS')
        .cutout('FILTER')
        .rename('SNPEFF_EFFECT', 'SNPEFF_EFFECT_%s' % snpeff_db)
        .rename('SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_FUNCTIONAL_CLASS_%s' % snpeff_db)
        .rename('SNPEFF_IMPACT', 'SNPEFF_IMPACT_%s' % snpeff_db)
    )

# <codecell>

tbl_snpeff_indel['Pf3D7july2015']

# <codecell>

tbl_snpeff['Pf3D7v91'].valuecounts('SNPEFF_EFFECT_Pf3D7v91').displayall()

# <codecell>

tbl_snpeff['Pf3D7v24'].valuecounts('SNPEFF_EFFECT_Pf3D7v24').displayall()

# <codecell>

tbl_snpeff['Pf3D7july2015'].valuecounts('SNPEFF_EFFECT_Pf3D7july2015').displayall()

# <codecell>

# Note these appear to be pseudogenes
tbl_snpeff['Pf3D7july2015'].selecteq('SNPEFF_EFFECT_Pf3D7july2015', 'TRANSCRIPT').displayall()

# <codecell>

# Note these appear to either rRNAs or else beyond last gene
tbl_snpeff['Pf3D7july2015'].selectnone('SNPEFF_EFFECT_Pf3D7july2015').displayall()

# <codecell>

tbl_snpeff_indel['Pf3D7v91'].valuecounts('SNPEFF_EFFECT_Pf3D7v91').displayall()

# <codecell>

(tbl_snpeff_indel['Pf3D7july2015']
    .addfield('Good', lambda rec: rec['VQSLOD'] > 5.5)
    .valuecounts('SNPEFF_EFFECT_Pf3D7july2015', 'Good', 'RegionType')
    .selectin('SNPEFF_EFFECT_Pf3D7july2015', ['FRAME_SHIFT', 'CODON_INSERTION', 'CODON_DELETION'])
    .displayall()
)

# <codecell>

# Note these appear to be multi-allelic (SNP and INDEL)
tbl_snpeff_indel['Pf3D7july2015'].selecteq('SNPEFF_EFFECT_Pf3D7july2015', 'NON_SYNONYMOUS_CODING').displayall()

# <codecell>

# Note these appear to either rRNAs or else beyond last gene
tbl_snpeff_indel['Pf3D7july2015'].selectnone('SNPEFF_EFFECT_Pf3D7july2015').displayall()

# <codecell>

tbl_snpeff['Pf3D7v24'].select(lambda rec: rec['POS'] >= 520989 and rec['POS'] <= 520989).displayall()

# <codecell>

tbl_snpeff['Pf3D7july2015']

# <codecell>

print(len(tbl_snpeff['Pf3D7v91']))
print(len(tbl_snpeff['Pf3D7july2015']))

# <codecell>

tbl_snpeffs = (tbl_snpeff['Pf3D7v91']
               .join(tbl_snpeff['Pf3D7july2015'], key=['CHROM', 'POS', 'REF', 'ALT']))

# <codecell>

print(len(tbl_snpeffs))
tbl_snpeffs

# <codecell>

tbl_snpeffs.select(lambda rec: rec['SNPEFF_EFFECT_Pf3D7v91'] != rec['SNPEFF_EFFECT_Pf3D7july2015']).displayall()

# <codecell>

(tbl_snpeffs
 .select(lambda rec: rec['SNPEFF_EFFECT_Pf3D7v91'] != rec['SNPEFF_EFFECT_Pf3D7july2015'])
 .valuecounts('SNPEFF_EFFECT_Pf3D7v91', 'SNPEFF_EFFECT_Pf3D7july2015')
).displayall()

# <codecell>

tbl_snpeffs.select(lambda rec: rec['POS'] >= 104704 and rec['POS'] <= 105209)

# <codecell>

tbl_snpeffs.select(lambda rec: rec['SNPEFF_EFFECT_Pf3D7v91'] == 'INTRON' and rec['SNPEFF_EFFECT_Pf3D7july2015'] == 'DOWNSTREAM')

# <codecell>

tbl_valuecounts = collections.OrderedDict()
for snpeff_db in ['Pf3D7v91', 'Pf3D7v24', 'Pf3D7july2015']:
    annotated_vcf_fn = "%s.%s.%s.%s.vcf.gz" % (
        os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered.annotated'),
        snpeff_db, 'SNP', 'Pf3D7_01_v3')
    tbl_temp = (
        etl
        .fromvcf(annotated_vcf_fn, samples=None)
        .vcfunpackinfo()
    )
    tbl_valuecounts[snpeff_db] = tbl_temp.valuecounts('SNPEFF_EFFECT')
    tbl_valuecounts[snpeff_db] = tbl_valuecounts[snpeff_db].replace('SNPEFF_EFFECT', None, 'None')

# <codecell>


# <codecell>

tbl_valuecounts['Pf3D7v91'].displayall()

# <codecell>

(tbl_valuecounts['Pf3D7v91']
    .outerjoin(tbl_valuecounts['Pf3D7v24'], key='SNPEFF_EFFECT')
    .outerjoin(tbl_valuecounts['Pf3D7july2015'], key='SNPEFF_EFFECT')
).displayall()
# .join(tbl_valuecounts['Pf3D7july2015'])

# <codecell>

tbl_temp = (
    etl
    .fromvcf(annotated_vcf_fn, samples=None)
    .vcfunpackinfo()
)

# <codecell>

tbl_temp.valuecounts('SNPEFF_EFFECT').displayall()

# <codecell>


