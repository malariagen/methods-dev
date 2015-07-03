# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'bams', 'bwa_mem')}
!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'gvcf', 'samples')}
!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'gvcf', 'lists')}
!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'gvcf', 'combined_gvcf')}
!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'recal')}

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
    .selectne('Sample name', 'ITA4?')
#     .addfield('bam_dir', lambda rec: os.path.dirname(rec['bam_fn']) if rec['bam_fn'] is not None else None)
    .addfield('remapped_sam_fn', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'bams',
        'bwa_mem',
        os.path.basename(rec['bam_fn']).replace('.bam', '.bwa_mem.sam')
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
)
tbl_samples_to_process.displayall(index_header=True)

# <codecell>

tbl_samples_to_process.values('bam_fn')[0]

# <headingcell level=1>

# Functions

# <codecell>

def remap_with_bwa_mem(input_fn = tbl_samples_to_process.values('bam_fn')[0],
                       bwa_mem_fn = tbl_samples_to_process.values('remapped_sam_fn')[0],
                       read_group_info = tbl_samples_to_process.values('RG')[0],
                       reference_fasta = REF_GENOME,
                       rewrite=False):
    if not os.path.isfile(bwa_mem_fn) or rewrite:
        !{samtools_exe} bamshuf -uOn 128 {input_fn} tmp | \
        {samtools_exe} bam2fq - | \
        {bwa_exe} mem -M -R '{read_group_info}' -p {reference_fasta} - > {bwa_mem_fn} 2> /dev/null

# <codecell>

remap_with_bwa_mem()

# <codecell>

def sort_mark_dup(bwa_mem_fn = tbl_samples_to_process.values('remapped_sam_fn')[0],
                  sorted_fn = tbl_samples_to_process.values('sorted_bam_fn')[0],
                  dedup_fn = tbl_samples_to_process.values('dedupped_bam_fn')[0],
                  dedup_metrics_fn = tbl_samples_to_process.values('dedup_metrics_fn')[0],
                  rewrite=False):
    
    if not os.path.isfile(sorted_fn) or rewrite:
        !{picard_exe} SortSam \
        INPUT={bwa_mem_fn} \
        OUTPUT={sorted_fn} \
        SORT_ORDER=coordinate 2> /dev/null
    
    if not os.path.isfile(dedup_fn) or rewrite:
        !{picard_exe} MarkDuplicates \
        INPUT={sorted_fn} \
        OUTPUT={dedup_fn} \
        METRICS_FILE={dedup_metrics_fn} 2> /dev/null

    if not os.path.isfile(dedup_fn.replace('.bam', '.bai')) or rewrite:
        !{picard_exe} BuildBamIndex \
        INPUT={dedup_fn} 2> /dev/null

    

# <codecell>

sort_mark_dup()

# <codecell>

def indel_realignment(dedup_fn = tbl_samples_to_process.values('dedupped_bam_fn')[0],
                      targets_fn = tbl_samples_to_process.values('targets_fn')[0],
                      realigned_fn = tbl_samples_to_process.values('realigned_bam_fn')[0],
                      reference_fasta = REF_GENOME,
                      rewrite=False):
    
    if not os.path.isfile(targets_fn) or rewrite:
        !{gatk_exe} \
        -T RealignerTargetCreator \
        -R {reference_fasta} \
        -I {dedup_fn} \
        -o {targets_fn}
    
    if not os.path.isfile(realigned_fn) or rewrite:
        !{gatk_exe} \
        -T IndelRealigner \
        -R {reference_fasta} \
        -I {dedup_fn} \
        -targetIntervals {targets_fn} \
        -o {realigned_fn}

# <codecell>

indel_realignment()

# <codecell>

def bqsr(realigned_fn = tbl_samples_to_process.values('realigned_bam_fn')[0],
         recal_table_fn = tbl_samples_to_process.values('recal_table_fn')[0],
         post_recal_table_fn = tbl_samples_to_process.values('post_recal_table_fn')[0],
         recal_plots_fn = tbl_samples_to_process.values('recal_plots_fn')[0],
         recal_bam_fn = tbl_samples_to_process.values('recal_bam_fn')[0],
         crosses_dir = crosses_dir,
         reference_fasta = REF_GENOME,
         rewrite=False):
    
    if not os.path.isfile(recal_table_fn) or rewrite:
        !{gatk_exe} \
        -T BaseRecalibrator \
        -R {reference_fasta} \
        -I {realigned_fn} \
        -knownSites {crosses_dir}/7g8_gb4.combined.final.vcf.gz \
        -knownSites {crosses_dir}/hb3_dd2.combined.final.vcf.gz \
        -knownSites {crosses_dir}/3d7_hb3.combined.final.vcf.gz \
        -o {recal_table_fn}

    if not os.path.isfile(post_recal_table_fn) or rewrite:
        !{gatk_exe} \
        -T BaseRecalibrator \
        -R {reference_fasta} \
        -I {realigned_fn} \
        -knownSites {crosses_dir}/7g8_gb4.combined.final.vcf.gz \
        -knownSites {crosses_dir}/hb3_dd2.combined.final.vcf.gz \
        -knownSites {crosses_dir}/3d7_hb3.combined.final.vcf.gz \
        -BQSR {recal_table_fn} \
        -o {post_recal_table_fn}

    if not os.path.isfile(recal_plots_fn) or rewrite:
        !{gatk_exe} \
        -T AnalyzeCovariates \
        -R {reference_fasta} \
        -before {recal_table_fn} \
        -after {post_recal_table_fn} \
        -plots {recal_plots_fn}

    if not os.path.isfile(recal_bam_fn) or rewrite:
        !{gatk_exe} \
        -T PrintReads \
        -R {reference_fasta} \
        -I {realigned_fn} \
        -BQSR {recal_table_fn} \
        -o {recal_bam_fn}

# <codecell>

bqsr()

# <codecell>

for rec in tbl_samples_to_process.data():
    print(rec)
    remap_with_bwa_mem(rec[9], rec[10], rec[20])
    sort_mark_dup(rec[10], rec[11], rec[12], rec[13])
    indel_realignment(rec[12], rec[14], rec[15])
    bqsr(rec[15], rec[16], rec[17], rec[18], rec[19])

# <codecell>

def haplotype_caller(recal_bam_fn = tbl_samples_to_process.values('recal_bam_fn')[0],
         gvcf_filestem = tbl_samples_to_process.values('gvcf_filestem')[0],
         reference_fasta = REF_GENOME,
         chromosomes = None,
         rewrite=False):
    
    if chromosomes is None:
        in_seq_handle = open(reference_fasta)
        chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")).keys()
        print(chromosomes)

    for chromosome in chromosomes:
        gvcf_fn = "%s.%s.g.vcf" % (gvcf_filestem, chromosome)
        if not os.path.isfile(gvcf_fn) or rewrite:
            !{gatk_exe} \
            -T HaplotypeCaller \
            -R {reference_fasta} \
            -I {recal_bam_fn} \
            -L {chromosome} \
            --emitRefConfidence GVCF \
            --variant_index_type LINEAR \
            --variant_index_parameter 128000 \
            -o {gvcf_fn} 2> /dev/null
        else:
            print("Skipping %s" % gvcf_fn)

# <codecell>

in_seq_handle = open(REF_GENOME)
chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")).keys()
for chromosome in chromosomes:
    print(chromosome)
    for rec in tbl_samples_to_process.data():
        print(rec[0])
        haplotype_caller(rec[19], rec[21], chromosomes=[chromosome])

# <codecell>

def combine_and_genotype_gvcfs(gvcf_filestems=tbl_samples_to_process.values('gvcf_filestem').array(),
                  list_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'gvcf', 'lists', 'gvcfs'),
                  combined_gvcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'gvcf', 'combined_gvcf', 'combined'),
                  genotyped_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'raw.snps.indels'),
                  chromosomes = None, max_number_per_list=200,
                  reference_fasta = REF_GENOME, rewrite=False):
    
    if chromosomes is None:
        in_seq_handle = open(REF_GENOME)
        chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")).keys()

    number_of_files = int(ceil(len(gvcf_filestems)/max_number_per_list))
    np.random.shuffle(gvcf_filestems)
    for chromosome in chromosomes:
        combined_gvcf_list_filename = "%s.%s.all.list" % (list_filestem, chromosome)
        fo = open(combined_gvcf_list_filename, 'w')
        for i in range(number_of_files):
            filestems_this_chunk = gvcf_filestems[(i*max_number_per_list):((i+1)*max_number_per_list)]
            gvcf_list_filename = "%s.%s.%d.list" % (list_filestem, chromosome, i)
            combined_gvcf_fn = "%s.%s.%d.g.vcf" % (combined_gvcf_filestem, chromosome, i)
            print(combined_gvcf_fn, file=fo)
            
            files_this_chrom = np.array(["%s.%s.g.vcf" % (gvcf_filestem, chromosome) for gvcf_filestem in filestems_this_chunk])
            np.savetxt(gvcf_list_filename, files_this_chrom, '%s')

            if not os.path.isfile(combined_gvcf_fn) or rewrite:
                !{gatk_exe} \
                -T CombineGVCFs \
                -R {reference_fasta} \
                --variant {gvcf_list_filename} \
                -o {combined_gvcf_fn}
            else:
                print("Skipping %s" % combined_gvcf_fn)
        fo.close()

        genotyped_vcf_fn = "%s.%s.vcf" % (genotyped_vcf_filestem, chromosome)
        
        if not os.path.isfile(genotyped_vcf_fn) or rewrite:
            !{gatk_exe} \
            -T GenotypeGVCFs \
            -R {reference_fasta} \
            --variant {combined_gvcf_list_filename} \
            -o {genotyped_vcf_fn}


# <codecell>

combine_and_genotype_gvcfs(chromosomes=['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_06_v3', 'Pf3D7_08_v3', 
                                        'Pf3D7_09_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3'],
                           max_number_per_list=2, rewrite=True)

# <codecell>

def train_vqsr(genotyped_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'raw.snps.indels'),
               recal_dir=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'recal'),
               annotations_line='-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP',
#                annotations_line='-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP',
               mode='SNP', maxGaussians=8, chromosomes=None,
               reference_fasta = REF_GENOME, rewrite=False):
    
#     raw_input_vcf_fn = "%s.%s.vcf" % (genotyped_vcf_filestem, chromosome)
    
    if chromosomes is None:
        in_seq_handle = open(reference_fasta)
        chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")).keys()
    
    raw_input_vcf_fns = ["%s.%s.vcf" % (genotyped_vcf_filestem, chromosome) for chromosome in chromosomes]
    input_line = "-input " + " -input ".join(raw_input_vcf_fns)
#     print(input_line)

    recal_fn = "%s/recalibrate_%s.recal" % (recal_dir, mode)
    tranches_fn = "%s/recalibrate_%s.tranches" % (recal_dir, mode)

    if not os.path.isfile(recal_fn) or rewrite:
        !{gatk_exe} \
        -T VariantRecalibrator \
        -R {reference_fasta} \
        {input_line} \
        -resource\:7g8_gb4,known=false,training=true,truth=true,prior=15.0 {crosses_dir}/7g8_gb4.combined.final.vcf.gz \
        -resource\:hb3_dd2,known=false,training=true,truth=true,prior=15.0 {crosses_dir}/hb3_dd2.combined.final.vcf \
        -resource\:3d7_hb3,known=false,training=true,truth=true,prior=15.0 {crosses_dir}/3d7_hb3.combined.final.vcf.gz \
        {annotations_line} \
        --maxGaussians {maxGaussians} \
        -mode {mode} \
        -recalFile {recal_fn} \
        -tranchesFile {tranches_fn}


# <codecell>

# Train on SNPs
train_vqsr(chromosomes=['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_06_v3', 'Pf3D7_08_v3', 
                        'Pf3D7_09_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3'],
           rewrite=True)

# <codecell>

# Train on indels
train_vqsr(chromosomes=['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_06_v3', 'Pf3D7_08_v3', 
                        'Pf3D7_09_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3'],
           mode='INDEL',
           annotations_line='-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum',
           maxGaussians=4,                  
           rewrite=True)

# <codecell>

def apply_vqsr(chromosome='Pf3D7_01_v3',
               genotyped_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'raw.snps.indels'),
               filtered_vcf_filestem=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'filtered'),
               recal_dir=os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'vcfs', 'vcf', 'recal'),
               ts_filter_level=99.5,
               mode='SNP',
               reference_fasta = REF_GENOME,
               rewrite=False):
    
    genotyped_vcf_fn = "%s.%s.vcf" % (genotyped_vcf_filestem, chromosome)
    filtered_vcf_fn = "%s.%s.%s.vcf" % (filtered_vcf_filestem, mode, chromosome)

    recal_fn = "%s/recalibrate_%s.recal" % (recal_dir, mode)
    tranches_fn = "%s/recalibrate_%s.tranches" % (recal_dir, mode)

    if not os.path.isfile(filtered_vcf_fn+'.gz') or rewrite:
        !{gatk_exe} \
        -T ApplyRecalibration \
        -R {reference_fasta} \
        -input {genotyped_vcf_fn} \
        -tranchesFile {tranches_fn} \
        -recalFile {recal_fn} \
        --ts_filter_level {ts_filter_level} \
        -mode {mode} \
        -o {filtered_vcf_fn}
        
        !bgzip -f {filtered_vcf_fn}
        !tabix -p vcf -f {filtered_vcf_fn}.gz

# <codecell>

apply_vqsr()

# <codecell>

for chromosome in ['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_06_v3', 'Pf3D7_08_v3', 
                    'Pf3D7_09_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3']:
    for mode in ['SNP', 'INDEL']:
        apply_vqsr(chromosome=chromosome, mode=mode)

# <codecell>


