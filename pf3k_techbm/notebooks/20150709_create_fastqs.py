# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'fastq')}

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
    .addfield('fastq_fn', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'fastq',
        os.path.basename(rec['bam_fn']).replace('.bam', '.fastq').replace('.cram', '.fastq')
    ))
#     .head(4)
)
tbl_samples_to_process.displayall(index_header=True)

# <headingcell level=1>

# Functions

# <codecell>

%%file first_last_100bp.py
import sys
from Bio import SeqIO

iterator = SeqIO.parse(sys.stdin, "fastq")

for firstRead in iterator:
    firstRead_first100 = firstRead[0:100]
    firstRead_first100.id = firstRead_first100.id.replace('#', '_first100#')
    firstRead_first100.description = ''
    firstRead_last100 = firstRead[150:250]
    firstRead_last100.id = firstRead_last100.id.replace('#', '_last100#')
    firstRead_last100.description = ''
    secondRead = next(iterator)
    secondRead_first100 = secondRead[0:100]
    secondRead_first100.id = secondRead_first100.id.replace('#', '_first100#')
    secondRead_first100.description = ''
    secondRead_last100 = secondRead[150:250]
    secondRead_last100.id = secondRead_last100.id.replace('#', '_last100#')
    secondRead_last100.description = ''
    SeqIO.write(firstRead_first100, sys.stdout, "fastq")
    SeqIO.write(secondRead_first100, sys.stdout, "fastq")
    SeqIO.write(firstRead_last100, sys.stdout, "fastq")
    SeqIO.write(secondRead_last100, sys.stdout, "fastq")


# <codecell>

def create_fastq(input_fn = tbl_samples_to_process.values('bam_fn')[0],
                   fastq_fn = tbl_samples_to_process.values('fastq_fn')[0],
#                    read_group_info = tbl_samples_to_process.values('RG')[0],
#                    reference_fasta = REF_GENOME,
                   rewrite=False):
    
    if not os.path.isfile(fastq_fn) or rewrite:
        if input_fn.endswith('.bam'):
            !{samtools_exe} bamshuf -uOn 128 {input_fn} tmp | \
            {samtools_exe} bam2fq - \
            > {fastq_fn} 2> /dev/null
        elif input_fn.endswith('.cram'):
            !{samtools_exe} view -b {input_fn} | \
            {samtools_exe} bamshuf -uOn 128 - tmp | \
            {samtools_exe} bam2fq - \
            > {fastq_fn} 2> /dev/null

# <codecell>

create_fastq()

# <codecell>

create_fastq(tbl_samples_to_process.values('bam_fn')[8], tbl_samples_to_process.values('fastq_fn')[8])

# <codecell>


