# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Plan
# - Save data in appropriate place
# - Parse sample spreadsheets (GATC and source bio)
# - Table of capillary samples with columns GATC, source bio, SQNM, illumina, Pf3k
# - Remove bams to make space
# - Download bams for all cappilary samples
# - Map capillary data with bwa mem and call with GATK haploid (SNPs and indels)
# - Check if any reads don't map Kelch13
# - Create vcf of all capillary samples combined
# - Pull sequenom data into table
# - Fix petlx bug with self.end and post pull request
# - Create vcf of appropriate samples/region from PGV4.0
# - Create table of all SNPs (chrom_pos_ref_alt), with columns for genotype in each assay
# - For each illumina discovered SNP, for each company
#     - number ILMN samples
#     - number capillary for each primer
#     - number rediscovered for each SNP
#     - % for each primer
# - For each capillary discovered variant 
#     - number capillary samples
#     - number illumina for each primer
#     - number rediscovered for each SNP
#     - % for each primer
# - Concordance (?)
# - Clean up notebook

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

GATC_EXCEL_FN = '/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/GATC-Biotech seq barcode sample list.xlsx'
GATC_SEQ_DIRs = [
    '/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150527 raw data',
    '/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150605 raw data'
]

SB_EXCEL_FN = '/data/plasmodium/pfalciparum/recon/data/raw/capillary/Source_Bioscience/source_bioscience_sample_list_20150513.xls'
SB_SEQ_DIR = '/data/plasmodium/pfalciparum/recon/data/raw/capillary/Source_Bioscience/run_2015_05_13'

CAPILLARY_SAMPLES_FN = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/initial_capillary_samples.xlsx'
!mkdir -p /data/plasmodium/pfalciparum/recon/data/processed/capillary

BAM_FILES_DIR = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/bams'
ALL_BAM_FILES_FOFN = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/bams/all_bams.list'
GATC_BAM_FILES_FOFN = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/bams/gatc_bams.list'
SB_BAM_FILES_FOFN = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/bams/sb_bams.list'
!mkdir -p {BAM_FILES_DIR}

VARIANT_FILES_DIR = '/data/plasmodium/pfalciparum/recon/data/processed/capillary/variants'
!mkdir -p {VARIANT_FILES_DIR}
ALL_CAPILLARY_VCF_FN = "%s/all_capillary.vcf" % VARIANT_FILES_DIR
GATC_CAPILLARY_VCF_FN = "%s/gatc_capillary.vcf" % VARIANT_FILES_DIR
SB_CAPILLARY_VCF_FN = "%s/sb_capillary.vcf" % VARIANT_FILES_DIR
MERGED_CALLS_EXCEL_FN = "%s/merged_calls.xlsx" % VARIANT_FILES_DIR

SQNM_VCF_FN = DATA_DIR + '/Sqnm_data_DK1066_W1378_20150603.vcf.gz'
SAMPLES_FN = DATA_DIR + '/%s_capillary_samples.txt'

# <headingcell level=1>

# Parse sample details

# <codecell>

def determine_gatc_seq_fn(rec):
    if rec['primer_name'].startswith('DGGE'):
        return("%s/%s.fas" % (GATC_SEQ_DIRs[1], rec['seqID']))
    else:
        return("%s/%s.fas" % (GATC_SEQ_DIRs[0], rec['seqID']))
        

# <codecell>

def which_primers(rec):
    if rec['primer_name'].startswith('Forward'):
        return('Forward')
    if rec['primer_name'].startswith('Reverse'):
        return('Reverse')
    if rec['primer_name'].endswith('seq1'):
        return('seq1')
    if rec['primer_name'].endswith('seq2'):
        return('seq2')
    else:
        return('unknown')
    

# <codecell>

from Bio import SeqIO

def convert_to_fastq(seq_fn='/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150527 raw data/00AJ08.fas',
                     fastq_fn=None, qual_good=60, qual_bad=0):
    if fastq_fn is None:
        fastq_fn = seq_fn.replace('.fas', '.fastq').replace('.seq', '.fastq')
        
    # make fastq
    with open(seq_fn, "r") as fasta, open(fastq_fn, "w") as fastq:
        for record in SeqIO.parse(fasta, "fasta"):
            record.seq = record.seq.tomutable()
            record.letter_annotations["phred_quality"] = [qual_good] * len(record)
            for i in range(len(record)):
                if i == 0:
                    if record.seq[i].islower() or record.seq[i+1].islower() or record.seq[i].upper()=='N':
                        record.letter_annotations["phred_quality"][i] = qual_bad
                elif i == (len(record) - 1):
                    if record.seq[i].islower() or record.seq[i-1].islower() or record.seq[i].upper()=='N':
                        record.letter_annotations["phred_quality"][i] = qual_bad
                else:
                    if record.seq[i].islower() or record.seq[i-1].islower() or record.seq[i+1].islower() or record.seq[i].upper()=='N':
                        record.letter_annotations["phred_quality"][i] = qual_bad
                # Change N's to a, as bwa seems to ignore these and treat as deletion. These should have qual=0 and hence not get called as SNPs
                if record.seq[i].upper()=='N':
                    record.seq[i] = 'a'
                        
            SeqIO.write(record, fastq, "fastq")
    
    return(fastq_fn)

convert_to_fastq()

# <codecell>

def sequence_length(seq_fn='/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150527 raw data/00AJ08.fas'):
    return(len(SeqIO.read(open(seq_fn, "r"), "fasta").seq))
sequence_length()

# <codecell>

tbl_gatc_metadata = (
    etl.fromxlsx(GATC_EXCEL_FN)
    .addfield('ox_code', lambda rec: rec['details'][0:8])
    .addfield('seq_fn', determine_gatc_seq_fn)
    .addfield('company', 'GATC')
    .addfield('primers', which_primers)
    .addfield('seq_length', lambda rec: sequence_length(rec['seq_fn']))
)
print(len(tbl_gatc_metadata))
tbl_gatc_metadata.displayall(index_header=True)

# <codecell>

# 424057201_PA0020C_For_12891_A02.ab1

tbl_sb_metadata = (
    etl.fromtsv(SB_EXCEL_FN)
    .selecteq('type', 'text')
    .rename('sample', 'ox_code')
    .addfield('seq_fn', lambda rec: "%s/%s" % (SB_SEQ_DIR, rec['file']))
    .addfield('company', 'Source_Bioscience')
    .addfield('seq_length', lambda rec: sequence_length(rec['seq_fn']))
)
print(len(tbl_sb_metadata))
tbl_sb_metadata.displayall(index_header=True)

# <codecell>

vcf_reader = vcf.Reader(filename=PGV4_VCF_FN)
pgv4_samples = etl.wrap([[x] for x in vcf_reader.samples]).pushheader(['ox_code']).addfield('is_in_pgv4_vcf', True)

# <codecell>

tbl_sample_sets

# <codecell>

tbl_pgv4_sample_manifest

# <codecell>

pgv4_samples

# <codecell>

tbl_capillary_samples = (
    tbl_gatc_metadata
    .cut(['ox_code', 'company'])
    .distinct('ox_code')
    .cat(tbl_sb_metadata.cut(['ox_code', 'company']).distinct('ox_code'))
    .leftjoin(tbl_sample_sets, lkey='ox_code', rkey='sample_code')
    .addfield('have_sequenom_data', lambda rec: rec['source_code'] is not None)
    .cut(['ox_code', 'company', 'have_sequenom_data'])
    .leftjoin(pgv4_samples)
    .leftjoin(tbl_pgv4_sample_manifest.selecteq('Exclude', 0), lkey='ox_code', rkey='Sample')
    .leftjoin(tbl_pgv4_sample_manifest, lkey='ox_code', rkey='Sample')
#     .cut(['ox_code', 'company', 'have_sequenom_data', 'is_in_pgv4_vcf', 'Pf3k v3'])
    .cut(['ox_code', 'company', 'have_sequenom_data', 'is_in_pgv4_vcf'])
    .replace('have_sequenom_data', False, None)
#     .convert('Pf3k v3', bool)
#     .replace('Pf3k v3', False, None)
    .leftjoin(tbl_gatc_metadata.selecteq('primer_name', 'Forward-DK_sqnm_12891').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'length_gatc')
    .leftjoin(tbl_sb_metadata.selecteq('primer_name', 'Forward - DK_sqnm_12891').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'length_sb')
    .addfield('forward', lambda rec: rec['length_gatc'] if rec['length_gatc'] is not None else rec['length_sb'])
    .cutout('length_gatc', 'length_sb')
    .leftjoin(tbl_gatc_metadata.selecteq('primer_name', 'Reverse-DK_sqnm_12825').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'length_gatc')
    .leftjoin(tbl_sb_metadata.selecteq('primer_name', 'Reverse - DK_sqnm_12825').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'length_sb')
    .addfield('reverse', lambda rec: rec['length_gatc'] if rec['length_gatc'] is not None else rec['length_sb'])
    .cutout('length_gatc', 'length_sb')
    .leftjoin(tbl_gatc_metadata.selecteq('primer_name', 'DGGE_12825_seq1').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'DGGE_seq1')
    .leftjoin(tbl_gatc_metadata.selecteq('primer_name', 'DGGE_12825_seq2').cut(['ox_code', 'seq_length']), lkey='ox_code', rkey='ox_code')
    .rename('seq_length', 'DGGE_seq2')
    .sort(['company', 'ox_code'])
)
print(len(tbl_capillary_samples.data()))
tbl_capillary_samples.displayall()

# <codecell>

tbl_capillary_samples.toxlsx(CAPILLARY_SAMPLES_FN)

# <codecell>

!{BWA_EXE} index {REF_GENOME}

# <codecell>

def call_variants(company='GATC', ox_code='PA0155-C', primers='Forward',
              seq_fn='/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150527 raw data/00AJ08.fas'):
    read_group_info = "@RG\tID:%s_%s_%s\tSM:%s_%s_%s\tPL:capillary" % (ox_code, company, primers, ox_code, primers, company)
    bam_fn = "%s/%s_%s.bam" % (VARIANT_FILES_DIR, ox_code, company, primers)
    vcf_fn = "%s/%s_%s.vcf" % (VARIANT_FILES_DIR, ox_code, company, primers)
    bwa_command = "%s mem -M -R '%s' %s \"%s\" 2> /dev/null | samtools view -bSh - > %s 2> /dev/null" % (
        BWA_EXE, read_group_info, REF_GENOME, seq_fn, bam_fn
    )
#     print(bwa_command)
    !{bwa_command}
    
    samtools_command = "samtools index %s" % bam_fn
#     print(samtools_command)
    !{samtools_command}
    
    gatk_command = "%s \
 -T UnifiedGenotyper \
 -R %s \
 -I %s \
 -ploidy 1 \
 --genotype_likelihoods_model BOTH \
 --defaultBaseQualities 60 \
 -stand_call_conf 0 \
 -stand_emit_conf 0 \
 -minIndelCnt 1 \
 -o %s > /dev/null 2>&1" % (GATK_EXE, REF_GENOME, bam_fn, vcf_fn)
#     print(gatk_command)
    !{gatk_command}

#      -minReadsPerAlignStart 0 \

#     !{BWA_EXE} mem -M -R {read_group_info} {REF_GENOME} "{seq_fn}" > {sam_fn}

# <codecell>

def run_bwa(ox_code='PA0155-C', company='GATC', primers='Forward',
            seq_fn='/data/plasmodium/pfalciparum/recon/data/raw/capillary/GATC/LIGHTrun 150527 raw data/00AJ08.fas'):
    read_group_info = "@RG\tID:%s_%s_%s\tSM:%s_%s_%s\tPL:capillary" % (ox_code, company, primers, ox_code, primers, company)
    bam_fn = "%s/%s_%s_%s.bam" % (BAM_FILES_DIR, ox_code, company, primers)
#     vcf_fn = "%s/%s_%s.vcf" % (VARIANT_FILES_DIR, ox_code, company, primers)
    
    fastq_fn = convert_to_fastq(seq_fn)
    
    bwa_command = "%s mem -M -R '%s' %s \"%s\" 2> /dev/null | samtools view -bSh - > %s 2> /dev/null" % (
        BWA_EXE, read_group_info, REF_GENOME, fastq_fn, bam_fn
    )
    !{bwa_command}
    
    samtools_command = "samtools sort %s %s" % (bam_fn, bam_fn.replace('.bam', ''))
    !{samtools_command}
    
    samtools_command = "samtools index %s" % bam_fn
    !{samtools_command}
    
    return(bam_fn)

# <codecell>

run_bwa()

# <codecell>

fo = open(ALL_BAM_FILES_FOFN, 'w')
for rec in tbl_gatc_metadata.data():
    bam_fn = run_bwa(rec[8], rec[10], rec[11], rec[9])
    print(bam_fn, file=fo)
for rec in tbl_sb_metadata.data():
    bam_fn = run_bwa(rec[3], rec[8], rec[4], rec[7])
    print(bam_fn, file=fo)
fo.close()

fo = open(GATC_BAM_FILES_FOFN, 'w')
for rec in tbl_gatc_metadata.data():
    bam_fn = run_bwa(rec[8], rec[10], rec[11], rec[9])
    print(bam_fn, file=fo)
fo.close()

fo = open(SB_BAM_FILES_FOFN, 'w')
for rec in tbl_sb_metadata.data():
    bam_fn = run_bwa(rec[3], rec[8], rec[4], rec[7])
    print(bam_fn, file=fo)
fo.close()

# <codecell>

def call_variants(bam_fn=BAM_FILES_FOFN, vcf_fn=CAPILLARY_VCF_FN):
    gatk_command = "%s \
 -T UnifiedGenotyper \
 -R %s \
 -I %s \
 -ploidy 1 \
 --genotype_likelihoods_model BOTH \
 -stand_call_conf 0 \
 -stand_emit_conf 0 \
 --min_base_quality_score 40 \
 -minIndelCnt 1 \
 -o %s > /dev/null 2>&1" % (GATK_EXE, REF_GENOME, bam_fn, vcf_fn)
    print(gatk_command)
    !{gatk_command}

#  --defaultBaseQualities 60 \

# <codecell>

calling_runs = {
    'all':[ALL_BAM_FILES_FOFN, ALL_CAPILLARY_VCF_FN],
    'GATC':[GATC_BAM_FILES_FOFN, GATC_CAPILLARY_VCF_FN],
    'SB':[SB_BAM_FILES_FOFN, SB_CAPILLARY_VCF_FN],
}
for calling_run in calling_runs:
    print(calling_run)
    call_variants(calling_runs[calling_run][0], calling_runs[calling_run][1])

# <headingcell level=1>

# Create table of capillary calls

# <codecell>

import petlx.bio

# <codecell>

tbl_capillary_calls = collections.OrderedDict()
for callset in ['all', 'GATC', 'SB']:
    tbl_capillary_calls[callset] = (etl
        .fromvcf(calling_runs[callset][1])
        .vcfmeltsamples()
        .vcfunpackcall()
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT'])
        .addfield('variant', lambda rec: "%s__%s__%s__%s" % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT']))
        .pivot('variant', 'SAMPLE', 'GT', max)
    )
    tbl_capillary_calls[callset].displayall()

# <codecell>

for callset in ['all', 'GATC', 'SB']:
    print(callset, len(tbl_capillary_calls[callset]))

# <headingcell level=1>

# Create table of merged calls from capillary, illumina and sequenom

# <codecell>

tbl_sequenom_calls = (etl
    .fromvcf(SQNM_VCF_FN)
).displayall()

# <codecell>

SQNM_VCF_FN

# <codecell>

tbl_sequenom_calls = collections.OrderedDict()
callsets = {'all':['GATC', 'Source_Bioscience'], 'GATC':['GATC'], 'SB':['Source_Bioscience']}
for callset in callsets:
    tbl_sequenom_calls[callset] = (etl
        .fromvcf(SQNM_VCF_FN)
        .select(lambda rec: rec['CHROM']=='Pf3D7_13_v3' and rec['POS'] >= 1724843 and rec['POS'] <= 1725958)
        .vcfmeltsamples()
        .vcfunpackcall()
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT'])
        .addfield('variant', lambda rec: "%s__%s__%s__%s" % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT']))
        .selectin('SAMPLE', tbl_capillary_samples.selectin('company', callsets[callset]).values('ox_code'))
        .convert('SAMPLE', lambda val: "%s__SQNM" % val)
        .pivot('variant', 'SAMPLE', 'GT', max)
    )
    print(callset)
    tbl_sequenom_calls[callset].displayall()

# <codecell>

tbl_capillary_calls.header()

# <codecell>

(tbl_capillary_calls
    .outerjoin(tbl_sequenom_calls)
    .cut(['variant', 'PH0022-C', 'PH0022-C_Forward_GATC', 'PH0022-C_Reverse_GATC', 'PH0022-C_seq1_GATC', 'PH0022-C_seq2_GATC'])
    .displayall()
)

# <codecell>

(tbl_capillary_calls
    .outerjoin(tbl_sequenom_calls)
    .cut(['variant', 'PA0020-C', 'PA0020-C_FOR_Source_Bioscience', 'PA0020-C_REV_Source_Bioscience'])
    .displayall()
)

# <codecell>

(tbl_capillary_calls
    .outerjoin(tbl_sequenom_calls)
    .cut(['variant', 'PH0024-C', 'PH0024-C_FOR_Source_Bioscience', 'PH0024-C_REV_Source_Bioscience'])
    .select(lambda rec: rec['PH0024-C'] is not None)
    .displayall()
)

# <codecell>

(tbl_capillary_calls
    .outerjoin(tbl_sequenom_calls)
    .cut(['variant', 'PM0090-C', 'PM0090-C_FOR_Source_Bioscience', 'PM0090-C_REV_Source_Bioscience'])
    .select(lambda rec: rec['PM0090-C'] is not None)
    .displayall()
)

# <codecell>

PGV4_VCF_FN

# <codecell>

# etl.fromvcf(PGV4_VCF_FN, chrom='Pf3D7_13_v3', start=1724843, stop=1725958)
etl.fromvcf(PGV4_VCF_FN, chrom='Pf3D7_13_v3')

# <codecell>

vcf_reader = vcf.Reader(filename=PGV4_VCF_FN)
pgv4_samples = np.array(vcf_reader.samples)
capillary_samples = tbl_capillary_samples.values('ox_code').array()
gatc_capillary_samples = tbl_capillary_samples.selecteq('company', 'GATC').values('ox_code').array()
sb_capillary_samples = tbl_capillary_samples.selecteq('company', 'Source_Bioscience').values('ox_code').array()

common_samples = np.array(np.intersect1d(capillary_samples, pgv4_samples), dtype=[('sample', '<U45')])
gatc_common_samples = np.array(np.intersect1d(gatc_capillary_samples, pgv4_samples), dtype=[('sample', '<U45')])
sb_common_samples = np.array(np.intersect1d(sb_capillary_samples, pgv4_samples), dtype=[('sample', '<U45')])
etl.fromarray(common_samples).data().totsv(SAMPLES_FN % 'all')
etl.fromarray(gatc_common_samples).data().totsv(SAMPLES_FN % 'GATC')
etl.fromarray(sb_common_samples).data().totsv(SAMPLES_FN % 'SB')

# <codecell>

for callset in ['all', 'GATC', 'SB']:
    output_fn = "%s/pgv4_gatk_subset_to_kelch13_%s.vcf.gz" % (DATA_DIR, callset)
    samples_fn = SAMPLES_FN % callset
    !{GATK_EXE} -T SelectVariants \
        -R {REF_GENOME} \
        -V {PGV4_VCF_FN} \
        -L Pf3D7_13_v3:1724843-1725958 \
        -o {output_fn} \
        --sample_file {samples_fn} \
        --unsafe LENIENT_VCF_PROCESSING

# <codecell>

etl.fromvcf("%s/pgv4_gatk_subset_to_kelch13_GATC.vcf.gz" % (DATA_DIR))

# <codecell>

tbl_illumina_calls = (etl
    .fromvcf("%s/pgv4_gatk_subset_to_kelch13.vcf.gz" % (DATA_DIR))
    .vcfmeltsamples()
    .vcfunpackcall()
    .cut(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'AD'])
    .addfield('variant', lambda rec: "%s__%s__%s__%s" % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT']))
    .pivot('variant', 'SAMPLE', 'AD', max)
)
tbl_illumina_calls.displayall()

# <codecell>

def determine_GT(AD):
    if sum(AD) < 5:
        return(None)
    if len(AD) == 3:
        if AD[0] > 1 and AD[1] > 1:
            return("0/1")
        if AD[0] > 1 and AD[2] > 1:
            return("0/2")
        if AD[1] > 1 and AD[2] > 1:
            return("1/2")
    if(AD[0] > 1 and AD[1] > 1):
        return("0/1")
    if(AD[0] <= 1 and AD[1] > 1):
        return("1/1")
    if(AD[0] > 1 and AD[1] <= 1):
        return("0/0")
    else:
        return("unknown")
    

# <codecell>

def is_variable(rec):
    if all([val in ['0/0', None] for val in np.array(rec)[1:]]):
        return False
    else:
        return True

# <codecell>

tbl_illumina_calls = collections.OrderedDict()
for callset in ['all', 'GATC', 'SB']:
    tbl_illumina_calls[callset] = (etl
        .fromvcf("%s/pgv4_gatk_subset_to_kelch13_%s.vcf.gz" % (DATA_DIR, callset))
        .vcfmeltsamples()
        .vcfunpackcall()
        .cut(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'AD'])
        .addfield('variant', lambda rec: "%s__%s__%s__%s" % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT']))
        .addfield('GT', lambda rec: determine_GT(rec['AD']))
        .convert('SAMPLE', lambda val: "%s__ILMN" % val)
        .selectnotnone('GT')
    #     .selectne('GT', '0/0')
        .pivot('variant', 'SAMPLE', 'GT', max)
        .select(is_variable)
    #     .replaceall(None, '0/0')
    )
    tbl_illumina_calls[callset].displayall()

# <codecell>

(etl
    .fromvcf("%s/pgv4_gatk_subset_to_kelch13.vcf.gz" % (DATA_DIR))
    .vcfmeltsamples()
    .vcfunpackcall()
    .cut(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'AD'])
    .addfield('variant', lambda rec: "%s__%s__%s__%s" % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT']))
    .addfield('GT', lambda rec: determine_GT(rec['AD']))
    .convert('SAMPLE', lambda val: "%s__ILMN" % val)
    .selectnone('GT')
).displayall()

# <codecell>

tbl_merged_calls = (tbl_capillary_calls
    .outerjoin(tbl_sequenom_calls)
    .outerjoin(tbl_illumina_calls)
    .rename('variant', 'Concat_variant')
)
tbl_merged_calls = tbl_merged_calls.cut(sort(tbl_merged_calls.header()).tolist())
tbl_merged_calls.toxlsx(MERGED_CALLS_EXCEL_FN)
tbl_merged_calls.displayall()

# <codecell>

tbl_merged_calls = collections.OrderedDict()
for callset in ['all', 'GATC', 'SB']:
    tbl_merged_calls[callset] = (tbl_capillary_calls[callset]
        .outerjoin(tbl_sequenom_calls[callset])
        .outerjoin(tbl_illumina_calls[callset])
        .rename('variant', 'Concat_variant')
    )
    tbl_merged_calls[callset] = tbl_merged_calls[callset].cut(sort(tbl_merged_calls[callset].header()).tolist())
    tbl_merged_calls[callset].toxlsx(MERGED_CALLS_EXCEL_FN, callset)
    print(callset)
    tbl_merged_calls[callset].displayall()

# <codecell>

sort(tbl_merged_calls.header()).tolist()

# <codecell>

1724943 + (15*68-8)

# <codecell>

1725954 - 1724943

# <codecell>

1726005 - 31

# <codecell>

1725986 - 31

# <codecell>

1724996 - 7

# <codecell>

1724989 

