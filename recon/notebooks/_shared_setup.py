# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Imports

# <codecell>

%run standard_imports.ipynb

# <headingcell level=1>

# Setup

# <codecell>

DATA_DIR = '/data/plasmodium/pfalciparum/recon/sqnm_assay_test_1'
CACHE_DIR = DATA_DIR + '/cache'
!mkdir -p {CACHE_DIR}

SQNM_EXCEL_FN = DATA_DIR + '/Sqnm_data_DK1066_W1378_20150603.xlsx'
SAMPLE_MANIFEST_FN = DATA_DIR + '/PGV4_mk5-RDP.xlsx'

REF_GENOME = '/data/plasmodium/pfalciparum/recon/roamato/Pf3D7_v3/3D7_sorted.fa'

PGV4_VCF_FN = '/nfs/team112_internal/production_files/Pf/4_0/pf_4_0_20140712_vfp1.newCoverageFilters_pass_5_99.5.HyperHet.vcf.gz'

BCFTOOLS_EXE = '/Users/rpearson/src/github/malariagen/methods-dev/pf3k_techbm/opt/bcftools/bcftools-1.2/bcftools'
GATK_EXE = 'java -jar /Users/rpearson/src/github/malariagen/methods-dev/pf3k_techbm/opt/gatk/GenomeAnalysisTK.jar'

# <headingcell level=1>

# Functions

# <codecell>

def etlcache(f):
    fn = os.path.join(CACHE_DIR, f.__name__)
    if not os.path.exists(fn):
        etl.topickle(f(), fn)
    return etl.frompickle(fn)
    
def nocache(f):
    fn = os.path.join(CACHE_DIR, f.__name__)
    if os.path.exists(fn):
        os.remove(fn)
    return f()

# <headingcell level=1>

# Load data

# <codecell>

@etlcache
def tbl_sample_sets():
    return (etl
        .fromxlsx(SQNM_EXCEL_FN, 'sample_sets')
#         .convertnumbers()
    )
print(len(tbl_sample_sets.data()))
tbl_sample_sets

# <codecell>

@etlcache
def tbl_pgv4_sample_manifest():
    return (etl
        .fromxlsx(SAMPLE_MANIFEST_FN, 'PGV4.0')
        .convertnumbers()
    )
print(len(tbl_pgv4_sample_manifest.data()))
tbl_pgv4_sample_manifest

# <codecell>

def unique_sample_code(rec):
    if rec[2] == 'WATER':
        return('WATER_%s' % rec[1])
    if rec[2] == 'WT3243-C':
        return('CEPH_human_pool_%s' % rec[1])
    if rec[2].startswith('WT'):
        return('%s_%s_%s' % (rec[2], rec[3], rec[1]))
    else:
        return(rec[2])

def extract_chrom(rec):
    vals = rec[0].split('_')
    chrom = "%s_%s_%s" % (vals[3], vals[4], vals[5])
    return(chrom)
    
def extract_pos(rec):
    vals = rec[0].split('_')
    pos = int(vals[6])
    return(pos)

# <codecell>

@etlcache
def tbl_pivoted_genotypes():
    return (etl
        .fromxlsx(SQNM_EXCEL_FN, 'pivoted_genotypes', data_only=True)
        .addfield('unique_sample_id', unique_sample_code)
#         .cutout(['plate_code', 'position-code', 'sample_code', 'source_code'])
#         .cutout(list(range(4)))
        .cut([14] + list(range(4, 14)))
        .sort('unique_sample_id')
        .transpose()
        .addfield('chrom', extract_chrom)
        .addfield('pos', extract_pos)
        .cut([385, 386] + list(range(1, 385)))
        .sort(['chrom', 'pos'])
        .rename('PH0833-C', 'PH0883-C')
    )
print(len(tbl_pivoted_genotypes.data()))
tbl_pivoted_genotypes

# <codecell>

len(tbl_pivoted_genotypes.toarray()[1])

# <codecell>


