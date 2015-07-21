# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Imports

# <codecell>

%run standard_imports.ipynb
%run plotting_setup.ipynb

# <headingcell level=1>

# Setup

# <codecell>

PGV_METADATA_FN = '../../meta/Pf/PGV4_mk5.xlsx'
PF3K_PANOPTES_FN = '../../meta/Pf/Pf3k_release3_panoptes_samples.txt'
PF_SOLARIS_FN = '../../meta/Pf/PF_samples_Richard_RDP.xlsx'
# ASSEMBLED_SAMPLES_FN = '../meta/PacBio samples draft 20150622.xlsx'
ASSEMBLED_SAMPLES_FN = '../meta/PacBio samples draft 26062015_RDP_20150706.xlsx'
PROCESSED_ASSEMBLED_SAMPLES_DIR = '/nfs/team112_internal/production_files/Pf3k/methods/assembled_samples'
!mkdir {PROCESSED_ASSEMBLED_SAMPLES_DIR}

# DATA_DIR = '../meta/20150617'
CACHE_DIR = '../meta/cache'
!mkdir -p {CACHE_DIR}

REF_GENOME = '/data/plasmodium/pfalciparum/recon/roamato/Pf3D7_v3/3D7_sorted.fa'

crosses_dir = '/data/plasmodium/pfalciparum/pf-crosses/data/public/1.0'
regions_fn = '../../../pf-crosses/meta/regions-20130225.bed.gz'

# <codecell>

install_dir = '../opt'
bwa_dir = install_dir + '/bwa'
bwa_exe = bwa_dir + '/bwa-0.7.12/bwa'

samtools_dir = install_dir + '/samtools'
samtools_exe = samtools_dir + '/samtools-1.2/samtools'

picard_dir = install_dir + '/picard'
picard_exe = 'java -jar ' + picard_dir + '/picard-tools-1.133/picard.jar'

bcftools_dir = install_dir + '/bcftools'
bcftools_exe = bcftools_dir + '/bcftools-1.2/bcftools'

gatk_dir = install_dir + '/gatk'
gatk_exe = 'java -jar ' + gatk_dir + '/GenomeAnalysisTK.jar'

snpeff_dir = install_dir + '/snpeff/snpEff'
snpeff_exe = 'java -jar ' + snpeff_dir + '/snpEff.jar'

# <codecell>

samtools_exe

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

def clean_year(year):
    if year is None:
        return None
    if year == '':
        return None
    else:
        return(int(year))

# <codecell>

@etlcache
def tbl_pgv_metadata():
    return (etl
        .fromxlsx(PGV_METADATA_FN, 'PGV4.0', data_only=True)
        .convert('Year', clean_year)
        .cut(['Sample', 'Pf3k v3', 'Exclude', 'AlfrescoNumber', 'Country', 'SolarisCountry', 'SolarisLocation', 'Code', 'Fws', 'Location', 'Region', 'Year'])
    )
print("tbl_pgv_metadata length = %d" % len(tbl_pgv_metadata.data()))

# <codecell>

@etlcache
def tbl_pgv_locations():
    return (etl
        .fromxlsx(PGV_METADATA_FN, 'Locations')
    )
print("tbl_pgv_locations length = %d" % len(tbl_pgv_locations.data()))

# <codecell>

@etlcache
def tbl_pf3k_metadata():
    return (etl
        .fromtsv(PF3K_PANOPTES_FN)
        .cut(['sample_id', 'study', 'country', 'samplingSite', 'collectionDate', 'basesof5xCoverage', 'averageReadLength', 'meanCoverage', 'PercentCallableGATK'])
    )
print("tbl_pf3k_metadata length = %d" % len(tbl_pf3k_metadata.data()))

# <codecell>

def determine_year(date):
    if date is None:
        return None
    if date == '':
        return None
    else:
        return date[0:4]

# <codecell>

@etlcache
def tbl_pf_solaris():
    return (etl
        .fromxlsx(PF_SOLARIS_FN, data_only=True)
        .cut(['ox_code', 'oxford_code', 'oxcode', 'creation_date', 'Dna_Extraction_date', 'sample_creation_date', 'contribution_arrival_date'])
    )
print("tbl_pf_solaris length = %d" % len(tbl_pf_solaris.data()))

# <codecell>

@etlcache
def tbl_assembled_samples():
    return (etl
        .fromxlsx(ASSEMBLED_SAMPLES_FN, data_only=True)
    )
print("tbl_assembled_samples length = %d" % len(tbl_assembled_samples.data()))

# <codecell>


