# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Imports

# <codecell>

%run standard_imports.ipynb

# <headingcell level=1>

# Setup

# <codecell>

PGV_METADATA_FN = '../../meta/Pf/PGV4_mk5.xlsx'
PF3K_PANOPTES_FN = '../../meta/Pf/Pf3k_release3_panoptes_samples.txt'
PF_SOLARIS_FN = '../../meta/Pf/PF_samples_Richard.xlsx'

DATA_DIR = '../meta/20150617'
CACHE_DIR = DATA_DIR + '/cache'
!mkdir -p {CACHE_DIR}

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

# etl.fromxlsx(PGV_METADATA_FN, 'PGV4.0', data_only=True)

# <codecell>

@etlcache
def tbl_pgv_metadata():
    return (etl
        .fromxlsx(PGV_METADATA_FN, 'PGV4.0', data_only=True)
#         .select(lambda rec: rec['Pf3k v3'] == True)
#         .select(lambda rec: rec['Exclude'] == False)
        .convert('Year', clean_year)
        .cut(['Sample', 'Pf3k v3', 'Exclude', 'AlfrescoNumber', 'Country', 'SolarisCountry', 'SolarisLocation', 'Code', 'Fws', 'Location', 'Region', 'Year'])
#         .convertnumbers()
    )
print(len(tbl_pgv_metadata.data()))
tbl_pgv_metadata
print(tbl_pgv_metadata.valuecounts('Location').displayall())

# <codecell>

tbl_pgv_metadata

# <codecell>

@etlcache
def tbl_pgv_locations():
    return (etl
        .fromxlsx(PGV_METADATA_FN, 'Locations')
#         .convertnumbers()
    )
print(len(tbl_pgv_locations.data()))
tbl_pgv_locations.selecteq('Location', 'Kassena').displayall()

# <codecell>

len(tbl_pgv_metadata.join(tbl_pgv_locations.distinct('Location'), 'Location').data())

# <codecell>

@etlcache
def tbl_pf3k_metadata():
    return (etl
        .fromtsv(PF3K_PANOPTES_FN)
        .cut(['sample_id', 'study', 'country', 'samplingSite', 'collectionDate', 'basesof5xCoverage', 'averageReadLength', 'meanCoverage', 'PercentCallableGATK'])
#         .convertnumbers()
    )
print(len(tbl_pf3k_metadata.data()))
tbl_pf3k_metadata
print(tbl_pf3k_metadata.valuecounts('samplingSite').displayall())

# <codecell>

def determine_year(date):
    if date is None:
        return None
    if date == '':
        return None
    else:
        return int(date[0:4])

# <codecell>

tbl_merged_samples = (tbl_pgv_metadata
    .outerjoin(tbl_pf3k_metadata, lkey='Sample', rkey='sample_id')
    .addfield('collectionYear', lambda rec: determine_year(rec['collectionDate']))
)
print(len(tbl_merged_samples.data()))
tbl_merged_samples

# <codecell>

tbl_merged_samples.valuecounts('Pf3k v3', 'Exclude').displayall()

# <codecell>

tbl_usable_samples = (tbl_merged_samples
    .select(lambda rec: (rec['Pf3k v3'] == True and rec['Exclude'] == False) or (rec['Pf3k v3'] is None))
)

# <codecell>

print(len(tbl_usable_samples.data()))

# <codecell>

tbl_usable_samples.valuecounts('Location', 'samplingSite').sort('samplingSite').displayall()

# <codecell>

np.sum(tbl_usable_samples.valuecounts('Location', 'samplingSite').values('count').array())

# <codecell>

2072+137

# <codecell>

tbl_usable_samples.valuecounts('Year', 'collectionYear').displayall()

# <codecell>

tbl_usable_samples.valuecounts('Year', 'collectionYear').select(lambda rec: rec['Year'] is not None and rec['collectionYear'] is not None and rec['Year'] != rec['collectionYear']).displayall()

# <codecell>

108+89+4+3

# <codecell>

tbl_usable_samples.select(lambda rec: rec['Year'] is not None and rec['collectionYear'] is None).valuecounts('AlfrescoNumber')

# <codecell>

tbl_usable_samples.select(lambda rec: rec['Year'] is None and rec['collectionYear'] is None).valuecounts('AlfrescoNumber').displayall()

# <codecell>

tbl_usable_samples.select(lambda rec: rec['Year'] is None and rec['collectionYear'] is None and rec['AlfrescoNumber'] == 1044)

# <codecell>

(etl
    .fromxlsx(PGV_METADATA_FN, 'PGV4.0', data_only=True)
    .select(lambda rec: rec['Year'] is None or rec['Year'] == '')
    .cut(['Sample', 'StudyCode', 'AlfrescoNumber', 'Study', 'Year'])
    .selecteq('AlfrescoNumber', 1044)
    .displayall()
)

# <codecell>

def determine_year_of_collection(rec):
    if rec['collectionYear'] is None:
        return(rec['Year'])
    else:
        return(rec['collectionYear'])

# <codecell>

(tbl_usable_samples
    .addfield('year_of_collection', determine_year_of_collection)
    .cut(['Sample', 'AlfrescoNumber'])
)

# <codecell>

@etlcache
def tbl_pf_solaris():
    return (etl
        .fromxlsx(PF_SOLARIS_FN, data_only=True)
#         .convert('Year', clean_year)
        .cut(['ox_code', 'oxford_code', 'oxcode', 'creation_date', 'Dna_Extraction_date', 'sample_creation_date', 'contribution_arrival_date'])
#         .cut(['Sample', 'Pf3k v3', 'Exclude', 'AlfrescoNumber', 'Country', 'SolarisCountry', 'SolarisLocation', 'Code', 'Fws', 'Location', 'Region', 'Year'])
#         .convertnumbers()
    )
print(len(tbl_pf_solaris.data()))
tbl_pf_solaris.display(100)
# print(tbl_pf_solaris.valuecounts('Location').displayall())

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['creation_date'] != rec['sample_creation_date']).displayall()

# <codecell>

print(len(tbl_pf_solaris.select(lambda rec: rec['Dna_Extraction_date'] is not None and rec['sample_creation_date'] is not None and rec['Dna_Extraction_date'] > rec['sample_creation_date'])))
print(len(tbl_pf_solaris.select(lambda rec: rec['Dna_Extraction_date'] is not None and rec['contribution_arrival_date'] is not None and rec['Dna_Extraction_date'] > rec['contribution_arrival_date'])))
print(len(tbl_pf_solaris.select(lambda rec: rec['sample_creation_date'] is not None and rec['contribution_arrival_date'] is not None and rec['sample_creation_date'] < rec['contribution_arrival_date'])))
 

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['Dna_Extraction_date'] is not None and rec['sample_creation_date'] is not None and rec['Dna_Extraction_date'] > rec['sample_creation_date']).display(10)

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['Dna_Extraction_date'] is not None and rec['contribution_arrival_date'] is not None and rec['Dna_Extraction_date'] > rec['contribution_arrival_date']).display(10)

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['sample_creation_date'] is not None and rec['contribution_arrival_date'] is not None and rec['sample_creation_date'] < rec['contribution_arrival_date']).displayall()
                      

# <codecell>

print('Dna_Extraction_date', len(tbl_pf_solaris.selectnone('Dna_Extraction_date')))
print('creation_date', len(tbl_pf_solaris.selectnone('creation_date')))
print('sample_creation_date', len(tbl_pf_solaris.selectnone('sample_creation_date')))
print('contribution_arrival_date', len(tbl_pf_solaris.selectnone('contribution_arrival_date')))

# <codecell>

tbl_pf_solaris.selectnone('sample_creation_date').displayall()

# <codecell>

len(tbl_pf_solaris.selectnone('ox_code'))

# <codecell>

len(tbl_pf_solaris.selectnone('oxcode'))

# <codecell>

len(tbl_pf_solaris.selectnone('oxford_code'))

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['ox_code'] is not None and rec['oxcode'] is not None and rec['ox_code'] != rec['oxcode']).displayall()

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['ox_code'] is not None and rec['oxford_code'] is not None and rec['ox_code'] != rec['oxford_code']).displayall()

# <codecell>

tbl_pf_solaris.select(lambda rec: rec['oxcode'] is not None and rec['oxford_code'] is not None and rec['oxcode'] != rec['oxford_code']).displayall()

# <codecell>

tbl_pf_solaris.duplicates('oxcode')

# <codecell>


