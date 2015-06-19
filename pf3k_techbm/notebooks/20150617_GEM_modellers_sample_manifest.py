# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <headingcell level=1>

# Setup

# <codecell>

SAMPLE_MANIFEST_FN = '../../meta/Pf/Pf3k_release3_modellers_sample_manifest.txt'
# LOCATIONS_FN = '../../meta/Pf/Pf3k_release3_modellers_locations.txt'
LOCATIONS_FN = '../../meta/Pf/Pf3k_release3_panoptes_locations.txt'

FIELDS_OF_INTEREST = ['Sample', 'YearOfCollection', 'samplingSite', 'Fws', 'PercentCallableGATK']

# <headingcell level=1>

# Create table

# <codecell>

tbl_merged_samples = (tbl_pgv_metadata
    .outerjoin(tbl_pf3k_metadata, lkey='Sample', rkey='sample_id')
    .addfield('collectionYear', lambda rec: determine_year(rec['collectionDate']))
    .leftjoin(tbl_pf_solaris, lkey='Sample', rkey='oxcode')
)
print(len(tbl_merged_samples.data()))
tbl_merged_samples

# <codecell>

tbl_merged_samples.valuecounts('Pf3k v3', 'Exclude').displayall()

# <codecell>

def determine_year_of_collection(rec):
    if rec['collectionYear'] is None:
        if rec['Year'] is None:
            return '<%s' % determine_year(str(rec['creation_date']))
        else:
            return(rec['Year'])
    else:
        return(rec['collectionYear'])

# <codecell>

tbl_usable_samples = (tbl_merged_samples
    .select(lambda rec: (rec['Pf3k v3'] == True and rec['Exclude'] == False) or (rec['Pf3k v3'] is None))
    .addfield('YearOfCollection', determine_year_of_collection)
    .cut(FIELDS_OF_INTEREST)
)

# <codecell>

print(len(tbl_usable_samples.data()))

# <codecell>

tbl_usable_samples

# <codecell>

print(np.sum(tbl_usable_samples.valuecounts('YearOfCollection').values('count').array()))
tbl_usable_samples.valuecounts('YearOfCollection').displayall()

# <codecell>

print(np.sum(tbl_usable_samples.valuecounts('samplingSite').values('count').array()))
tbl_usable_samples.valuecounts('samplingSite').displayall()

# <codecell>

tbl_panoptes_locations = etl.fromtsv(LOCATIONS_FN)
tbl_panoptes_locations.displayall()

# <codecell>

len(tbl_usable_samples.join(tbl_panoptes_locations, lkey='samplingSite', rkey='location').data())

# <codecell>

tbl_usable_samples.totsv(SAMPLE_MANIFEST_FN)

# <codecell>


