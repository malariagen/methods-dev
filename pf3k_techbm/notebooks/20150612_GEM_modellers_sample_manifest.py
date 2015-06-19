# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run standard_imports.ipynb

PGV_METADATA_FN = '../../meta/Pf/PGV4_mk5.xlsx'
PF3K_PANOPTES_FN = '../../meta/Pf/Pf3k_release3_panoptes_samples.txt'

# <headingcell level=1>

# Load data

# <codecell>

tbl_pgv_metadata = etl.fromxlsx(PGV_METADATA_FN)
tbl_pgv_metadata

