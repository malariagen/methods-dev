# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

crosses_dir = '/data/plasmodium/pfalciparum/pf-crosses/data/public/1.0'
!mkdir -p {crosses_dir}

# <headingcell level=1>

# Transfer data

# <codecell>

!rsync -avL oscar.well.ox.ac.uk:{crosses_dir}/*combined* {crosses_dir}

# <codecell>


