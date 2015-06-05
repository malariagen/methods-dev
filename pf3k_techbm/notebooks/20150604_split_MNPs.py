# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

from collections import OrderedDict

# <codecell>

methods = ['freebayes', 'HC', 'mpileup', 'pilon', 'platypus', 'UGT']

vcf_files = OrderedDict()
vcf_files['feb_2015'] = OrderedDict()
vcf_files['nov_2014'] = OrderedDict()

for method in methods:
    

