# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# standard libraries
from __future__ import division, print_function
import sys
import os
import operator
import itertools
import collections
import functools
import glob
import csv
import datetime
import bisect
import sqlite3
import subprocess
import traceback
import random
import gc

# third party libraries
# N.B., don't import matplotlib because notebooks and scripts may want to import it in different ways (e.g., if using different backends)
import numpy as np
import scipy
import scipy.stats as stats
import pandas
import numexpr
import sklearn
import sklearn.decomposition
import sklearn.manifold
import pyfasta
import pysam
# import pysamstats
import petl
import petl as etl
etl.repr_index_header = True
import petl.io
import petlx
import vcf
# import vcfnp
import h5py
import tables
import vcfplt
import sh
import Bio
from Bio import SeqIO

# <codecell>

def print_library_versions():
    print('python', sys.version)
    print('numpy', np.__version__)
    print('scipy', scipy.__version__)
    print('pandas', pandas.__version__)
    print('numexpr', numexpr.__version__)
    print('pysam', pysam.__version__)
#     print('pysamstats', pysamstats.__version__)
    print('petl', petl.__version__)
    etl.repr_index_header = True
    print('petlx', petlx.__version__)
    print('vcf', vcf.VERSION)
#     print('vcfnp', vcfnp.__version__)
    print('h5py', h5py.__version__)
    print('tables', tables.__version__)
    print('vcfplt', vcfplt.__version__)

# <codecell>

print_library_versions()

# <codecell>


