#!/usr/bin/env python

import sys
import os
import logging
from logging import info, debug
import optparse
from glob import glob
from subprocess import call, check_call
import vcfnp
import numpy as np


def build_arrays(vcf_fn, region, samples=None, force=False, ploidy=1, fields=['AD'], arities=dict(AD=2), progress=10000, calldata_file_suffix='.calldata.npy'):
    
    variants_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.variants.npy'
    if force or not os.path.exists(variants_array_fn):
        print >>sys.stderr, 'building', variants_array_fn
        V = vcfnp.variants(vcf_fn, region=region, progress=progress,
                           fields=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'num_alleles', 'is_snp', 'svlen', 'CODING', 'NS', 'UQ', 'DP'],
                           dtypes=[('CHROM', 'S12'), ('POS', '<i4'), ('ID', 'S12'), ('REF', 'S12'), ('ALT', 'S12'), ('QUAL', '<f4'), ('FILTER', [('Biallelic', '?'), ('CodingType', '?'), ('HetUniq', '?'), ('MinAlt', '?'), ('MonoAllelic', '?'), ('NoAltAllele', '?'), ('Region', '?'), ('triallelic', '?'), ('PASS', '?')]), ('num_alleles', 'u1'), ('is_snp', '?'), ('svlen', '<i4'), ('CODING', '?'), ('NS', '<i4'), ('UQ', '<i4'), ('DP', '<i4')]
                           # dtypes=['S12', '<i4', 'S12', 'S12', 'S12', '<f4', [('Biallelic', '?'), ('CodingType', '?'), ('HetUniq', '?'), ('MinAlt', '?'), ('MonoAllelic', '?'), ('NoAltAllele', '?'), ('Region', '?'), ('triallelic', '?'), ('PASS', '?')], 'u1', '?', '<i4', '<i4', '<i4']
                           )
        np.save(variants_array_fn, V)
    else:
        print >>sys.stderr, 'skipping', variants_array_fn 
        
    # info_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.info.npy'
    # if force or not os.path.exists(info_array_fn):
    #     print >>sys.stderr, 'building', info_array_fn
    #     I = vcfnp.info(vcf_fn, region=region, progress=progress, fields=['NS', 'UQ', 'CODING', 'DP', 'AD'], vcf_types={'DP': 'Integer', 'AD': 'Integer'}, arities={'AD': 2})
    #     np.save(info_array_fn, I)
    # else:
    #     print >>sys.stderr, 'skipping', info_array_fn
    #
    calldata_array_fn = vcf_fn + '.' + region.replace('-:', '_') + calldata_file_suffix
    if force or not os.path.exists(calldata_array_fn):
        print >>sys.stderr, 'building', calldata_array_fn
        C = vcfnp.calldata(vcf_fn, region=region, progress=progress, ploidy=ploidy, fields=fields, arities=arities, samples=samples)
        np.save(calldata_array_fn, C)
    else:
        print >>sys.stderr, 'skipping', calldata_array_fn

def main(options, args):
    
    # init task index
    if options.task_index is not None:
        task = options.task_index
    elif 'SGE_TASK_ID' in os.environ:
        task = int(os.environ['SGE_TASK_ID'])
    else:
        raise Exception('could not determine task index; %s', options)
    
    info('begin')
    
    # chromosomes
    chromosomes = ['Pf3D7_%02d_v3' % i for i in range(1, 15)]

    # dereference task index to chromosome
    chrom = chromosomes[task-1]
    info(chrom)
    
    # load samples
    samples = np.loadtxt(options.sample_fn, dtype='a15')
    
    # load array
    build_arrays(options.vcf_fn, region=chrom, samples=samples, calldata_file_suffix=options.calldata_file_suffix)
    
    # done
    info('end')
    
    
if __name__ == '__main__':

    # handle command line options
    parser = optparse.OptionParser()
    # N.B., use either -I or -i -s depending on whether input BAM has been split by chromosome or not
    parser.add_option('-v', '--vcf-file', dest='vcf_fn', metavar='FILE', default=None,
                      help='input vcf file')
    parser.add_option('-s', '--sample-file', dest='sample_fn', metavar='FILE',
                      help='file of samples to be loaded, text, one sample per line')
    parser.add_option('-c', '--calldata-file-suffix', dest='calldata_file_suffix', metavar='FILE',
                      default='.calldata.npy',
                      help='suffix for calldata numpy files')
    parser.add_option('-n', '--task-index', dest='task_index', metavar='INTEGER',
                      default=None, type='int',
                      help='override task index (otherwise obtained via environment variable)')
    parser.add_option('-l', '--log-level', dest='log_level', metavar='LEVEL',
                      default='DEBUG',
                      help='logging level')
    options, args = parser.parse_args()
    
    # init logging
    logging.basicConfig(level=getattr(logging, options.log_level.upper()),
                        format='[%(levelname)s] [%(filename)s] %(asctime)s %(message)s')

    main(options, args)
    
