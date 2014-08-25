#!/usr/bin/env python

import sys
from optparse import OptionParser
import errno

import vcf
import numpy as np
from scipy import stats
from vcf.parser import _Info as Info
from vcf.parser import _Filter as Filter
from collections import OrderedDict

__version__ = "0.0.1"

if __name__ == '__main__':
	
	usage = 'usage: %prog [options]'
	description = "Populate the FILTER field of a vcf file based on thresholds supplied. Will remove any existing MaxCoverage and MinCoverage filters."
	epilog = """

Examples:
	
	coverageFilters.py --input_vcf_fn input.vcf --output_vcf_fn output.vcf --min_coding_coverage 215121 --max_coding_coverage 325206

Version: {version} (pyvcf {vcfversion})

""".format(version=__version__, vcfversion=vcf.VERSION)

	OptionParser.format_epilog = lambda self, formatter: self.epilog
	parser = OptionParser(usage=usage, description=description, epilog=epilog)
	parser.add_option('-i', '--input_vcf_fn', dest='input_vcf_fn', help='input vcf filename', default=None)
	parser.add_option('-o', '--output_vcf_fn', dest='output_vcf_fn', help='output vcf filename', default=None)
	parser.add_option('-a', '--min_coding_coverage', dest='min_coding_coverage', type='int', help='Threshold for MinCoverage filter in coding SNPs. Currently defined as 25th percentile of DP of core region coding SNPs.', default=None)
	parser.add_option('-z', '--max_coding_coverage', dest='max_coding_coverage', type='int', help='Threshold for MaxCoverage filter in coding SNPs. Currently defined as 95th percentile of DP of core region coding SNPs.', default=None)
    # parser.add_option('-b', '--min_noncoding_coverage', dest='min_noncoding_coverage', type='int', help='Threshold for MinCoverage filter in noncoding SNPs. Currently defined as 40th percentile of DP of core region noncoding SNPs.', default=None)
    # parser.add_option('-y', '--max_noncoding_coverage', dest='max_noncoding_coverage', type='int', help='Threshold for MaxCoverage filter in noncoding SNPs. Currently defined as 95th percentile of DP of core region noncoding SNPs.', default=None)
	parser.add_option('-c', '--chromosome', dest='chromosome', help='only this chromosome will be output in the vcf. Specifying this can make the code quicker if, for example, the input file contains the whole genome and the output is only required for a single chromosome or region', default=None)
	parser.add_option('-s', '--start', dest='start', type='int', help='only the region starting at this position will be output in the vcf. Specifying this can make the code quicker if, for example, the input file contains the whole genome and the output is only required for a single chromosome or region', default=None)
	parser.add_option('-e', '--end', dest='end', type='int', help='only the region ending at this position will be output in the vcf. Specifying this can make the code quicker if, for example, the input file contains the whole genome and the output is only required for a single chromosome or region', default=None)
	parser.add_option('-p', '--progress', dest='progress', type='int', help='report progress every N rows [1000]', default=1000)
	options, args = parser.parse_args()
	
	try:
		
        # def change_coverage_filters(input_vcf_fn, output_vcf_fn, min_coding_coverage, max_coding_coverage, min_noncoding_coverage, max_noncoding_coverage, chromosome, start, end, progress=1000, **kwargs):
		def change_coverage_filters(input_vcf_fn, output_vcf_fn, min_coding_coverage, max_coding_coverage, chromosome, start, end, progress=1000, **kwargs):
			with open(input_vcf_fn, 'rb') as input_file, open(output_vcf_fn, 'wb') as output_file:
				
				# set up VCF reader and writer
				if chromosome != None and start != None and end != None:
					inputReader = vcf.Reader(input_file).fetch(chromosome, start, end)
				else:
					inputReader = vcf.Reader(input_file)
				
				# add new info to header for DP and AD as not currently included
				inputReader.infos['DP'] = Info(id='DP', num=1, type='Integer', desc='Total post-SNP-o-matic depth')
				inputReader.infos['AD'] = Info(id='AD', num=-2, type='Integer', desc='Total post-SNP-o-matic depth at each allele')
				inputReader.infos['TRI2BI'] = Info(id='TRI2BI', num=0, type='Flag', desc='Originally tri-allelic but third allele removed')
				
                # inputReader.filters['MinCoverage'] = Filter('MinCoverage', 'coding:%dx, noncoding:%dx. These are based on DP from INFO.' % (min_coding_coverage, min_noncoding_coverage))
                # inputReader.filters['MaxCoverage'] = Filter('MaxCoverage', 'coding:%dx, noncoding:%dx. These are based on DP from INFO.' % (max_coding_coverage, max_noncoding_coverage))
				inputReader.filters['MinCoverage'] = Filter('MinCoverage', 'coding:%dx. This is based on DP from INFO.' % (min_coding_coverage))
				inputReader.filters['MaxCoverage'] = Filter('MaxCoverage', 'coding:%dx. This is based on DP from INFO.' % (max_coding_coverage))
				
				writer = vcf.Writer(output_file, template=inputReader)
				 
				# counters
				i, passed = 0, 0
				
				try:
					
					for cur_var in inputReader:
						
						if i > 0 and i % progress == 0:
							print i, passed, cur_var.CHROM, cur_var.POS
						
						cur_var.FILTER = list(set(cur_var.FILTER)-set(['MinCoverage', 'MaxCoverage']))
						if cur_var.INFO['DP'] < min_coding_coverage and 'CODING' in cur_var.INFO.keys():
							cur_var.add_filter('MinCoverage')
						if cur_var.INFO['DP'] > max_coding_coverage and 'CODING' in cur_var.INFO.keys():
							cur_var.add_filter('MaxCoverage')
                        # if cur_var.INFO['DP'] < min_noncoding_coverage and not ('CODING' in cur_var.INFO.keys()):
                        #     cur_var.add_filter('MinCoverage')
                        # if cur_var.INFO['DP'] > max_noncoding_coverage and not ('CODING' in cur_var.INFO.keys()):
                        #     cur_var.add_filter('MaxCoverage')
							
                        # if len(cur_var.FILTER) == 0: # PASS
						if cur_var.FILTER is None or len(cur_var.FILTER) == 0: # PASS
							cur_var.FILTER = 'PASS'
							passed += 1

						writer.write_record(cur_var)

						# move on
						i += 1

				except StopIteration:
					print 'done (%d variants, %d passed)' % (i, passed)
			
		print(options.input_vcf_fn)
		
		change_coverage_filters(input_vcf_fn=options.input_vcf_fn,
			output_vcf_fn=options.output_vcf_fn,
			min_coding_coverage=options.min_coding_coverage,
			max_coding_coverage=options.max_coding_coverage,
            # min_noncoding_coverage=options.min_noncoding_coverage,
            # max_noncoding_coverage=options.max_noncoding_coverage,
			chromosome=options.chromosome,
			start=options.start,
			end=options.end,
			progress=options.progress)
	
	except IOError as e:
		if e.errno == errno.EPIPE:
			pass # ignore broken pipe
		else:
			raise