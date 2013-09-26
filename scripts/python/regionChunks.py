#!/usr/bin/env python

import sys
from optparse import OptionParser
import errno


__version__ = "0.0.1"

if __name__ == '__main__':
	
	usage = 'usage: %prog [options]'
	description = "Create genome chunks manifest file and associated vcf fofn"
	epilog = """

Examples:
	
	regionChunks --genome data/genome/PvivaxGenomic_PlasmoDB-6.0.shortnames.lowercase.reordered.fasta --chunk_size 10000 --chromosome_name_format gb_CM%06d --first_chromosome_number 442 --last_chromosome_number 455 --chunks_manifest_fn ./meta/pv_1_0_regions_chunks_10000bp.tab

Version: {version}

""".format(version=__version__)

	OptionParser.format_epilog = lambda self, formatter: self.epilog
	parser = OptionParser(usage=usage, description=description, epilog=epilog)
	parser.add_option('-g', '--genome', dest='genome', help='reference genome in FASTA format. This MUST have an associated .fai file', default=None)
	parser.add_option('-s', '--chunk_size', dest='chunk_size', type='int', help='size of each genome chunk in bp', default=10000)
	parser.add_option('-q', '--chromosome_name_format', dest='chromosome_name_format', help='format string for names of chromosomes', default=None)
	parser.add_option('-f', '--first_chromosome_number', dest='first_chromosome_number', type='int', help='integer number of first chromosome', default=None)
	parser.add_option('-l', '--last_chromosome_number', dest='last_chromosome_number', type='int', help='integer number of last chromosome', default=None)
	parser.add_option('-c', '--chunks_manifest_fn', dest='chunks_manifest_fn', help='output chunks manifest filename', default=None)
	parser.add_option('-t', '--vcf_fn_format', dest='vcf_fn_format', help='vcf filename format', default=None)
	parser.add_option('-v', '--vcf_fofn', dest='vcf_fofn', help='output vcf file of filenames (fofn)', default=None)
	options, args = parser.parse_args()
	
	try:
		
		def write_region_chunks(chunk_size, chromosomes, chromosome_lengths, chunks_manifest_fn, vcf_fn_format, vcf_fofn):
			fo = open(chunks_manifest_fn, "w")
			if vcf_fn_format != None and vcf_fofn != None:
			    fvo = open(vcf_fofn, "w")
			for chromosome in chromosomes:
				chunk_start = 1
				while chunk_start <= chromosome_lengths[chromosome]:
					chunk_end = min((chunk_start + chunk_size) - 1, chromosome_lengths[chromosome])
					print >>fo, '%s:%d-%d\t%s:%d..%d\t%s\t%d\t%d' % (chromosome, chunk_start, chunk_end, chromosome, chunk_start, chunk_end, chromosome, chunk_start, chunk_end)
					if vcf_fn_format != None and vcf_fofn != None:
					    print >>fvo, vcf_fn_format % '%s:%d-%d' % (chromosome, chunk_start, chunk_end)
					chunk_start += chunk_size
					chunk_end += chunk_size
			fo.close()
			fvo.close()
			
		genome_fai = options.genome + '.fai'
		chromosome_lengths = dict(zip([line.split()[0] for line in open(genome_fai)], [int(line.split()[1]) for line in open(genome_fai)]))
		if options.chromosome_name_format == None: # here we use every chromosome in the reference genome
			chromosomes = chromosome_lengths.keys()
			chromosomes.sort()
		else: # here we only use chromosomes specified using supplied parameters - this is important in Pv where there are many short contigs
			chromosomes = [options.chromosome_name_format % n for n in range(options.first_chromosome_number, options.last_chromosome_number + 1)]
		if options.chunks_manifest_fn == None:
			chunks_manifest_fn = 'meta/pv_1_0_regions_chunks_%dbp.tab' % options.chunk_size
		else:
			chunks_manifest_fn = options.chunks_manifest_fn

		write_region_chunks(chunk_size=options.chunk_size,
			chromosomes=chromosomes,
			chromosome_lengths=chromosome_lengths,
			chunks_manifest_fn=chunks_manifest_fn,
			vcf_fn_format=options.vcf_fn_format,
			vcf_fofn=options.vcf_fofn)
	
	except IOError as e:
		if e.errno == errno.EPIPE:
			pass # ignore broken pipe
		else:
			raise
