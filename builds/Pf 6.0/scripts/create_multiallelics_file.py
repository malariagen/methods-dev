#!/usr/bin/env python

import sys
from optparse import OptionParser
import errno

import vcf

__version__ = "0.0.1"

if __name__ == '__main__':
    
    usage = 'usage: %prog [options]'
    description = "Determine whether variants are biallelic, biallelic with spanning deletions ot multi-allelic"
    epilog = """
Examples:
    
    create_multiallelics_file.py --input_vcf_fn /nfs/team112_internal/production_files/Pf/6_0/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz --output_txt_fn /lustre/scratch109/malaria/rp7/data/methods-dev/builds/Pf6.0/20161122_Pf60_final_vcfs/SNP_INDEL_Pf3D7_01_v3.combined.filtered.multiallelics.txt
Version: {version}
""".format(version=__version__)

    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage=usage, description=description, epilog=epilog)
    parser.add_option('-i', '--input_vcf_fn', dest='input_vcf_fn', help='input vcf filename')
    parser.add_option('-o', '--output_txt_fn', dest='output_txt_fn', help='output txt filename', default=None)
    parser.add_option('-p', '--progress', dest='progress', type='int', help='report progress every N rows [10000]', default=10000)
    options, args = parser.parse_args()
    
    try:
        
        def create_multiallelics_file(input_vcf_fn, multiallelics_fn=None, progress=10000):
            if multiallelics_fn is None:
                multiallelics_fn = input_vcf_fn.replace('.vcf.gz', '.txt')
            fo = open(multiallelics_fn, 'w')
            input_vcf_reader = vcf.Reader(filename=input_vcf_fn)
            i = 0
            for record in input_vcf_reader:
                if i > 0 and i % progress == 0:
                    print(i, record.CHROM, record.POS)
                if len(record.ALT) == 1:
                    multiallelic = 'BI'
                elif len(record.ALT) == 2 and (record.ALT[0] == '*' or record.ALT[1] == '*'):
                    multiallelic = 'SD'
                else:
                    multiallelic = 'MU'
                print("%s\t%s\t%s\t%s\t%s" % (
                        record.CHROM,
                        record.POS,
                        record.REF,
                        ",".join([str(x) for x in record.ALT]),
                        multiallelic), file=fo)
                i += 1
            fo.close()

        
        create_multiallelics_file(
            input_vcf_fn=options.input_vcf_fn,
            multiallelics_fn=options.output_txt_fn,
            progress=options.progress
        )
    
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass # ignore broken pipe
        else:
            raise
