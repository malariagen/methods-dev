#!/usr/bin/env python

import sys
from optparse import OptionParser
import errno

import allel
import h5py
import collections
import numpy as np

__version__ = "0.0.1"

if __name__ == '__main__':
    
    usage = 'usage: %prog [options]'
    description = "Create summary of genotypes for a single sample"
    epilog = """
Examples:
    
    create_sample_summary.py --hdf5_fn /lustre/scratch111/malaria/rp7/data/methods-dev/builds/Pf6.0/20161128_HDF5_build/hdf5/Pf_60_npy_PID_a12.h5 --index 0
Version: {version}
""".format(version=__version__)

    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage=usage, description=description, epilog=epilog)
    parser.add_option('-f', '--hdf5_fn', dest='hdf5_fn', help='input HDF5 filename', default='/lustre/scratch111/malaria/rp7/data/methods-dev/builds/Pf6.0/20161128_HDF5_build/hdf5/Pf_60_npy_PID_a12.h5')
    parser.add_option('-i', '--index', dest='index', type='int', help='sample index (0-based)', default=0)
    options, args = parser.parse_args()
    
    try:

        def create_sample_summary(index=0, hdf5_fn='/lustre/scratch111/malaria/rp7/data/methods-dev/builds/Pf6.0/20161128_HDF5_build/hdf5/Pf_60_npy_PID_a12.h5'):
            results = collections.OrderedDict()

            hdf = h5py.File(hdf5_fn, 'r')
            samples = hdf['samples'][:]
            results['sample_id'] = samples[index].decode('ascii')
        #     output_fn = "%s/%s.txt" % (output_filestem, results['sample_id'])
        #     fo = open(output_fn, 'w')
            is_pass = hdf['variants']['FILTER_PASS'][:]

            genotypes = allel.GenotypeArray(hdf['calldata']['genotype'][:, [index], :])
        #     genotypes_pass = genotypes[is_pass]

            svlen = hdf['variants']['svlen'][:]
            svlen1 = svlen[np.arange(svlen.shape[0]), genotypes[:, 0, 0] - 1]
            svlen1[np.in1d(genotypes[:, 0, 0], [-1, 0])] = 0 # Ref and missing are considered non-SV
            svlen2 = svlen[np.arange(svlen.shape[0]), genotypes[:, 0, 1] - 1]
            svlen2[np.in1d(genotypes[:, 0, 1], [-1, 0])] = 0 # Ref and missing are considered non-SV
            indel_len = svlen1
            het_indels = (svlen1 != svlen2)
            indel_len[het_indels] = svlen1[het_indels] + svlen2[het_indels]

            is_indel = (indel_len != 0)
            is_inframe = ((indel_len != 0) & (indel_len%3 == 0))
            is_frameshift = ((indel_len != 0) & (indel_len%3 != 0))

            ac = hdf['variants']['AC'][:]
            ac1 = ac[np.arange(ac.shape[0]), genotypes[:, 0, 0] - 1]
            ac1[np.in1d(genotypes[:, 0, 0], [-1, 0])] = 0 # Ref and missing are considered non-singleton
            ac2 = ac[np.arange(ac.shape[0]), genotypes[:, 0, 1] - 1]
            ac2[np.in1d(genotypes[:, 0, 1], [-1, 0])] = 0 # Ref and missing are considered non-singleton

            is_snp = (hdf['variants']['VARIANT_TYPE'][:] == b'SNP')
            is_bi = (hdf['variants']['MULTIALLELIC'][:] == b'BI')
            is_sd = (hdf['variants']['MULTIALLELIC'][:] == b'SD')
            is_mu = (hdf['variants']['MULTIALLELIC'][:] == b'MU')
            is_ins = ((svlen1 > 0) | (svlen2 > 0))
            is_del = ((svlen1 < 0) | (svlen2 < 0))
            is_ins_del_het = (((svlen1 > 0) & (svlen2 < 0)) | ((svlen1 < 0) & (svlen2 > 0))) # These are hets where one allele is insertion and the other allele is deletion (these are probably rare)
            is_coding = (hdf['variants']['CDS'][:])
            is_vqslod6 = (hdf['variants']['VQSLOD'][:] >= 6.0)
            is_hq_snp = (is_pass & is_snp & is_bi & is_coding)
            is_vhq_snp = (is_pass & is_vqslod6 & is_snp & is_bi & is_coding)
            is_nonsynonymous = (hdf['variants']['SNPEFF_EFFECT'][:] == b'NON_SYNONYMOUS_CODING')
            is_synonymous = (hdf['variants']['SNPEFF_EFFECT'][:] == b'SYNONYMOUS_CODING')
            is_frameshift_snpeff = (hdf['variants']['SNPEFF_EFFECT'][:] == b'FRAME_SHIFT')
            is_inframe_snpeff = np.in1d(hdf['variants']['SNPEFF_EFFECT'][:], [b'CODON_INSERTION', b'CODON_DELETION', b'CODON_CHANGE_PLUS_CODON_DELETION', b'CODON_CHANGE_PLUS_CODON_INSERTION'])

            is_singleton = (
                ((ac1 == 1) & (genotypes[:, 0, 0] > 0)) |
                ((ac2 == 1) & (genotypes[:, 0, 1] > 0)) |
                ((ac1 == 2) & (genotypes[:, 0, 0] > 0) & (genotypes[:, 0, 1] > 0))
            )

            is_hom_ref = ((genotypes[:, 0, 0] == 0) & (genotypes[:, 0, 1] == 0))
            is_het = ((genotypes[:, 0, 0] != genotypes[:, 0, 1]))
            is_hom_alt = ((genotypes[:, 0, 0] > 0) & (genotypes[:, 0, 0] == genotypes[:, 0, 1]))
            is_non_ref = ((genotypes[:, 0, 0] > 0) | (genotypes[:, 0, 1] > 0))
            is_missing = ((genotypes[:, 0, 0] == -1))
            is_called = ((genotypes[:, 0, 0] >= 0))

            GQ = hdf['calldata']['GQ'][:, index]
            is_GQ_30 = (GQ >= 30)
            is_GQ_99 = (GQ >= 99)
            DP = hdf['calldata']['DP'][:, index]
            PGT = hdf['calldata']['PGT'][:, index]
            is_phased = np.in1d(PGT, [b'.', b''], invert=True)

            mutations = np.char.add(hdf['variants']['REF'][:][(is_pass & is_snp & is_bi & is_non_ref)], hdf['variants']['ALT'][:, 0][(is_pass & is_snp & is_bi & is_non_ref)])
            is_transition = np.in1d(mutations, [b'AG', b'GA', b'CT', b'TC'])
            is_transversion = np.in1d(mutations, [b'AC', b'AT', b'GC', b'GT', b'CA', b'CG', b'TA', b'TG'])
            is_AT_to_AT = np.in1d(mutations, [b'AT', b'TA'])
            is_CG_to_CG = np.in1d(mutations, [b'CG', b'GC'])
            is_AT_to_CG = np.in1d(mutations, [b'AC', b'AG', b'TC', b'TG'])
            is_CG_to_AT = np.in1d(mutations, [b'CA', b'GA', b'CT', b'GT'])

            results['num_variants']             = genotypes.shape[0]
            results['num_pass_variants']        = np.count_nonzero(is_pass)
            results['num_missing']              = np.count_nonzero(is_missing)
            results['num_pass_missing']         = np.count_nonzero(is_pass & is_missing)
            results['num_called']               = np.count_nonzero(~is_missing)
        #     results['num_called']               = (results['num_variants'] - results['num_missing'])
            results['num_pass_called']          = np.count_nonzero(is_pass & is_called)
        #     results['num_pass_called_2']        = np.count_nonzero(is_pass & ~is_missing)
        #     results['num_pass_called']          = (results['num_pass_variants'] - results['num_pass_missing'])

            results['num_hom_ref']              = np.count_nonzero(is_hom_ref)
            results['num_het']                  = np.count_nonzero(is_het)
            results['num_pass_het']             = np.count_nonzero(is_pass & is_het)
            results['num_hom_alt']              = np.count_nonzero(is_hom_alt)
            results['num_pass_hom_alt']         = np.count_nonzero(is_pass & is_hom_alt)
        #     results['num_pass_non_ref']         = (results['num_pass_het'] + results['num_pass_hom_alt'])
            results['num_pass_non_ref']         = np.count_nonzero(is_pass & is_non_ref)
        #     results['num_variants_2']           = results['num_hom_ref'] + results['num_het'] + results['num_hom_alt'] + results['num_missing']

            results['num_biallelic_het']        = np.count_nonzero(is_pass & is_bi & is_het)
            results['num_biallelic_hom_alt']    = np.count_nonzero(is_pass & is_bi & is_hom_alt)
            results['num_spanning_del_het']     = np.count_nonzero(is_pass & is_sd & is_het)
            results['num_spanning_del_hom_alt'] = np.count_nonzero(is_pass & is_sd & is_hom_alt)
            results['num_multiallelic_het']     = np.count_nonzero(is_pass & is_mu & is_het)
            results['num_multiallelic_hom_alt'] = np.count_nonzero(is_pass & is_mu & is_hom_alt)

        #     results['num_snp_hom_ref']          = np.count_nonzero(is_pass & is_snp & is_het)
            results['num_snp_het']              = np.count_nonzero(is_pass & is_snp & is_het)
            results['num_snp_hom_alt']          = np.count_nonzero(is_pass & is_snp & is_hom_alt)
            results['num_snp']                  = (results['num_snp_het'] + results['num_snp_hom_alt'])
        #     results['num_indel_hom_ref']        = genotypes_pass.subset(~is_snp).count_hom_ref(axis=0)[0]
            results['num_indel_het']            = np.count_nonzero(is_pass & ~is_snp & is_het)
            results['num_indel_hom_alt']        = np.count_nonzero(is_pass & ~is_snp & is_het)
            results['num_indel']                  = (results['num_indel_het'] + results['num_indel_hom_alt'])

            results['num_ins_het']              = np.count_nonzero(is_pass & is_ins & is_het)
            results['num_ins_hom_alt']          = np.count_nonzero(is_pass & is_ins & is_hom_alt)
            results['num_ins']                  = (results['num_ins_hom_alt'] + results['num_ins_het'])
            results['num_del_het']              = np.count_nonzero(is_pass & is_del & is_het)
            results['num_del_hom_alt']          = np.count_nonzero(is_pass & is_del & is_hom_alt)
            results['num_del']                  = (results['num_del_hom_alt'] + results['num_del_het'])

            results['num_coding_het']           = np.count_nonzero(is_pass & is_coding & is_het)
            results['num_coding_hom_alt']       = np.count_nonzero(is_pass & is_coding & is_hom_alt)
            results['num_coding']               = (results['num_coding_het'] + results['num_coding_hom_alt'])

            results['num_hq_snp_called']        = np.count_nonzero(is_hq_snp & ~is_missing)
            results['num_hq_snp_hom_ref']       = np.count_nonzero(is_hq_snp & is_hom_ref)
            results['num_hq_snp_het']           = np.count_nonzero(is_hq_snp & is_het)
            results['num_hq_snp_hom_alt']       = np.count_nonzero(is_hq_snp & is_hom_alt)
            results['num_vhq_snp_called']       = np.count_nonzero(is_vhq_snp & ~is_missing)
            results['num_vhq_snp_hom_ref']      = np.count_nonzero(is_vhq_snp & is_hom_ref)
            results['num_vhq_snp_het']          = np.count_nonzero(is_vhq_snp & is_het)
            results['num_vhq_snp_hom_alt']      = np.count_nonzero(is_vhq_snp & is_hom_alt)

            results['num_singleton']            = np.count_nonzero(is_pass & is_singleton)
            results['num_biallelic_singleton']  = np.count_nonzero(is_pass & is_bi & is_singleton)
            results['num_hq_snp_singleton']     = np.count_nonzero(is_hq_snp & is_singleton)
            results['num_vhq_snp_singleton']    = np.count_nonzero(is_vhq_snp & is_singleton)

            results['num_bi_nonsynonymous']     = np.count_nonzero(is_pass & is_bi & is_snp & is_non_ref & is_nonsynonymous)
            results['num_bi_synonymous']        = np.count_nonzero(is_pass & is_bi & is_snp & is_non_ref & is_synonymous)
        #     results['num_frameshift']           = np.count_nonzero(is_pass & is_indel & is_non_ref & is_coding & is_frameshift)
        #     results['num_inframe']              = np.count_nonzero(is_pass & is_indel & is_non_ref & is_coding & is_inframe)
            results['num_frameshift']           = np.count_nonzero(is_pass & is_indel & is_coding & is_frameshift)
            results['num_inframe']              = np.count_nonzero(is_pass & is_indel & is_coding & is_inframe)
            results['num_bi_frameshift']        = np.count_nonzero(is_pass & is_bi & is_indel & is_coding & is_non_ref & is_frameshift)
            results['num_bi_inframe']           = np.count_nonzero(is_pass & is_bi & is_indel & is_coding & is_non_ref & is_inframe)
            results['num_hq_frameshift']        = np.count_nonzero(is_pass & is_vqslod6 & is_bi & is_indel & is_coding & is_non_ref & is_frameshift)
            results['num_hq_inframe']           = np.count_nonzero(is_pass & is_vqslod6 & is_bi & is_indel & is_coding & is_non_ref & is_inframe)
            results['num_bi_frameshift_snpeff'] = np.count_nonzero(is_pass & is_bi & ~is_snp & is_non_ref & is_frameshift_snpeff)
            results['num_bi_inframe_snpeff']    = np.count_nonzero(is_pass & is_bi & ~is_snp & is_non_ref & is_inframe_snpeff)

            results['num_bi_transition']        = np.count_nonzero(is_transition)
            results['num_bi_transversion']      = np.count_nonzero(is_transversion)
            results['num_bi_AT_to_AT']          = np.count_nonzero(is_AT_to_AT)
            results['num_bi_CG_to_CG']          = np.count_nonzero(is_CG_to_CG)
            results['num_bi_AT_to_CG']          = np.count_nonzero(is_AT_to_CG)
            results['num_bi_CG_to_AT']          = np.count_nonzero(is_CG_to_AT)

            results['num_phased']               = np.count_nonzero(is_pass & is_phased)
            results['num_phased_non_ref']       = np.count_nonzero(is_pass & is_phased & is_non_ref)
            results['num_phased_hom_ref']       = np.count_nonzero(is_pass & is_phased & is_hom_ref)
            results['num_phased_missing']       = np.count_nonzero(is_pass & is_phased & is_missing)

            results['num_GQ_30']                = np.count_nonzero(is_pass & is_called & is_GQ_30)
            results['num_het_GQ_30']            = np.count_nonzero(is_pass & is_het & is_GQ_30)
            results['num_hom_alt_GQ_30']        = np.count_nonzero(is_pass & is_hom_alt & is_GQ_30)
            results['num_GQ_99']                = np.count_nonzero(is_pass & is_called & is_GQ_99)
            results['num_het_GQ_99']            = np.count_nonzero(is_pass & is_het & is_GQ_99)
            results['num_hom_alt_GQ_99']        = np.count_nonzero(is_pass & is_hom_alt & is_GQ_99)

            results['pc_pass']                  = 0.0 if results['num_called'] == 0 else \
                results['num_pass_called'] / results['num_called']
            results['pc_missing']               = 0.0 if results['num_variants'] == 0 else \
                results['num_missing'] / results['num_variants']
            results['pc_pass_missing']          = 0.0 if results['num_pass_variants'] == 0 else \
                results['num_pass_missing'] / results['num_pass_variants']
            results['pc_het']                   = 0.0 if results['num_called'] == 0 else \
                results['num_het'] / results['num_called']
            results['pc_pass_het']              = 0.0 if results['num_pass_called'] == 0 else \
                results['num_pass_het'] / results['num_pass_called']
            results['pc_hq_snp_het']            = 0.0 if results['num_hq_snp_called'] == 0 else \
                results['num_hq_snp_het'] / results['num_hq_snp_called']
            results['pc_vhq_snp_het']           = 0.0 if results['num_vhq_snp_called'] == 0 else \
                results['num_vhq_snp_het'] / results['num_vhq_snp_called']
            results['pc_hom_alt']               = 0.0 if results['num_called'] == 0 else \
                results['num_hom_alt'] / results['num_called']
            results['pc_pass_hom_alt']          = 0.0 if results['num_pass_called'] == 0 else \
                results['num_pass_hom_alt'] / results['num_pass_called']
            results['pc_snp']                   = 0.0 if results['num_pass_non_ref'] == 0 else \
                (results['num_snp_het'] + results['num_snp_hom_alt']) / results['num_pass_non_ref']
        #     results['pc_snp_v2']                = (results['num_snp_het'] + results['num_snp_hom_alt']) / (results['num_snp_het'] + results['num_snp_hom_alt'] + results['num_indel_het'] + results['num_indel_hom_alt'])
            results['pc_biallelic']             = 0.0 if results['num_pass_non_ref'] == 0 else \
                (results['num_biallelic_het'] + results['num_biallelic_hom_alt']) / results['num_pass_non_ref']
            results['pc_spanning_del']          = 0.0 if results['num_pass_non_ref'] == 0 else \
                (results['num_spanning_del_het'] + results['num_spanning_del_hom_alt']) / results['num_pass_non_ref']
            results['pc_mutliallelic']          = 0.0 if results['num_pass_non_ref'] == 0 else \
                (results['num_multiallelic_het'] + results['num_multiallelic_hom_alt']) / results['num_pass_non_ref']
            results['pc_ins']                   = 0.0 if (results['num_ins'] + results['num_del']) == 0 else \
                (results['num_ins'] / (results['num_ins'] + results['num_del']))
            results['pc_coding']                = 0.0 if results['num_pass_non_ref'] == 0 else \
                results['num_coding'] / results['num_pass_non_ref']
        #     results['pc_bi_nonsynonymous']      = results['num_bi_nonsynonymous'] / (results['num_biallelic_het'] + results['num_biallelic_hom_alt'])
            results['pc_bi_nonsynonymous']      = 0.0 if (results['num_bi_nonsynonymous'] + results['num_bi_synonymous']) == 0 else \
                results['num_bi_nonsynonymous'] / (results['num_bi_nonsynonymous'] + results['num_bi_synonymous'])
            results['pc_frameshift']            = 0.0 if (results['num_frameshift'] + results['num_inframe']) == 0 else \
                results['num_frameshift'] / (results['num_frameshift'] + results['num_inframe'])
            results['pc_bi_frameshift']         = 0.0 if (results['num_bi_frameshift'] + results['num_bi_inframe']) == 0 else \
                results['num_bi_frameshift'] / (results['num_bi_frameshift'] + results['num_bi_inframe'])
            results['pc_hq_frameshift']         = 0.0 if (results['num_hq_frameshift'] + results['num_hq_inframe']) == 0 else \
                results['num_hq_frameshift'] / (results['num_hq_frameshift'] + results['num_hq_inframe'])
            results['pc_bi_frameshift_snpeff']  = 0.0 if (results['num_bi_frameshift_snpeff'] + results['num_bi_inframe_snpeff']) == 0 else \
                results['num_bi_frameshift_snpeff'] / (results['num_bi_frameshift_snpeff'] + results['num_bi_inframe_snpeff'])
            results['pc_bi_transition']         = 0.0 if (results['num_bi_transition'] + results['num_bi_transversion']) == 0 else \
                results['num_bi_transition'] / (results['num_bi_transition'] + results['num_bi_transversion'])
            results['pc_bi_AT_to_AT']           = 0.0 if (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT']) == 0 else \
                results['num_bi_AT_to_AT'] / (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT'])
            results['pc_bi_CG_to_CG']           = 0.0 if (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT']) == 0 else \
                results['num_bi_CG_to_CG'] / (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT'])
            results['pc_bi_AT_to_CG']           = 0.0 if (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT']) == 0 else \
                results['num_bi_AT_to_CG'] / (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT'])
            results['pc_bi_CG_to_AT']           = 0.0 if (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT']) == 0 else \
                results['num_bi_CG_to_AT'] / (results['num_bi_AT_to_AT'] + results['num_bi_CG_to_CG'] + results['num_bi_AT_to_CG'] + results['num_bi_CG_to_AT'])
            results['pc_phased']                = 0.0 if results['num_pass_non_ref'] == 0 else \
                results['num_phased_non_ref'] / results['num_pass_non_ref']
            results['pc_phased_hom_ref']        = 0.0 if results['num_phased'] == 0 else \
                results['num_phased_hom_ref'] / results['num_phased']
            results['pc_phased_missing']        = 0.0 if results['num_phased'] == 0 else \
                results['num_phased_missing'] / results['num_phased']
            results['pc_GQ_30']                 = 0.0 if results['num_pass_called'] == 0 else \
                results['num_GQ_30'] / results['num_pass_called']
            results['pc_het_GQ_30']             = 0.0 if results['num_pass_het'] == 0 else \
                results['num_het_GQ_30'] / results['num_pass_het']
            results['pc_hom_alt_GQ_30']         = 0.0 if results['num_pass_hom_alt'] == 0 else \
                results['num_hom_alt_GQ_30'] / results['num_pass_hom_alt']
            results['pc_GQ_99']                 = 0.0 if results['num_pass_called'] == 0 else \
                results['num_GQ_99'] / results['num_pass_called']
            results['pc_het_GQ_99']             = 0.0 if results['num_pass_het'] == 0 else \
                results['num_het_GQ_99'] / results['num_pass_het']
            results['pc_hom_alt_GQ_99']         = 0.0 if results['num_pass_hom_alt'] == 0 else \
                results['num_hom_alt_GQ_99'] / results['num_pass_hom_alt']

            results['mean_GQ']                  = np.mean(GQ[is_pass])
            results['mean_GQ_hom_ref']          = np.mean(GQ[is_pass & is_hom_ref])
            results['mean_GQ_het']              = np.mean(GQ[is_pass & is_het])
            results['mean_GQ_hom_alt']          = np.mean(GQ[is_pass & is_hom_alt])
            results['mean_DP']                  = np.mean(DP[is_pass])
            results['mean_DP_hom_ref']          = np.mean(DP[is_pass & is_hom_ref])
            results['mean_DP_het']              = np.mean(DP[is_pass & is_het])
            results['mean_DP_hom_alt']          = np.mean(DP[is_pass & is_hom_alt])
        #     results['mean_GQ_2']                = np.nanmean(GQ[is_pass])
        #     results['mean_DP_2']                = np.nanmean(DP[is_pass])
        #     results['mean_DP']                  = np.mean(DP)
        #     results['mean_DP_2']                = np.nanmean(DP)

            results['mean_indel_len']           = np.mean(indel_len[is_pass])
            results['total_indel_len']          = np.sum(indel_len[is_pass])

            print('\t'.join([str(x) for x in list(results.keys())]))
            print('\t'.join([str(x) for x in list(results.values())]))

        
        create_sample_summary(
            hdf5_fn=options.hdf5_fn,
            index=options.index
        )
    
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass # ignore broken pipe
        else:
            raise
