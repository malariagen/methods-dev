
# coding: utf-8

# # Setup

# In[3]:

get_ipython().magic('run _shared_setup.ipynb')


# In[2]:

release4_vcfnp_dir = collections.OrderedDict()
release4_vcfnp_dir['7G8'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['GB4'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.GB4.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['KH02'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740936.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['KE01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740937.vcf.gz.vcfnp_cache'
release4_vcfnp_dir['GN01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.ERS740940.vcf.gz.vcfnp_cache'

validation5_vcfnp_dir = collections.OrderedDict()
validation5_vcfnp_dir['7G8'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.7G8.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['GB4'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.GB4.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['KH02'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.KH02.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['KE01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.KE01.vcf.gz.vcfnp_cache'
validation5_vcfnp_dir['GN01'] = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/validation_WG.GN01.vcf.gz.vcfnp_cache'


# In[ ]:




# # Analysis

# In[71]:

def create_variants_array(
    variants_npy_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache/variants.npy'
    calls_npy_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/vcfnp/WG.7G8.vcf.gz.vcfnp_cache/calldata_2d.npy'
):
    variants_array = np.load(variants_npy_fn)
    calls_array = np.load(calls_npy_fn)
    flt_spanning_del = np.array([x[0] != b'*' for x in variants_array['ALT']])
    print(np.unique(flt_spanning_del, return_counts=True))
    alt_allele_num = np.array([int(x.astype(str)[0]) for x in calls_array['GT'][:,0]])
    is_het = np.array([x[0] != x[1] for x in calls_array['GT'][:,0]])
    highest_cov_allele = np.argmax(calls_array['AD'][:,0,:], axis=1)
    alt_allele_num[is_het] = highest_cov_allele[is_het]
    flt_minor_het = (alt_allele_num > 0)
    print(np.unique(flt_minor_het, return_counts=True))
    alt_allele_num = alt_allele_num-1
    alt_allele_num[alt_allele_num == -1] = 0
    alts = variants_array['ALT'][np.arange(len(alt_allele_num)), alt_allele_num]
    variants_array = np.lib.recfunctions.append_fields(
        variants_array,
        'alt',
        alts
    ).data
    flt_all = (
        flt_spanning_del &
        flt_minor_het
    )
    return(variants_array[flt_all])


# In[72]:

variants_array = create_variants_array()


# In[74]:

variants_array.shape


# In[81]:

np.unique(variants_array['SNPEFF_EFFECT'], return_counts=True)


# In[78]:

flt_pass = variants_array['VQSLOD'] > 0
print(np.unique(flt_pass, return_counts=True))
flt_core = variants_array['RegionType'] == b'Core'
print(np.unique(flt_core, return_counts=True))
flt_core_pass = (flt_pass & flt_core)
print(np.unique(flt_core_pass, return_counts=True))


# In[76]:

variants_array['RegionType']


# In[88]:

np.unique((variants_array[['is_snp', 'SNPEFF_EFFECT']]), return_counts=True)
# , np.in1d(variants_array['SNPEFF_EFFECT'], [b'INTERGENIC', b'INTRAGENIC'])))


# In[86]:

variants_array[['is_snp', 'SNPEFF_EFFECT']]


# In[90]:

etl.fromarray(variants_array).display(index_header=True)


# In[91]:

etl.fromarray(variants_array[['CHROM', 'POS', 'REF', 'alt']]).display(index_header=True)


# In[122]:

import copy
def create_consensus(
    variants_array,
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    pos_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/pos/7G8.pos',
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
):
    ref_dict=SeqIO.to_dict(SeqIO.parse(open(ref_genome_fn), "fasta"))
    pos_dir = os.path.dirname(pos_fn)
    if not os.path.exists(pos_dir):
        os.makedirs(pos_dir)
    fasta_dir = os.path.dirname(consensus_fasta_fn)
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)
    pos_fo = open(pos_fn, 'w')
    fasta_fo = open(consensus_fasta_fn, 'w')
    mutable_ref_dict = copy.deepcopy(ref_dict)
    for chrom in mutable_ref_dict:
        mutable_ref_dict[chrom].seq = mutable_ref_dict[chrom].seq.tomutable()
    previous_new_pos = 1
    current_offset = 0
    current_chrom = ''
    for i, rec in enumerate(etl.fromarray(variants_array).data()):
        if i > 0 and current_chrom != rec[0].decode("utf-8"):
            print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  
            current_chrom = rec[0].decode("utf-8")
            current_offset = 0
            previous_new_pos = 1
        reflen = len(rec[2])
        altlen = len(rec[3])
        pos = int(rec[1])
        new_pos = pos + current_offset
        print("%s\t%d\t%d\t%d\t%d" % (rec[0].decode("utf-8"), pos, new_pos, previous_new_pos, current_offset), file=pos_fo)  
        previous_new_pos = new_pos
        startpos = pos + current_offset - 1
        endpos = pos + current_offset + reflen - 1   
        mutable_ref_dict[rec[0].decode("utf-8")].seq[startpos:endpos] = rec[3]
        current_offset = current_offset + altlen - reflen
#         if i%1000 == 0:
#             print(current_chrom, current_offset)

    print("%s\t%d\t%d\t%d\t%d" % (current_chrom, 0, 9999999, previous_new_pos, current_offset), file=pos_fo)  
    
    for chrom in mutable_ref_dict:
        SeqIO.write(mutable_ref_dict[chrom], fasta_fo, "fasta")
    


# In[ ]:

def run_nucmer(
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
    assembly_fasta_fn='/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf7G8.Jul2015.fasta',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
    delta_filter_i=99,
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23'
):
    
    get_ipython().system('{nucmer_dir}/nucmer -p {nucmer_filestem} {consensus_fasta_fn} {assembly_fasta_fn}')
    get_ipython().system('{nucmer_dir}/delta-filter -m -i {delta_filter_i} -l 1000 {nucmer_filestem}.delta > {nucmer_filestem}.filter.delta')
    get_ipython().system('{nucmer_dir}/show-coords -THqcl -o {nucmer_filestem}.filter.delta > {nucmer_filestem}.filter.coords')
    get_ipython().system('{nucmer_dir}/show-snps -CHTr {nucmer_filestem}.filter.delta > {nucmer_filestem}.Csnp')


# In[ ]:

def convert_Csnp(
    sample='7G8',
    consensus_fasta_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/fasta/7G8.fasta',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
):

    Csnp_fn=nucmer_filestem + ".Csnp"
    nucmer_vcf_fn=nucmer_filestem + ".vcf"
    ref_dict=SeqIO.to_dict(SeqIO.parse(open(consensus_fasta_fn), "fasta"))

#     Write VCF header
    fo = open(nucmer_vcf_fn, 'w')
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("##description=This file created with convert_Csnp function in 20150929_consensus_alignment.ipynb\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
#     Write data
    ref_chroms = ['M76611', 'PFC10_API_IRAB'] + ["Pf3D7_%02d_v3" % i for i in range(1, 15)]
    query_chroms = ["Pf%s_MT" % sample, "Pf%s_API" % sample] + ["Pf%s_%02d" % (sample, i) for i in range(1, 15)]
    chrom_dict = dict(zip(ref_chroms, query_chroms))

    tbl_csnp = (etl.fromtsv(Csnp_fn).select(lambda rec: rec[9] == chrom_dict[rec[8]]))
    
    current_state = "SNP"
    ins_sequence = ''
    previous_chrom = ''
    for rec in tbl_csnp:
        (pos, ref, alt, pos2, buff, dist, frm, frm2, chrom, query_chrom) = rec
        pos = int(pos)
        pos2 = int(pos2)
        if (
            frm != '1' or
            not(ref in ['A', 'C', 'T', 'G', '.']) or
            not(alt in ['A', 'C', 'T', 'G', '.'])
        ):
            return("error")
        if alt == '.':
            variant_type = 'Del'
        elif ref == '.':
            variant_type = 'Ins'
        else:
            variant_type = 'SNP'
        if variant_type == 'SNP':
            if current_state == 'Del':
                del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                del_alt = ref_dict[previous_chrom][variant_start-1]
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
            if current_state == 'Ins':
                ins_ref = ref_dict[previous_chrom][variant_start-1]
                ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                ins_sequence = ''
            fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, pos, ref, alt))
            current_state = "SNP"
        elif variant_type == 'Del':
            if current_state == 'Del':
                if pos > (variant_end+1) or chrom != previous_chrom: # i.e. we have moved into a different del so need to print out previous
                    if variant_start < 1 or variant_start > len(ref_dict[chrom]):
                        print(chrom, variant_start, ref_dict[chrom], len(ref_dict[chrom]))
                    del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                    del_alt = ref_dict[previous_chrom][variant_start-1]
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
                    variant_start = pos-1
            if current_state == 'Ins':
                ins_ref = ref_dict[previous_chrom][variant_start-1]
                ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                variant_start = pos-1
            if current_state == 'SNP':
                variant_start = pos-1
            variant_end = pos
            current_state = "Del"
        elif variant_type == 'Ins':
            if current_state == 'Del':
                del_ref = str(ref_dict[previous_chrom][variant_start-1:variant_end].seq)
                del_alt = ref_dict[previous_chrom][variant_start-1]
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, del_ref, del_alt))
                variant_start = pos
            if current_state == 'Ins':
                if pos > (variant_end+1) or chrom != previous_chrom: # i.e. we have moved into a different del so need to print out previous
                    ins_ref = ref_dict[previous_chrom][variant_start-1]
                    ins_alt = ref_dict[previous_chrom][variant_start-1] + ins_sequence
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (previous_chrom, variant_start, ins_ref, ins_alt))
                    variant_start = pos
                    ins_sequence = ''
            if current_state == 'SNP':
                variant_start = pos
            ins_sequence = ins_sequence + alt
            variant_end = pos
            current_state = "Ins"
        else:
            return("error")
        previous_chrom = chrom
    fo.close()
    
    return(nucmer_vcf_fn)


# In[ ]:

def convert_coords(
    pos_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/pos/7G8.pos',
    nucmer_filestem='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8',
):
    
    nucmer_vcf_fn = nucmer_filestem + ".vcf"
    converted_vcf_fn = nucmer_filestem + ".3d7coordinates.vcf"
    
    tbl_vcf = etl.fromvcf(nucmer_vcf_fn).convert('POS', int).addfield('POS1', lambda rec: rec['POS']+1)
    tbl_pos = (etl.fromtsv(pos_fn)
        .setheader(('CHROM', 'var_pos', 'end', 'start', 'diff'))
        .convert('var_pos', int)
        .convert('start', int)
        .convert('end', int)
        .convert('diff', int)
        .select(lambda rec: rec['end'] > rec['start'])
    )
    
    tbl_vcf_converted = (tbl_vcf
        .intervalleftjoin(tbl_pos, lstart='POS', lstop='POS1', lkey='CHROM',
                          rstart='start', rstop='end', rkey='CHROM')
        .addfield('NEW_POS', lambda rec: rec['POS'] if rec['diff'] is None else rec['POS'] - rec['diff'])
    )

#     Write VCF header
    fo = open(converted_vcf_fn, 'w')
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("##description=This file created with tbl_converted_coords function in 20150929_consensus_alignment.ipynb\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

#     Write VCF rows
    for rec in tbl_vcf_converted.data():
        fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (rec[0], rec[13], rec[3], rec[4][0]))
    fo.close()
        
    return tbl_vcf_converted
    
    


# In[ ]:

def annotate_vcf(
    unannotated_vcf_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.3d7coordinates.vcf',
    rewrite=False
)
    
    sorted_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".sorted.vcf")
    leftaligned_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".leftaligned.vcf")
    snpeff_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".snpeff.vcf")
    snpeff_annotated_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".snpeff_annotated.vcf")
    annotated_vcf_fn = unannotated_vcf_fn.replace(".vcf", ".annotated.vcf")
    
    get_ipython().system('vcfstreamsort < {unannotated_vcf_fn} > {sorted_vcf_fn}')

    if (not os.path.isfile(leftaligned_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{gatk_exe} -T LeftAlignAndTrimVariants         -R {REF_GENOME}         -V {sorted_vcf_fn}         -o {leftaligned_vcf_fn} #         2> /dev/null')

    if (not os.path.isfile(snpeff_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{snpeff_exe}         -v -o gatk Pf3D7july2015         {leftaligned_vcf_fn}         -no-downstream         -no-upstream         > {snpeff_vcf_fn} #         2> /dev/null')

    if (not os.path.isfile(snpeff_annotated_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
        get_ipython().system('{gatk_exe}         -T VariantAnnotator         -R {REF_GENOME}         -A TandemRepeatAnnotator         -A SnpEff         --variant {leftaligned_vcf_fn}         --snpEffFile {snpeff_vcf_fn}         -o {snpeff_annotated_vcf_fn} #         2> /dev/null')

    if not os.path.isfile(annotated_vcf_fn+'.gz') or rewrite:
        get_ipython().system("cat {snpeff_annotated_vcf_fn}         | vcf-annotate -a {regions_fn}            -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.'            -c CHROM,FROM,TO,INFO/RegionType         > {annotated_vcf_fn}")

        get_ipython().system('bgzip -f {annotated_vcf_fn}')
        get_ipython().system('tabix -p vcf -f {annotated_vcf_fn}.gz')

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.SNP.vcf')+'.gz') or rewrite:
        get_ipython().system("{gatk_exe}         -T SelectVariants         -R {REF_GENOME}         -V {annotated_vcf_fn}.gz         -selectType SNP         -o {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')} 2> /dev/null")

        get_ipython().system("bgzip -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}")
        get_ipython().system("tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}.gz")

    if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')+'.gz') or rewrite:
        get_ipython().system("{gatk_exe}         -T SelectVariants         -R {REF_GENOME}         -V {annotated_vcf_fn}.gz         -xlSelectType SNP         -o {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')} 2> /dev/null")

        get_ipython().system("bgzip -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}")
        get_ipython().system("tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}.gz")


# In[ ]:

import vcfnp

def run_vcfnp(
    annotated_vcf_fn = '/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.3d7coordinates.annoatated.vcf'
):
    v = vcfnp.variants(
        vcf_fn=annotated_vcf_fn,
        progress=10000,
        arities={
            'ALT': 2,
            'AF': 2,
            'AC': 2,
            'MLEAF': 2,
            'MLEAC': 2,
            'RPA': 3
        },
        dtypes={
            'REF': 'a400', 
            'ALT': 'a400',
            'RegionType': 'a25', 
            'VariantType': 'a40',
            'RU': 'a40',
            'set': 'a40',
            'SNPEFF_AMINO_ACID_CHANGE':'a20',
            'SNPEFF_CODON_CHANGE':'a20',
            'SNPEFF_EFFECT':'a33',
            'SNPEFF_EXON_ID':'a2',
            'SNPEFF_FUNCTIONAL_CLASS':'a8',
            'SNPEFF_GENE_BIOTYPE':'a14',
            'SNPEFF_GENE_NAME':'a20',
            'SNPEFF_IMPACT':'a8',
            'SNPEFF_TRANSCRIPT_ID':'a20',
            'VariantType':'a60',
            'culprit':'a14',
        },
        cache=True
    )
    return(v)


# In[5]:

def find_nearest(a, b):
    nearest_before = abs(a - b[np.searchsorted(b, a)-1])
    nearest_after = abs(a - b[np.searchsorted(b, a) % len(b)])
    return(np.minimum(nearest_before, nearest_after))
find_nearest(
    np.array([-1, 1, 2, 5, 9, 12, 22]),
    np.array([0, 1, 4, 6, 12])
)


# In[ ]:

def distance_to_nearest(
    variants_array,
    nucmer_array
):
    distance_array = np.zeros(len(variants_array), dtype='u4')
    for chrom in np.unique(variants_array['CHROM']):
        print(chrom)
        variants_pos_this_chrom = variants_array['POS'][variants_array['CHROM'] == chrom]
        nucmer_pos_this_chrom = nucmer_array['POS'][nucmer_array['CHROM'] == chrom]
        nearest_this_chrom = find_nearest(variants_pos_this_chrom, nucmer_pos_this_chrom)
        distance_array[variants_array['CHROM'] == chrom] = nearest_this_chrom
    new_array = variants_array.addfield('nearest_nucmer', distance_array)
    return(new_array)
        


# In[ ]:

def within_accessible(
    variants_array,
    assembly_sample='7G8',
    filter_coords_fn='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment/nucmer/7G8.filter.coords'
):
#     Figure out how to do this in numpy. Numpy intervals? Or just use in1d?
    is_accessible = np.zeros(len(variants_array), dtype='bool')
    array_coords = (
        etl.fromtsv(truth_coords_fn)
        .pushheader(['ref_start', 'ref_end', 'assembly_start', 'assembly_end', 'ref_length', 'assembly_length', 'identity',
             'ref_chrom_length', 'assembly_chrom_length', 'ref_pct', 'assembly_pct', 'ref_chrom', 'assembly_chrom',
             'ignore'])
        .cut(['ref_chrom', 'assembly_chrom', 'ref_start', 'ref_end'])
        .convertnumbers()
        .toarray()
    )
    for chrom in np.unique(variants_array['CHROM']):
        assembly_chrom = "Pf%s_%s" % (assembly_sample, chrom.decode("utf-8")[6:8])
        starts = array_coords['ref_start'][
            (array_coords['ref_chrom'] == chrom.decode("utf-8")) &
            (array_coords['assembly_chrom'] == assembly_chrom)
        ]

    is_accessible = 
    new_array = variants_array.addfield('is_accessible', is_accessible)
    return(new_array)
        


# In[ ]:




# In[ ]:

def determine_assembly_fasta(sample='7G8'):
    if sample in ('7G8', 'GB4'):
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/ReleaseJul/Pf%s.Jul2015.fasta" % isolate_code)
    else:
        return("/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/Pf%s.Jun2015.fasta" % isolate_code)


# In[125]:

def run_fdr_pipeline(
    sample='7G8',
    ref_genome_fn='/lustre/scratch110/malaria/rp7/Pf3k/GATKbuild/Pfalciparum_GeneDB_Aug2015/Pfalciparum.genome.fasta',
    consensus_dir='/lustre/scratch110/malaria/rp7/Pf3k/release4_candidate/consensus_alignment',
    delta_filter_i=99,
    nucmer_dir='~/src/github/malariagen/methods-dev/pf3k_techbm/opt_4/mummer/MUMmer3.23'
):

    variants_npy_fn = release4_vcfnp_dir[sample] + '/variants.npy'
    calls_npy_fn = release4_vcfnp_dir[sample] + '/calldata_2d.npy'
    pos_fn = '%s/pos/%s.pos' % (consensus_dir, sample)
    consensus_fasta_fn = '%s/fasta/%s.fasta' % (consensus_dir, sample)
    assembly_fasta_fn = determine_assembly_fasta(sample)
    nucmer_filestem = '%s/nucmer/%s' % (consensus_dir, sample)
    unannotated_vcf_fn = "%s.3d7coordinates.vcf" % nucmer_filestem
    annotated_vcf_fn = "%s.3d7coordinates.annotated.vcf" % nucmer_filestem

    print("Creating array %s" % sample)
    variants_array = create_variants_array(**kwargs)
    flt_pass = variants_array['VQSLOD'] > 0
    flt_core = variants_array['RegionType'] == b'Core'
    flt_core_pass = (flt_pass & flt_core)

    print("Creating consensus %s" % sample)
    create_consensus(variants_array[flt_core_pass], **kwargs)
    
    print("Running nucmer %s" % sample)
    run_nucmer(**kwargs)

    print("Coverting nucmer output to vcf %s" % sample)
    convert_Csnp(**kwargs)

    print("Converting nucmer vcf coordinates %s" % sample)
    _ = convert_coords(**kwargs)

    print("Annotate nucmer vcf %s" % sample)
    annotate_vcf(**kwargs)
    
    print("Create array of nucmer vcf %s" % sample)
    nucmer_array = run_vcfnp(**kwargs)
    
    print("Determining closest nucmer calls %s" % sample)
    variants_nearest_array = distance_to_nearest(variants_array[flt_core_pass], nucmer_array)

    print("Plotting histogram of distances %s" % sample)
    print("Determing accessibility %s" % sample)
    
    print("Plotting histogram of distances %s" % sample)
    print("Plotting histogram of distances %s" % sample)
    print("Plotting histogram of distances %s" % sample)
    
        


# In[124]:

run_fdr_pipeline()


# In[ ]:



