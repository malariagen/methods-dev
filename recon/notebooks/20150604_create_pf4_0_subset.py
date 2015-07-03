# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

SQNM_VCF_FN = DATA_DIR + '/Sqnm_data_DK1066_W1378_20150603.vcf'
SAMPLES_FN = DATA_DIR + '/sequenom_samples.txt'

# <headingcell level=1>

# Look at data

# <codecell>

tbl_pivoted_genotypes.displayall(index_header=True)

# <headingcell level=1>

# Functions

# <headingcell level=1>

# Determine ref and alt

# <codecell>

all_calls_array = tbl_pivoted_genotypes.cutout('chrom', 'pos').toarray()
all_calls_array_2d = np.vstack([all_calls_array[item] for item in all_calls_array.dtype.names]).T

# <codecell>

all_calls_array_2d.dtype

# <codecell>

in_seq_handle = open(REF_GENOME)
ref_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

# <codecell>

number_of_columns = len(tbl_pivoted_genotypes.header())

ref_bases = np.zeros(10, '<U4')
alt_bases = np.zeros(10, '<U4')
for i, rec in enumerate(tbl_pivoted_genotypes.data()):
    ref_bases[i] = ref_dict[rec[0]][rec[1]-1]
    alt_bases_this_snp = []
    for j in range(2, number_of_columns):
        if len(rec[j]) == 1 and rec[j] != 'X' and rec[j] != ref_bases[i] and not(rec[j] in alt_bases_this_snp):
            alt_bases_this_snp.append(rec[j])
    alt_bases[i] = ",".join(alt_bases_this_snp)
    
print(ref_bases)
print(alt_bases)

# <codecell>

ref_calls_array_2d = np.vstack([ref_bases for item in all_calls_array.dtype.names]).T
print(shape(ref_calls_array_2d))
ref_calls_array_2d

# <codecell>

alt_calls_array_2d = np.vstack([alt_bases for item in all_calls_array.dtype.names]).T
print(shape(alt_calls_array_2d))
alt_calls_array_2d

# <codecell>

def determine_vcf_call(genotype, ref, alt):
    if genotype == 'X':
        return("./.")
    if genotype == ref:
        return("0/0")
    if genotype == alt:
        return("1/1")
    if (genotype == ref + alt) or (genotype == alt + ref):
        return("0/1")
    else:
        return("?")
vec_determine_vcf_call = np.vectorize(determine_vcf_call)

# <codecell>

vcf_calls = vec_determine_vcf_call(all_calls_array_2d, ref_calls_array_2d, alt_calls_array_2d)

# <codecell>

len(vcf_calls)

# <codecell>

for temp in vcf_calls:
    print("\t".join(temp))

# <codecell>

vcf_calls.view(np.recarray)

# <codecell>

etl.fromarray(vcf_calls.view(np.recarray))

# <headingcell level=1>

# Write VCF file

# <codecell>

samples_string = "\t".join(np.array(tbl_pivoted_genotypes.header())[2:])

# <codecell>

fo = open(SQNM_VCF_FN, 'w')
fo.write("##fileformat=VCFv4.1\n")
fo.write("##description=This file created with 20150604_create_pf4_0_subset.ipynb\n")
fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % samples_string)
for i, rec in enumerate(tbl_pivoted_genotypes.data()):
    variant_info_string = "%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT" % (rec[0], rec[1], ref_bases[i], alt_bases[i])
    genotypes_string = "\t".join(vcf_calls[i, :])
    fo.write(variant_info_string + "\t" + genotypes_string + "\n")
fo.close()

# <codecell>

!bgzip -f {SQNM_VCF_FN}
!tabix -p vcf -f {SQNM_VCF_FN}.gz

# <headingcell level=1>

# Write samples files

# <codecell>

vcf_reader = vcf.Reader(filename=PGV4_VCF_FN)
pgv4_samples = np.array(vcf_reader.samples)

# <codecell>

sqnm_samples = tbl_pivoted_genotypes.transpose().data().values('pos').array()

# <codecell>

np.intersect1d(sqnm_samples, pgv4_samples)

# <codecell>

common_samples = np.array(np.intersect1d(sqnm_samples, pgv4_samples), dtype=[('sample', '<U45')])
etl.fromarray(common_samples).data().totsv(SAMPLES_FN)

# <codecell>

common_samples[0]

# <codecell>

tbl_pivoted_genotypes.transpose().data().cut('pos').select(lambda rec: rec[0][0]=='P' or rec[0][0]=='Q').data().totsv(SAMPLES_FN)

# <headingcell level=1>

# Create subset of PGV4.0 vcf

# <codecell>

!{BCFTOOLS_EXE} isec -p {DATA_DIR} -o pgv4_subset_to_sqnm_assay_test.vcf.gz -O z {PGV4_VCF_FN} {SQNM_VCF_FN}.gz

# <codecell>

/software/java/bin/java -jar /nfs/team112_internal/production/tools/bin/gatk/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T SelectVariants -R /lustre/scratch109/malaria/pfalciparum/resources/3D7_V3.fasta -V /nfs/team112_internal/production_files/Pf/4_0/pf_4_0_20140712_vfp1.newCoverageFilters_pass_5_99.5.HyperHet.vcf.gz -L ~/test.bed -o ~/test_out.vcf

# <codecell>

!java -version

# <codecell>

!{GATK_EXE} -T SelectVariants \
    -R {REF_GENOME} \
    -V {PGV4_VCF_FN} \
    -L {SQNM_VCF_FN}.gz \
    -o {DATA_DIR}/pgv4_gatk_subset_to_sqnm_assay_test.vcf.gz \
    --sample_file {SAMPLES_FN} \
    --unsafe LENIENT_VCF_PROCESSING

# <codecell>


