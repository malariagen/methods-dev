# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import vcf

# <codecell>

VCF_FN = '/nfs/team112_internal/production_files/Pf/4_0/pf_4_0_20140712_vfp1.newCoverageFilters_pass_5_99.5.HyperHet.vcf.gz'

# <codecell>

inputReader = vcf.Reader(open(VCF_FN, 'rb'))
samples = inputReader.samples
# i = 0
# for record in inputReader:
#     i = i + 1
#     if i >= 5:
#         break
#     print(record)
cur_var = next(inputReader)

# <codecell>

print(cur_var)

# <codecell>

cur_var.samples[1]['AD']

# <codecell>

i = 0
while True:
    cur_var = next(inputReader)
    i = i + 1
    if i 
    if cur_var.FILTER == ['PASS']:
        print(i)
        break

# <codecell>

def replaceNoneWithZeros(object, numOfElements):
    if object == None:
        return([0 for i in range(numOfElements)])
    elif isinstance(object, int):
        object = [object]
    elif len(object) < numOfElements:
        for i in range(numOfElements-len(object)):
            object.append(0)
    return object

# <codecell>

sample_calls = [cur_var.genotype(s) for s in vcf_samples]
num_alleles = np.size(cur_var.alleles)
num_alts = len(cur_var.ALT)
if num_alleles != num_alts + 1:
    print("error: num_alleles=%d, num_alts=%d" % (num_alleles, num_alts))
ADs = np.array([replaceNoneWithZeros(cur_var.genotype(s)['AD'], num_alleles) for s in samples])

# <codecell>

shape(ADs)

# <codecell>


