# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=2>

# Prepare VCF

# <codecell>

import os

# <codecell>

vcf_fn = '../data/production_files/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.gz'
!ls -l {vcf_fn}*

# <codecell>

!zcat {vcf_fn} | head -n19

# <codecell>

!tabix -p vcf {vcf_fn}

# <codecell>

vcf_bgz_fn = '../data/plasmodium/pfalciparum/methods_dev/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.bgz'

# <codecell>

if not os.path.exists(vcf_bgz_fn):
    !zcat {vcf_fn} | bgzip -c > {vcf_bgz_fn}

# <codecell>

if not os.path.exists(vcf_bgz_fn + '.tbi'):
    !tabix -p vcf {vcf_bgz_fn}

# <headingcell level=2>

# Build arrays

# <codecell>

import vcfnp
import numpy as np
import sys

# <codecell>

def build_arrays(vcf_fn, region, samples=None, force=False):
    
    variants_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.variants.npy'
    if force or not os.path.exists(variants_array_fn):
        print >>sys.stderr, 'building', variants_array_fn
        V = vcfnp.variants(vcf_fn, region=region, progress=20000)
        np.save(variants_array_fn, V)
    else:
        print >>sys.stderr, 'skipping', variants_array_fn 
        
    info_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.info.npy'
    if force or not os.path.exists(info_array_fn):
        print >>sys.stderr, 'building', info_array_fn
        I = vcfnp.info(vcf_fn, region=region, progress=20000, fields=['NS', 'UQ', 'CODING', 'DP', 'AD'], vcf_types={'DP': 'Integer', 'AD': 'Integer'}, arities={'AD': 2})
        np.save(info_array_fn, I)
    else:
        print >>sys.stderr, 'skipping', info_array_fn 
        
    

# <codecell>

for chrom in ['Pf3D7_%02d_v3' % i for i in range(1, 15)]:
    build_arrays(vcf_bgz_fn, region=chrom)

# <codecell>

def load_arrays(vcf_fn, region):
    variants_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.variants.npy'
    info_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.info.npy'
    return np.load(variants_array_fn).view(np.recarray), np.load(info_array_fn).view(np.recarray)

# <headingcell level=2>

# Investigate variants

# <codecell>

chrom = 'Pf3D7_01_v3'

# <codecell>

V, I = load_arrays(vcf_bgz_fn, region=chrom)

# <codecell>

V

# <codecell>

I

# <codecell>

figure(figsize=(16, 4))
hist(I.DP, bins=linspace(0, 700000, 100), color='w')
xlabel('DP')
ylabel('frequency')
title('all SNPs (%s)' % chrom)

# <codecell>

figure(figsize=(16, 4))
hist(I.DP[I.CODING], bins=linspace(0, 600000, 150), histtype='step', color='b', linewidth=2, label='all coding')
hist(I.DP[~I.CODING], bins=linspace(0, 600000, 150), histtype='step', color='r', linewidth=2, label='all non-coding')
xlabel('DP')
ylabel('frequency')
title('all SNPs (%s)' % chrom)
legend()

# <codecell>

figure(figsize=(16, 4))
all_coding = I.DP[I.CODING]
hist(all_coding, bins=linspace(0, 600000, 150), histtype='step', color='b', linewidth=2, label='all coding')
axvline(percentile(all_coding, 15), color='k', linestyle=':')
axvline(percentile(all_coding, 85), color='k', linestyle=':')
xlabel('DP')
ylabel('frequency')
title('all SNPs (%s)' % chrom)
legend()

# <codecell>

from petl.interactive import etl
import petlx.array
import petlx.interval

# <codecell>

tbl_regions = (etl
    .fromtsv('../data/plasmodium/pfalciparum/methods_dev/Pf/misc/regions-20130225.onebased.txt')
    .pushheader(['chrom', 'start', 'stop', 'region'])
    .convertnumbers()
)
tbl_regions.head()

# <codecell>

tbl_variants = etl.fromarray(V).cut('CHROM', 'POS')
tbl_variants.head()

# <codecell>

tbl_variants_with_regions = tbl_variants.intervalleftjoin(tbl_regions, lfacet='CHROM', rfacet='chrom', lstart='POS', lstop='POS', rstart='start', rstop='stop', proximity=1)
tbl_variants_with_regions.head()

# <codecell>

region = tbl_variants_with_regions.cut('region').torecarray(dtype={'region': 'a25'}).region

# <codecell>

region

# <codecell>

set(region)

# <codecell>

figure(figsize=(16, 4))
bins = linspace(0, 600000, 150)

core_coding = I.DP[(region == 'Core') & I.CODING]
f, _, _ = hist(core_coding, bins=bins, histtype='step', linewidth=2, label='core coding')
axvline(percentile(core_coding, 15), color='k', linestyle=':')
axvline(percentile(core_coding, 85), color='k', linestyle=':')

hypervariable_coding = I.DP[((region == 'SubtelomericHypervariable') | (region == 'InternalHypervariable')) & I.CODING]
hist(hypervariable_coding, bins=linspace(0, 600000, 150), histtype='step', linewidth=2, label='hypervariable coding')

xlabel('DP')
ylabel('frequency')
title('all SNPs (%s)' % chrom)
legend()

# <codecell>

figure(figsize=(16, 4))
hist(I.DP[(region == 'Core') & I.CODING], bins=linspace(0, 600000, 150), histtype='step', linewidth=2, label='core coding')
hist(I.DP[(region == 'SubtelomericHypervariable') & I.CODING], bins=linspace(0, 600000, 150), histtype='step', linewidth=2, label='hypervariable coding')
hist(I.DP[(region == 'Core') & ~I.CODING], bins=linspace(0, 600000, 150), histtype='step', linewidth=2, label='core non-coding')
hist(I.DP[(region == 'SubtelomericHypervariable') & ~I.CODING], bins=linspace(0, 600000, 150), histtype='step', linewidth=2, label='hypervariable non-coding')
xlabel('DP')
ylabel('frequency')
title('all SNPs (%s)' % chrom)
legend()

# <codecell>

colors = {'Core': 'w', 'SubtelomericHypervariable': 'r', 'InternalHypervariable': 'r', 'SubtelomericRepeat': 'orange', 'Centromere': 'k'}
figure(figsize=(16, 4))

subplot2grid((4, 1), (0, 0), rowspan=3)
window = 10000
bins = arange(0, 650000, window)
X = bins[:-1] + window/2
Y, _ = histogram(V.POS, bins=bins)
plot(X, Y*1./window, linewidth=2, label='all SNPs')
Y, _ = histogram(V.POS[I.CODING], bins=bins)
plot(X, Y*1./window, linewidth=2, label='coding SNPs')
Y, _ = histogram(V.POS[~I.CODING], bins=bins)
plot(X, Y*1./window, linewidth=2, label='non-coding SNPs')

title('SNP density')
xlim(0, max(V.POS))
xticks([])
ylabel('SNP density')
legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

subplot2grid((4, 1), (3, 0))
xranges = [(r.start, r.stop-r.start) for r in tbl_regions.eq('chrom', chrom).records()]
broken_barh(xranges, (0, 1), color=[colors[r.region] for r in tbl_regions.eq('chrom', chrom).records()])
yticks([])
xlim(0, max(V.POS))
xlabel('%s position (bp)' % chrom)

# <codecell>

import scipy.stats as stats

colors = {'Core': 'w', 'SubtelomericHypervariable': 'r', 'InternalHypervariable': 'r', 'SubtelomericRepeat': 'orange', 'Centromere': 'k'}
figure(figsize=(16, 4))

subplot2grid((4, 1), (0, 0), rowspan=3)
window = 10000
bins = arange(0, 650000, window)
X = bins[:-1] + window/2
pos = V.POS[I.CODING]
values = I.DP[I.CODING]

# outer percentiles
Y1, _, _ = stats.binned_statistic(pos, values, statistic=lambda a: percentile(a, 5), bins=bins)
Y2, _, _ = stats.binned_statistic(pos, values, statistic=lambda a: percentile(a, 95), bins=bins)
check = (Y2 > 0) & (Y2 > Y1)
fill_between(X[check], Y1[check], Y2[check], color='b', alpha=.2)
# quartiles
Y1, _, _ = stats.binned_statistic(pos, values, statistic=lambda a: percentile(a, 25), bins=bins)
Y2, _, _ = stats.binned_statistic(pos, values, statistic=lambda a: percentile(a, 75), bins=bins)
check = (Y2 > 0) & (Y2 > Y1)
fill_between(X[check], Y1[check], Y2[check], color='b', alpha=.4)
# median
M, _, _ = stats.binned_statistic(pos, values, statistic='median', bins=bins)
plot(X, M, color='b', linewidth=1)
        
title('depth at coding SNPs')
xlim(0, max(V.POS))
xticks([])
ylabel('DP')
legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

subplot2grid((4, 1), (3, 0))
xranges = [(r.start, r.stop-r.start) for r in tbl_regions.eq('chrom', chrom).records()]
broken_barh(xranges, (0, 1), color=[colors[r.region] for r in tbl_regions.eq('chrom', chrom).records()])
yticks([])
xlim(0, max(V.POS))
xlabel('%s position (bp)' % chrom)

# <codecell>


