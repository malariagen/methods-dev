# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=2>

# Determine duplicate, clonal and 3D7 samples

# <codecell>

from petlx.xlsx import fromxlsx
import petlx.array
from petl.interactive import etl
from petlx.interval import intervaljoin
from petl import totext, tocsv
from scipy import stats
from scipy.spatial.distance import pdist, cdist, squareform
import vcf
import vcfnp
import numpy as np
import sys
import os
import gc
from collections import OrderedDict

# <codecell>

def tabulate(x):
    u, indices = np.unique(x, return_inverse=True)
    return dict(zip(u, np.bincount(indices)))

# <codecell>

PLOT_DIR = '/data/plasmodium/pfalciparum/methods_dev/Pf/3_0/20140114_pf_community_3_0_coverage_filter_v2'
!mkdir -p {PLOT_DIR}

# <codecell>

metadata_2_2_fn = '/data/plasmodium/pfalciparum/3_0/methods-dev/meta/metadata2.2.xlsx'
vcf_gz_fn = '/data/mirror/nfs/team112_internal/production_files/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.gz'
vcf_bgz_fn = '/data/plasmodium/pfalciparum/methods_dev/Pf/3_0/merged_hetuniq_newbiallelic_20130809.vcf.bgz'
evaluation_samples_fmt = '/data/plasmodium/pfalciparum/3_0/methods-dev/meta/evaluation_samples_%d.txt'

numberOfClonalSamplesToUse = 100

# <codecell>

tbl_regions = (etl
    .fromtsv('/data/plasmodium/pfalciparum/methods_dev/Pf/misc/regions-20130225.onebased.txt')
    .pushheader(['chrom', 'start', 'stop', 'region'])
    .convertnumbers()
)
tbl_regions.head()

# <codecell>

metadata_2_2 = etl.fromxlsx(metadata_2_2_fn, 'Sheet1')

# <codecell>

requiredColumns = {'Sample': 'a25', 'Study': 'a25', 'Study Code': 'a25', 'IsDuplicate': 'b', 'Fws': 'f4', 'Notes': 'a50'}

# <codecell>

sampleMetadata = metadata_2_2.cut(requiredColumns.keys()).convert('Fws', float).toarray(dtype=requiredColumns)
Fws = ma.array(sampleMetadata['Fws'], mask=np.isnan(sampleMetadata['Fws']))
print(max(Fws[~Fws.mask]))
print(min(Fws[~Fws.mask]))

# <codecell>

FwsThresholdForClonal = 0.987
FwsThresholdForNonClonal = 0.95

# <codecell>

figure(figsize=(16, 4))
hist(Fws, bins=linspace(0, 1.0, 101), color='w')
xlabel('Fws')
ylabel('Frequency (number of samples)')
axvline(FwsThresholdForClonal, color='red')
axvline(FwsThresholdForNonClonal, color='blue')
savefig(PLOT_DIR + '/Fws_distribution_full.pdf', bbox_inches='tight')

# <codecell>

figure(figsize=(16, 4))
hist(Fws, bins=linspace(0.90, 1.0, 201), color='w')
xlabel('Fws')
ylabel('Frequency (number of samples)')
axvline(FwsThresholdForClonal, color='red')
axvline(FwsThresholdForNonClonal, color='blue')
savefig(PLOT_DIR + '/Fws_distribution_highFws.pdf', bbox_inches='tight')

# <codecell>

samplesInVcf = np.array(vcf.Reader(open(vcf_gz_fn, 'rb')).samples)
len(samplesInVcf)

# <headingcell level=4>

# Identify duplicates in different ways

# <codecell>

duplicates_2_2_1 = sampleMetadata['Sample'][(np.char.find(sampleMetadata['Notes'], 'Duplicate of') > -1) & (np.char.find(sampleMetadata['Notes'], '?') == -1)]
duplicates_2_2_2 = np.char.replace(sampleMetadata['Notes'][(np.char.find(sampleMetadata['Notes'], 'Duplicate of') > -1) & (np.char.find(sampleMetadata['Notes'], '?') == -1)], 'Duplicate of ', '')
duplicates_x_1 = samplesInVcf[np.char.endswith(samplesInVcf, '-Cx')]
duplicates_x_2 = np.char.replace(duplicates_x_1, '-Cx', '-C')
duplicates_W_1 = samplesInVcf[np.char.endswith(samplesInVcf, '-CW')]
duplicates_W_2 = np.char.replace(duplicates_W_1, '-CW', '-C')

# <codecell>

np.in1d(duplicates_W_2, samplesInVcf)

# <codecell>

np.random.seed(12345)
samplesForAnalysis = dict()
samplesForAnalysis['3D7'] = sampleMetadata['Sample'][np.char.find(sampleMetadata['Study Code'], '3D7') > -1]
samplesForAnalysis['3D7'] = samplesForAnalysis['3D7'][np.in1d(samplesForAnalysis['3D7'], samplesInVcf)]
samplesForAnalysis['clonal'] = np.random.choice(sampleMetadata['Sample'][sampleMetadata['Fws'] > FwsThresholdForClonal], numberOfClonalSamplesToUse)
samplesForAnalysis['clonal'] = samplesForAnalysis['clonal'][np.in1d(samplesForAnalysis['clonal'], samplesInVcf)]
samplesForAnalysis['MoI'] = np.random.choice(sampleMetadata['Sample'][(sampleMetadata['Fws'] > 0.0) & (sampleMetadata['Fws'] < FwsThresholdForNonClonal)], numberOfClonalSamplesToUse)
samplesForAnalysis['MoI'] = samplesForAnalysis['MoI'][np.in1d(samplesForAnalysis['MoI'], samplesInVcf)]
samplesForAnalysis['crosses'] = sampleMetadata['Sample'][np.char.find(sampleMetadata['Study'], 'PFproj') > -1]
samplesForAnalysis['crosses'] = samplesForAnalysis['crosses'][np.in1d(samplesForAnalysis['crosses'], samplesInVcf)]
samplesForAnalysis['duplicates_2_2_1'] = duplicates_2_2_1[(np.in1d(duplicates_2_2_1, samplesInVcf)) & (np.in1d(duplicates_2_2_2, samplesInVcf))]
samplesForAnalysis['duplicates_2_2_2'] = duplicates_2_2_2[(np.in1d(duplicates_2_2_1, samplesInVcf)) & (np.in1d(duplicates_2_2_2, samplesInVcf))]
samplesForAnalysis['duplicates_x_1'] = duplicates_x_1[(np.in1d(duplicates_x_1, samplesInVcf)) & (np.in1d(duplicates_x_2, samplesInVcf))]
samplesForAnalysis['duplicates_x_2'] = duplicates_x_2[(np.in1d(duplicates_x_1, samplesInVcf)) & (np.in1d(duplicates_x_2, samplesInVcf))]
samplesForAnalysis['duplicates_W_1'] = duplicates_W_1[(np.in1d(duplicates_W_1, samplesInVcf)) & (np.in1d(duplicates_W_2, samplesInVcf))]
samplesForAnalysis['duplicates_W_2'] = duplicates_W_2[(np.in1d(duplicates_W_1, samplesInVcf)) & (np.in1d(duplicates_W_2, samplesInVcf))]

# <codecell>

for sampleSetName in samplesForAnalysis:
    print(sampleSetName, len(samplesForAnalysis[sampleSetName]))

# <codecell>

allSamplesForAnalysis = np.empty(0, dtype='a15')
for sampleSetName in samplesForAnalysis:
    allSamplesForAnalysis = np.append(allSamplesForAnalysis, samplesForAnalysis[sampleSetName])
allUniqueSamplesForAnalysis = np.unique(allSamplesForAnalysis)
allUniqueSamplesOrderedAsInVcf = np.array(samplesInVcf)[np.in1d(np.array(samplesInVcf), allSamplesForAnalysis)]

# <codecell>

# sanity checks on numbers of samples
print(len(allSamplesForAnalysis))
print(len(np.unique(allSamplesForAnalysis)))
print(len(allUniqueSamplesForAnalysis))
print(len(allUniqueSamplesOrderedAsInVcf))

# <codecell>

sampleIndexes = dict()
for sampleSetName in samplesForAnalysis:
    sampleIndexes[sampleSetName] = np.array(map(lambda x: np.where(allUniqueSamplesOrderedAsInVcf==x)[0][0], samplesForAnalysis[sampleSetName]))

# <codecell>

for sampleIndex in sampleIndexes:
    print(sampleIndex, len(sampleIndexes[sampleIndex]))

# <codecell>

np.savetxt(evaluation_samples_fmt % len(allUniqueSamplesOrderedAsInVcf), allUniqueSamplesOrderedAsInVcf, fmt='%s')

# <codecell>


# <headingcell level=2>

# Build arrays

# <codecell>

# ./scripts/cluster/setups/pf_3_0_build_numpy_arrays_663_samples.sh

# <codecell>

def load_arrays(vcf_fn, region):
    variants_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.variants.npy'
    info_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.info.npy'
    calldata_array_fn = vcf_fn + '.' + region.replace('-:', '_') + '.calldata_663samples.npy'
    return np.load(variants_array_fn).view(np.recarray), np.load(info_array_fn).view(np.recarray), vcfnp.view2d(np.load(calldata_array_fn))

# <headingcell level=2>

# Combine all chromosomes

# <codecell>

V_fn = PLOT_DIR + '/V.npy'
I_fn = PLOT_DIR + '/I.npy'
C_fn = PLOT_DIR + '/C.npy'

if not os.path.exists(V_fn) or not os.path.exists(I_fn) or not os.path.exists(C_fn):
    Vl=range(14)
    Il=range(14)
    Cl=range(14)
    chromosomes = ['Pf3D7_%02d_v3' % i for i in range(1, 15)]
    for i, chromosome in enumerate(chromosomes):
        print i, chromosome
        Vl[i], Il[i], Cl[i] = load_arrays(vcf_bgz_fn, region=chromosome)
    V = np.concatenate(Vl)
    del Vl
    gc.collect()
    I = np.concatenate(Il)
    del Il
    gc.collect()
    C = np.concatenate(Cl)
    del Cl
    gc.collect()
    np.save(V_fn, V)
    np.save(I_fn, I)
    np.save(C_fn, C)
else:
    V = np.load(V_fn)
    I = np.load(I_fn)
    C = np.load(C_fn)

# <codecell>

print(shape(I))
shape(C['AD'])

# <codecell>

tbl_variants = etl.fromarray(V).cut('CHROM', 'POS')

region_fn = PLOT_DIR + '/region.npy'
gc300_fn = PLOT_DIR + '/gc300.npy'
preSoMcoverage_fn = PLOT_DIR + '/preSoMcoverage.npy'

if not os.path.exists(region_fn):
    tbl_variants_with_regions = tbl_variants.intervalleftjoin(tbl_regions, lfacet='CHROM', rfacet='chrom', lstart='POS', lstop='POS', rstart='start', rstop='stop', proximity=1)
    region = tbl_variants_with_regions.cut('region').torecarray(dtype={'region': 'a25'}).region
    np.save(region_fn, region)
else:
    region = np.load(region_fn)    

if not os.path.exists(gc300_fn):
    tbl_gc = (etl
        .fromtsv('/data/plasmodium/pfalciparum/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.gc.300.txt.gz')
        .pushheader(['chrom', 'pos', 'gc300'])
        .convertnumbers()
    )
    tbl_variants_with_gc = (
        tbl_variants
        .rename({'CHROM': 'chrom', 'POS': 'pos'})
        .leftjoin(tbl_gc, key=('chrom', 'pos'), presorted=True, missing=-1)
    )
    gc300 = tbl_variants_with_gc.cut('gc300').torecarray(dtype={'gc300': 'f4'}).gc300
    np.save(gc300_fn, gc300)
else:
    gc300 = np.load(gc300_fn)    

if not os.path.exists(preSoMcoverage_fn):
    tbl_preSoMcoverage = (etl
        .fromtsv('/data/mirror/lustre/scratch108/malaria/pf3_depth/joined.depth')
        .pushheader(['chrom', 'pos', 'preSoMcoverage'])
        .convertnumbers()
    )
    tbl_variants_with_preSoMcoverage = (
        tbl_variants
        .rename({'CHROM': 'chrom', 'POS': 'pos'})
        .leftjoin(tbl_preSoMcoverage, key=('chrom', 'pos'), presorted=True, missing=-1)
    )
    preSoMcoverage = tbl_variants_with_preSoMcoverage.cut('preSoMcoverage').torecarray(dtype={'preSoMcoverage': 'i4'}).preSoMcoverage
    np.save(preSoMcoverage_fn, preSoMcoverage)
else:
    preSoMcoverage = np.load(preSoMcoverage_fn)    

# <headingcell level=2>

# Calls-based analyses

# <codecell>

print(shape(C))
print(shape(C)[0] * shape(C)[1])

# <codecell>

GT_fn = PLOT_DIR + '/GT.npy'
if not os.path.exists(GT_fn):
    GT = np.array((
        (((C['AD'][:,:,0] >= 2) & (C['AD'][:,:,1] >= 2) & (C['AD'][:,:,0] + C['AD'][:,:,1] >= 5)) * 2) +   # het calls
        (((C['AD'][:,:,0] >= 5) & (C['AD'][:,:,1] < 2)) * 1) +                                             # hom ref calls
        (((C['AD'][:,:,0] < 2) & (C['AD'][:,:,1] >= 5)) * 3) +                                             # hom alt calls
        (((C['AD'][:,:,0] <= 2) & (C['AD'][:,:,1] <= 2)) * 0)                                              # "missing" calls
        ), dtype=int8
    )
    np.save(GT_fn, GT)
else:
    GT = np.load(GT_fn)

# <codecell>

GT.dtype

# <headingcell level=4>

# Compute number of homozygote call discordances between all pairwise comparisons of samples

# <codecell>

SNPsets = {'all': np.ones(len(V), dtype=np.bool), 'PASS': V['FILTER']['PASS']}

# <codecell>

cleanCallDiscordanceMatrices_fn = PLOT_DIR + '/cleanCallDiscordanceMatrices.npy'
cleanCallsMatrices_fn = PLOT_DIR + '/cleanCallsMatrices.npy'
cleanCallDiscordanceProportionMatrices_fn = PLOT_DIR + '/cleanCallDiscordanceProportionMatrices.npy'
if not os.path.exists(cleanCallDiscordanceMatrices_fn) or not os.path.exists(cleanCallsMatrices_fn) or not os.path.exists(cleanCallDiscordanceProportionMatrices_fn):
    cleanCallDiscordanceMatrices = dict()
    cleanCallsMatrices = dict()
    cleanCallDiscordanceProportionMatrices = dict()
    for SNPset in SNPsets:
        cleanCallDiscordanceMatrices[SNPset] = dict()
        cleanCallsMatrices[SNPset] = dict()
        cleanCallDiscordanceProportionMatrices[SNPset] = dict()
        for sampleIndex in sampleIndexes:
            print SNPset, sampleIndex, "creating call discordance matrix...",
            cleanCallDiscordanceMatrices[SNPset][sampleIndex] = squareform(pdist(GT[SNPsets[SNPset], :][:, sampleIndexes[sampleIndex]].transpose(), lambda u, v: count_nonzero(((u == 1) & (v == 3)) | ((u == 3) & (v == 1)))))
            print "creating calls matrix...",
            cleanCallsMatrices[SNPset][sampleIndex] = squareform(pdist(GT[SNPsets[SNPset], :][:, sampleIndexes[sampleIndex]].transpose(), lambda u, v: count_nonzero(((u == 1) | (v == 3)) & ((u == 3) | (v == 1)))))
            print "creating call discordance proportion matrix..."
            cleanCallDiscordanceProportionMatrices[SNPset][sampleIndex] = numpy.ma.array(np.where(cleanCallsMatrices[SNPset][sampleIndex]==0, np.nan, cleanCallDiscordanceMatrices[SNPset][sampleIndex] * 1.0 / cleanCallsMatrices[SNPset][sampleIndex]), mask=(cleanCallsMatrices[SNPset][sampleIndex]==0))
    np.save(cleanCallDiscordanceMatrices_fn, cleanCallDiscordanceMatrices)
    np.save(cleanCallsMatrices_fn, cleanCallsMatrices)
    np.save(cleanCallDiscordanceProportionMatrices_fn, cleanCallDiscordanceProportionMatrices)
else:
    cleanCallDiscordanceMatrices = np.load(cleanCallDiscordanceMatrices_fn).item()
    cleanCallsMatrices = np.load(cleanCallsMatrices_fn).item()
    cleanCallDiscordanceProportionMatrices = np.load(cleanCallDiscordanceProportionMatrices_fn).item()

# <codecell>

def cleanCalls(gt, samples1, samples2):
    calls = (
        ((gt[:,samples1] == 1) | (gt[:,samples1] == 3)) &
        ((gt[:,samples2] == 1) | (gt[:,samples2] == 3))
    )
    return calls

def cleanCallDiscordances(gt, samples1, samples2):
    discordantCalls = (
        ((gt[:,samples1] == 1) & (gt[:,samples2] == 3)) |
        ((gt[:,samples1] == 3) & (gt[:,samples2] == 1))
    )
    return discordantCalls

def calledCalls(gt, samples1, samples2):
    calls = (
        ((gt[:,samples1] == 1) | (gt[:,samples1] == 2) | (gt[:,samples1] == 3)) &
        ((gt[:,samples2] == 1) | (gt[:,samples2] == 2) | (gt[:,samples2] == 3))
    )
    return calls

def calledCallDiscordances(gt, samples1, samples2):
    discordantCalls = (
        ((gt[:,samples1] == 1) & ((gt[:,samples2] == 2) | (gt[:,samples2] == 3))) |
        ((gt[:,samples1] == 2) & ((gt[:,samples2] == 1) | (gt[:,samples2] == 3))) |
        ((gt[:,samples1] == 3) & ((gt[:,samples2] == 1) | (gt[:,samples2] == 2)))
    )
    return discordantCalls

def allCallDiscordances(gt, samples1, samples2):
    discordantCalls = gt[:,samples1] != gt[:,samples2]
    return discordantCalls

# <codecell>

duplicateCallAnalyses_fn = PLOT_DIR + '/duplicateCallAnalyses.npy'
duplicateSets = ['2_2', 'x', 'W']
if not os.path.exists(duplicateCallAnalyses_fn):
    duplicateCallAnalyses = dict()
    for SNPset in SNPsets:
        duplicateCallAnalyses[SNPset] = dict()
        for duplicateSet in duplicateSets:
            print SNPset, duplicateSet,
            duplicateCallAnalyses[SNPset][duplicateSet] = dict()
            samples1 = sampleIndexes['duplicates_' + duplicateSet + '_1']
            samples2 = sampleIndexes['duplicates_' + duplicateSet + '_2']
            print 'duplicateCleanCalls',
            duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCalls'] = np.sum(cleanCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCallsPerSamplePair'] = np.sum(cleanCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantCleanCalls',
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCalls'] = np.sum(cleanCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCallsPerSamplePair'] = np.sum(cleanCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            duplicateCallAnalyses[SNPset][duplicateSet]['proportionDiscordantCleanCallsPerSamplePair'] = duplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCallsPerSamplePair']
            print 'duplicateCalledCalls',
            duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCalledCalls'] = np.sum(calledCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCalledCallsPerSamplePair'] = np.sum(calledCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantCalledCalls',
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantCalledCalls'] = np.sum(calledCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantCalledCallsPerSamplePair'] = np.sum(calledCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantAllCalls'
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantAllCalls'] = np.sum(allCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            duplicateCallAnalyses[SNPset][duplicateSet]['discordantAllCallsPerSamplePair'] = np.sum(allCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
    np.save(duplicateCallAnalyses_fn, duplicateCallAnalyses)
else:
    duplicateCallAnalyses = np.load(duplicateCallAnalyses_fn).item() #saved as numpy array, so need item() to get back dict

# <codecell>

duplicateCallAnalyses['all']['2_2']['discordantCleanCallsPerSamplePair']

# <codecell>

trueDuplicateThreshold = 2e-4
sampleSetsToPlot = ['3D7', 'clonal', 'MoI', 'crosses']
numberOfPlots = len(sampleSetsToPlot) + len(duplicateSets)
fig = figure(figsize=(8.27, 11.69))
ax_full = fig.add_subplot(1, 1, 1)    # The big subplot

# Turn off axis lines and ticks of the big subplot
ax_full.spines['top'].set_color('none')
ax_full.spines['bottom'].set_color('none')
ax_full.spines['left'].set_color('none')
ax_full.spines['right'].set_color('none')
ax_full.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

plotRowNumber = 0
for sampleIndex in sampleSetsToPlot:
    plotRowNumber += 1
    ax = fig.add_subplot(numberOfPlots, 1, plotRowNumber)
    hist(np.reshape(cleanCallDiscordanceProportionMatrices['all'][sampleIndex], (len(sampleIndexes[sampleIndex])**2, 1)), bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='all SNPs')
    hist(np.reshape(cleanCallDiscordanceProportionMatrices['PASS'][sampleIndex], (len(sampleIndexes[sampleIndex])**2, 1)), bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='PASS SNPs')
    xlim(0, 0.015)
    axvline(trueDuplicateThreshold, color='red')
    ax.set_xticks([])
    title('Pairwise clean call discordances between samples - ' + sampleIndex)
    if(plotRowNumber==1):
        legend()
for duplicateSet in duplicateSets:
    plotRowNumber += 1
    ax = fig.add_subplot(numberOfPlots, 1, plotRowNumber)
    cleanCallDiscordanceProportionsAll = duplicateCallAnalyses['all'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses['all'][duplicateSet]['duplicateCleanCallsPerSamplePair']
    hist(cleanCallDiscordanceProportionsAll, bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='all SNPs')
    cleanCallDiscordanceProportionsPASS = duplicateCallAnalyses['PASS'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses['PASS'][duplicateSet]['duplicateCleanCallsPerSamplePair']
    hist(cleanCallDiscordanceProportionsPASS, bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='PASS SNPs')
    xlim(0, 0.015)
    axvline(trueDuplicateThreshold, color='red')
    if(plotRowNumber < numberOfPlots):
        ax.set_xticks([])
    title('Pairwise clean call discordances between duplicate samples - ' + duplicateSet)

# Set common labels
ax_full.legend()
ax_full.set_xlabel('Proportion of discordant clean (homozygote) calls')
ax_full.set_ylabel('Frequency (number of sample pairs)')

savefig(PLOT_DIR + '/Pairwise_clean_call_discordances.pdf', bbox_inches='tight')

# <codecell>

sampleIndexesTrueDuplicates = dict()
for duplicateSet in duplicateSets:
    sampleIndexesTrueDuplicates['duplicates_' + duplicateSet + '_1'] = sampleIndexes['duplicates_' + duplicateSet + '_1'][(duplicateCallAnalyses['PASS'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses['PASS'][duplicateSet]['duplicateCleanCallsPerSamplePair']) < trueDuplicateThreshold]
    sampleIndexesTrueDuplicates['duplicates_' + duplicateSet + '_2'] = sampleIndexes['duplicates_' + duplicateSet + '_2'][(duplicateCallAnalyses['PASS'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses['PASS'][duplicateSet]['duplicateCleanCallsPerSamplePair']) < trueDuplicateThreshold]
        

# <codecell>

for sampleIndex in sampleIndexesTrueDuplicates:
    print(sampleIndex, len(sampleIndexesTrueDuplicates[sampleIndex]))

# <codecell>

trueDuplicateCallAnalyses_fn = PLOT_DIR + '/trueDuplicateCallAnalyses.npy'
trueDuplicateThreshold = 2e-4
if not os.path.exists(trueDuplicateCallAnalyses_fn):
    trueDuplicateCallAnalyses = dict()
    for SNPset in SNPsets:
        trueDuplicateCallAnalyses[SNPset] = dict()
        for duplicateSet in duplicateSets:
            print SNPset, duplicateSet,
            trueDuplicateCallAnalyses[SNPset][duplicateSet] = dict()
            samples1 = sampleIndexesTrueDuplicates['duplicates_' + duplicateSet + '_1']
            samples2 = sampleIndexesTrueDuplicates['duplicates_' + duplicateSet + '_2']
            print 'duplicateCleanCalls',
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCalls'] = np.sum(cleanCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCallsPerSamplePair'] = np.sum(cleanCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantCleanCalls',
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCalls'] = np.sum(cleanCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCallsPerSamplePair'] = np.sum(cleanCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['proportionDiscordantCleanCallsPerSamplePair'] = duplicateCallAnalyses[SNPset][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / duplicateCallAnalyses[SNPset][duplicateSet]['duplicateCleanCallsPerSamplePair']
            print 'duplicateCalledCalls',
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['duplicateCalledCalls'] = np.sum(calledCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['duplicateCalledCallsPerSamplePair'] = np.sum(calledCalls(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantCalledCalls',
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantCalledCalls'] = np.sum(calledCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantCalledCallsPerSamplePair'] = np.sum(calledCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
            print 'discordantAllCalls'
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantAllCalls'] = np.sum(allCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=1)
            trueDuplicateCallAnalyses[SNPset][duplicateSet]['discordantAllCallsPerSamplePair'] = np.sum(allCallDiscordances(GT[SNPsets[SNPset], :], samples1, samples2), axis=0)
    np.save(trueDuplicateCallAnalyses_fn, trueDuplicateCallAnalyses)
else:
    trueDuplicateCallAnalyses = np.load(trueDuplicateCallAnalyses_fn).item()

# <headingcell level=4>

# Sanity check that correct duplicates have been removed

# <codecell>

trueDuplicateThreshold = 2e-4
for duplicateSet in duplicateSets:
    figure(figsize=(16, 4))
    cleanCallDiscordanceProportionsAll = trueDuplicateCallAnalyses['all'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / trueDuplicateCallAnalyses['all'][duplicateSet]['duplicateCleanCallsPerSamplePair']
    hist(cleanCallDiscordanceProportionsAll, bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='all SNPs')
    cleanCallDiscordanceProportionsPASS = trueDuplicateCallAnalyses['PASS'][duplicateSet]['discordantCleanCallsPerSamplePair'] * 1.0 / trueDuplicateCallAnalyses['PASS'][duplicateSet]['duplicateCleanCallsPerSamplePair']
    hist(cleanCallDiscordanceProportionsPASS, bins=linspace(0.0, 0.015, 300), histtype='step', linewidth=2, label='PASS SNPs')
    xlim(0, 0.015)
    xlabel('Proportion of discordant clean (homozygote) calls')
    ylabel('Frequency (number of sample pairs)')
    axvline(trueDuplicateThreshold, color='red')
    title('Pairwise clean call discordances between duplicate samples - ' + duplicateSet)
    legend()

# <codecell>

missingness_threshold=5
het_threshold=2
AFs_fn = PLOT_DIR + '/AFs.npy'
if not os.path.exists(AFs_fn):
    AFs = dict()
    for sampleSetName in ['clonal', 'MoI']:
        print sampleSetName
        AFs[sampleSetName] = dict()
        ADsThisSampleSet = C['AD'][:, sampleIndexes[sampleSetName], :]
        aaf = (ADsThisSampleSet[:, :, 1]) * 1.0 / ((ADsThisSampleSet[:, :, 0]) + (ADsThisSampleSet[:, :, 1]))
        aaf[ADsThisSampleSet[:, :, 1] < het_threshold] = 0.0
        aaf[ADsThisSampleSet[:, :, 0] < het_threshold] = 1.0
        aaf[(ADsThisSampleSet[:, :, 0] + ADsThisSampleSet[:, :, 1]) < missingness_threshold] = np.nan
        AFs[sampleSetName]['nraf'] = np.apply_along_axis(stats.stats.nanmean, 1, aaf)
        AFs[sampleSetName]['maf'] = np.where(AFs[sampleSetName]['nraf'] <= 0.5, AFs[sampleSetName]['nraf'], 1.0-AFs[sampleSetName]['nraf'])
        AFs[sampleSetName]['nraf'] = ma.array(AFs[sampleSetName]['nraf'], mask=np.isnan(AFs[sampleSetName]['nraf']))
        AFs[sampleSetName]['maf'] = ma.array(AFs[sampleSetName]['maf'], mask=np.isnan(AFs[sampleSetName]['maf']))
    np.save(AFs_fn, AFs)
else:
    AFs = np.load(AFs_fn).item()

# <codecell>

tabulate(np.isnan(AFs['clonal']['nraf']))

# <codecell>

# sanity check
tabulate(AFs['MoI']['maf'] > 0.5)

# <headingcell level=4>

# Sanity check of AFs

# <codecell>

figure(figsize=(16, 4))
hist(AFs['clonal']['nraf'], bins=linspace(0, 1, 100), log=True, histtype='step', linewidth=2, label='clonal NRAF')
hist(AFs['clonal']['maf'], bins=linspace(0, 1, 100), log=True, histtype='step', linewidth=2, label='clonal MAF')
hist(AFs['MoI']['nraf'], bins=linspace(0, 1, 100), log=True, histtype='step', linewidth=2, label='MoI NRAF')
hist(AFs['MoI']['maf'], bins=linspace(0, 1, 100), log=True, histtype='step', linewidth=2, label='MoI MAF')
legend()

# <codecell>

np.sum(GT[:, sampleIndexes['3D7']] == 3, axis=0)

# <codecell>

noCalls3D7_fn = PLOT_DIR + '/noCalls3D7.npy'
if not os.path.exists(noCalls3D7_fn):
    noCalls3D7 = np.sum(GT[:, sampleIndexes['3D7']] == 0, axis=1)
    np.save(noCalls3D7_fn, noCalls3D7)
else:
    noCalls3D7 = np.load(noCalls3D7_fn)

homAlts3D7_fn = PLOT_DIR + '/homAlts3D7.npy'
if not os.path.exists(homAlts3D7_fn):
    homAlts3D7 = np.sum(GT[:, sampleIndexes['3D7']] == 3, axis=1)
    np.save(homAlts3D7_fn, homAlts3D7)
else:
    homAlts3D7 = np.load(homAlts3D7_fn)

hets3D7_fn = PLOT_DIR + '/hets3D7.npy'
if not os.path.exists(hets3D7_fn):
    hets3D7 = np.sum(GT[:, sampleIndexes['3D7']] == 2, axis=1)
    np.save(hets3D7_fn, hets3D7)
else:
    hets3D7 = np.load(hets3D7_fn)

numHetClonal_fn = PLOT_DIR + '/numHetClonal.npy'
if not os.path.exists(numHetClonal_fn):
    numHetClonal = np.sum(GT[:, sampleIndexes['clonal']] == 2, axis=1)
    np.save(numHetClonal_fn, numHetClonal)
else:
    numHetClonal = np.load(numHetClonal_fn)

numHetMoI_fn = PLOT_DIR + '/numHetMoI.npy'
if not os.path.exists(numHetMoI_fn):
    numHetMoI = np.sum(GT[:, sampleIndexes['MoI']] == 2, axis=1)
    np.save(numHetMoI_fn, numHetMoI)
else:
    numHetMoI = np.load(numHetMoI_fn)

# <codecell>

tabulate(numpy.core.defchararray.add(np.where(region == 'Core', 'core', 'non-core'), np.where(V['FILTER']['PASS'], ' PASS', ' non-PASS')))

# <headingcell level=4>

# Create recarray to be used in plots

# <codecell>

def is_transition(v):
    return (v['is_snp'] & (v['num_alleles'] == 2) & 
            ((v['REF'] == 'A') & (v['ALT'] == 'G')) 
             | ((v['REF'] == 'G') & (v['ALT'] == 'A')) 
             | ((v['REF'] == 'C') & (v['ALT'] == 'T')) 
             | ((v['REF'] == 'T') & (v['ALT'] == 'C')))


def is_transversion(v):
    return (v['is_snp'] & (v['num_alleles'] == 2) & ~is_transition(v))

# <codecell>

variants_fn = PLOT_DIR + '/variants.npy'
force=True
if force==True or not os.path.exists(variants_fn):
    variants = np.core.records.fromarrays([
        V['CHROM'],
        V['POS'],
        I['NS'],
        I['UQ'],
        I['CODING'],
        I['DP'],
        preSoMcoverage,
        preSoMcoverage-I['DP'],
        (preSoMcoverage-I['DP'])*1.0/preSoMcoverage,
        gc300,
        is_transition(V),
        is_transversion(V),
        AFs['clonal']['nraf'],
        AFs['clonal']['maf'],
        AFs['MoI']['nraf'],
        AFs['MoI']['maf'],
        noCalls3D7,
        homAlts3D7,
        hets3D7,
        numHetClonal,
        numHetMoI,
        trueDuplicateCallAnalyses['all']['2_2']['duplicateCleanCalls'],
        trueDuplicateCallAnalyses['all']['x']['duplicateCleanCalls'],
        trueDuplicateCallAnalyses['all']['W']['duplicateCleanCalls'],
        where(trueDuplicateCallAnalyses['all']['2_2']['duplicateCleanCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['2_2']['discordantCleanCalls']/trueDuplicateCallAnalyses['all']['2_2']['duplicateCleanCalls']),
        where(trueDuplicateCallAnalyses['all']['x']['duplicateCleanCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['x']['discordantCleanCalls']/trueDuplicateCallAnalyses['all']['x']['duplicateCleanCalls']),
        where(trueDuplicateCallAnalyses['all']['W']['duplicateCleanCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['W']['discordantCleanCalls']/trueDuplicateCallAnalyses['all']['W']['duplicateCleanCalls']),
        trueDuplicateCallAnalyses['all']['2_2']['duplicateCalledCalls'],
        trueDuplicateCallAnalyses['all']['x']['duplicateCalledCalls'],
        trueDuplicateCallAnalyses['all']['W']['duplicateCalledCalls'],
        where(trueDuplicateCallAnalyses['all']['2_2']['duplicateCalledCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['2_2']['discordantCalledCalls']/trueDuplicateCallAnalyses['all']['2_2']['duplicateCalledCalls']),
        where(trueDuplicateCallAnalyses['all']['x']['duplicateCalledCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['x']['discordantCalledCalls']/trueDuplicateCallAnalyses['all']['x']['duplicateCalledCalls']),
        where(trueDuplicateCallAnalyses['all']['W']['duplicateCalledCalls']==0, 0.0, 1.0*trueDuplicateCallAnalyses['all']['W']['discordantCalledCalls']/trueDuplicateCallAnalyses['all']['W']['duplicateCalledCalls']),
        1.0 * trueDuplicateCallAnalyses['all']['2_2']['discordantAllCalls'] / len(sampleIndexesTrueDuplicates['duplicates_2_2_1']),
        1.0 * trueDuplicateCallAnalyses['all']['x']['discordantAllCalls'] / len(sampleIndexesTrueDuplicates['duplicates_x_1']),
        1.0 * trueDuplicateCallAnalyses['all']['W']['discordantAllCalls'] / len(sampleIndexesTrueDuplicates['duplicates_W_1']),
        where(
            (
                trueDuplicateCallAnalyses['all']['2_2']['duplicateCleanCalls'] +
                trueDuplicateCallAnalyses['all']['x']['duplicateCleanCalls'] +
                trueDuplicateCallAnalyses['all']['W']['duplicateCleanCalls']
            )==0,
            0.0,
            1.0*(
                trueDuplicateCallAnalyses['all']['2_2']['discordantCleanCalls'] +
                trueDuplicateCallAnalyses['all']['x']['discordantCleanCalls'] +
                trueDuplicateCallAnalyses['all']['W']['discordantCleanCalls']
            )/(
                trueDuplicateCallAnalyses['all']['2_2']['duplicateCleanCalls'] +
                trueDuplicateCallAnalyses['all']['x']['duplicateCleanCalls'] +
                trueDuplicateCallAnalyses['all']['W']['duplicateCleanCalls']
            )
        ),
        where(
            (
                trueDuplicateCallAnalyses['all']['2_2']['duplicateCalledCalls'] +
                trueDuplicateCallAnalyses['all']['x']['duplicateCalledCalls'] +
                trueDuplicateCallAnalyses['all']['W']['duplicateCalledCalls']
            )==0,
            0.0,
            1.0*(
                trueDuplicateCallAnalyses['all']['2_2']['discordantCalledCalls'] +
                trueDuplicateCallAnalyses['all']['x']['discordantCalledCalls'] +
                trueDuplicateCallAnalyses['all']['W']['discordantCalledCalls']
            )/(
                trueDuplicateCallAnalyses['all']['2_2']['duplicateCalledCalls'] +
                trueDuplicateCallAnalyses['all']['x']['duplicateCalledCalls'] +
                trueDuplicateCallAnalyses['all']['W']['duplicateCalledCalls']
            )
        ),
        np.where(region == 'Core', 1, 0),
        np.where(V['FILTER']['PASS'], 1, 0),
        np.where((~V['FILTER']['Biallelic'])&(~V['FILTER']['CodingType'])&(~V['FILTER']['HetUniq'])&(~V['FILTER']['MinAlt'])&(~V['FILTER']['MonoAllelic'])&(~V['FILTER']['NoAltAllele'])&(~V['FILTER']['Region'])&(~V['FILTER']['triallelic']), 1, 0),
        np.where(I['CODING'], 1, 0)
    ]
    ,names='CHROM,POS,NS,UQ,CODING,DP,preSoMcoverage,DPdropDueToSoM,propDPdropDueToSoM,gc300,is_transition,is_transversion,clonalNRAF,clonalMAF,MoINRAF,MoIMAF,' +
        'noCalls3D7,homAlts3D7,hets3D7,numHetClonal,numHetMoI,' + 
        'duplicateCleanCalls_2_2,duplicateCleanCalls_x,duplicateCleanCalls_W,' + 
        'proportionDiscordantCleanCalls_2_2,proportionDiscordantCleanCalls_x,proportionDiscordantCleanCalls_W,' + 
        'duplicateCalledCalls_2_2,duplicateCalledCalls_x,duplicateCalledCalls_W,' + 
        'proportionDiscordantCalledCalls_2_2,proportionDiscordantCalledCalls_x,proportionDiscordantCalledCalls_W,' +
        'proportionDiscordantAllCalls_2_2,proportionDiscordantAllCalls_x,proportionDiscordantAllCalls_W,' + 
        'proportionDiscordantCleanCalls_all,proportionDiscordantCalledCalls_all,' +
        'isInCore,PASS,PASSexceptCoverage,CodingNumeric'
    )
    np.save(variants_fn, variants)
else:
    variants = np.load(variants_fn)
force=False

# <codecell>

tabulate(variants['CHROM'])

# <codecell>

attributeDescriptions = OrderedDict([
    ('clonalNRAF',                          {'description': 'Mean non-ref allele frequency of 100 clonal samples',            'statistic': 'mean',  'colour': '#FF0000'}),
    ('clonalMAF',                           {'description': 'Mean minor allele frequency of 100 clonal samples',              'statistic': 'mean',  'colour': '#FF5500'}),
    ('MoINRAF',                             {'description': 'Mean non-ref allele frequency of 100 MoI samples',               'statistic': 'mean',  'colour': '#FFAA00'}),
    ('MoIMAF',                              {'description': 'Mean minor allele frequency of 100 MoI samples',                 'statistic': 'mean',  'colour': '#FFFF00'}),
    ('noCalls3D7',                          {'description': 'Number of missing calls in 4 samples of 3D7',                    'statistic': 'sum',   'colour': '#FF99FF'}),
    ('homAlts3D7',                          {'description': 'Number of hom alt calls in 4 samples of 3D7',                    'statistic': 'sum',   'colour': '#FF99BB'}),
    ('hets3D7',                             {'description': 'Number of het calls in 4 samples of 3D7',                        'statistic': 'sum',   'colour': '#FF9966'}),
    ('numHetClonal',                        {'description': 'Mean heterozygosity of 100 clonal samples',                      'statistic': 'mean',  'colour': '#CC00FF'}),
    ('numHetMoI',                           {'description': 'Mean heterozygosity of 100 MoI samples',                         'statistic': 'mean',  'colour': '#CC00BB'}),
    ('duplicateCleanCalls_2_2',             {'description': 'Mean number of duplicate clean calls in 2.2 duplicates',         'statistic': 'mean',  'colour': '#FF22FF'}),
    ('duplicateCleanCalls_x',               {'description': 'Mean number of duplicate clean calls in \'x\' duplicates',       'statistic': 'mean',  'colour': '#FF22BB'}),
    ('duplicateCleanCalls_W',               {'description': 'Mean number of duplicate clean calls in \'W\' duplicates',       'statistic': 'mean',  'colour': '#FF2266'}),
    ('proportionDiscordantCleanCalls_2_2',  {'description': 'Mean proportion of discordant clean calls in 2.2 duplicates',    'statistic': 'mean',  'colour': '#FF66FF'}),
    ('proportionDiscordantCleanCalls_x',    {'description': 'Mean proportion of discordant clean calls in \'x\' duplicates',  'statistic': 'mean',  'colour': '#FF66BB'}),
    ('proportionDiscordantCleanCalls_W',    {'description': 'Mean proportion of discordant clean calls in \'W\' duplicates',  'statistic': 'mean',  'colour': '#FF6666'}),
    ('proportionDiscordantCalledCalls_2_2', {'description': 'Mean proportion of discordant called calls in 2.2 duplicates',   'statistic': 'mean',  'colour': '#FFCCFF'}),
    ('proportionDiscordantCalledCalls_x',   {'description': 'Mean proportion of discordant called calls in \'x\' duplicates', 'statistic': 'mean',  'colour': '#FFCCBB'}),
    ('proportionDiscordantCalledCalls_W',   {'description': 'Mean proportion of discordant called calls in \'W\' duplicates', 'statistic': 'mean',  'colour': '#FFCC66'}),
    ('isInCore',                            {'description': 'Proportion of variants in core genome',                          'statistic': 'mean',  'colour': '#CCFFFF'}),
    ('PASS',                                {'description': 'Proportion of PASS variants',                                    'statistic': 'mean',  'colour': '#CCFFAA'}),
    ('PASSexceptCoverage',                  {'description': 'Proportion of PASS variants excluding coverage filter',          'statistic': 'mean',  'colour': '#CCFF55'}),
    ('CodingNumeric',                       {'description': 'Proportion of coding variants',                                  'statistic': 'mean',  'colour': '#CCFF00'}),
    ('NS',                                  {'description': 'Mean number of samples called',                                  'statistic': 'mean',  'colour': '#66FFFF'}),
    ('UQ',                                  {'description': 'Mean uniqueness of variants',                                    'statistic': 'mean',  'colour': '#66FFBB'}),
    ('DP',                                  {'description': 'Mean depth at variants from VCF file',                           'statistic': 'mean',  'colour': '#66FF77'}),
    ('preSoMcoverage',                      {'description': 'Mean coverage from pre-SoM bam files',                           'statistic': 'mean',  'colour': '#66FF33'}),
    ('gc300',                               {'description': 'Mean GC content',                                                'statistic': 'mean',  'colour': '#66FF00'})
])

# <headingcell level=3>

# Create filter decision plots

# <codecell>

binwidth=2000

VCF_DP_core_coding = variants['DP'][(variants['isInCore']==1) & (variants['CodingNumeric']==1)]
VCF_DP_noncore_coding = variants['DP'][(variants['isInCore']==0) & (variants['CodingNumeric']==1)]
PreSoM_core_coding = variants['preSoMcoverage'][(variants['isInCore']==1) & (variants['CodingNumeric']==1)]
PreSoM_noncore_coding = variants['preSoMcoverage'][(variants['isInCore']==0) & (variants['CodingNumeric']==1)]
PreSoM_core_coding_rescaled = PreSoM_core_coding*np.median(VCF_DP_core_coding)/np.median(PreSoM_core_coding)

VCF_DP_core_noncoding = variants['DP'][(variants['isInCore']==1) & (variants['CodingNumeric']==0)]

#lowerBound = np.median(VCF_DP_core_coding) * 0.8
#upperBound = np.median(VCF_DP_core_coding) * 1.2
lowerBoundDP = percentile(variants['DP'][(variants['CodingNumeric']==1)], 15)
upperBoundDP = percentile(variants['DP'][(variants['CodingNumeric']==1)], 85)
lowerBoundPreSoM = percentile(variants['preSoMcoverage'][(variants['CodingNumeric']==1)], 15)
upperBoundPreSoM = percentile(variants['preSoMcoverage'][(variants['CodingNumeric']==1)], 85)
print(lowerBoundDP)
print(upperBoundDP)
print(lowerBoundPreSoM)
print(upperBoundPreSoM)

bins=arange(0, 600001, binwidth)

figure(figsize=(16, 4))
hist(VCF_DP_core_coding, bins=bins, histtype='step', linewidth=2, color='blue', label='VCF DP core coding')
hist(VCF_DP_noncore_coding, bins=bins, histtype='step', linewidth=2, color='green', label='VCF DP non-core coding')
hist(PreSoM_core_coding, bins=bins, histtype='step', linewidth=2, color='red', label='Pre-SoM coverage core coding')
hist(PreSoM_noncore_coding, bins=bins, histtype='step', linewidth=2, color='orange', label='Pre-SoM coverage non-core coding')
hist(PreSoM_core_coding_rescaled, bins=bins, histtype='step', linewidth=2, color='purple', label='Pre-SoM coverage core coding rescaled')
axvline(lowerBoundDP, color='blue')
axvline(upperBoundDP, color='blue')
legend()
xlabel('Coverage')
ylabel('Frequency (number of SNPs)')
savefig(PLOT_DIR + '/coverageThresholdsCoding.pdf', bbox_inches='tight')

# <codecell>

def characterisation_plot(v, xattr, yattrDescriptions=attributeDescriptions, yattrToIgnore=['clonalNRAF', 'MoINRAF', 'numHetMoI', 'PASS', 'PASSexceptCoverage', 'CodingNumeric'],
        verticalLineValues=None, verticalLineColours=None, thresholds=None, bins=200, shouldAddJitter=True, plotFilestem=PLOT_DIR + '/characterisation',
        extendedXlabels=False, distributionToUseForLabels=None, medianToUseForLabels=None, SDtoUseForLabels=None,
        plotHeight=24, plotWidth=24):
    def diff(a, b):
        b = set(b)
        return [aa for aa in a if aa not in b]
    if shouldAddJitter:
        values = 1.0*v[xattr] + np.random.normal(loc=0.0, scale=1e-10, size=size(v[xattr]))
    else:
        values = v[xattr]
    if thresholds==None:
        thresholds = list()
        for i in (np.array(range(bins + 1)) * 100.0)/bins:
            thresholds.append(stats.scoreatpercentile(values, i))
    if verticalLineValues != None:
        if issubclass(v[xattr].dtype.type, numpy.integer):
            verticalLinePos = [stats.percentileofscore(v[xattr], x, kind='strict') * bins / 100.0 for x in verticalLineValues]
        else:
            verticalLinePos = [stats.percentileofscore(values, x) * bins / 100.0 for x in verticalLineValues]
    N = len(thresholds) - 1
    metrics = OrderedDict()
    for attributeDescription in diff(yattrDescriptions.keys(), yattrToIgnore):
        metrics[attributeDescription] = list()
    metrics['TiTv'] = list()
    for i in range(N):
        flt = (values >= thresholds[i]) & (values < thresholds[i+1])
        is_transition = v[flt]['is_transition']
        transitions = count_nonzero(is_transition)
        transversions = count_nonzero(~is_transition)
        titv =  (transitions * 1. / transversions) if transversions > 0 else 0  
        metrics['TiTv'].append(titv)
        for attributeDescription in diff(yattrDescriptions.keys(), yattrToIgnore):
            if yattrDescriptions[attributeDescription]['statistic'] == 'mean':
                metrics[attributeDescription].append(stats.stats.nanmean(v[flt][attributeDescription]))
            if yattrDescriptions[attributeDescription]['statistic'] == 'sum':
                metrics[attributeDescription].append(sum(v[flt][attributeDescription]))
            if yattrDescriptions[attributeDescription]['statistic'] == 'count':
                metrics[attributeDescription].append(count_nonzero(v[flt][attributeDescription]))
    ind = np.arange(N)  # the x locations for the groups
    width = 1       # the width of the bars
    numberOfAttributes = len(diff(yattrDescriptions.keys(), yattrToIgnore))
    fig = plt.figure(figsize=(plotWidth, plotHeight))
    rows, cols = numberOfAttributes+1, 1
    
    for (i, attributeDescription) in enumerate(diff(yattrDescriptions.keys(), yattrToIgnore)):
        ax = fig.add_subplot(rows, cols, i+1)
        ax.bar(ind, metrics[attributeDescription], width, color=yattrDescriptions[attributeDescription]['colour'])
        ax.set_ylabel(yattrDescriptions[attributeDescription]['description'], rotation='horizontal', horizontalalignment='right')
        ax.set_xticks([])
        ax.set_xlim(0, N)
        ax.tick_params(axis='y', which='major', labelsize=6)
        if verticalLineValues != None:
            for i, x in enumerate(verticalLinePos):
                if verticalLineColours == None:
                    ax.axvline(x)
                else:
                    ax.axvline(x, color=verticalLineColours[i])
            #ax.axvline(verticalLinePos)

    ax = fig.add_subplot(rows, cols, numberOfAttributes+1)
    ax.bar(ind, metrics['TiTv'], width, color='b')
    ax.set_xlabel(xattr)
    ax.set_ylabel('Mean Ti/Tv ratio', rotation='horizontal', horizontalalignment='right')
    ax.set_xticks(ind)
    if issubclass(v[xattr].dtype.type, numpy.integer):
        xLabels = np.array(map(str, [int(round(threshold)) for threshold in thresholds]))[arange(bins)]
    else:
        xLabels = np.array(map(str, [round(threshold, 6) for threshold in thresholds]))[arange(bins)]
    if extendedXlabels:
        if distributionToUseForLabels==None:
            distributionToUseForLabels=v[xattr]
        thresholdsAsProportionOfMedian = np.array(thresholds)[arange(bins)] / np.median(distributionToUseForLabels)
        thresholdsAsSDsFromMedian = (np.array(thresholds)[arange(bins)] - np.median(distributionToUseForLabels)) / np.std(distributionToUseForLabels)
        percentiles = np.array([stats.percentileofscore(distributionToUseForLabels, x) for x in np.array(thresholds)[arange(bins)]])
        labelExtenstions = np.core.defchararray.add(np.core.defchararray.add(np.array([", %.2f" % number for number in thresholdsAsProportionOfMedian]), np.array([", %.2f" % number for number in thresholdsAsSDsFromMedian])), np.array([", %.1f" % number for number in percentiles]))
        xLabels = np.core.defchararray.add(xLabels, labelExtenstions)
    ax.set_xticklabels(xLabels)
    xticks(rotation='vertical')
    ax.set_xlim(0, N)
    if verticalLineValues != None:
        for i, x in enumerate(verticalLinePos):
            if verticalLineColours == None:
                ax.axvline(x)
            else:
                ax.axvline(x, color=verticalLineColours[i])
    
    ax.tick_params(axis='both', which='major', labelsize=6)
    #ax.tick_params(axis='both', which='minor', labelsize=6)
    
    # save
    if plotFilestem is not None:
        fig.savefig(plotFilestem + '.' + xattr + '.pdf', bbox_inches='tight')

    return True

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'preSoMcoverage',
    bins=200, extendedXlabels=True,
    plotFilestem=PLOT_DIR + '/AllCoding')

# <codecell>

characterisation_plot(variants[(variants['isInCore']==1) & (variants['CodingNumeric']==1)], 'preSoMcoverage',
    bins=200, extendedXlabels=True,
    plotFilestem=PLOT_DIR + '/CoreCoding')

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'DP',
    bins=200, extendedXlabels=True,
    verticalLineValues= [percentile(variants['DP'][(variants['CodingNumeric']==1)], 35), percentile(variants['DP'][(variants['CodingNumeric']==1)], 95)],
    verticalLineColours=['blue', 'blue'],
    plotFilestem=PLOT_DIR + '/AllCoding_35_95')

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'DP',
    bins=200, extendedXlabels=True,
    verticalLineValues= [percentile(variants['DP'][(variants['CodingNumeric']==1)], 25), percentile(variants['DP'][(variants['CodingNumeric']==1)], 99)],
    verticalLineColours=['blue', 'blue'],
    plotFilestem=PLOT_DIR + '/AllCoding_25_99')

# <codecell>

characterisation_plot(variants[(variants['isInCore']==1) & (variants['CodingNumeric']==1)], 'DP',
    bins=200, extendedXlabels=True,
    plotFilestem=PLOT_DIR + '/CoreCoding')

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'propDPdropDueToSoM',
    bins=200, extendedXlabels=True,
    plotFilestem=PLOT_DIR + '/AllCoding')

# <codecell>

DPcoreCoding = variants['DP'][(variants['isInCore']==1) & (variants['CodingNumeric']==1)]
DPallCoding = variants['DP'][(variants['CodingNumeric']==1)]
PreSoMcoreCoding = variants['preSoMcoverage'][(variants['isInCore']==1) & (variants['CodingNumeric']==1)]
DPcoreNoncoding = variants['DP'][(variants['isInCore']==1) & (variants['CodingNumeric']==0)]
DPallNoncoding = variants['DP'][(variants['CodingNumeric']==0)]

# <codecell>

binwidth=2000

lowerBoundCoding = percentile(DPallCoding, 25)
upperBoundCoding = percentile(DPallCoding, 99)
lowerBoundNoncoding = percentile(DPallNoncoding, 40)
upperBoundNoncoding = percentile(DPallNoncoding, 95)
print(lowerBoundCoding)
print(upperBoundCoding)
print(lowerBoundNoncoding)
print(upperBoundNoncoding)

bins=arange(0, 600001, binwidth)

figure(figsize=(16, 4))
hist(DPallCoding, bins=bins, histtype='step', linewidth=2, color='blue', label='VCF DP all coding')
hist(DPallNoncoding, bins=bins, histtype='step', linewidth=2, color='pink', label='VCF DP all non-coding')
axvline(lowerBoundCoding, color='blue')
axvline(upperBoundCoding, color='blue')
axvline(lowerBoundNoncoding, color='pink')
axvline(upperBoundNoncoding, color='pink')
legend()
xlabel('Coverage')
ylabel('Frequency (number of SNPs)')
savefig(PLOT_DIR + '/coverageThresholdsFinal.pdf', bbox_inches='tight')

# <codecell>

coverageThresholds = OrderedDict([
    ('Current',                         {'variable': 'DP',             'lower': 81641,                                                'upper':  300910}),
    ('Median core post-SoM +/- 1 20%',  {'variable': 'DP',             'lower': np.median(DPcoreCoding)*0.8,                          'upper':  np.median(DPcoreCoding)*1.2}),
    ('post-SoM core 25-95 percentiles', {'variable': 'DP',             'lower': percentile(DPcoreCoding, 25),                         'upper':  percentile(DPcoreCoding, 95)}),
    ('post-SoM all 35-95 percentiles',  {'variable': 'DP',             'lower': percentile(DPallCoding, 35),                          'upper':  percentile(DPallCoding, 95)}),
    ('post-SoM all 25-99 percentiles',  {'variable': 'DP',             'lower': percentile(DPallCoding, 25),                          'upper':  percentile(DPallCoding, 99)}),
    ('Median core post-Som +/- 1 SD',   {'variable': 'DP',             'lower': np.median(DPcoreCoding)-np.std(DPcoreCoding),         'upper':  np.median(DPcoreCoding)+np.std(DPcoreCoding)}),
    ('Median core pre-SoM +/- 1 SD',    {'variable': 'preSoMcoverage', 'lower': np.median(PreSoMcoreCoding)-np.std(PreSoMcoreCoding), 'upper':  np.median(PreSoMcoreCoding)+np.std(PreSoMcoreCoding)})
])

# <codecell>

thresholdSummaryDtypes = [
    ('lowCoverageThreshold', int),
    ('highCoverageThreshold', int),
    ('numberRemovedLowCoverage', int),
    ('numberRemovedHighCoverage', int),
    ('numberRetainedAfterCoverageFilters', int),
    ('numberRemovedLowCoverageExclusively', int),
    ('numberRemovedHighCoverageExclusively', int),
    ('numberRetainedAfterAllFilters', int),
    ('mean3D7missingRetainedVariants', float),
    ('mean3D7homAltsRetainedVariants', float),
    ('meanNumberOfHetsInClonalRetainedVariants', float),
    ('meanNumberOfDiscordantCleanCallsRetainedVariants', float),
    ('meanNumberOfDiscordantCalledCallsRetainedVariants', float)
]

def thresholdEffectsSummary(allPotentialVariants, lower, upper, variable, passOtherFilters=None):
    if passOtherFilters==None:
        passOtherFilters = (allPotentialVariants['PASSexceptCoverage']==1)
    variantsRemovedByLowerThreshold = allPotentialVariants[(allPotentialVariants[variable] < lower)]
    variantsRemovedByUpperThreshold = allPotentialVariants[(allPotentialVariants[variable] > upper)]
    variantsRetained = allPotentialVariants[(allPotentialVariants[variable] >= lower) & (allPotentialVariants[variable] <= upper)]
    return(
        np.array(
            [(
                lower,
                upper,
                len(variantsRemovedByLowerThreshold),
                len(variantsRemovedByUpperThreshold),
                len(variantsRetained),
                len(allPotentialVariants[(allPotentialVariants[variable] < lower) & passOtherFilters]),
                len(allPotentialVariants[(allPotentialVariants[variable] > upper) & passOtherFilters]),
                len(allPotentialVariants[(allPotentialVariants[variable] >= lower) & (allPotentialVariants[variable] <= upper) & passOtherFilters]),
                np.sum(variantsRetained['noCalls3D7']),
                np.sum(variantsRetained['homAlts3D7']),
                np.mean(variantsRetained['numHetClonal']),
                np.mean((variantsRetained['proportionDiscordantCleanCalls_2_2'] + variantsRetained['proportionDiscordantCleanCalls_x'] + variantsRetained['proportionDiscordantCleanCalls_W']) / 3.0),
                np.mean((variantsRetained['proportionDiscordantCalledCalls_2_2'] + variantsRetained['proportionDiscordantCalledCalls_x'] + variantsRetained['proportionDiscordantCalledCalls_W']) / 3.0)
            )],
            dtype = thresholdSummaryDtypes
        )
    )
            

# <codecell>

max(variants['DP'][(variants['CodingNumeric']==1) & (variants['noCalls3D7'] > 0)])

# <codecell>

thresholdsSummaryTable = np.empty([0], dtype = thresholdSummaryDtypes)
for coverageThreshold in coverageThresholds:
    print(coverageThreshold)
    thresholdsSummaryTable = np.append(
        thresholdsSummaryTable,
        thresholdEffectsSummary(
            variants[(variants['CodingNumeric']==1)],
            coverageThresholds[coverageThreshold]['lower'],
            coverageThresholds[coverageThreshold]['upper'],
            coverageThresholds[coverageThreshold]['variable']
        )
    )

# <codecell>

etl(petlx.array.fromarray(thresholdsSummaryTable)).transpose().pushheader([''] + coverageThresholds.keys())

# <codecell>

tocsv(etl(petlx.array.fromarray(thresholdsSummaryTable)).transpose().pushheader([''] + coverageThresholds.keys()), PLOT_DIR + '/thresholdsSummaryTable.csv')

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==0)], 'DP',
    bins=200, extendedXlabels=True,
    verticalLineValues= [percentile(variants['DP'][(variants['CodingNumeric']==0)], 40), percentile(variants['DP'][(variants['CodingNumeric']==0)], 95)],
    verticalLineColours=['red', 'red'],
    plotFilestem=PLOT_DIR + '/AllNoncoding_40_95')

# <headingcell level=3>

# Suggested plot for supplementary material of future manuscript

# <codecell>

manuscriptAttributeDescriptions = OrderedDict([
    ('homAlts3D7',                          {'description': 'Number of hom alt calls in 4 samples of 3D7',                    'statistic': 'sum',   'colour': '#FF99BB'}),
    ('numHetClonal',                        {'description': 'Mean heterozygosity of 100 clonal samples',                      'statistic': 'mean',  'colour': '#CC00FF'}),
    ('proportionDiscordantCleanCalls_all',  {'description': 'Mean proportion of discordant clean calls in all duplicates',    'statistic': 'mean',  'colour': '#FF66FF'}),
    ('proportionDiscordantCalledCalls_all', {'description': 'Mean proportion of discordant called calls in all duplicates',   'statistic': 'mean',  'colour': '#FFCCFF'})
])

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'DP', yattrDescriptions=manuscriptAttributeDescriptions,
    bins=200, extendedXlabels=True, plotHeight=12,
    verticalLineValues= [percentile(variants['DP'][(variants['CodingNumeric']==1)], 35), percentile(variants['DP'][(variants['CodingNumeric']==1)], 95)],
    verticalLineColours=['blue', 'blue'],
    plotFilestem=PLOT_DIR + '/AllCoding_35_95_manuscript')

# <codecell>

characterisation_plot(variants[(variants['CodingNumeric']==1)], 'DP', yattrDescriptions=manuscriptAttributeDescriptions,
    bins=200, extendedXlabels=True, plotHeight=12,
    verticalLineValues= [percentile(variants['DP'][(variants['CodingNumeric']==1)], 25), percentile(variants['DP'][(variants['CodingNumeric']==1)], 99)],
    verticalLineColours=['blue', 'blue'],
    plotFilestem=PLOT_DIR + '/AllCoding_25_99_manuscript')

# <codecell>


