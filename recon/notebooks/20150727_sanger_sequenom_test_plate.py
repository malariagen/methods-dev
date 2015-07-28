# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

sanger_sequenom_fn = '/nfs/team112_internal/rp7/recon/sanger_sequenom/processed_data/DK_KR_iPLEX_DK1079_W1381_heights.xls'

# <codecell>

# import brewer2mpl
# set1 = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
genotype_colours = ['lightgray', 'black', 'red', 'blue', 'green']

# <headingcell level=1>

# Functions

# <codecell>

def determine_GT(AD):
    if sum(AD) < 5:
        return(None)
    if len(AD) == 3:
        if AD[0] > 1 and AD[1] > 1:
            return("0/1")
        if AD[0] > 1 and AD[2] > 1:
            return("0/2")
        if AD[1] > 1 and AD[2] > 1:
            return("1/2")
    if(AD[0] > 1 and AD[1] > 1):
        return("0/1")
    if(AD[0] <= 1 and AD[1] > 1):
        return("1/1")
    if(AD[0] > 1 and AD[1] <= 1):
        return("0/0")
    else:
        return("unknown")
    

# <headingcell level=1>

# Load data

# <codecell>

tbl_sanger_sequenom = etl.fromtsv(sanger_sequenom_fn)
tbl_sanger_sequenom

# <headingcell level=1>

# Scatter plots

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_sanger_sequenom
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)

    alleles = (tbl_sanger_sequenom
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('ALLELE')
        .values('ALLELE')
        .array()
    )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_sanger_sequenom
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('GENOTYPE_ID')
        .addfield('genotype_length', lambda rec: len(rec['GENOTYPE_ID']))
        .sort('genotype_length')
        .values('GENOTYPE_ID')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    for j, genotype in enumerate(insert(genotypes, 0, None)):
        x_heights = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('GENOTYPE_ID', genotype)
            .selecteq('ALLELE', alleles[0])
            .values('HEIGHT')
            .array()
        )
        y_heights = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('GENOTYPE_ID', genotype)
            .selecteq('ALLELE', alleles[1])
            .values('HEIGHT')
            .array()
        )
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)

fig.tight_layout()

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_sanger_sequenom
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)

    alleles = (tbl_sanger_sequenom
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('ALLELE')
        .values('ALLELE')
        .array()
    )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_sanger_sequenom
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('GENOTYPE_ID')
        .addfield('genotype_length', lambda rec: len(rec['GENOTYPE_ID']))
        .sort('genotype_length')
        .values('GENOTYPE_ID')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    for j, genotype in enumerate(insert(genotypes, 0, None)):
        x_heights = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('GENOTYPE_ID', genotype)
            .selecteq('ALLELE', alleles[0])
            .values('HEIGHT')
            .array()
        )
        y_heights = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('GENOTYPE_ID', genotype)
            .selecteq('ALLELE', alleles[1])
            .values('HEIGHT')
            .array()
        )
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)

fig.tight_layout()

# <codecell>

(tbl_sanger_sequenom
.convertnumbers()
.selecteq('ASSAY_ID', 'gb_CM000450_1697797')
.selectgt('HEIGHT', 50)
).displayall()

# <codecell>

(tbl_sanger_sequenom
.convertnumbers()
.selecteq('ASSAY_ID', 'gb_CM000452_238200')
.selectgt('HEIGHT', 50)
).displayall()

# <codecell>

(etl.fromvcf(PGV4_VCF_FN, chrom='Pf3D7_13_v3', start=1725259, stop=1725259)
    .vcfmeltsamples()
    .vcfunpackcall()
    .cutout('GT')
    .addfield('GT', lambda rec: determine_GT(rec['AD']))
    .valuecounts('GT')
)

# <codecell>


# <codecell>

fig = figure(figsize=(12, 5))
for i, ASSAY_ID in enumerate(tbl_sanger_sequenom
    .distinct('ASSAY_ID')
    .select(lambda rec: rec['ASSAY_ID'].startswith('Pf'))
    .values('ASSAY_ID')
    .array()
):
    assay_split=ASSAY_ID.split('_')
    if(assay_split[0]=='Pf'):
        print(ASSAY_ID)
        assay_chrom = '_'.join((assay_split[3], assay_split[4], assay_split[5]))
        assay_pos = int(assay_split[6])

        alleles = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .distinct('ALLELE')
            .values('ALLELE')
            .array()
        )

        # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
        genotypes = (tbl_sanger_sequenom
            .selecteq('ASSAY_ID', ASSAY_ID)
            .distinct('GENOTYPE_ID')
            .addfield('genotype_length', lambda rec: len(rec['GENOTYPE_ID']))
            .sort('genotype_length')
            .values('GENOTYPE_ID')
            .array()
        )

        tbl_pgv4 = (etl.fromvcf(PGV4_VCF_FN, chrom=assay_chrom, start=assay_pos, stop=assay_pos)
            .vcfmeltsamples()
            .vcfunpackcall()
            .cutout('GT')
            .addfield('GT', lambda rec: determine_GT(rec['AD']))
            .cut(['SAMPLE', 'GT'])
        )

        ax = fig.add_subplot(2, 5, i+1)
        for j, genotype in enumerate([None, 'WATER', '0/0', '1/1', '0/1']):
            x_heights = (tbl_sanger_sequenom
                .selecteq('ASSAY_ID', ASSAY_ID)
                .selecteq('ALLELE', alleles[0])
                .leftjoin(tbl_pgv4, lkey='SAMPLE_ID', rkey='SAMPLE')
                .update('GT', 'WATER', where=lambda rec: rec['SAMPLE_ID'] == 'WATER')
                .selecteq('GT', genotype)
                .values('HEIGHT')
                .array()
            )
            y_heights = (tbl_sanger_sequenom
                .selecteq('ASSAY_ID', ASSAY_ID)
                .selecteq('ALLELE', alleles[1])
                .leftjoin(tbl_pgv4, lkey='SAMPLE_ID', rkey='SAMPLE')
                .update('GT', 'WATER', where=lambda rec: rec['SAMPLE_ID'] == 'WATER')
                .selecteq('GT', genotype)
                .values('HEIGHT')
                .array()
            )
#             y_heights = (tbl_sanger_sequenom
#                 .selecteq('ASSAY_ID', ASSAY_ID)
#                 .selecteq('GENOTYPE_ID', genotype)
#                 .selecteq('ALLELE', alleles[1])
#                 .values('HEIGHT')
#                 .array()
#             )
            ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
            ax.set_title("%s %d" % (assay_chrom, assay_pos))

fig.tight_layout()

# <codecell>

(tbl_sanger_sequenom
    .selecteq('ASSAY_ID', ASSAY_ID)
    .selecteq('ALLELE', alleles[0])
    .leftjoin(tbl_pgv4, lkey='SAMPLE_ID', rkey='SAMPLE')
    .update('GT', 'WATER', where=lambda rec: rec['SAMPLE_ID'] == 'WATER')
    .update('GT', 'LAB', where=lambda rec: rec['SAMPLE_ID'].startswith('WT'))
    .valuecounts('GT')
    .displayall()
)

# <codecell>

150+125+70+18+16+5

# <codecell>

print("Number of samples in PGV4 = %d" % (150+70+5))

# <codecell>


