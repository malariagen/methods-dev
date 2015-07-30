# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

sanger_sequenom_fn = '/nfs/team112_internal/rp7/recon/sanger_sequenom/processed_data/DK_KR_iPLEX_DK1079_W1381.txt'
calling_parameters_fn = '/nfs/team112_internal/rp7/recon/sanger_sequenom/processed_data/DK1079_W1381_sanger_test/data_DK1079_W1381_parameters.txt'

# <codecell>


# <codecell>

# import brewer2mpl
# set1 = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
genotype_colours = ['lightgray', 'black', 'red', 'blue', 'green']
# genotype_colours = ['lightgray', 'black', 'red', 'blue', 'green']

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

tbl_sanger_sequenom = etl.fromtsv(sanger_sequenom_fn).convertnumbers().sort(['PLATE', 'WELL_POSITION', 'ASSAY_ID', 'ALLELE'])
tbl_sanger_sequenom

# <codecell>

tbl_calling_parameters_sep = etl.fromtsv(calling_parameters_fn).selectne('assay_code', 'amelogenin_XY_SNP1_E').convertnumbers()
# tbl_calling_parameters_sep
tbl_calling_parameters_1 = (tbl_calling_parameters_sep
    .addrownumbers()
    .select(lambda rec: rec['row'] % 2 == 1)
    .rename('allele', 'allele_1')
    .rename('lift_degrees', 'lift_degrees_1')
    .rename('genotype_threshold_degrees', 'genotype_threshold_degrees_1')
    .rename('low_intensity_cutoff', 'low_intensity_cutoff_1')
    .cut(['plate_code', 'assay_code', 'alleles_pair', 'allele_1', 'lift_degrees_1', 'genotype_threshold_degrees_1', 'low_intensity_cutoff_1'])
)
tbl_calling_parameters_2 = (tbl_calling_parameters_sep
    .addrownumbers()
    .select(lambda rec: rec['row'] % 2 == 0)
    .rename('allele', 'allele_2')
    .rename('lift_degrees', 'lift_degrees_2')
    .rename('genotype_threshold_degrees', 'genotype_threshold_degrees_2')
    .rename('low_intensity_cutoff', 'low_intensity_cutoff_2')
    .cut(['allele_2', 'lift_degrees_2', 'genotype_threshold_degrees_2', 'low_intensity_cutoff_2'])
)
tbl_calling_parameters = (
    tbl_calling_parameters_1
    .annex(tbl_calling_parameters_2)
)
tbl_calling_parameters.displayall()

# <codecell>

tbl_sanger_sequenom.valuecounts('STATUS').displayall()

# <codecell>

tbl_sanger_sequenom.addrownumbers().select(lambda rec: rec['row'] % 2 == 0).valuecounts('STATUS').displayall()

# <codecell>

def calc_low_intensity(x=3, y=2, theta=45):
    tan_sq_theta = math.tan(math.radians(theta))**2
    return(math.sqrt((tan_sq_theta + 1) * (((x**2)*(y**2))/(((x**2)*tan_sq_theta) + y**2))))

# <codecell>

def determine_genotype(rec):
    if rec['INTENSITY'] < rec['low_intensity_cutoff']:
        return('./.')
    elif rec['THETA'] < rec['genotype_threshold_degrees_1']:
        return('0/0')
    elif rec['THETA'] > rec['genotype_threshold_degrees_2']:
        return('1/1')
    else:
        return('0/1')

# <codecell>

tbl_first_allele = (tbl_sanger_sequenom
    .addrownumbers()
    .select(lambda rec: rec['row'] % 2 == 1)
    .rename('ALLELE', 'ALLELE_1')
    .rename('HEIGHT', 'HEIGHT_1')
    .rename('MASS', 'MASS_1')
    .rename('GENOTYPE_ID', 'SEQUENOM_GENOTYPE')
    .cut(['ASSAY_ID', 'SEQUENOM_GENOTYPE', 'SAMPLE_ID', 'STATUS', 'WELL_POSITION', 'ALLELE_1', 'HEIGHT_1', 'MASS_1'])
)

tbl_second_allele = (tbl_sanger_sequenom
    .addrownumbers()
    .select(lambda rec: rec['row'] % 2 == 0)
    .rename('ALLELE', 'ALLELE_2')
    .rename('HEIGHT', 'HEIGHT_2')
    .rename('MASS', 'MASS_2')
    .rename('GENOTYPE_ID', 'SEQUENOM_GENOTYPE')
    .cut(['ALLELE_2', 'HEIGHT_2', 'MASS_2'])
)

tbl_results = (
    tbl_first_allele
    .annex(tbl_second_allele)
    .addfield('THETA', lambda rec: math.degrees(math.atan(np.float64(rec['HEIGHT_2'])/np.float64(rec['HEIGHT_1']))))
    .addfield('INTENSITY', lambda rec: math.sqrt(math.pow(np.float64(rec['HEIGHT_1']), 2) + math.pow(np.float64(rec['HEIGHT_2']), 2)))
    .join(tbl_calling_parameters, lkey='ASSAY_ID', rkey='assay_code')
#     .convert('low_intensity_cutoff_2', lambda val: val+1)
    .addfield('low_intensity_cutoff', lambda rec: calc_low_intensity(rec['low_intensity_cutoff_1'], rec['low_intensity_cutoff_2'], rec['THETA']))
    .addfield('genotype', determine_genotype)
)
tbl_results

# <codecell>

tbl_results.valuecounts('SEQUENOM_GENOTYPE', 'genotype', 'ALLELE_1', 'ALLELE_2').displayall()
# tbl_results.valuecounts('SEQUENOM_GENOTYPE').displayall()

# <headingcell level=1>

# Scatter plots

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_results
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)
    fields_of_interest = ['genotype_threshold_degrees_1', 'genotype_threshold_degrees_2',
                          'low_intensity_cutoff_1', 'low_intensity_cutoff_2']
    (gtd1, gtd2, lic1, lic2) = tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).cut(fields_of_interest).data()[0]

#     alleles = (tbl_results
#         .selecteq('ASSAY_ID', ASSAY_ID)
#         .distinct('ALLELE')
#         .values('ALLELE')
#         .array()
#     )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_results
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('genotype')
#         .addfield('genotype_length', lambda rec: len(rec['SEQUENOM_GENOTYPE']))
#         .sort('genotype_length')
        .values('genotype')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    max_x = 0
    max_y = 0
    for j, genotype in enumerate(insert(genotypes, 0, None)):
#         heights = (tbl_results
#             .selecteq('ASSAY_ID', ASSAY_ID)
#             .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .cut('HEIGHT_1', 'HEIGHT_2')
#             .toarray()
#         )
        x_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[0])
            .values('HEIGHT_1')
            .array()
        )
        y_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[1])
            .values('HEIGHT_2')
            .array()
        )
#         ax.scatter(heights['HEIGHT_1'], heights['HEIGHT_2'], color=genotype_colours[j], alpha=0.5, label=genotype)
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
    
#         print(x_heights)
        
        if len(x_heights) > 0 and np.nanmax(x_heights) > max_x:
            max_x = np.nanmax(x_heights)
        if len(y_heights) > 0 and np.nanmax(y_heights) > max_y:
            max_y = np.nanmax(y_heights)
        
        
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)
    
    ax.plot([0, max_x], [0, max_x*math.tan(math.radians(gtd1))], color='blue', linestyle='-', linewidth=1)
    ax.plot([0, max_y*math.tan(math.radians(90.0-gtd2))], [0, max_y], color='red', linestyle='-', linewidth=1)
    
    li_x = np.linspace(0, lic1, 100)
    li_y = np.sqrt( (lic2**2) * (1 - ( (li_x**2) / (lic1**2) ) ) )
    ax.plot(li_x, li_y, color='black', linestyle='-', linewidth=1)
    
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)

fig.tight_layout()

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_results
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)
    print(tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).data()[0])
    fields_of_interest = ['genotype_threshold_degrees_1', 'genotype_threshold_degrees_2',
                          'low_intensity_cutoff_1', 'low_intensity_cutoff_2']
    (gtd1, gtd2, lic1, lic2) = tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).cut(fields_of_interest).data()[0]

#     alleles = (tbl_results
#         .selecteq('ASSAY_ID', ASSAY_ID)
#         .distinct('ALLELE')
#         .values('ALLELE')
#         .array()
#     )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_results
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('SEQUENOM_GENOTYPE')
        .addfield('genotype_length', lambda rec: len(rec['SEQUENOM_GENOTYPE']))
        .sort('genotype_length')
        .values('SEQUENOM_GENOTYPE')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    for j, genotype in enumerate(insert(genotypes, 0, None)):
#         heights = (tbl_results
#             .selecteq('ASSAY_ID', ASSAY_ID)
#             .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .cut('HEIGHT_1', 'HEIGHT_2')
#             .toarray()
#         )
        x_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .selecteq('ALLELE', alleles[0])
            .values('HEIGHT_1')
            .array()
        )
        y_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .selecteq('ALLELE', alleles[1])
            .values('HEIGHT_2')
            .array()
        )
#         ax.scatter(heights['HEIGHT_1'], heights['HEIGHT_2'], color=genotype_colours[j], alpha=0.5, label=genotype)
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)

fig.tight_layout()

# <codecell>

70*math.tan(math.radians(10))

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_results
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)
    fields_of_interest = ['genotype_threshold_degrees_1', 'genotype_threshold_degrees_2',
                          'low_intensity_cutoff_1', 'low_intensity_cutoff_2']
    (gtd1, gtd2, lic1, lic2) = tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).cut(fields_of_interest).data()[0]

#     alleles = (tbl_results
#         .selecteq('ASSAY_ID', ASSAY_ID)
#         .distinct('ALLELE')
#         .values('ALLELE')
#         .array()
#     )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_results
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('genotype')
#         .addfield('genotype_length', lambda rec: len(rec['SEQUENOM_GENOTYPE']))
#         .sort('genotype_length')
        .values('genotype')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    max_x = 0
    max_y = 0
    for j, genotype in enumerate(insert(genotypes, 0, None)):
#         heights = (tbl_results
#             .selecteq('ASSAY_ID', ASSAY_ID)
#             .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .cut('HEIGHT_1', 'HEIGHT_2')
#             .toarray()
#         )
        x_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[0])
            .values('HEIGHT_1')
            .array()
        )
        y_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[1])
            .values('HEIGHT_2')
            .array()
        )
#         ax.scatter(heights['HEIGHT_1'], heights['HEIGHT_2'], color=genotype_colours[j], alpha=0.5, label=genotype)
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
    
#         print(x_heights)
        
        if len(x_heights) > 0 and np.nanmax(x_heights) > max_x:
            max_x = np.nanmax(x_heights)
        if len(y_heights) > 0 and np.nanmax(y_heights) > max_y:
            max_y = np.nanmax(y_heights)
        
        
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)
    
    ax.plot([0, max_x], [0, max_x*math.tan(math.radians(gtd1))], color='blue', linestyle='-', linewidth=1)
    ax.plot([0, max_y*math.tan(math.radians(90.0-gtd2))], [0, max_y], color='red', linestyle='-', linewidth=1)
    
    li_x = np.linspace(0, lic1, 100)
    li_y = np.sqrt( (lic2**2) * (1 - ( (li_x**2) / (lic1**2) ) ) )
    ax.plot(li_x, li_y, color='black', linestyle='-', linewidth=1)
    
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)

fig.tight_layout()

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_results
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)
    fields_of_interest = ['genotype_threshold_degrees_1', 'genotype_threshold_degrees_2',
                          'low_intensity_cutoff_1', 'low_intensity_cutoff_2']
    (gtd1, gtd2, lic1, lic2) = tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).cut(fields_of_interest).data()[0]

    genotypes = (tbl_results
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('genotype')
#         .addfield('genotype_length', lambda rec: len(rec['SEQUENOM_GENOTYPE']))
#         .sort('genotype_length')
        .values('genotype')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    max_x = 0
    max_y = 0
    for j, genotype in enumerate(insert(genotypes, 0, None)):
        x_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[0])
            .values('THETA')
            .array()
        )
        y_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[1])
            .values('INTENSITY')
            .array()
        )
#         ax.scatter(heights['HEIGHT_1'], heights['HEIGHT_2'], color=genotype_colours[j], alpha=0.5, label=genotype)
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
    
#         print(x_heights)
        
#         if len(x_heights) > 0 and np.nanmax(x_heights) > max_x:
#             max_x = np.nanmax(x_heights)
#         if len(y_heights) > 0 and np.nanmax(y_heights) > max_y:
#             max_y = np.nanmax(y_heights)
        
        
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)
    
    ax.axvline(gtd1, color='blue')
    ax.axvline(gtd2, color='red')
    
#     ax.plot([0, max_x], [0, max_x*math.tan(math.radians(gtd1))], color='blue', linestyle='-', linewidth=1)
#     ax.plot([0, max_y*math.tan(math.radians(90.0-gtd2))], [0, max_y], color='red', linestyle='-', linewidth=1)
    
    li_x = np.linspace(0, 90, 100)
    li_y = np.array([calc_low_intensity(lic1, lic2, theta) for theta in li_x])
#     li_y = np.sqrt( (lic2**2) * (1 - ( (li_x**2) / (lic1**2) ) ) )
    ax.plot(li_x, li_y, color='black', linestyle='-', linewidth=1)
    ax.set_xlim(0, 90)
    ax.set_ylim(bottom=0)

fig.tight_layout()

# <codecell>

fig = figure(figsize=(12, 7.5))
for i, ASSAY_ID in enumerate(tbl_results
    .distinct('ASSAY_ID')
    .values('ASSAY_ID')
    .array()
):
    print(ASSAY_ID)
    fields_of_interest = ['genotype_threshold_degrees_1', 'genotype_threshold_degrees_2',
                          'low_intensity_cutoff_1', 'low_intensity_cutoff_2']
    (gtd1, gtd2, lic1, lic2) = tbl_calling_parameters.selecteq('assay_code', ASSAY_ID).cut(fields_of_interest).data()[0]

#     alleles = (tbl_results
#         .selecteq('ASSAY_ID', ASSAY_ID)
#         .distinct('ALLELE')
#         .values('ALLELE')
#         .array()
#     )

    # genotypes = [alleles[0], alleles[1], alleles[0]+alleles[1], '']
    genotypes = (tbl_results
        .selecteq('ASSAY_ID', ASSAY_ID)
        .distinct('genotype')
#         .addfield('genotype_length', lambda rec: len(rec['SEQUENOM_GENOTYPE']))
#         .sort('genotype_length')
        .values('genotype')
        .array()
    )

    ax = fig.add_subplot(3, 5, i+1)
    max_x = 0
    max_y = 0
    for j, genotype in enumerate(insert(genotypes, 0, None)):
#         heights = (tbl_results
#             .selecteq('ASSAY_ID', ASSAY_ID)
#             .selecteq('SEQUENOM_GENOTYPE', genotype)
#             .cut('HEIGHT_1', 'HEIGHT_2')
#             .toarray()
#         )
        x_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[0])
            .values('HEIGHT_1')
            .array()
        )
        y_heights = (tbl_results
            .selecteq('ASSAY_ID', ASSAY_ID)
            .selecteq('genotype', genotype)
#             .selecteq('ALLELE', alleles[1])
            .values('HEIGHT_2')
            .array()
        )
#         ax.scatter(heights['HEIGHT_1'], heights['HEIGHT_2'], color=genotype_colours[j], alpha=0.5, label=genotype)
        ax.scatter(x_heights, y_heights, color=genotype_colours[j], alpha=0.5, label=genotype)
    
#         print(x_heights)
        
        if len(x_heights) > 0 and np.nanmax(x_heights) > max_x:
            max_x = np.nanmax(x_heights)
        if len(y_heights) > 0 and np.nanmax(y_heights) > max_y:
            max_y = np.nanmax(y_heights)
        
        
        assay_split=ASSAY_ID.split('_')
        assay_title = ASSAY_ID if not ASSAY_ID.startswith('Pf') else "%s %s" % (
            '_'.join((assay_split[3], assay_split[4], assay_split[5])),
            assay_split[6]
        )
        ax.set_title(assay_title)
    
    ax.plot([0, max_x], [0, max_x*math.tan(math.radians(gtd1))], color='blue', linestyle='-', linewidth=1)
    ax.plot([0, max_y*math.tan(math.radians(90.0-gtd2))], [0, max_y], color='red', linestyle='-', linewidth=1)
    
    li_x = np.linspace(0, lic1, 100)
    li_y = np.sqrt( (lic2**2) * (1 - ( (li_x**2) / (lic1**2) ) ) )
    ax.plot(li_x, li_y, color='black', linestyle='-', linewidth=1)
    
    ax.set_xlim(0, 2*lic1)
    ax.set_ylim(0, 2*lic2)

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


