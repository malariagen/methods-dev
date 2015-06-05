# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

SANGER_BAM_DIR = '/nfs/team112_internal/production_files/Pf/4_0'

# <headingcell level=1>

# Load data

# <codecell>

tbl_sample_sets

# <codecell>

tbl_sample_bam_info = (
    tbl_sample_sets
    .join(tbl_pgv4_sample_manifest, lkey='sample_code', rkey='Sample')
    .cut(['sample_code', 'Study'])
    .addfield('sanger_bam_fn', lambda rec: "%s/%s/%s/%s.bam" % (
        SANGER_BAM_DIR,
        rec['Study'],
        rec['sample_code'].replace('-', '_'),
        rec['sample_code'].replace('-', '_')
    ))
    .addfield('local_dir', lambda rec: "%s/%s/%s/" % (
        SANGER_BAM_DIR,
        rec['Study'],
        rec['sample_code'].replace('-', '_')
    ))
)
print(len(tbl_sample_bam_info.data()))
tbl_sample_bam_info.display(index_header=True)

# <codecell>

for rec in tbl_sample_bam_info.data():
    print(rec[2])
    !mkdir -p {rec[3]}
    !rsync -avL malsrv2:{rec[2]} {rec[3]}
    !samtools index {rec[2]}

# <codecell>


