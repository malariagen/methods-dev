# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <headingcell level=1>

# Transfer bams

# <codecell>

(tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
    .addfield('bam_dir', lambda rec: os.path.dirname(rec['bam_fn']) if rec['bam_fn'] is not None else None)
).displayall(index_header=True)

# <codecell>

for rec in (tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
    .addfield('bam_dir', lambda rec: os.path.dirname(rec['bam_fn']) if rec['bam_fn'] is not None else None)
).data():
    print(rec[0], rec[4])
    !mkdir -p {rec[10]}
    !rsync -avL malsrv2:{rec[9]} {rec[10]}
    !samtools index {rec[9]}

# <codecell>


