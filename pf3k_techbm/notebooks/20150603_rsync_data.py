# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

team112_dir = '/nfs/team112_internal/production_files/Pf3k/methods'
tdo_dir = '/nfs/pathogen003/tdo/Pf3K/SNPset/Pf3D7vsPfIT'

# <codecell>

!mkdir -p {team112_dir}
!mkdir -p {tdo_dir}
!mkdir ../log

# <codecell>

!rsync -avL rp7@malsrv2:{team112_dir} {team112_dir} > ../log/20150603_team112_dir.log
!rsync -avL rp7@malsrv2:{tdo_dir} {tdo_dir} > ../log/20150603_tdo_dir.log

# <markdowncell>

# Abandoned the above and ran from the command line

# <codecell>


