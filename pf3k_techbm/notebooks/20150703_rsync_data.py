# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

team112_dir = '/nfs/team112_internal/production_files/Pf3k/methods'
tdo_dir = '/nfs/pathogen003/tdo/Pf3K/SNPset/Pf3D7vsPfIT'
tdo2_dir = '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/embl_V1/Fasta_inChromosomes'
tdo3_dir = '/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011'
tdo4_dir = '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram'
IT_dir = '/nfs/team112_internal/production_files/Pf3k/methods/GATKbuild/IT'

# <codecell>

!mkdir -p {team112_dir}
!mkdir -p {tdo_dir}
!mkdir -p {tdo2_dir}
!mkdir -p {tdo3_dir}
!mkdir -p {tdo4_dir}
!mkdir -p {IT_dir}
!mkdir ../log

# <codecell>

# !rsync -avL rp7@malsrv2:{team112_dir}/ {team112_dir} > ../log/20150626_team112_dir.log
# !rsync -avL rp7@malsrv2:{tdo_dir}/ {tdo_dir} > ../log/20150626_tdo_dir.log
# !rsync -avL rp7@malsrv2:{tdo2_dir}/ {tdo2_dir} > ../log/20150626_tdo2_dir.log
# !rsync -avL rp7@malsrv2:{tdo3_dir}/ {tdo3_dir} > ../log/20150703_tdo3_dir.log
# !rsync -avL rp7@malsrv2:{tdo3_dir}/ {tdo3_dir} > ../log/20150703_tdo3_dir.log
!rsync -avL rp7@malsrv2:{tdo4_dir}/ {tdo4_dir} > ../log/20150706_tdo4_dir.log

# <codecell>

!rsync -avL rp7@malsrv2:{IT_dir}/ {IT_dir} > ../log/20150705_IT_dir.log

# <codecell>


