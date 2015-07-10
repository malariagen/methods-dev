# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Downloads
# First download latest version (3.4-0) from http://downloads.sourceforge.net/project/cortexassembler/cortex_var/latest/CORTEX_release_v1.0.5.21.tgz

# <codecell>

profile_fn = '~/.profile'

install_dir = '../opt'

cortex_dir = install_dir + '/cortex'
!mkdir -p {cortex_dir}
cortex_exe = cortex_dir + '/???'

# <codecell>

current_dir = !pwd
current_dir = current_dir[0]

# <codecell>

current_dir

# <headingcell level=1>

# Install Cortex

# <codecell>

!cp ~/Downloads/CORTEX_release_v1.0.5.21.tgz {cortex_dir}

# <codecell>

%cd {cortex_dir}
!tar -xvf CORTEX_release_v1.0.5.21.tgz 2> /dev/null

# <codecell>

%cd CORTEX_release_v1.0.5.21
!bash install.sh
!make cortex_var

# <codecell>

%

