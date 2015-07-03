# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Downloads
# First download latest version (3.4-0) from https://www.broadinstitute.org/gatk/download
# 
# Download bwa from http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2
# 
# Download samtools from http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2
# 
# Download picard from https://github.com/broadinstitute/picard/releases/download/1.133/picard-tools-1.133.zip
# 
# Download bcftools from https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 (not required by GATK, but thought would be good)
# 
# IGV installed using mac version at https://www.broadinstitute.org/software/igv/download

# <headingcell level=1>

# Setup

# <codecell>

profile_fn = '~/.profile'

install_dir = '../opt'
!mkdir -p {install_dir}

bwa_dir = install_dir + '/bwa'
!mkdir -p {bwa_dir}
bwa_exe = bwa_dir + '/bwa-0.7.12/bwa'

samtools_dir = install_dir + '/samtools'
!mkdir -p {samtools_dir}
samtools_exe = samtools_dir + '/samtools-1.2/samtools'

picard_dir = install_dir + '/picard'
!mkdir -p {picard_dir}
picard_exe = 'java -jar ' + picard_dir + '/picard-tools-1.133/picard.jar'

bcftools_dir = install_dir + '/bcftools'
!mkdir -p {bcftools_dir}
bcftools_exe = bcftools_dir + '/bcftools-1.2/bcftools'

gatk_dir = install_dir + '/gatk'
!mkdir -p {gatk_dir}
gatk_exe = 'java -jar ' + gatk_dir + '/GenomeAnalysisTK.jar'

# <codecell>

bwa_exe

# <headingcell level=1>

# Install bwa

# <codecell>

!cp ~/Downloads/bwa-0.7.12.tar.bz2 {bwa_dir}

# <codecell>

%cd {bwa_dir}
!tar -xvf bwa-0.7.12.tar.bz2

# <codecell>

%cd bwa-0.7.12
!make

# <codecell>

%cd ../../../notebooks/

# <codecell>

# Add path to .profile
!echo "" >> {profile_fn}
!echo "# added by Ipython notebook 20150603_install_GATK.ipynb" >> {profile_fn}
!echo "export PATH=\"/Users/rpearson/src/github/malariagen/methods-dev/recon/opt/bwa/bwa-0.7.12:\$$PATH\"" >> {profile_fn}



# <codecell>

# Check bwa works
!{bwa_exe}

# <headingcell level=1>

# Install samtools

# <codecell>

!cp /Users/rpearson/Downloads/samtools-1.2.tar.bz2 {samtools_dir}
%cd {samtools_dir}
!tar xvzf samtools-1.2.tar.bz2 
%cd samtools-1.2 
!make
%cd ../../../notebooks/

# <codecell>

# Add path to .profile
!echo "" >> {profile_fn}
!echo "# added by Ipython notebook 20150603_install_GATK.ipynb" >> {profile_fn}
!echo "export PATH=\"/Users/rpearson/src/github/malariagen/methods-dev/recon/opt/samtools/samtools-1.2:\$$PATH\"" >> {profile_fn}

# <codecell>

# Check samtools works
!{samtools_exe}

# <codecell>


# <headingcell level=1>

# Install picard

# <codecell>

!cp /Users/rpearson/Downloads/picard-tools-1.133.zip {picard_dir}
%cd {picard_dir}
!tar xjf picard-tools-1.133.zip
# Add environment vairable to .profile
!echo "" >> {profile_fn}
!echo "# added by Ipython notebook 20150603_install_GATK.ipynb" >> {profile_fn}
!echo "export PICARD=\"/Users/rpearson/src/github/malariagen/methods-dev/recon/opt/picard/picard-tools-1.133/picard.jar\"" >> {profile_fn}
%cd ../../notebooks/

# <codecell>

# Test picard
!{picard_exe} GenotypeConcordance -h

# <headingcell level=1>

# Install bcftools

# <codecell>

!cp /Users/rpearson/Downloads/bcftools-1.2.tar.bz2 {bcftools_dir}
%cd {bcftools_dir}
!tar xvzf bcftools-1.2.tar.bz2 
%cd bcftools-1.2 
!make
%cd ../../../notebooks/

# <codecell>

# Add path to .profile
!echo "" >> {profile_fn}
!echo "# added by Ipython notebook 20150603_install_GATK.ipynb" >> {profile_fn}
!echo "export PATH=\"/Users/rpearson/src/github/malariagen/methods-dev/recon/opt/bcftools/bcftools-1.2:\$$PATH\"" >> {profile_fn}

# <codecell>

# Check bcftools works
!{bcftools_exe}

# <headingcell level=1>

# Install GATK

# <codecell>

!cp /Users/rpearson/Downloads/GenomeAnalysisTK-3.4-0.tar.bz2 {gatk_dir}
%cd {gatk_dir}
!tar xjf GenomeAnalysisTK-3.4-0.tar.bz2
# Add environment vairable to .profile
!echo "" >> {profile_fn}
!echo "# added by Ipython notebook 20150603_install_GATK.ipynb" >> {profile_fn}
!echo "export GATK=\"/Users/rpearson/src/github/malariagen/methods-dev/recon/opt/gatk/GenomeAnalysisTK.jar\"" >> {profile_fn}
%cd ../../notebooks/

# <codecell>

# Test GATK
!{gatk_exe} -T VariantsToAllelicPrimitives -h

# <codecell>


