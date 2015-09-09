# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run standard_imports.ipynb

# <codecell>

latest_genedb_fasta_url = 'ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/2015-08/Pfalciparum.genome.fasta.gz'
latest_plasmodb_fasta_url = 'http://www.plasmodb.org/common/downloads/release-25/Pfalciparum3D7/fasta/data/PlasmoDB-25_Pfalciparum3D7_Genome.fasta'
my_fasta = '/lustre/scratch110/malaria/dj6/pf3k/resources/Pfalciparum.genome.fasta'
old_fasta = '/lustre/scratch109/malaria/pfalciparum/resources/3D7_V3.fasta'

# <codecell>

!wget {latest_genedb_fasta_url}

# <codecell>

!wget {latest_plasmodb_fasta_url}

# <codecell>

!gunzip Pfalciparum.genome.fasta.gz

# <codecell>

my_ref = SeqIO.to_dict(SeqIO.parse(open(my_fasta), "fasta"))

# <codecell>

for chromosome in my_ref:
    print(chromosome, len(my_ref[chromosome]))

# <codecell>

genedb_aug2015_ref = SeqIO.to_dict(SeqIO.parse(open('Pfalciparum.genome.fasta'), "fasta"))

# <codecell>

for chromosome in genedb_aug2015_ref:
    print(chromosome, len(genedb_aug2015_ref[chromosome]))

# <codecell>

plasmodb_25_ref = SeqIO.to_dict(SeqIO.parse(open('PlasmoDB-25_Pfalciparum3D7_Genome.fasta'), "fasta"))
for chromosome in plasmodb_25_ref:
    print(chromosome, len(plasmodb_25_ref[chromosome]))

# <codecell>

old_ref = SeqIO.to_dict(SeqIO.parse(open(old_fasta), "fasta"))
for chromosome in old_ref:
    print(chromosome, len(old_ref[chromosome]))

# <codecell>

str(old_ref['PFC10_API_IRAB'].seq)

# <codecell>

str(genedb_aug2015_ref['Pf3D7_API_v3'].seq)

# <codecell>

str(my_ref['PF_apicoplast_genome_1'].seq)

# <codecell>


