# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run standard_imports.ipynb

# <codecell>

thomas_genome_fn = '/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011/Pf3D7_v3.fasta'
thomas_api_fn = '/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011/PFC10_API_IRAB_RDP.embl'
plasmodb_24_fn = 'PlasmoDB-24_Pfalciparum3D7_Genome.fasta'

# <codecell>

!wget http://www.plasmodb.org/common/downloads/release-24/Pfalciparum3D7/fasta/data/PlasmoDB-24_Pfalciparum3D7_Genome.fasta

# <codecell>

in_seq_handle = open(thomas_genome_fn)
chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

# <codecell>

in_seq_handle = open(plasmodb_24_fn)
plasmodb_24 = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

# <codecell>

len(chromosomes['PF_apicoplast_genome_1'])

# <codecell>

len(plasmodb_24['PFC10_API_IRAB'])

# <codecell>

in_seq_handle = open(thomas_api_fn)
api = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "embl"))

# <codecell>

len(api['XXX.4'])

# <codecell>


