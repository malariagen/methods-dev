# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run _shared_setup.ipynb

# <codecell>

plasmodb24_fn = '/Users/rpearson/Downloads/PlasmoDB-24_Pfalciparum3D7.gff' # http://www.plasmodb.org/common/downloads/release-24/Pfalciparum3D7/gff/data/PlasmoDB-24_Pfalciparum3D7.gff
genedb_gff_fn = '/Users/rpearson/Downloads/Pfalciparum.gff3' # ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz 
genedb_fasta_fn = '/Users/rpearson/Downloads/Pfalciparum.genome.fasta' # ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.genome.fasta.gz 
snpeff_config_fn = snpeff_dir + '/snpEff.config'

# <codecell>

current_dir = !pwd

# <codecell>

print(current_dir[0])

# <codecell>

# Setup directory for PlasmoDB snpeff database
%cd {current_dir[0]}
!mkdir {snpeff_dir}/data/Pf3D7v24
%cd {snpeff_dir}/data/Pf3D7v24
!cp {plasmodb24_fn} genes.gff

# Setup directory for GeneDB snpeff database
%cd {current_dir[0]}
!mkdir {snpeff_dir}/data/Pf3D7july2015
%cd {snpeff_dir}/data/Pf3D7july2015
!cp {genedb_gff_fn} genes.gff
!cp {genedb_fasta_fn} sequences.fa

# <codecell>

%cd {current_dir[0]}

fo = open(snpeff_config_fn, "a")

print('''
# Pf3D7 PlasmoDB v24. Added by Richard Pearson. Downloaded from http://www.plasmodb.org/common/downloads/release-24/Pfalciparum3D7/gff/data/PlasmoDB-24_Pfalciparum3D7.gff
Pf3D7v24.genome: Plasmodium_falciparum
        Pf3D7v24.M76611.codonTable: Protozoan_Mitochondrial
        Pf3D7v24.PFC10_API_IRAB.codonTable: Bacterial_and_Plant_Plastid
        Pf3D7v24.reference: http://www.plasmodb.org/common/downloads/release-24/Pfalciparum3D7/gff/data/PlasmoDB-24_Pfalciparum3D7.gff
        Pf3D7v24.cds_mal_mito_1-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.cds_mal_mito_2-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.cds_mal_mito_3-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.rna_mal_mito_1-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.rna_mal_mito_2-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.rna_mal_mito_3-1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.mal_mito_1.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.mal_mito_2.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.mal_mito_3.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_10\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_15\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_16\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_17\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.mal_mito_RNA19\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_1\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_20\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_21\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_9\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_LSUC\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_LSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_LSUG\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA11\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA12\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA14\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA18\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA22\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA4\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA5\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA6\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_RNA7\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_rna_SSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7v24.malmito_SSUB\:rRNA.codonTable : Vertebrate_Mitochondrial

# Pf3D7 GeneDB July 2015. Added by Richard Pearson. Downloaded July 8, 2015 from ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz 
Pf3D7july2015.genome: Plasmodium_falciparum
        Pf3D7july2015.Pf_M76611.codonTable: Protozoan_Mitochondrial
        Pf3D7july2015.reference: ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.gff3.gz
        Pf3D7july2015.mal_mito_1.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_2.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_3.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_10\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_15\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_16\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_17\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.mal_mito_RNA19\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_1\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_20\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_21\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_9\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUC\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_LSUG\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA11\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA12\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA14\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA18\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA22\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA4\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA5\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA6\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_RNA7\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_rna_SSUF\:rRNA.codonTable : Vertebrate_Mitochondrial
        Pf3D7july2015.malmito_SSUB\:rRNA.codonTable : Vertebrate_Mitochondrial
''', file=fo)

fo.close()

# <codecell>

fo.close()

# <codecell>

%cd {snpeff_dir}
!java -jar snpEff.jar build -gff3 -v Pf3D7v24 -upDownStreamLen 1 > /dev/null
%cd {current_dir[0]}

# <codecell>

%cd {snpeff_dir}
!java -jar snpEff.jar build -gff3 -v Pf3D7july2015
%cd {current_dir[0]}

# <codecell>


