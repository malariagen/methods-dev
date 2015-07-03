# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

thomas_example_gff_fn = '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/embl_V1/Fasta_inChromosomes/Reference/3D7_PfGui.Pf3D7_04_v3.Mutations.gff'
output_dir = '/nfs/team112_internal/production_files/Pf3k/methods/truth_sets'
!mkdir {output_dir}
out_vcf_fn = "%s/3D7_PfGui.Pf3D7_04_v3.Mutations.vcf" % output_dir
left_aligned_vcf_fn = "%s/3D7_PfGui.Pf3D7_04_v3.Mutations.leftAligned.vcf" % output_dir

REF_GENOME = '/data/plasmodium/pfalciparum/recon/roamato/Pf3D7_v3/3D7_sorted.fa'

gatk_exe = 'java -jar ' + '../../opt/gatk/GenomeAnalysisTK.jar'
# gatk_exe = '../../opt/gatk/GenomeAnalysisTK.jar'

# <codecell>

in_seq_handle = open(REF_GENOME)
ref_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

# <codecell>

print(str(ref_dict['Pf3D7_04_v3'][116602:116603].seq))
print(str(ref_dict['Pf3D7_04_v3'][116601]))
print(str(ref_dict['Pf3D7_04_v3'][116602]))
print(str(ref_dict['Pf3D7_04_v3'][116603]))
print()
print(str(ref_dict['Pf3D7_04_v3'][117276]))
print(str(ref_dict['Pf3D7_04_v3'][117277]))
print(str(ref_dict['Pf3D7_04_v3'][117278]))
print(str(ref_dict['Pf3D7_04_v3'][117279]))
117278

# <codecell>

def convert_gff(gff_fn=thomas_example_gff_fn, vcf_fn=out_vcf_fn):
#     Write VCF header
    fo = open(vcf_fn, 'w')
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("##description=This file created with 20150629_convert_gff_to_vcf.ipynb\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
#     Write data
    tbl_gff = (etl.fromtsv(gff_fn))
    current_state = "SNP"
    ins_sequence = ''
    for rec in tbl_gff:
        (chrom, bba, variant_type, pos, pos2, zero, strand, dot, note) = rec
        pos = int(pos)
        pos2 = int(pos2)
        zero = int(zero)
        if (
            bba != 'BBA' or
            not(variant_type in ['SNP', 'Del', 'Ins', 'Synteny']) or
            (variant_type != 'Synteny' and pos != pos2) or
            zero != 0 or
            dot != '.' or
            not(note.startswith('note'))
        ):
            return("error")
        if variant_type == 'SNP':
            if current_state == 'Del':
                ref = str(ref_dict[chrom][variant_start-1:variant_end].seq)
                alt = ref_dict[chrom][variant_start-1]
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
            if current_state == 'Ins':
                ref = ref_dict[chrom][variant_start-1]
                alt = ref_dict[chrom][variant_start-1] + ins_sequence
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                ins_sequence = ''
            ref = ref_dict[chrom][pos-1]
            alt = note[20]
            fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, pos, ref, alt))
            current_state = "SNP"
        elif variant_type == 'Del':
            if current_state == 'Del':
                if pos > (variant_end+1): # i.e. we have moved into a different del so need to print out previous
                    ref = str(ref_dict[chrom][variant_start-1:variant_end].seq)
                    alt = ref_dict[chrom][variant_start-1]
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                    variant_start = pos-1
            if current_state == 'Ins':
                ref = ref_dict[chrom][variant_start-1]
                alt = ref_dict[chrom][variant_start-1] + ins_sequence
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                variant_start = pos-1
            if current_state == 'SNP':
                variant_start = pos-1
            variant_end = pos
            current_state = "Del"
        elif variant_type == 'Ins':
            if current_state == 'Del':
                ref = str(ref_dict[chrom][variant_start-1:variant_end].seq)
                alt = ref_dict[chrom][variant_start-1]
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                variant_start = pos
            if current_state == 'Ins':
                if pos > (variant_end+1): # i.e. we have moved into a different del so need to print out previous
                    ref = ref_dict[chrom][variant_start-1]
                    alt = ref_dict[chrom][variant_start-1] + ins_sequence
                    fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                    variant_start = pos
                    ins_sequence = ''
            if current_state == 'SNP':
                variant_start = pos
            ins_sequence = ins_sequence + note[26]
            variant_end = pos
            current_state = "Ins"
        elif variant_type == 'Synteny':
            if current_state == 'Del':
                ref = str(ref_dict[chrom][variant_start-1:variant_end].seq)
                alt = ref_dict[chrom][variant_start]
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
            if current_state == 'Ins':
                ref = ref_dict[chrom][variant_start-1]
                alt = ref_dict[chrom][variant_start-1] + ins_sequence
                fo.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\n" % (chrom, variant_start, ref, alt))
                ins_sequence = ''
            current_state = "Synteny"
        else:
            return("error")
    fo.close()
    
    return(vcf_fn)
            
                
        
    

# <codecell>

convert_gff()

# <codecell>

thomas_example_gff_fn

# <codecell>

!{gatk_exe} -T LeftAlignAndTrimVariants \
    -R {REF_GENOME} \
    -V {out_vcf_fn} \
    -o {left_aligned_vcf_fn}

# <codecell>

#     --unsafe LENIENT_VCF_PROCESSING

# <codecell>

!echo "{gatk_exe}"

# <codecell>


