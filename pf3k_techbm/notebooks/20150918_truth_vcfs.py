# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

%run _shared_setup.ipynb

# <codecell>

!mkdir -p {os.path.join(PROCESSED_ASSEMBLED_SAMPLES_DIR, 'truth_vcfs_2')}

# <headingcell level=1>

# Check data

# <codecell>

tbl_samples_to_process = (tbl_assembled_samples
    .cutout('Notes')
    .selectnotnone('bam_fn')
    .selecteq('To be used for', 'Validation')
    .addfield('thomas_gff_filestem', lambda rec: os.path.join(
        '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/SNP/Reference',
        "Pf%s" % rec['Isolate code']
    ))
    .addfield('truth_vcf_filestem', lambda rec: os.path.join(
        PROCESSED_ASSEMBLED_SAMPLES_DIR,
        'truth_vcfs_2',
        "truth_%s" % rec['Isolate code']
    ))
    .cut([0, 1, 14, 15])
#     .head(4)
)
tbl_samples_to_process.displayall(index_header=True)

# <codecell>

tbl_samples_to_process.values('thomas_gff_filestem')[0]

# <headingcell level=1>

# Functions

# <codecell>

def convert_gff(gff_fn="%s.Pf3D7_01_v3.Mutations.gff" % tbl_samples_to_process.values('thomas_gff_filestem')[0],
                vcf_fn="%s.Pf3D7_01_v3.vcf" % tbl_samples_to_process.values('truth_vcf_filestem')[0],
                ref_dict=SeqIO.to_dict(SeqIO.parse(open(REF_GENOME), "fasta"))):
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

rewrite=True

in_seq_handle = open(REF_GENOME)
chromosomes = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")).keys()

for chromosome in sort(list(chromosomes)):
    if chromosome == 'PFC10_API_IRAB' or chromosome == 'M76611':
        continue
#     if chromosome == 'PFC10_API_IRAB':
#         gff_chromosome='PF_apicoplast_genome_1' 
#     elif chromosome == 'M76611':
#         gff_chromosome='Pf_M76611' 
#     else:
#         gff_chromosome=chromosome
    print(chromosome)
    for rec in tbl_samples_to_process.data():
        truth_gff_gn = "%s.%s.Mutations.gff" % (rec[2], chromosome)
        unannotated_vcf_fn = "%s.%s.unannotated.vcf" % (rec[3], chromosome)
        left_aligned_vcf_fn = "%s.%s.leftaligned.vcf" % (rec[3], chromosome)
        snpeff_vcf_fn = "%s.%s.snpeff.vcf" % (rec[3], chromosome)
        snpeff_annotated_vcf_fn = "%s.%s.snpeff_annotated.vcf" % (rec[3], chromosome)
        annotated_vcf_fn = "%s.%s.annotated.vcf" % (rec[3], chromosome)
        
#         if os.path.exists("%s.%s.Mutations.gff" % (rec[2], chromosome)):
        if (not os.path.isfile(unannotated_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
            convert_gff(
                gff_fn=truth_gff_gn,
                vcf_fn=unannotated_vcf_fn
            )

        if (not os.path.isfile(left_aligned_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
            !{gatk_exe} -T LeftAlignAndTrimVariants \
            -R {REF_GENOME} \
            -V {unannotated_vcf_fn} \
            -o {left_aligned_vcf_fn} \
            2> /dev/null

        if (not os.path.isfile(snpeff_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
            !{snpeff_exe} \
            -v -o gatk Pf3D7july2015 \
            {left_aligned_vcf_fn} \
            -no-downstream \
            -no-upstream \
            > {snpeff_vcf_fn} \
            2> /dev/null

        if (not os.path.isfile(snpeff_annotated_vcf_fn) and not os.path.isfile(annotated_vcf_fn+'.gz')) or rewrite:
            !{gatk_exe} \
            -T VariantAnnotator \
            -R {REF_GENOME} \
            -A TandemRepeatAnnotator \
            -A SnpEff \
            --variant {left_aligned_vcf_fn} \
            --snpEffFile {snpeff_vcf_fn} \
            -o {snpeff_annotated_vcf_fn} \
            2> /dev/null

        if not os.path.isfile(annotated_vcf_fn+'.gz') or rewrite:
            !cat {snpeff_annotated_vcf_fn} \
            | vcf-annotate -a {regions_fn} \
               -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. SubtelomericRepeat: repetitive regions at the ends of the chromosomes. SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
               -c CHROM,FROM,TO,INFO/RegionType \
            > {annotated_vcf_fn}

            !bgzip -f {annotated_vcf_fn}
            !tabix -p vcf -f {annotated_vcf_fn}.gz

        if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.SNP.vcf')+'.gz') or rewrite:
            !{gatk_exe} \
            -T SelectVariants \
            -R {REF_GENOME} \
            -V {annotated_vcf_fn}.gz \
            -selectType SNP \
            -o {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')} 2> /dev/null

            !bgzip -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}
            !tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.SNP.vcf')}.gz
            
        if not os.path.isfile(annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')+'.gz') or rewrite:
            !{gatk_exe} \
            -T SelectVariants \
            -R {REF_GENOME} \
            -V {annotated_vcf_fn}.gz \
            -xlSelectType SNP \
            -o {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')} 2> /dev/null

            !bgzip -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}
            !tabix -p vcf -f {annotated_vcf_fn.replace('.vcf', '.INDEL.vcf')}.gz

#         !rm {unannotated_vcf_fn}
#         !rm {left_aligned_vcf_fn}
#         !rm {snpeff_vcf_fn}
#         !rm {snpeff_annotated_vcf_fn}
#         !rm {unannotated_vcf_fn}.idx
#         !rm {left_aligned_vcf_fn}.idx
#         !rm {snpeff_vcf_fn}.idx
#         !rm {snpeff_annotated_vcf_fn}.idx

# <codecell>

sort(list(chromosomes))

# <codecell>


