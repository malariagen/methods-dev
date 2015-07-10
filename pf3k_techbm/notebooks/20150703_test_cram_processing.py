# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run _shared_setup.ipynb

# <codecell>

cram_fn = '/nfs/pathogen003/tdo/Pfalciparum/PF3K/Reference12Genomes/Sharing_25Jun2015/cram/16503_1#1.cram'
fastq_fn = '/tmp/16503_1#1.fastq'
_100bp_fastq_fn = '/tmp/16503_1#1_100bp.fastq'

# <codecell>

# !{samtools_exe} bamshuf -uOn 128 {cram_fn} tmp | {samtools_exe} bam2fq - > {fastq_fn}

# <codecell>

!{samtools_exe} view -b {cram_fn} | {samtools_exe} bamshuf -uOn 128 - tmp | {samtools_exe} bam2fq - > {fastq_fn}

# <codecell>

record = next(SeqIO.parse(fastq_fn, "fastq"))

# <codecell>

len(record)

# <codecell>

print(record[0:100].format("fastq"))

# <codecell>

record.name.replace('#', '_first100#')

# <codecell>

type(record)

# <codecell>

record_first100 = record[0:100]
# record_first100.id = record_first100.id
record_first100.id = record_first100.id.replace('#', '_first100#')
print(record_first100.id)
print(record_first100.name)
print(record_first100.description)
record_first100.name = ''
record_first100.description = ''
# print(record_first100.name.replace('#', '_first100#'))
print(record_first100.format("fastq"))

# <codecell>

fo = open(_100bp_fastq_fn, "w")

iterator = SeqIO.parse(fastq_fn, "fastq")

for firstRead in iterator:
    firstRead_first100 = firstRead[0:100]
    firstRead_first100.id = firstRead_first100.id.replace('#', '_first100#')
    firstRead_first100.description = ''
    firstRead_last100 = firstRead[150:250]
    firstRead_last100.id = firstRead_last100.id.replace('#', '_last100#')
    firstRead_last100.description = ''
    secondRead = next(iterator)
    secondRead_first100 = secondRead[0:100]
    secondRead_first100.id = secondRead_first100.id.replace('#', '_first100#')
    secondRead_first100.description = ''
    secondRead_last100 = secondRead[150:250]
    secondRead_last100.id = secondRead_last100.id.replace('#', '_last100#')
    secondRead_last100.description = ''
    SeqIO.write(firstRead_first100, fo, "fastq")
    SeqIO.write(secondRead_first100, fo, "fastq")
    SeqIO.write(firstRead_last100, fo, "fastq")
    SeqIO.write(secondRead_last100, fo, "fastq")

fo.close()

# <codecell>

%%file /tmp/first_last_100bp.py
import sys
from Bio import SeqIO

iterator = SeqIO.parse(sys.stdin, "fastq")

for firstRead in iterator:
    firstRead_first100 = firstRead[0:100]
    firstRead_first100.id = firstRead_first100.id.replace('#', '_first100#')
    firstRead_first100.description = ''
    firstRead_last100 = firstRead[150:250]
    firstRead_last100.id = firstRead_last100.id.replace('#', '_last100#')
    firstRead_last100.description = ''
    secondRead = next(iterator)
    secondRead_first100 = secondRead[0:100]
    secondRead_first100.id = secondRead_first100.id.replace('#', '_first100#')
    secondRead_first100.description = ''
    secondRead_last100 = secondRead[150:250]
    secondRead_last100.id = secondRead_last100.id.replace('#', '_last100#')
    secondRead_last100.description = ''
    SeqIO.write(firstRead_first100, sys.stdout, "fastq")
    SeqIO.write(secondRead_first100, sys.stdout, "fastq")
    SeqIO.write(firstRead_last100, sys.stdout, "fastq")
    SeqIO.write(secondRead_last100, sys.stdout, "fastq")


# <codecell>

!cat /tmp/16503_1#1.fastq | python /tmp/first_last_100bp.py > /tmp/16503_1#1_100bp_with_script.fastq

# <codecell>


