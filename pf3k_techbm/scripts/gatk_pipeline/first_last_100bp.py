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
