#!/usr/bin/env python
import sys,random
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import sample
from subprocess import call 


input_genome = sys.argv[1]
random_int = int(sys.argv[2])  
output_fasta1 = sys.argv[3]
output_fasta2= sys.argv[4]
output_fastq = sys.argv[5]
input_genome2= sys.argv[6]

seq = list(SeqIO.parse(input_genome,"fasta"))             
total_length = len(seq)                                   

outlist1 = []
outlist2 = []
outlist3 = []
outlist4 = []

for i in range(random_int):
  position1 = random.randint(1,total_length-1)               
  for j in range(len(outlist1)):                  
    if position1 == outlist1[j]:
       print "Same sequence ERROR"
       position1= random.randint(1,total_length-1)           
  output_seq = seq[position1]
  outlist2.append(output_seq) 
SeqIO.write(outlist2, output_fasta1, "fasta")

for reads in SeqIO.parse(output_fasta1, "fasta"):	
	id1=reads.id
	seq=reads.seq[1:50]
	seq1=str(seq)
    output_seq1 = id1, seq1
    outlist3.append(output_seq1) 
for (id1,seq1) in outlist3:
	read1 = SeqRecord(Seq(seq1),id=id1,description="")
    outlist4.append(read1)
SeqIO.write(outlist4, output_fasta2, "fasta")

with open(output_fasta2, "r") as fasta, open(output_fastq, "w") as fastq:
    for rec in SeqIO.parse(fasta, "fasta"):
        rec.letter_annotations["phred_quality"] = [40] * len(rec)
        SeqIO.write(rec, fastq, "fastq")


process1 = subprocess.call(['bowtie2-build',  '-q', '-f', input_genome2, input_genome2])
process2 = subprocess.call(['bowtie2', '-x', input_genome2, '-U', output_fastq, 'S', 'Output_sam'])
samfile = pysam.AlignmentFile("Output_sam", "rb")
process3= subprocess.call(['samtools idxstats', 'samfile'])



print "Input Genomic sequence in fasta format =", input_genome
print "Number of random sequences formed =", rand_int
print "Fasta Sequences with random reads = ", output_fasta1
print "Fasta Sequence with random reads of 550 bp length = ", output_fasta2
print "Fastq Sequence formed from Fasta adding random quality values =", output_fastq

















 