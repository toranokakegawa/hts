#!/usr/bin/python -tt


import os
import re
import sys
from detect import detect
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP

def viral() #recherche par mot(c)-clé(s)
print "tape your key words separated by an underscore"
key_word = raw_input() 
for seq_record in SeqIO.parse("gbvrl_nucleic.fst", "fasta"): 
    seq_id=seq_record.id
    string_id=seq_id.tostring()
    if string_id == key_word
     	output = open("papillo.fasta", "w")
    	output.write(str(seq_record.id) + "\n")
    	output.write(str(seq_record.seq) + "\n")

output.close()

def trim_reads() #parser un fichier fastq et trimmer en fonction  d'un argument passé en paramètre
records = (rec[:sys.argv[2] for rec in SeqIO.parse(sys.argv[1], "fastq"))
handle = open("trimmed_data.fastq", "w")
count = SeqIO.write =(records, handle, "fastq")
print( "saved %i reads" % count)
handle.close()

def phred()#parser un fastq file et retenir les séquences au dessus d'un seuil 
records = (rec for rec in \
	   SeqIO.parse(sys.argv[1], "fastq") \
	   if min(rec.letter_annotations[ "phred_quality"]) >= 20)
count = SeqIO.write(records, "quality_seq.fastq", "fastq")
print( "%i reads with a phred score above 20" % count)


def convertTofasta()#convertir un fichier fastq en fasta

def bowtie2()
print("tape the file name of the reference genome and its extension")
file_name = raw_input()
os.system("bowtie2-build file_name index")
os.system("bowtie2 index trimmed_data.fastq > output.sam")#portabilité le choix duc fichier à mapper
os.system("samtools faidx file_name")
os.system("samtools import file_name.fai output.sam output.bam")
os.system("samtools sorted output.bam output.sorted")
os.system("samtools sorted index output.sorted.bam")
os.system("samtools tview output.sorted.bam file_name")

def main():
viral()

    

