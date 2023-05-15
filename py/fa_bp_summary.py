#!/usr/bin/env python
# conda activate fa_stats
# python py/fa_bp_summary.py [file.fa(.gz)] [report.yml]
# runtime on gzipped seal_genome ~ 15 - 20min
import sys
import gzip
import io
import os
import time
from Bio import SeqIO
from progress.bar import Bar

filename = sys.argv[1]
yml_out = sys.argv[2]
# prefix = filename.split("/")[-1].split(".")[0]
# yml_out = prefix + ".yml"

if filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
    fh = io.TextIOWrapper(io.BufferedReader(gzip.open(filename)))
elif filename.endswith(".fa") or filename.endswith(".fasta"):
    fh = open(filename,'r')
else:
    raise Exception("Unsupported File Type")

# determine number of seqs in fasta
seq_len = 0
for record in SeqIO.parse(fh, "fasta"): 
    seq_len+=1

c=0
a=0
g=0
t=0
C=0
A=0
G=0
T=0
n=0
total_bases = 0
n_seq=0

if filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
    fh = io.TextIOWrapper(io.BufferedReader(gzip.open(filename)))
elif filename.endswith(".fa") or filename.endswith(".fasta"):
    fh = open(filename,'r')

parser = SeqIO.parse(fh, "fasta")
with Bar('Processing fasta:', max = seq_len, fill = '█') as bar:
    for record in parser:
        n_seq+=1
        # print( "%.0f st seq" % (n_seq)) 
        for x in str(record.seq):
            total_bases+=1
            time.sleep(.01)
            match x:
                case "C":
                    C+=1
                case "G":
                    G+=1
                case "A":
                    A+=1
                case "T":
                    T+=1
                case "N":
                    n+=1
                case "c":
                    c+=1
                case "g":
                    g+=1
                case "a":
                    a+=1
                case "t":
                    t+=1
                case "n":
                    n+=1
        bar.next()

non_n_bp=a+t+g+c+A+T+G+C
unmasked_bp=A+T+G+C
soft_masked_bp=a+t+g+c
gc_content=(g+c+G+C)/non_n_bp
gc_content_unmasked=(G+C)/(A+T+G+C)

with open(yml_out, "a") as f:
    print("genome:", file = f)
    print( "  file_name: '%s'" % (filename) , file = f)
    print( "  n_scaffolds: %.0f" % (n_seq) , file = f)
    print( "  GC_content: %.3f" % (gc_content) , file = f)
    print( "  GC_content_unmasked: %.3f" % (gc_content_unmasked) , file = f)
    print( "  bp_stats:", file = f)
    print( "    'total':", file = f)
    print( "      'total_bases': %.0f" % (total_bases) , file = f)
    print( "      'A': %.0f" % (a+A) , file = f)
    print( "      'T': %.0f" % (t+T) , file = f)
    print( "      'G': %.0f" % (g+G) , file = f)
    print( "      'C': %.0f" % (c+C) , file = f)
    print( "    'masked':", file = f)
    print( "      'a_masked': %.0f" % (a) , file = f)
    print( "      't_masked': %.0f" % (t) , file = f)
    print( "      'g_masked': %.0f" % (g) , file = f)
    print( "      'c_masked': %.0f" % (c) , file = f)
    print( "    'n': %.0f" % (n) , file = f)

