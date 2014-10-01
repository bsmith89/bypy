#!/usr/bin/env python

# Convert FASTQ formatted files to FASTA (hopefully) quickly

import sys
from Bio import SeqIO

def main():
    for record in SeqIO.parse(sys.argv[1], 'fastq'):
        SeqIO.write(record, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
