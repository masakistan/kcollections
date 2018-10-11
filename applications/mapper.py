#!/usr/bin/env python
import argparse, pdb

from kcollections import Kdict

from Bio import SeqIO

KMER_SIZE = 19

def get_kmers(seq, k=KMER_SIZE):
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i:i + k])
        yield kmer
    
def build_reference(reference, k=KMER_SIZE):
    ref_dict = Kdict(k)
    with open(reference, 'rU') as ref_fh:
        for chromosome in SeqIO.parse(ref_fh, 'fasta'):
            build_chromosome(chromosome, ref_dict)


def build_chromosome(chromosome, ref_dict, k=KMER_SIZE):
    chr_id = str(chromosome.id)
    for kmer in get_kmers(chromosome.seq):
        if kmer not in ref_dict:
            ref_dict[kmer] = []
        ref_dict[kmer] += [chr_id + ':::' + str(i)]

        
def map_reads(reads, ref_dict):
    with open(reads, 'rU') as reads_fh:
        loci = []
        for read in SeqIO.parse(reads_fh, 'fastq'):
            for kmer in get_kmers(read.seq):
                if kmer in ref_dict:
                    loci += [ref_dict[kmer]]

                    
def score_loci(loci):
    for locus in loci:
        if len(locus) > 1:
            pass
        else:
            pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Map reads using kcollections.')
    parser.add_argument('reference', type=str, 
                        help='Path to the reference genome in FASTA format')
    parser.add_argument('reads', type=str, nargs='+', 
                        help='One or more paths to reads in FASTQ format')
    
    args = parser.parse_args()
    
    ref_dict = build_reference(args.reference)
    
