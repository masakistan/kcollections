#!/usr/bin/env python
import argparse, pdb
from collections import Counter

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
            print('@SQ SN:{} LN:{}'.format(chromosome.id.split()[0], len(chromosome.seq)))
            build_chromosome(chromosome, ref_dict)
    return ref_dict


def build_chromosome(chromosome, ref_dict, k=KMER_SIZE):
    locus = 1
    for kmer in get_kmers(chromosome.seq):
        if kmer not in ref_dict:
            ref_dict[kmer] = []
        ref_dict[kmer] += [locus]
        locus += 1


def map_reads(reads, ref_dict):
    with open(reads, 'rU') as reads_fh:
        for read in SeqIO.parse(reads_fh, 'fastq'):
            loci = Counter()
            for offset, kmer in enumerate(get_kmers(read.seq)):
                if kmer in ref_dict:
                    for locus in ref_dict[kmer]:
                        loci[locus - offset] += 1
            score, locus = score_loci(loci)
            print('*\t0\tNC_007779.1\t{}\t0\t*\t*\t0\t{}\t{}\t*'.format(locus, len(read.seq), read.seq))


def score_loci(loci):
    max_score, max_locus = 0, 0
    for locus, score in loci.iteritems():
        if score > max_score:
            max_score, max_locus = score, locus
    return max_score, max_locus


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Map reads using kcollections.')
    parser.add_argument('reference', type=str, default='../scripts/source.fa',
                        help='Path to the reference genome in FASTA format (default: ../scripts/source.fa)')
    parser.add_argument('reads', type=str, nargs='+', default='./reads.fq',
                        help='One or more paths to reads in FASTQ format')

    args = parser.parse_args()

    ref_dict = build_reference(args.reference)

    for reads in args.reads:
        map_reads(reads, ref_dict)
