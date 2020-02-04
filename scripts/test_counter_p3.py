#!/usr/bin/env python
import argparse

from kcollections import Kcounter, Kset
from collections import Counter

from Bio import SeqIO

KMER_SIZE = 23

def get_kmers(seq, k=KMER_SIZE):
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i:i + k])
        yield kmer


def test_kcounter(sequence):
    kcounter = Kcounter(KMER_SIZE)
    kset = Kset(KMER_SIZE)
    counter = Counter()
    tk = 'TGCCGGATGCGCTTTGCTTATCC'
    with open(sequence, 'rU') as seq_fh:
        for record in SeqIO.parse(seq_fh, 'fasta'):
            for kmer in get_kmers(record.seq):
                kset.add(kmer)
                kcounter[kmer] += 1
                counter[kmer] += 1
                assert kcounter[tk] == counter[tk], "{} {} {}".format(kcounter[tk], counter[tk], kmer)
    
    for kmer, count in counter.items():
        if kmer not in kcounter:
            print('ERROR: kmer present in Counter, but not in Kcounter')
            continue
        if kcounter[kmer] != count:
            print('ERROR: count mismatch, Counter count:', count, 'Kcounter count:', kcounter[kmer], kmer)

    for kmer, count in kcounter.items():
        if kmer not in counter:
            print('ERROR: kmer present in Kcounter, but not in Counter')
            continue
        if counter[kmer] != count:
            print('ERROR: count mismatch, Kcounter count:', count, 'Counter count:', counter[kmer], kmer)


if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Test Kcounter.')
   parser.add_argument('sequence', type=str, default='./source.fa',
                       help='Path to a sequence file to test on')

   args = parser.parse_args()

   test_kcounter(args.sequence)
