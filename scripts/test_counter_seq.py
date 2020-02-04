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
    kc = Kcounter(KMER_SIZE)
    counter = Counter()
    tk = 'TGCCGGATGCGCTTTGCTTATCC'
    with open(sequence, 'rU') as seq_fh:
        for record in SeqIO.parse(seq_fh, 'fasta'):
            kcounter.add_seq(str(record.seq))
            for idx, kmer in enumerate(get_kmers(record.seq)):
                counter[kmer] += 1
                #print(idx, idx + KMER_SIZE, kmer)

    print(len(kcounter), len(counter))
    
    for kmer, count in counter.items():
        try:
            if kmer not in kcounter:
                print('ERROR: kmer present in Counter, but not in Kcounter {}'.format(kmer))
                continue
            if kcounter[kmer] != count:
                print('ERROR: count mismatch, Counter count:', count, 'Kcounter count:', kcounter[kmer], kmer)
        except Exception as e:
            print(e)


    for kmer, count in kcounter.items():
        if kmer not in counter:
            print('ERROR: kmer present in Kcounter, but not in Counter {}'.format(kmer))
            continue
        if counter[kmer] != count:
            print('ERROR: count mismatch, Kcounter count:', count, 'Counter count:', counter[kmer], kmer)
    print("done")


if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Test Kcounter.')
   parser.add_argument('sequence', type=str, default='./source.fa',
                       help='Path to a sequence file to test on')

   args = parser.parse_args()

   test_kcounter(args.sequence)
