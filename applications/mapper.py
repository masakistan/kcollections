#!/usr/bin/env python3
"""Minimal k-mer read mapping example using Kdict."""

import argparse
from collections import Counter

from kcollections import Kdict

try:
    from Bio import SeqIO
except ImportError as exc:
    raise SystemExit(
        "BioPython is required for this example. Install with: pip install biopython"
    ) from exc

KMER_SIZE = 19


def get_kmers(seq, k=KMER_SIZE):
    for i in range(len(seq) - k + 1):
        yield str(seq[i : i + k])


def rev_comp(kmer):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[bp] for bp in kmer[::-1])


def build_reference(reference, k=KMER_SIZE):
    ref_dict = Kdict((list, str), k)
    with open(reference) as ref_fh:
        for chromosome in SeqIO.parse(ref_fh, "fasta"):
            chromosome_id = chromosome.id.split()[0]
            print(f"@HD\tVN:1.6\tSO:unsorted")
            print(f"@SQ\tSN:{chromosome_id}\tLN:{len(chromosome.seq)}")
            build_chromosome(chromosome, chromosome_id, ref_dict)
    return ref_dict


def build_chromosome(chromosome, chromosome_id, ref_dict, k=KMER_SIZE):
    locus = 1
    for kmer in get_kmers(chromosome.seq):
        kmer = kmer.upper()
        if kmer not in ref_dict:
            ref_dict[kmer] = []
        ref_dict[kmer].append(f"{chromosome_id}:{locus}")
        locus += 1


def main():
    parser = argparse.ArgumentParser(description="Simple k-mer mapper demo")
    parser.add_argument("reference", help="Reference FASTA")
    parser.add_argument("reads", help="Reads FASTA/FASTQ")
    args = parser.parse_args()
    ref = build_reference(args.reference)
    print(f"Indexed {len(ref)} unique {KMER_SIZE}-mers")


if __name__ == "__main__":
    main()
