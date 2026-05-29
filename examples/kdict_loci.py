#!/usr/bin/env python3
"""Map each k-mer to start positions in a sequence (research pattern)."""

from kcollections import Kdict

K = 15
DNA = "AAACTGTCTTCAAACTGTCTTT"


kd = Kdict((list, int), K)
for i in range(len(DNA) - K + 1):
    kmer = DNA[i : i + K]
    if kmer in kd:
        kd[kmer].append(i)
    else:
        kd[kmer] = [i]

for kmer, positions in kd.items():
    print(f"{kmer}: {positions}")
