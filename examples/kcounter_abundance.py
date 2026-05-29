#!/usr/bin/env python3
"""Count k-mer abundances in a sequence."""

from kcollections import Kcounter

K = 21
DNA = "AAACTGTCTTCAAACTGTCTTTAAACTGTCTTC"

kc = Kcounter(K)
kc.add_seq(DNA)

print(f"Distinct {K}-mers: {len(kc)}")
for kmer, count in kc.most_common(5):
    print(f"  {kmer}: {count}")
