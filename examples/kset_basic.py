#!/usr/bin/env python3
"""Index k-mers from a DNA string and save to disk."""

from kcollections import Kset

K = 21
DNA = (
    "AAACTGTCTTCCTTTATTTGTTCAGGGATCGTGTCAGTA"
    "AAACTGTCTTCCTTTATTTGTTCAGGG"
)

ks = Kset(K)
ks.add_seq(DNA)
print(f"Indexed {len(ks)} unique {K}-mers")

path = "example.kset"
ks.save(path)

ks2 = Kset.from_file(path)
print(f"Reloaded {len(ks2)} k-mers; match={len(ks) == len(ks2)}")
