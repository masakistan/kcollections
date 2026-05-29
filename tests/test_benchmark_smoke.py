"""Lightweight performance smoke tests (not strict benchmarks)."""

import time

from kcollections import Kcounter, Kset

K = 11
SEQ = "AAACTGTCTTCAAACTGTCTTT"


def test_kset_add_seq_throughput():
    ks = Kset(K)
    t0 = time.perf_counter()
    ks.add_seq(SEQ)
    elapsed = time.perf_counter() - t0
    assert len(ks) == len(SEQ) - K + 1
    assert elapsed < 2.0, f"add_seq too slow: {elapsed:.2f}s"


def test_kcounter_bulk_insert():
    kc = Kcounter(K)
    kmers = [SEQ[i : i + K] for i in range(len(SEQ) - K + 1)]
    t0 = time.perf_counter()
    for kmer in kmers:
        kc[kmer] = kc.get(kmer, 0) + 1
    elapsed = time.perf_counter() - t0
    assert len(kc) == len(kmers)
    assert elapsed < 2.0, f"bulk insert too slow: {elapsed:.2f}s"
