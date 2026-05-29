"""Core tests for kcollections (no BioPython required)."""

import os
import tempfile

import pytest

import kcollections
from kcollections import Kcounter, Kdict, Kset


K = 11
KMER_A = "AAACTGTCTTC"
KMER_B = "AAACTGTCTTT"
SEQ = "AAACTGTCTTCAAACTGTCTTT"


@pytest.fixture
def kset():
    return Kset(K)


@pytest.fixture
def kdict():
    return Kdict(str, K)


@pytest.fixture
def kcounter():
    return Kcounter(K)


class TestKset:
    def test_add_and_contains(self, kset):
        kset.add(KMER_A)
        assert KMER_A in kset
        assert KMER_B not in kset

    def test_add_seq(self, kset):
        kset.add_seq(SEQ)
        assert len(kset) == len(SEQ) - K + 1

    def test_iteration(self, kset):
        kset.add(KMER_A)
        kset.add(KMER_B)
        assert set(kset) == {KMER_A, KMER_B}

    def test_remove(self, kset):
        kset.add(KMER_A)
        del kset[KMER_A]
        assert KMER_A not in kset

    def test_serialization_roundtrip(self, kset):
        kset.add_seq(SEQ)
        with tempfile.NamedTemporaryFile(delete=False) as fh:
            path = fh.name
        try:
            kset.save(path)
            other = Kset.from_file(path)
            assert len(other) == len(kset)
            for kmer in kset:
                assert kmer in other
        finally:
            os.unlink(path)

    def test_set_algebra(self, kset):
        kset.add(KMER_A)
        other = Kset(K)
        other.add(KMER_B)
        assert kset.isdisjoint(other)
        union = kset | other
        assert KMER_A in union and KMER_B in union
        inter = kset & other
        assert len(inter) == 0

    def test_intersection_update(self, kset):
        kset.add(KMER_A)
        kset.add(KMER_B)
        other = {KMER_A}
        kset.intersection_update(other)
        assert list(kset) == [KMER_A]

    def test_symmetric_difference_update(self, kset):
        kset.add(KMER_A)
        kset.symmetric_difference_update([KMER_B])
        assert KMER_A in kset and KMER_B in kset
        kset.symmetric_difference_update([KMER_A])
        assert KMER_A not in kset and KMER_B in kset

    @pytest.mark.skip(reason="parallel_add consumer sync regressed after pthread→std::thread migration")
    def test_parallel_add(self, kset):
        kset.parallel_add_init(4)
        kset.parallel_add_seq(SEQ)
        kset.parallel_add_join()
        assert len(kset) == len(SEQ) - K + 1


class TestKdict:
    def test_set_get(self, kdict):
        kdict[KMER_A] = "alpha"
        assert kdict[KMER_A] == "alpha"

    def test_missing_key(self, kdict):
        with pytest.raises(KeyError):
            _ = kdict[KMER_B]

    def test_pop_default(self, kdict):
        assert kdict.pop(KMER_B, "missing") == "missing"

    def test_items(self, kdict):
        kdict[KMER_A] = "x"
        assert list(kdict.items()) == [(KMER_A, "x")]


class TestKcounter:
    def test_increment(self, kcounter):
        kcounter[KMER_A] += 1
        kcounter[KMER_A] += 2
        assert kcounter[KMER_A] == 3

    def test_missing_returns_zero(self, kcounter):
        assert kcounter[KMER_B] == 0

    def test_add_seq(self, kcounter):
        kcounter.add_seq(SEQ)
        assert kcounter[KMER_A] >= 1

    def test_most_common(self, kcounter):
        kcounter[KMER_A] = 5
        kcounter[KMER_B] = 1
        assert kcounter.most_common(1)[0] == (KMER_A, 5)


class TestTextIO:
    def test_export_import_kset(self, kset, tmp_path):
        kset.add_seq(SEQ)
        path = tmp_path / "kmers.txt"
        n = kcollections.export_kmers(kset, str(path))
        assert n == len(kset)
        other = Kset(K)
        assert kcollections.import_kmers(other, str(path)) == n
        assert set(other) == set(kset)

    def test_export_import_kcounter(self, kcounter, tmp_path):
        kcounter[KMER_A] = 3
        kcounter[KMER_B] = 1
        path = tmp_path / "counts.txt"
        kcollections.export_kmers(kcounter, str(path))
        other = Kcounter(K)
        kcollections.import_kmers(other, str(path))
        assert other[KMER_A] == 3
        assert other[KMER_B] == 1


class TestImports:
    def test_version(self):
        assert kcollections.__version__ == "3.2.0"

    def test_serialization_format_constant(self):
        assert kcollections.SERIALIZATION_FORMAT == "kcollections-v2"

    def test_kdict_rejects_nested_types(self):
        with pytest.raises(TypeError):
            Kdict((list, list), 10)

    def test_kdict_from_file(self, tmp_path):
        kd = Kdict(int, K)
        kd[KMER_A] = 7
        path = tmp_path / "d.kc"
        kd.save(str(path))
        kd2 = kcollections.kdict_from_file(int, str(path))
        assert kd2[KMER_A] == 7
