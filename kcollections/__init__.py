"""Memory-efficient k-mer collections built on a Bloom Filter Trie."""

from __future__ import annotations

from typing import Any, Callable, Iterable, Iterator, Optional, Union

from . import _kcollections as _kc
from ._compat import SERIALIZATION_FORMAT

_kdict_mod = _kc
KsetParent = _kc.Kset
KcounterParent = _kc.Kcounter

__version__ = "3.3.1"

__all__ = [
    "Kset",
    "Kdict",
    "Kcounter",
    "kdict_from_file",
    "export_kmers",
    "import_kmers",
    "SERIALIZATION_FORMAT",
    "__version__",
]

_SCALAR_TYPES: dict[type, str] = {
    int: "int",
    float: "float",
    bool: "bool",
    str: "str",
}


class _Persistent:
    def save(self, path: str) -> None:
        self.write(path)  # type: ignore[misc]

    def load(self, path: str) -> None:
        self.read(path)  # type: ignore[misc]

    @classmethod
    def from_file(cls, path: str):
        inst = cls()
        inst.load(path)
        return inst


def _resolve_kdict_class(val_type: Union[type, tuple]) -> type:
    if isinstance(val_type, type):
        if val_type not in _SCALAR_TYPES:
            raise TypeError(
                f"unsupported Kdict value type {val_type!r}; "
                "use int, float, bool, str, or (list, T)"
            )
        name = _SCALAR_TYPES[val_type]
        return getattr(_kdict_mod, f"Kdict_{name}")
    if not (isinstance(val_type, tuple) and len(val_type) == 2 and val_type[0] is list):
        raise TypeError(
            "Kdict collection values must be (list, T) with T in int, float, bool, str"
        )
    elem = val_type[1]
    if elem not in _SCALAR_TYPES:
        raise TypeError(f"unsupported list element type {elem!r}")
    return getattr(_kdict_mod, "Kdict_vector_" + _SCALAR_TYPES[elem])


def _resolve_casters(
    val_type: Union[type, tuple],
) -> tuple[Optional[Callable], Optional[Callable], Optional[Callable]]:
    if isinstance(val_type, type):
        return None, None, None
    elem = val_type[1]
    parts = ["vector", _SCALAR_TYPES[elem]]
    caster_name = "o" + "_".join(parts)
    seq_caster_name = "ovector_" + "_".join(parts)
    caster = getattr(_kdict_mod, caster_name, None)
    seq_caster = getattr(_kdict_mod, seq_caster_name, None)
    return caster, seq_caster, elem


def create_kdict(base: type) -> type:
    class tkdict(_Persistent, base):
        def __init__(
            self,
            k: int = 0,
            caster: Optional[Callable] = None,
            seq_caster: Optional[Callable] = None,
            rcaster: Optional[Callable] = None,
        ):
            super().__init__(k)
            self.caster = caster
            self.rcaster = rcaster
            self.seq_caster = seq_caster

        def __str__(self) -> str:
            return "{" + ",".join(f"{key}:{val}" for key, val in self.items()) + "}"

        def __repr__(self) -> str:
            return self.__str__()

        def _cast_value(self, val: Any) -> Any:
            if self.rcaster is not None:
                try:
                    return self.rcaster(val)
                except (TypeError, ValueError, AttributeError):
                    pass
            return val

        def __setitem__(self, key: str, val: Any) -> None:
            if self.caster is not None:
                val = self.caster(val)
            super().__setitem__(key, val)

        def items(self) -> Iterator[tuple[str, Any]]:
            for kmer, val in self.__iter__():
                yield kmer, self._cast_value(val)

        def keys(self) -> Iterator[str]:
            for kmer, _ in self.__iter__():
                yield kmer

        def values(self) -> Iterator[Any]:
            for _, val in self.__iter__():
                yield self._cast_value(val)

        def copy(self) -> "tkdict":
            new_kdict = type(self)(self.k, self.caster, self.seq_caster, self.rcaster)
            for kmer in self:
                new_kdict[kmer] = self[kmer]
            return new_kdict

        def get(self, key: str, default: Any = None) -> Any:
            if key in self:
                return self._cast_value(super().__getitem__(key))
            return default

        def popitem(self) -> tuple[str, Any]:
            kmer, item = next(self.items())
            del self[kmer]
            return kmer, item

        def setdefault(self, key: str, value: Any = None) -> Any:
            if key not in self:
                self[key] = value
                return value
            return self._cast_value(self[key])

        def pop(self, key: str, *default: Any) -> Any:
            if key in self:
                value = self._cast_value(self[key])
                del self[key]
                return value
            if default:
                return default[0]
            raise KeyError(key)

        def update(self, *others: Iterable) -> None:
            for other in others:
                for item in other:
                    try:
                        key, val = item
                    except (TypeError, ValueError):
                        key = item
                        val = other[key]
                    try:
                        self[key] = self.caster(val) if self.caster else val
                    except (TypeError, ValueError):
                        self[key] = val

        def parallel_add_seq(self, seq: str, values: Iterable) -> None:
            if self.seq_caster and self.caster:
                values = self.seq_caster(map(self.caster, values))
            super().parallel_add_seq(seq, values)

        def add_seq(self, seq: str, values: Iterable) -> None:
            if self.seq_caster and self.caster:
                values = self.seq_caster(map(self.caster, values))
            super().add_seq(seq, values)

    return tkdict


def Kdict(val_type: Union[type, tuple], k: int) -> Any:
    """Build a k-mer dictionary (scalar or ``(list, T)`` values)."""
    base_cls = _resolve_kdict_class(val_type)
    caster, seq_caster, rcaster = _resolve_casters(val_type)
    wrapper = create_kdict(base_cls)
    return wrapper(k, caster, seq_caster, rcaster)


def kdict_from_file(val_type: Union[type, tuple], path: str) -> Any:
    kd = Kdict(val_type, 0)
    kd.load(path)
    return kd


class Kset(_Persistent, KsetParent):
    def __init__(self, k: int = 0):
        super().__init__(k)

    def __str__(self) -> str:
        return "{" + ",".join(self) + "}"

    def __repr__(self) -> str:
        return self.__str__()

    def copy(self) -> "Kset":
        new_set = Kset(self.k)
        for kmer in self:
            new_set.add(kmer)
        return new_set

    def update(self, *iters: Iterable[str]) -> "Kset":
        for _iter in iters:
            for item in _iter:
                self.add(item)
        return self

    def discard(self, item: str) -> None:
        if item in self:
            del self[item]

    def pop(self) -> str:
        kmer = next(iter(self))
        del self[kmer]
        return kmer

    def isdisjoint(self, other: Iterable[str]) -> bool:
        return not any(kmer in self for kmer in other)

    def issubset(self, other: Iterable[str]) -> bool:
        return all(kmer in other for kmer in self)

    def issuperset(self, other: Iterable[str]) -> bool:
        return all(kmer in self for kmer in other)

    def intersection(self, *other_sets: Iterable[str]) -> "Kset":
        new_set = Kset(self.k)
        for kmer in self:
            if all(kmer in other for other in other_sets):
                new_set.add(kmer)
        return new_set

    def intersection_update(self, *other_sets: Iterable[str]) -> "Kset":
        for kmer in list(self):
            if not all(kmer in other for other in other_sets):
                del self[kmer]
        return self

    def difference(self, other: Iterable[str]) -> "Kset":
        new_set = Kset(self.k)
        for kmer in self:
            if kmer not in other:
                new_set.add(kmer)
        return new_set

    def difference_update(self, *other_sets: Iterable[str]) -> "Kset":
        for other_set in other_sets:
            for kmer in other_set:
                self.discard(kmer)
        return self

    def symmetric_difference(self, other_set: Iterable[str]) -> "Kset":
        new_set = Kset(self.k)
        for kmer in self:
            if kmer not in other_set:
                new_set.add(kmer)
        for kmer in other_set:
            if kmer not in self:
                new_set.add(kmer)
        return new_set

    def symmetric_difference_update(self, other_set: Iterable[str]) -> "Kset":
        for kmer in other_set:
            if kmer in self:
                del self[kmer]
            else:
                self.add(kmer)
        return self

    def union(self, *other_sets: Iterable[str]) -> "Kset":
        new_set = self.copy()
        for other_set in other_sets:
            new_set.update(other_set)
        return new_set

    def __and__(self, other: Iterable[str]) -> "Kset":
        return self.intersection(other)

    def __or__(self, other: Iterable[str]) -> "Kset":
        return self.union(other)

    def __xor__(self, other: Iterable[str]) -> "Kset":
        return self.symmetric_difference(other)

    def __sub__(self, other: Iterable[str]) -> "Kset":
        return self.difference(other)

    def __iand__(self, other: Iterable[str]) -> "Kset":
        return self.intersection_update(other)

    def __ior__(self, other: Iterable[str]) -> "Kset":
        return self.update(other)

    def __ixor__(self, other: Iterable[str]) -> "Kset":
        return self.symmetric_difference_update(other)

    def __isub__(self, other: Iterable[str]) -> "Kset":
        return self.difference_update(other)


class Kcounter(_Persistent, KcounterParent):
    def __init__(self, k: int = 0):
        super().__init__(k)

    def __getitem__(self, key: str) -> int:
        if key in self:
            return super().__getitem__(key)
        return 0

    def __str__(self) -> str:
        return "{" + ",".join(f"{key}:{val}" for key, val in self.items()) + "}"

    def __repr__(self) -> str:
        return self.__str__()

    def items(self) -> Iterator[tuple[str, int]]:
        for kmer, val in self.__iter__():
            yield kmer, val

    def keys(self) -> Iterator[str]:
        for kmer, _ in self.__iter__():
            yield kmer

    def values(self) -> Iterator[int]:
        for _, val in self.__iter__():
            yield val

    def copy(self) -> "Kcounter":
        new_kcounter = Kcounter(self.k)
        for kmer in self:
            new_kcounter[kmer] = self[kmer]
        return new_kcounter

    def get(self, key: str, default: int = 0) -> int:
        if key in self:
            return self[key]
        return default

    def popitem(self) -> tuple[str, int]:
        kmer, item = next(self.items())
        del self[kmer]
        return kmer, item

    def setdefault(self, key: str, value: int = 0) -> int:
        if key not in self:
            self[key] = value
            return value
        return self[key]

    def pop(self, key: str, *default: Any) -> int:
        if key in self:
            value = self[key]
            del self[key]
            return value
        if default:
            return default[0]
        raise KeyError(key)

    def update(self, *others: Iterable) -> None:
        for other in others:
            for item in other:
                try:
                    key, val = item
                except (TypeError, ValueError):
                    key = item
                    val = other[key]
                self[key] = val

    def most_common(self, n: Optional[int] = None) -> list[tuple[str, int]]:
        items = sorted(self.items(), key=lambda x: x[1], reverse=True)
        if n is None:
            return items
        return items[:n]


def _export_kmers(container: Any, path: str, *, include_counts: bool = False) -> int:
    """Write one k-mer per line (optional count column for Kcounter)."""
    n = 0
    with open(path, "w", encoding="ascii") as fh:
        if include_counts:
            for kmer, count in container.items():
                fh.write(f"{kmer}\t{count}\n")
                n += 1
        elif hasattr(container, "items"):
            for kmer, _ in container.items():
                fh.write(f"{kmer}\n")
                n += 1
        else:
            for kmer in container:
                fh.write(f"{kmer}\n")
                n += 1
    return n


def export_kmers(container: Any, path: str) -> int:
    return _export_kmers(container, path, include_counts=isinstance(container, Kcounter))


def import_kmers(container: Any, path: str, *, clear: bool = False) -> int:
    """Load k-mers from a text file (one per line; optional tab count)."""
    if clear:
        container.clear()
    n = 0
    with open(path, encoding="ascii") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                kmer, val = line.split("\t", 1)
                if isinstance(container, Kcounter):
                    container[kmer] = int(val)
                else:
                    container[kmer] = val
            elif isinstance(container, Kset):
                container.add(line)
            elif isinstance(container, Kcounter):
                container[line] = container.get(line, 0) + 1
            else:
                container[line] = True
            n += 1
    return n
