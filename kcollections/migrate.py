"""Migrate on-disk kcollections archives (v1 → v2)."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Tuple, Union

from . import Kcounter, Kdict, Kset

ContainerKind = Literal["kset", "kdict", "kcounter"]
_KIND_BYTE = {1: "kset", 2: "kdict", 3: "kcounter"}


def probe_archive(path: Union[str, Path]) -> Tuple[int, ContainerKind]:
    """Return ``(version, kind)`` for a kcollections binary archive."""
    data = Path(path).read_bytes()
    if len(data) < 7 or data[:4] != b"KCOL":
        raise ValueError(f"not a kcollections archive: {path}")
    version = int.from_bytes(data[4:6], "little")
    kind_byte = data[6]
    if kind_byte not in _KIND_BYTE:
        raise ValueError(f"unknown container kind byte {kind_byte} in {path}")
    if version not in (1, 2):
        raise ValueError(f"unsupported archive version {version} in {path}")
    return version, _KIND_BYTE[kind_byte]  # type: ignore[return-value]


def migrate_archive(
    src: Union[str, Path],
    dst: Union[str, Path],
    *,
    kdict_value_type: type = int,
) -> ContainerKind:
    """Load ``src`` (v1 or v2) and write ``dst`` as serialization v2.

    For ``Kdict`` archives, pass ``kdict_value_type`` (e.g. ``int``, ``str``, or ``(list, int)``).
    """
    src_p = Path(src)
    dst_p = Path(dst)
    version, kind = probe_archive(src_p)

    if version == 2 and src_p.resolve() == dst_p.resolve():
        return kind

    if kind == "kset":
        obj = Kset.from_file(str(src_p))
    elif kind == "kcounter":
        obj = Kcounter.from_file(str(src_p))
    else:
        obj = Kdict(kdict_value_type)
        obj.load(str(src_p))

    obj.save(str(dst_p))
    return kind
