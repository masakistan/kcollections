"""``python -m kcollections migrate SRC DST`` — rewrite v1 archives as v2."""

from __future__ import annotations

import argparse
import sys

from .migrate import migrate_archive, probe_archive


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m kcollections")
    sub = parser.add_subparsers(dest="command", required=True)

    migrate_p = sub.add_parser("migrate", help="convert a v1 (2.2) archive to v2")
    migrate_p.add_argument("src", help="input .kc archive")
    migrate_p.add_argument("dst", help="output path (v2)")
    migrate_p.add_argument(
        "--kdict-type",
        default="int",
        choices=["int", "float", "bool", "str"],
        help="Kdict value type when migrating Kdict archives",
    )

    probe_p = sub.add_parser("probe", help="print archive version and container kind")
    probe_p.add_argument("path")

    args = parser.parse_args(argv)

    if args.command == "probe":
        version, kind = probe_archive(args.path)
        print(f"version={version} kind={kind}")
        return 0

    type_map = {
        "int": int,
        "float": float,
        "bool": bool,
        "str": str,
    }
    kind = migrate_archive(
        args.src,
        args.dst,
        kdict_value_type=type_map[args.kdict_type],
    )
    print(f"migrated {kind} -> {args.dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
