"""Internal trie introspection for research and debugging (unstable API)."""

from __future__ import annotations

from typing import Any


def inspect(container: Any) -> dict[str, Any]:
    """Return low-level trie statistics. API may change without notice."""
    stats = dict(container._trie_stats())
    stats["k"] = container.k
    stats["size"] = len(container)
    return stats
