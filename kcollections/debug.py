"""Internal trie introspection for research and debugging (unstable API)."""

from __future__ import annotations

from typing import Any


def inspect(container: Any) -> dict[str, Any]:
    """Return low-level trie statistics. API may change without notice."""
    root = container.get_root()
    return {
        "k": container.k,
        "size": len(container),
        "root_vs_size": container.get_vs_size(root),
        "root_uc_size": container.get_uc_size(root),
    }
