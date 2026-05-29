"""Compatibility helpers and deprecation policy."""

from __future__ import annotations

import warnings

# Documented on-disk format generation (Boost binary archive; same kcollections major)
SERIALIZATION_FORMAT = "kcollections-v1"

# Next major release where deprecated Python APIs are removed
DEPRECATION_REMOVAL_VERSION = "3.0"


def deprecate(message: str, *, stacklevel: int = 3) -> None:
    warnings.warn(
        f"{message} Will be removed in kcollections {DEPRECATION_REMOVAL_VERSION}.",
        DeprecationWarning,
        stacklevel=stacklevel,
    )
