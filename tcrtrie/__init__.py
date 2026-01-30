from __future__ import annotations

import pathlib

from ._tcrtrie import Trie, AIRREntity
from .lazy import LazyObject
from .vdjdb_cache import (
    cache_root,
    read_cached_latest_tag,
    warmup_vdjdb_cache,
    get_cached_airr_path,
)

__all__ = ["Trie", "AIRREntity", "VDJdb"]


def _build_vdjdb_trie() -> Trie:
    root = cache_root()

    try:
        airr = warmup_vdjdb_cache(root=root)
        return Trie(str(airr))
    except Exception:
        pass

    tag = read_cached_latest_tag(root)
    if tag:
        airr = get_cached_airr_path(root, tag)
        if airr.exists():
            return Trie(str(airr))

    raise RuntimeError(
        "VDJdb cache is empty and update check/download failed. "
        "Connect to network and call VDJdb again (or run tcrtrie-vdjdb-update)."
    )


VDJdb = LazyObject(_build_vdjdb_trie)
