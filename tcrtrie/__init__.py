from __future__ import annotations

import pathlib

from ._tcrtrie import Trie, AIRREntity
from .lazy import LazyObject
from .vdjdb_cache import (
    cache_root,
    read_cached_latest_tag,
    warmup_vdjdb_cache,
    get_cached_airr_path,
    get_cached_sqlite_path,
)
from .vdjdb_client import VDJdbClient

__all__ = ["Trie", "AIRREntity", "VDJdb"]


def _build_vdjdb() -> VDJdbClient:
    root = cache_root()

    try:
        warmup_vdjdb_cache(root=root)
        tag = read_cached_latest_tag(root)
        if not tag:
            raise RuntimeError("VDJdb warmup succeeded but LATEST marker missing.")
        airr = get_cached_airr_path(root, tag)
        sqlite = get_cached_sqlite_path(root, tag)
        if airr.exists() and sqlite.exists():
            return VDJdbClient(trie=Trie(str(airr)), sqlite_path=sqlite)
    except Exception:
        pass

    tag = read_cached_latest_tag(root)
    if tag:
        airr = get_cached_airr_path(root, tag)
        sqlite = get_cached_sqlite_path(root, tag)
        if airr.exists() and sqlite.exists():
            return VDJdbClient(trie=Trie(str(airr)), sqlite_path=sqlite)

    raise RuntimeError(
        "VDJdb cache is empty or incomplete (missing AIRR/SQLite), and update failed. "
        "Connect to network and call VDJdb again (or run tcrtrie-vdjdb-update)."
    )


VDJdb = LazyObject(_build_vdjdb)
