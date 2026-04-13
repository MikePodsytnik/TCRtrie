from __future__ import annotations

from typing import TYPE_CHECKING, cast

from ._tcrtrie import Trie, AIRREntity
from .lazy import LazyObject
from .vdjdb_cache import (
    cache_root,
    read_cached_active_tag,
    get_cached_airr_path,
    get_cached_sqlite_path,
)
from .vdjdb_client import VDJdbClient

__all__ = ["Trie", "AIRREntity", "VDJdbClient", "VDJdb"]


def _build_vdjdb() -> VDJdbClient:
    root = cache_root()
    tag = read_cached_active_tag(root)
    if not tag:
        raise RuntimeError("VDJdb cache is empty. Run: tcrtrie-vdjdb-update --tag <tag> (or --list).")

    airr = get_cached_airr_path(root, tag)
    sqlite = get_cached_sqlite_path(root, tag)
    if not (airr.exists() and sqlite.exists()):
        raise RuntimeError(
            f"VDJdb cache for tag '{tag}' is incomplete. Run: tcrtrie-vdjdb-update --tag {tag}"
        )

    return VDJdbClient(trie=Trie(str(airr)), sqlitePath=sqlite)

if TYPE_CHECKING:
    VDJdb: VDJdbClient
else:
    VDJdb = cast(VDJdbClient, LazyObject(_build_vdjdb))
