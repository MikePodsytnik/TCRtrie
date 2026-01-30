from __future__ import annotations

import pathlib
from typing import Optional

from .vdjdb_loader import (
    fetch_latest_vdjdb_tag,
    fetch_latest_vdjdb_txt,
    vdjdb_txt_to_airr_minimal,
)

import sys
import time

def _make_progress_printer(prefix: str = "Downloading"):
    last_print = 0.0

    def printer(done: int, total: int | None):
        nonlocal last_print
        now = time.time()
        if now - last_print < 0.1:
            return
        last_print = now

        if total:
            pct = done / total * 100
            bar_len = 30
            filled = int(bar_len * done / total)
            bar = "#" * filled + "-" * (bar_len - filled)
            sys.stderr.write(f"\r{prefix}: [{bar}] {pct:5.1f}% ({done/1e6:.1f}/{total/1e6:.1f} MB)")
        else:
            sys.stderr.write(f"\r{prefix}: {done/1e6:.1f} MB")
        sys.stderr.flush()

    return printer


def cache_root() -> pathlib.Path:
    return pathlib.Path.home() / ".cache" / "tcrtrie" / "vdjdb"


def _latest_marker(root: pathlib.Path) -> pathlib.Path:
    return root / "LATEST"


def read_cached_latest_tag(root: pathlib.Path) -> Optional[str]:
    try:
        return _latest_marker(root).read_text(encoding="utf-8").strip() or None
    except Exception:
        return None


def write_cached_latest_tag(root: pathlib.Path, tag: str) -> None:
    root.mkdir(parents=True, exist_ok=True)
    _latest_marker(root).write_text(tag + "\n", encoding="utf-8")


def get_cached_airr_path(root: pathlib.Path, tag: str) -> pathlib.Path:
    return root / tag / "vdjdb_airr.tsv"


def warmup_vdjdb_cache(*, root: pathlib.Path | None = None) -> pathlib.Path:
    """
    Ensures cache contains the latest VDJdb release (by GitHub latest tag).
    Returns path to cached AIRR TSV.
    """
    root = root or cache_root()

    tag = fetch_latest_vdjdb_tag()
    airr = get_cached_airr_path(root, tag)
    if airr.exists():
        write_cached_latest_tag(root, tag)
        return airr

    progress = _make_progress_printer(prefix=f"Downloading VDJdb {tag}")
    vdjdb_txt, _ = fetch_latest_vdjdb_txt(cache_dir=root, version=tag, on_progress=progress)
    sys.stderr.write("\n")
    airr.parent.mkdir(parents=True, exist_ok=True)
    vdjdb_txt_to_airr_minimal(vdjdb_txt=vdjdb_txt, out_airr_tsv=airr)

    write_cached_latest_tag(root, tag)
    return airr
