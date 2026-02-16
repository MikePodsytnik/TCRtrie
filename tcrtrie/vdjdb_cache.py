from __future__ import annotations

import pathlib
from typing import Optional
import sys
import time

from .vdjdb_loader import (
    fetch_vdjdb_releases,
    fetch_latest_vdjdb_tag,
    fetch_latest_vdjdb_txt,
    fetch_vdjdb_txt_from_asset_url,
    _http_get_text,
    vdjdb_txt_to_airr_and_sqlite,
)


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
            sys.stderr.write(
                f"\r{prefix}: [{bar}] {pct:5.1f}% ({done/1e6:.1f}/{total/1e6:.1f} MB)"
            )
        else:
            sys.stderr.write(f"\r{prefix}: {done/1e6:.1f} MB")
        sys.stderr.flush()

    return printer


def cache_root() -> pathlib.Path:
    return pathlib.Path.home() / ".cache" / "tcrtrie" / "vdjdb"


def _active_marker(root: pathlib.Path) -> pathlib.Path:
    return root / "ACTIVE"


def read_cached_active_tag(root: pathlib.Path) -> Optional[str]:
    try:
        return _active_marker(root).read_text(encoding="utf-8").strip() or None
    except Exception:
        return None


def write_cached_active_tag(root: pathlib.Path, tag: str) -> None:
    root.mkdir(parents=True, exist_ok=True)
    _active_marker(root).write_text(tag + "\n", encoding="utf-8")


def get_cached_airr_path(root: pathlib.Path, tag: str) -> pathlib.Path:
    return root / tag / "vdjdb_airr.tsv"


def get_cached_sqlite_path(root: pathlib.Path, tag: str) -> pathlib.Path:
    return root / tag / "vdjdb.sqlite"


def list_vdjdb_releases() -> list[tuple[str, str]]:
    rels = fetch_vdjdb_releases()
    out: list[tuple[str, str]] = []
    for r in rels:
        tag = r.get("tag_name")
        published = r.get("published_at") or ""
        if tag:
            out.append((tag, published))
    return out


def install_vdjdb_tag(*, tag: str, root: pathlib.Path | None = None) -> pathlib.Path:
    root = root or cache_root()

    airr = get_cached_airr_path(root, tag)
    sqlite = get_cached_sqlite_path(root, tag)
    if airr.exists() and sqlite.exists():
        write_cached_active_tag(root, tag)
        return airr

    progress = _make_progress_printer(prefix=f"Downloading VDJdb {tag}")
    vdjdb_txt, _ = fetch_latest_vdjdb_txt(cache_dir=root, version=tag, on_progress=progress)
    sys.stderr.write("\n")

    vdjdb_txt_to_airr_and_sqlite(vdjdb_txt=vdjdb_txt, out_airr_tsv=airr, out_sqlite=sqlite)
    write_cached_active_tag(root, tag)
    return airr


def install_vdjdb_latest(*, root: pathlib.Path | None = None) -> pathlib.Path:
    root = root or cache_root()
    tag = fetch_latest_vdjdb_tag()
    return install_vdjdb_tag(tag=tag, root=root)


def _web_dir(root: pathlib.Path) -> pathlib.Path:
    return root / "web"


def _clear_dir(path: pathlib.Path) -> None:
    if not path.exists():
        return
    for p in sorted(path.rglob("*"), reverse=True):
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
            elif p.is_dir():
                p.rmdir()
        except Exception:
            pass
    try:
        path.rmdir()
    except Exception:
        pass


def _read_web_latest_url() -> str:
    url = "https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt"
    headers = {"User-Agent": "tcrtriepy-vdjdb-loader"}
    text = _http_get_text(url, headers=headers)

    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            continue
        return s

    raise RuntimeError("latest-version.txt is empty or contains no valid URLs.")


def install_vdjdb_web_latest(*, root: pathlib.Path | None = None) -> pathlib.Path:
    root = root or cache_root()

    web_dir = _web_dir(root)
    _clear_dir(web_dir)
    web_dir.mkdir(parents=True, exist_ok=True)

    asset_url = _read_web_latest_url()

    progress = _make_progress_printer(prefix="Downloading VDJdb (web)")
    vdjdb_txt = fetch_vdjdb_txt_from_asset_url(cache_dir=web_dir, asset_url=asset_url, on_progress=progress)
    sys.stderr.write("\n")

    airr = web_dir / "vdjdb_airr.tsv"
    sqlite = web_dir / "vdjdb.sqlite"
    vdjdb_txt_to_airr_and_sqlite(vdjdb_txt=vdjdb_txt, out_airr_tsv=airr, out_sqlite=sqlite)

    write_cached_active_tag(root, "web")
    return airr

