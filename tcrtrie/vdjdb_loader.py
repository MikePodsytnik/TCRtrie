from __future__ import annotations

import csv
import json
import os
import pathlib
import tarfile
import urllib.request
import zipfile
import sqlite3
from dataclasses import dataclass


@dataclass(frozen=True)
class VDJdbRelease:
    tag: str
    asset_name: str
    asset_url: str


class VDJdbDownloadError(RuntimeError):
    pass


class VDJdbFormatError(RuntimeError):
    pass


def _http_get_json(url: str, headers: dict[str, str]) -> object:
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _download(
    url: str,
    out_path: pathlib.Path,
    headers: dict[str, str],
    *,
    on_progress=None,  # callable(downloaded_bytes: int, total_bytes: int | None)
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp, open(tmp_path, "wb") as f:
        total = resp.headers.get("Content-Length")
        total_bytes = int(total) if total and total.isdigit() else None

        downloaded = 0
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if on_progress:
                on_progress(downloaded, total_bytes)

    tmp_path.replace(out_path)


def _pick_release_asset(assets: list[dict]) -> dict:
    if not assets:
        raise VDJdbDownloadError("VDJdb GitHub release has no assets.")

    def score(a: dict) -> int:
        name = (a.get("name") or "").lower()
        if name.endswith(".zip"):
            return 3
        if name.endswith((".tar.gz", ".tgz")):
            return 2
        return 1

    return max(assets, key=score)


def _extract_archive(archive_path: pathlib.Path, dest_dir: pathlib.Path) -> None:
    name = archive_path.name.lower()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(dest_dir)
        return

    if name.endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive_path, "r:gz") as tf:
            tf.extractall(dest_dir)
        return

    raise VDJdbDownloadError(f"Unknown archive format: {archive_path.name}")


def _find_vdjdb_txt(root: pathlib.Path) -> pathlib.Path:
    matches = list(root.rglob("vdjdb.txt"))
    if not matches:
        raise VDJdbDownloadError("vdjdb.txt not found inside extracted release.")
    matches.sort(key=lambda p: len(p.parts))
    return matches[0]


def _github_headers(github_token: str | None = None) -> dict[str, str]:
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "tcrtriepy-vdjdb-loader",
    }
    token = github_token or os.getenv("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def fetch_vdjdb_releases(*, github_token: str | None = None, per_page: int = 100) -> list[dict]:
    headers = _github_headers(github_token)
    api = f"https://api.github.com/repos/antigenomics/vdjdb-db/releases?per_page={per_page}"
    rels = _http_get_json(api, headers=headers)
    if not isinstance(rels, list):
        raise VDJdbDownloadError("Unexpected response from GitHub releases endpoint.")
    return rels


def fetch_latest_vdjdb_tag(*, github_token: str | None = None) -> str:
    headers = _github_headers(github_token)
    api = "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest"
    rel = _http_get_json(api, headers=headers)
    if not isinstance(rel, dict):
        raise VDJdbDownloadError("Unexpected response from GitHub latest release endpoint.")
    tag = rel.get("tag_name")
    if not tag:
        raise VDJdbDownloadError("Cannot read tag_name from GitHub latest release response.")
    return tag


def fetch_latest_vdjdb_txt(
    *,
    cache_dir: pathlib.Path,
    version: str,
    github_token: str | None = None,
    on_progress=None,
) -> tuple[pathlib.Path, VDJdbRelease]:
    headers = _github_headers(github_token)

    api = f"https://api.github.com/repos/antigenomics/vdjdb-db/releases/tags/{version}"
    rel = _http_get_json(api, headers=headers)
    if not isinstance(rel, dict):
        raise VDJdbDownloadError("Unexpected response from GitHub tag endpoint.")

    tag = rel.get("tag_name")
    if not tag:
        raise VDJdbDownloadError("Cannot read tag_name from GitHub release response.")

    asset = _pick_release_asset(rel.get("assets") or [])
    asset_name = asset.get("name")
    asset_url = asset.get("browser_download_url")
    if not asset_name or not asset_url:
        raise VDJdbDownloadError("Chosen asset missing name or browser_download_url.")

    meta = VDJdbRelease(tag=tag, asset_name=asset_name, asset_url=asset_url)

    rel_dir = cache_dir / tag
    archive_path = rel_dir / asset_name
    extract_dir = rel_dir / "extracted"
    canonical_txt = extract_dir / "vdjdb.txt"

    if canonical_txt.exists():
        return canonical_txt, meta

    rel_dir.mkdir(parents=True, exist_ok=True)

    if not archive_path.exists():
        _download(asset_url, archive_path, headers=headers, on_progress=on_progress)

    _extract_archive(archive_path, extract_dir)
    real_txt = _find_vdjdb_txt(extract_dir)

    if real_txt.resolve() != canonical_txt.resolve():
        canonical_txt.write_bytes(real_txt.read_bytes())

    return canonical_txt, meta


def _init_vdjdb_sqlite(sqlite_path: pathlib.Path, cols: list[str]) -> sqlite3.Connection:
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(sqlite_path)
    con.execute("PRAGMA journal_mode=WAL;")
    con.execute("DROP TABLE IF EXISTS vdjdb;")
    col_defs = ", ".join([f'"{c}" TEXT' for c in cols])
    con.execute(f'CREATE TABLE vdjdb (idx INTEGER PRIMARY KEY, {col_defs});')
    return con


def _extract_airr_triplet(row: dict[str, str]) -> tuple[str, str, str] | None:
    cdr3 = (row.get("cdr3") or "").strip()
    v = (row.get("v.segm") or "").strip()
    j = (row.get("j.segm") or "").strip()

    v = v.replace('\t', '').replace('\n', '').replace('\r', '')
    j = j.replace('\t', '').replace('\n', '').replace('\r', '')

    if not (cdr3 and v and j):
        return None
    return cdr3, v, j


def vdjdb_txt_to_airr_and_sqlite(
    *,
    vdjdb_txt: pathlib.Path,
    out_airr_tsv: pathlib.Path,
    out_sqlite: pathlib.Path,
) -> None:
    out_airr_tsv.parent.mkdir(parents=True, exist_ok=True)

    with open(vdjdb_txt, "r", newline="", encoding="utf-8") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise VDJdbFormatError("vdjdb.txt has no header")

        needed = {"cdr3", "v.segm", "j.segm"}
        missing = needed - set(reader.fieldnames)
        if missing:
            raise VDJdbFormatError(f"vdjdb.txt missing columns: {sorted(missing)}")

        cols = list(reader.fieldnames)

        con = _init_vdjdb_sqlite(out_sqlite, cols)
        try:
            airr_fields = ["junction_aa", "v_call", "j_call"]
            with open(out_airr_tsv, "w", newline="", encoding="utf-8") as fout:
                w = csv.DictWriter(fout, fieldnames=airr_fields, delimiter="\t")
                w.writeheader()

                idx = 0
                col_list = ", ".join([f'"{c}"' for c in cols])
                placeholders = ",".join(["?"] * (len(cols) + 1))
                insert_sql = f'INSERT INTO vdjdb (idx, {col_list}) VALUES ({placeholders});'

                for row in reader:
                    triplet = _extract_airr_triplet(row)
                    if triplet is None:
                        continue
                    cdr3, v, j = triplet

                    w.writerow({"junction_aa": cdr3, "v_call": v, "j_call": j})
                    values = [row.get(c, "") for c in cols]
                    con.execute(insert_sql, [idx] + values)
                    idx += 1

            con.commit()
        finally:
            con.close()
