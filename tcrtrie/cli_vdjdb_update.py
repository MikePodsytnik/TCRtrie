from __future__ import annotations

import argparse
import sys

from .vdjdb_cache import (
    install_vdjdb_latest,
    install_vdjdb_tag,
    install_vdjdb_web_latest,
    list_vdjdb_releases,
)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="tcrtrie-vdjdb-update")

    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument("--list", action="store_true", help="List available VDJdb release tags")
    g.add_argument("--tag", type=str, help="Install specific VDJdb release tag (e.g. 2025-12-29)")
    g.add_argument(
        "--web",
        action="store_true",
        help="Install VDJdb version from latest-version.txt (web channel) into cache/.../vdjdb/web",
    )

    args = p.parse_args(argv)

    if args.list:
        rels = list_vdjdb_releases()
        if not rels:
            print("No releases found.", file=sys.stderr)
            return 1

        tag_w = max(len("tag"), max(len(tag) for tag, _ in rels))
        date_w = max(len("published_at"), max(len(dt) for _, dt in rels))

        print(f"{'tag':<{tag_w}}  {'published_at':<{date_w}}")
        print(f"{'-' * tag_w}  {'-' * date_w}")

        for tag, published in rels:
            print(f"{tag:<{tag_w}}  {published:<{date_w}}")

        return 0

    if args.web:
        airr = install_vdjdb_web_latest()
        print(str(airr))
        return 0

    if args.tag:
        airr = install_vdjdb_tag(tag=args.tag)
        print(str(airr))
        return 0

    airr = install_vdjdb_latest()
    print(str(airr))
    return 0
